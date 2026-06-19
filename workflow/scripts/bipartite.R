# =============================================================================
# bipartite.R
# Build a bipartite taxon-function network from stratified data and compute
# network-level + species-level specialization metrics.
# =============================================================================

suppressPackageStartupMessages({
    library(tidyverse)
    library(data.table)
    library(bipartite)
    library(pheatmap)
})

log_con <- file(snakemake@log[[1]], open = "wt")
sink(log_con, type = "output"); sink(log_con, type = "message")

strat_path <- snakemake@input[["strat"]]
out_mat    <- snakemake@output[["matrix"]]
out_met    <- snakemake@output[["metrics"]]
out_plot   <- snakemake@output[["plot"]]
top_t      <- as.integer(snakemake@params[["top_taxa"]])
top_f      <- as.integer(snakemake@params[["top_func"]])

dir.create(dirname(out_mat), recursive = TRUE, showWarnings = FALSE)

strat <- fread(strat_path)
sample_cols <- setdiff(colnames(strat),
                       c("feature", "function_name", "taxon"))

# Total abundance per (taxon, function) summed across samples
agg <- strat[, .(total = rowSums(.SD, na.rm = TRUE)),
             .SDcols = sample_cols,
             by = .(taxon, function_name)]

# Pivot to taxon x function matrix
bip_dt <- dcast(agg, taxon ~ function_name,
                value.var = "total", fill = 0)
mat <- as.matrix(bip_dt[, -1])
rownames(mat) <- bip_dt$taxon

cat("[bipartite] Full bipartite matrix:", dim(mat), "(taxa x functions)\n")

# Trim to top N taxa and top N functions by marginal abundance
taxa_totals <- rowSums(mat)
func_totals <- colSums(mat)
keep_taxa <- names(sort(taxa_totals, decreasing = TRUE))[1:min(top_t, nrow(mat))]
keep_func <- names(sort(func_totals, decreasing = TRUE))[1:min(top_f, ncol(mat))]
mat_top <- mat[keep_taxa, keep_func, drop = FALSE]
cat("[bipartite] Trimmed to:", dim(mat_top), "for visualization\n")

fwrite(as.data.table(mat, keep.rownames = "taxon"), out_mat, sep = "\t")

# ---- Specialization metrics ------------------------------------------------
# Compute on the trimmed top matrix to keep runtime tractable. H2fun on a
# full thousands x hundreds matrix is combinatorially expensive when values
# are not integers.
mat_safe <- mat_top[rowSums(mat_top) > 0, colSums(mat_top) > 0, drop = FALSE]
# Scale to integers for speed and to satisfy H2fun's integer requirement.
mat_int <- round(mat_safe * 1e6)
h2 <- tryCatch(H2fun(mat_int),
               error = function(e) { warning("H2fun failed: ", conditionMessage(e)); c(NA) })

sl_higher <- tryCatch(specieslevel(mat_int, level = "higher"),
                      error = function(e) NULL)
sl_lower  <- tryCatch(specieslevel(mat_int, level = "lower"),
                      error = function(e) NULL)

if (!is.null(sl_higher) && !is.null(sl_lower)) {
    metrics_dt <- rbind(
        cbind(level = "function", node = rownames(sl_higher), as.data.table(sl_higher)),
        cbind(level = "taxon", node = rownames(sl_lower), as.data.table(sl_lower)),
        fill = TRUE
    )
    fwrite(metrics_dt, out_met, sep = "\t")
}
cat("[bipartite] H2 specialization:", h2[1], "\n")

# ---- Plot bipartite network ------------------------------------------------
# plotweb's argument signature varies across bipartite versions. Use a defensive
# call that tries the full-featured form first, falls back progressively, and
# uses a heatmap if all plotweb calls fail.
pdf(out_plot, width = 12, height = 8)
plot_done <- FALSE

# Attempt 1: full-featured (modern bipartite)
plot_done <- tryCatch({
    bipartite::plotweb(mat_top,
                        method = "cca",
                        text.rot = 90,
                        col.high = "steelblue",
                        col.low = "tomato",
                        bor.col.interaction = NA,
                        labsize = 0.7,
                        ybig = 1.2)
    TRUE
}, error = function(e) {
    cat("[bipartite] full plotweb failed:", conditionMessage(e), "\n")
    FALSE
})

# Attempt 2: minimal plotweb call (just the matrix)
if (!plot_done) {
    plot_done <- tryCatch({
        bipartite::plotweb(mat_top)
        TRUE
    }, error = function(e) {
        cat("[bipartite] minimal plotweb failed:", conditionMessage(e), "\n")
        FALSE
    })
}

# Attempt 3: heatmap fallback (always works)
if (!plot_done) {
    cat("[bipartite] falling back to heatmap visualization\n")
    pheatmap::pheatmap(
        log10(mat_top + 1e-12),
        cluster_rows = TRUE,
        cluster_cols = TRUE,
        main = "Bipartite taxon-function heatmap (log10 abundance)",
        fontsize_row = 7, fontsize_col = 7,
        color = colorRampPalette(c("white", "tomato", "darkred"))(50)
    )
} else {
    title(main = "Top taxa-function bipartite network",
          sub = paste("Top", top_t, "taxa x top", top_f, "functions; H2 =",
                      signif(h2[1], 3)),
          cex.main = 1.1)
}
dev.off()

cat("[bipartite] Done.\n")
sink(type = "message"); sink(type = "output"); close(log_con)
