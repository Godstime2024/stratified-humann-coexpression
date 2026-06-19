# =============================================================================
# procrustes.R
# Procrustes test of concordance between taxonomic and functional sample
# distance structures (compositional / Aitchison).
# =============================================================================

suppressPackageStartupMessages({
    library(tidyverse)
    library(data.table)
    library(vegan)
})

log_con <- file(snakemake@log[[1]], open = "wt")
sink(log_con, type = "output"); sink(log_con, type = "message")

taxa_path  <- snakemake@input[["taxa_clr"]]
unstrat_path <- snakemake@input[["unstrat_clr"]]
meta_path  <- snakemake@input[["meta"]]
out_res    <- snakemake@output[["result"]]
out_plot   <- snakemake@output[["plot"]]
nperm      <- as.integer(snakemake@params[["permutations"]])
group_col  <- snakemake@params[["group_col"]]

dir.create(dirname(out_res), recursive = TRUE, showWarnings = FALSE)

taxa <- fread(taxa_path)
func <- fread(unstrat_path)
meta <- fread(meta_path)

sample_cols <- intersect(colnames(taxa), colnames(func))
sample_cols <- setdiff(sample_cols, c("taxon", "feature"))

taxa_mat <- t(as.matrix(taxa[, ..sample_cols]))
func_mat <- t(as.matrix(func[, ..sample_cols]))

# Aitchison distance = Euclidean on CLR data
taxa_d <- dist(taxa_mat, method = "euclidean")
func_d <- dist(func_mat, method = "euclidean")

mds_t <- cmdscale(taxa_d, k = 2)
mds_f <- cmdscale(func_d, k = 2)

proc <- procrustes(mds_t, mds_f)
ptest <- protest(mds_t, mds_f, permutations = nperm)

# Tabulate result
res <- data.table(
    metric = c("m12_squared", "correlation", "significance",
               "n_samples", "permutations"),
    value  = c(ptest$ss, ptest$t0, ptest$signif, nrow(taxa_mat), nperm)
)
fwrite(res, out_res, sep = "\t")

# Plot
pdf(out_plot, width = 7, height = 7)
meta_align <- meta[match(rownames(taxa_mat), meta[[1]]), ]
grp <- meta_align[[group_col]]
cols <- as.numeric(factor(grp))
plot(proc, main = sprintf("Procrustes: taxa vs functions (m2=%.3f, p=%.3f)",
                          ptest$ss, ptest$signif))
points(proc$X, col = cols, pch = 19, cex = 1.2)
points(proc$Yrot, col = cols, pch = 1, cex = 1.2)
legend("topright", legend = levels(factor(grp)),
       col = seq_along(levels(factor(grp))), pch = 19, bty = "n")
dev.off()

cat("[procrustes] m2 =", ptest$ss, "p =", ptest$signif, "\n")
cat("[procrustes] Done.\n")
sink(type = "message"); sink(type = "output"); close(log_con)
