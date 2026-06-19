# =============================================================================
# kegg_enrichment.R
# KEGG pathway enrichment per WGCNA module using KO IDs.
# Skips gracefully if input features are not KO IDs (e.g., COG categories).
# =============================================================================

suppressPackageStartupMessages({
    library(tidyverse)
    library(data.table)
    library(KEGGREST)
    library(clusterProfiler)
    library(pheatmap)
})

log_con <- file(snakemake@log[[1]], open = "wt")
sink(log_con, type = "output"); sink(log_con, type = "message")

modules_path <- snakemake@input[["modules"]]
out_table    <- snakemake@output[["table"]]
out_summary  <- snakemake@output[["summary"]]
out_heatmap  <- snakemake@output[["heatmap"]]
cache_dir    <- snakemake@params[["cache_dir"]]
top_n_paths  <- as.integer(snakemake@params[["top_n_paths"]])
pval_cutoff  <- as.numeric(snakemake@params[["pval_cutoff"]])

dir.create(dirname(out_table), recursive = TRUE, showWarnings = FALSE)
dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)

# ---- Load module assignments -----------------------------------------------
modules <- fread(modules_path)
cat("[kegg] Modules loaded:", nrow(modules), "features\n")

# ---- Detect if features are KO IDs -----------------------------------------
ko_pattern <- "^K\\d{5}$"
is_ko <- grepl(ko_pattern, modules$feature)
n_ko <- sum(is_ko)
cat("[kegg] Features matching KO pattern (K#####):", n_ko, "of", nrow(modules), "\n")

if (n_ko < 10) {
    cat("[kegg] Not enough KO features (<10). Skipping KEGG enrichment.\n")
    cat("[kegg] (Input is likely COG or non-KO data.)\n")
    fwrite(data.table(message = "KEGG enrichment skipped: features are not KO IDs"),
           out_table, sep = "\t")
    fwrite(data.table(message = "KEGG enrichment skipped: features are not KO IDs"),
           out_summary, sep = "\t")
    pdf(out_heatmap, width = 6, height = 4)
    plot.new(); title("KEGG enrichment skipped: features are not KO IDs")
    dev.off()
    sink(type = "message"); sink(type = "output"); close(log_con)
    quit(status = 0)
}

# ---- Download/load KO-to-pathway mapping -----------------------------------
ko2path_cache <- file.path(cache_dir, "ko_to_pathway.rds")
pathnames_cache <- file.path(cache_dir, "pathway_names.rds")

konames_cache <- file.path(cache_dir, "ko_names.rds")

if (file.exists(ko2path_cache) && file.exists(pathnames_cache) && file.exists(konames_cache)) {
    cat("[kegg] Loading cached KEGG mapping...\n")
    ko2path <- readRDS(ko2path_cache)
    path_names <- readRDS(pathnames_cache)
    ko_names <- readRDS(konames_cache)
} else {
    cat("[kegg] Downloading KEGG mappings from REST API...\n")
    cat("[kegg] (This may take 2-3 minutes; cached for future runs.)\n")

    # KO -> pathway links
    raw_link <- keggLink("pathway", "ko")
    ko2path <- data.table(
        ko = sub("ko:", "", names(raw_link)),
        pathway = sub("path:", "", raw_link)
    )
    ko2path <- ko2path[grepl("^ko\\d+$", pathway)]
    saveRDS(ko2path, ko2path_cache)

    # Pathway names
    raw_names <- keggList("pathway", "ko")
    path_names <- data.table(
        pathway = sub("path:", "", names(raw_names)),
        name = unname(raw_names)
    )
    saveRDS(path_names, pathnames_cache)

    # KO names (human-readable)
    raw_ko <- keggList("ko")
    ko_names <- data.table(
        ko = sub("ko:", "", names(raw_ko)),
        ko_name = unname(raw_ko)
    )
    # Also write a TSV companion so downstream scripts can read it without RDS
    saveRDS(ko_names, konames_cache)
    fwrite(ko_names, file.path(cache_dir, "ko_names.tsv"), sep = "\t")

    cat("[kegg] Cached", nrow(ko2path), "KO-pathway links and",
        nrow(ko_names), "KO names\n")
}

# Build TERM2GENE table for enricher
t2g <- data.frame(term = ko2path$pathway, gene = ko2path$ko)
universe <- unique(modules$feature[is_ko])

# ---- Enrichment per module -------------------------------------------------
mod_levels <- unique(modules$module)
mod_levels <- mod_levels[mod_levels != "grey"]  # grey = unassigned
cat("[kegg] Running enrichment on", length(mod_levels), "modules\n")

all_results <- list()
for (m in mod_levels) {
    mod_kos <- modules$feature[modules$module == m & is_ko]
    if (length(mod_kos) < 5) next

    enr <- tryCatch(
        enricher(gene = mod_kos, TERM2GENE = t2g, universe = universe,
                 pvalueCutoff = 1, qvalueCutoff = 1,
                 minGSSize = 3, maxGSSize = 500),
        error = function(e) { cat("[kegg]", m, "failed:", conditionMessage(e), "\n"); NULL }
    )
    if (is.null(enr) || nrow(enr@result) == 0) next

    res <- as.data.table(enr@result)
    res[, module := m]
    res[, n_module_kos := length(mod_kos)]
    all_results[[m]] <- res
}

if (length(all_results) == 0) {
    cat("[kegg] No enrichment results returned. Writing empty outputs.\n")
    fwrite(data.table(message = "No significant enrichment results"),
           out_table, sep = "\t")
    fwrite(data.table(message = "No significant enrichment results"),
           out_summary, sep = "\t")
    pdf(out_heatmap, width = 6, height = 4)
    plot.new(); title("No KEGG enrichment results")
    dev.off()
    sink(type = "message"); sink(type = "output"); close(log_con)
    quit(status = 0)
}

results <- rbindlist(all_results, fill = TRUE)

# Add pathway names
results <- merge(results, path_names, by.x = "ID", by.y = "pathway", all.x = TRUE)
setcolorder(results, c("module", "ID", "name", "Description",
                       "GeneRatio", "BgRatio", "pvalue", "p.adjust",
                       "qvalue", "Count", "n_module_kos", "geneID"))

fwrite(results, out_table, sep = "\t")
cat("[kegg] Full results written:", nrow(results), "rows\n")

# ---- Summary: top pathways per module --------------------------------------
sig <- results[p.adjust < pval_cutoff]
top_per_mod <- sig[order(p.adjust),
                   .SD[1:min(.N, top_n_paths)],
                   by = module]
fwrite(top_per_mod, out_summary, sep = "\t")
cat("[kegg] Significant (FDR <", pval_cutoff, "):", nrow(sig), "module-pathway pairs\n")

# ---- Heatmap: -log10(p.adjust) for top pathways x modules -------------------
if (nrow(sig) > 0) {
    sig[, log_padj := -log10(p.adjust + 1e-10)]
    # Pick top pathways across all modules
    path_max <- sig[, .(max_log = max(log_padj)), by = .(ID, name)][
        order(-max_log)][1:min(.N, 30)]
    heat_dt <- sig[ID %in% path_max$ID,
                   .(module, ID, name, log_padj)]

    # Pivot to matrix
    heat_wide <- dcast(heat_dt, name ~ module,
                       value.var = "log_padj", fill = 0)
    heat_mat <- as.matrix(heat_wide[, -1])
    rownames(heat_mat) <- heat_wide$name

    pdf(out_heatmap, width = max(8, ncol(heat_mat) * 0.6 + 4),
        height = max(6, nrow(heat_mat) * 0.25 + 2))
    pheatmap(heat_mat,
             color = colorRampPalette(c("white", "#fd8d3c", "#bd0026"))(50),
             cluster_rows = TRUE, cluster_cols = TRUE,
             main = "KEGG pathway enrichment per WGCNA module (-log10 FDR)",
             fontsize_row = 7, fontsize_col = 8,
             angle_col = 45)
    dev.off()
    cat("[kegg] Heatmap saved\n")
} else {
    pdf(out_heatmap, width = 6, height = 4)
    plot.new(); title("No pathways pass FDR cutoff")
    dev.off()
}

cat("[kegg] Done.\n")
sink(type = "message"); sink(type = "output"); close(log_con)
