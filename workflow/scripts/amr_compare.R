# =============================================================================
# amr_compare.R
# Compare AMR gene burden (AMRFinderPlus) between CON and SCFP using corrected
# Wisconsin metadata (n=13 with RUM_KSU_42 included).
# =============================================================================

suppressPackageStartupMessages({
    library(data.table)
})

DC_FILE  <- "/mnt/d/Wisconsin_data1/amr_bagel_results_hpc/merged/amr_drugclass_matrix.tsv"
HIT_FILE <- "/mnt/d/Wisconsin_data1/amr_bagel_results_hpc/merged/all_amrfinder_hits.tsv"
OUT_DIR  <- "results_ko/amr_compare"
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

# Corrected metadata (n=13, RUM_KSU_42 in CON)
meta <- data.table(
    sample_id = c("RUM_KSU_17","RUM_KSU_18","RUM_KSU_19","RUM_KSU_20","RUM_KSU_21",
                  "RUM_KSU_23","RUM_KSU_42","RUM_KSU_43","RUM_KSU_44","RUM_KSU_45",
                  "RUM_KSU_46","RUM_KSU_47","RUM_KSU_48"),
    treatment = c("CON","SCFP","SCFP","CON","CON",
                  "SCFP","CON","CON","SCFP","SCFP",
                  "CON","CON","SCFP")
)
cat("Groups:\n"); print(table(meta$treatment))

# ---- 1. Drug class matrix --------------------------------------------------
dc <- fread(DC_FILE)
dc <- merge(meta, dc, by.x = "sample_id", by.y = "Sample")

cat("\n=== Total AMR hits per sample ===\n")
class_cols <- setdiff(colnames(dc), c("sample_id", "treatment"))
dc[, total_hits := rowSums(.SD), .SDcols = class_cols]
print(dc[, .(sample_id, treatment, total_hits)])

cat("\n=== Mean total hits by group ===\n")
print(dc[, .(mean = mean(total_hits), sd = sd(total_hits), n = .N), by = treatment])

# Wilcoxon test on total hits
w_total <- wilcox.test(total_hits ~ treatment, data = dc)
cat("\nWilcoxon total AMR hits CON vs SCFP: W =", w_total$statistic,
    " p =", signif(w_total$p.value, 4), "\n")

# ---- 2. Per drug class test ------------------------------------------------
cat("\n=== Per drug class Wilcoxon (CON vs SCFP) ===\n")
class_stats <- rbindlist(lapply(class_cols, function(cl) {
    con  <- dc[treatment == "CON",  get(cl)]
    scfp <- dc[treatment == "SCFP", get(cl)]
    if (sum(con) + sum(scfp) == 0) return(NULL)
    w <- suppressWarnings(wilcox.test(con, scfp))
    data.table(
        drug_class = cl,
        CON_mean   = mean(con),
        SCFP_mean  = mean(scfp),
        diff       = mean(con) - mean(scfp),
        W          = unname(w$statistic),
        p          = w$p.value
    )
}))
class_stats[, q := p.adjust(p, method = "BH")]
setorder(class_stats, p)
print(class_stats)
fwrite(class_stats, file.path(OUT_DIR, "drugclass_compare.tsv"), sep = "\t")

# ---- 3. Gene-level AMRFinderPlus ------------------------------------------
hits <- fread(HIT_FILE)
setnames(hits, names(hits), make.names(names(hits)))
cat("\n=== AMRFinderPlus hits columns ===\n")
print(colnames(hits))
cat("Rows:", nrow(hits), "\n")

# Find the gene-symbol column and sample column robustly
gene_col <- intersect(c("Gene.symbol", "gene_symbol", "Element.symbol"), colnames(hits))[1]
samp_col <- intersect(c("Sample", "sample", "Name"), colnames(hits))[1]
cat("Using gene column:", gene_col, " sample column:", samp_col, "\n")

# Per sample per gene count
gene_long <- hits[, .N, by = c(samp_col, gene_col)]
setnames(gene_long, c("sample_id", "gene", "n_hits"))
gene_wide <- dcast(gene_long, sample_id ~ gene, value.var = "n_hits", fill = 0)
gene_wide <- merge(meta, gene_wide, by = "sample_id")
gene_cols <- setdiff(colnames(gene_wide), c("sample_id", "treatment"))

cat("\n=== Per gene Wilcoxon (top 20 by p) ===\n")
gene_stats <- rbindlist(lapply(gene_cols, function(g) {
    con  <- gene_wide[treatment == "CON",  get(g)]
    scfp <- gene_wide[treatment == "SCFP", get(g)]
    if (sum(con) + sum(scfp) < 2) return(NULL)
    w <- suppressWarnings(wilcox.test(con, scfp))
    data.table(
        gene      = g,
        CON_mean  = mean(con),
        SCFP_mean = mean(scfp),
        diff      = mean(con) - mean(scfp),
        n_samples_with_hit = sum(con > 0) + sum(scfp > 0),
        W         = unname(w$statistic),
        p         = w$p.value
    )
}))
gene_stats[, q := p.adjust(p, method = "BH")]
setorder(gene_stats, p)
print(head(gene_stats, 20))
fwrite(gene_stats, file.path(OUT_DIR, "gene_compare.tsv"), sep = "\t")

# ---- 4. Aminoglycoside subset (cross-ref MEsalmon finding) ----------------
cat("\n=== Aminoglycoside genes (cross-reference MEsalmon aacA-aphD) ===\n")
amg <- gene_stats[grepl("aac|aph|ant|str|aad|kan", gene, ignore.case = TRUE)]
print(amg)
fwrite(amg, file.path(OUT_DIR, "aminoglycoside_genes.tsv"), sep = "\t")

cat("\nOutputs:\n"); print(list.files(OUT_DIR, full.names = TRUE))
