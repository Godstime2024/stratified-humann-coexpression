# =============================================================================
# metacyc_compare.R
# CON vs SCFP test on HUMAnN3 MetaCyc pathway relative abundance.
# Independent ontology from KEGG.
# =============================================================================
suppressPackageStartupMessages({ library(data.table) })

PA_FILE <- "/mnt/d/Wisconsin_data1/merged_results/subset/humann3_pathabundance_relab_subset.tsv"
OUT_DIR <- "results_ko/metacyc"
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

meta <- data.table(
    sample_id = c("RUM_KSU_17","RUM_KSU_18","RUM_KSU_19","RUM_KSU_20","RUM_KSU_21",
                  "RUM_KSU_23","RUM_KSU_42","RUM_KSU_43","RUM_KSU_44","RUM_KSU_45",
                  "RUM_KSU_46","RUM_KSU_47","RUM_KSU_48"),
    treatment = c("CON","SCFP","SCFP","CON","CON",
                  "SCFP","CON","CON","SCFP","SCFP",
                  "CON","CON","SCFP")
)

pa <- fread(PA_FILE)
setnames(pa, "# Pathway", "pathway")
# Trim "_genes_Abundance" suffix from sample columns
samp_cols <- grep("_Abundance$", colnames(pa), value = TRUE)
new_names <- sub("_genes_Abundance$", "", samp_cols)
setnames(pa, samp_cols, new_names)
cat("Pathways:", nrow(pa), "  Samples:", length(new_names), "\n")

# Drop UNMAPPED/UNINTEGRATED rows
pa <- pa[!pathway %in% c("UNMAPPED", "UNINTEGRATED")]
# Drop stratified rows (pipe in name)
pa_top <- pa[!grepl("\\|", pathway)]
cat("Non-stratified MetaCyc pathways:", nrow(pa_top), "\n")

# Match metadata
samp_use <- intersect(new_names, meta$sample_id)
treat <- meta$treatment[match(samp_use, meta$sample_id)]
cat("Samples in test:", length(samp_use), " (",
    sum(treat == "CON"), "CON,", sum(treat == "SCFP"), "SCFP)\n")

pa_mat <- as.matrix(pa_top[, ..samp_use])
rownames(pa_mat) <- pa_top$pathway

# Per-pathway Wilcoxon test
results <- rbindlist(lapply(seq_len(nrow(pa_mat)), function(i) {
    v <- pa_mat[i, ]
    con  <- v[treat == "CON"]
    scfp <- v[treat == "SCFP"]
    if ((sum(con > 0) + sum(scfp > 0)) < 4) return(NULL)
    wt <- suppressWarnings(wilcox.test(con, scfp))
    data.table(
        pathway = rownames(pa_mat)[i],
        CON_mean = signif(mean(con), 3),
        SCFP_mean = signif(mean(scfp), 3),
        log2_fc = signif(log2((mean(scfp) + 1e-9) / (mean(con) + 1e-9)), 3),
        W = unname(wt$statistic),
        p = wt$p.value
    )
}))
results[, q := p.adjust(p, method = "BH")]
setorder(results, p)
fwrite(results, file.path(OUT_DIR, "metacyc_all_pathways.tsv"), sep = "\t")
cat("\n=== Top 25 MetaCyc pathways by raw p ===\n")
print(head(results, 25))

# Methanogenesis-related MetaCyc pathways
patterns <- c("methan", "PWY-5198", "PWY-5247", "PWY-5266", "PWY-5305", "PWY-5469",
              "METHANOGENESIS", "PWY-5677", "acetyl-CoA", "acetate", "PWY-1042",
              "WOOD", "HYDROGEN", "propionate", "PROPIONATE", "PWY-7117",
              "PWY-5494", "P162-PWY", "P122-PWY", "PWY-5747")
methano <- results[grepl(paste(patterns, collapse = "|"), pathway, ignore.case = TRUE)]
setorder(methano, p)
fwrite(methano, file.path(OUT_DIR, "metacyc_methanogenesis_pathways.tsv"), sep = "\t")
cat("\n=== Methanogenesis / acetogenesis / propionate MetaCyc pathways ===\n")
print(methano)

cat("\nDone. Outputs in:", OUT_DIR, "\n")
