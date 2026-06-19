# =============================================================================
# metacyc_lda.R
# LEfSe-style LDA effect size biomarker discovery on MetaCyc pathway abundance.
#
# Algorithm (Segata et al. 2011):
#   1. Per pathway Kruskal-Wallis between CON and SCFP
#   2. For pathways with K-W p < threshold, compute LDA effect size:
#      LDA effect = sign(mean_SCFP - mean_CON) * log10(1 + 1e6 * |mean diff|)
#      (Scaling factor 1e6 brings relative-abundance scale into a usable range,
#       per the original LEfSe formulation for relative-abundance data.)
#   3. Signed LDA scores: negative = biomarker for CON, positive = biomarker for SCFP
# =============================================================================
suppressPackageStartupMessages({
    library(data.table)
    library(MASS)
})

PA_FILE <- "/mnt/d/Wisconsin_data1/merged_results/subset/humann3_pathabundance_relab_subset.tsv"
OUT_DIR <- "results_ko/metacyc_lda"
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

# Thresholds:
# - LEfSe default: KW p < 0.05 AND |LDA score| > 2.0 (on count data)
# - For our small-n compositional data we relax to KW p < 0.10 and ignore
#   the 2.0 LDA cutoff (which is calibrated for OTU count abundances).
# - Top 15 biomarkers by |LDA score| are reported.
KW_PVAL_THRESHOLD  <- 0.10
TOP_N_BIOMARKERS   <- 15
SCALE_FACTOR       <- 1e6    # relative-abundance to ppm scaling

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
samp_cols <- grep("_Abundance$", colnames(pa), value = TRUE)
new_names <- sub("_genes_Abundance$", "", samp_cols)
setnames(pa, samp_cols, new_names)
pa <- pa[!pathway %in% c("UNMAPPED", "UNINTEGRATED")]
pa_top <- pa[!grepl("\\|", pathway)]

samp_use <- intersect(new_names, meta$sample_id)
treat <- meta$treatment[match(samp_use, meta$sample_id)]

pa_mat <- as.matrix(pa_top[, ..samp_use])
rownames(pa_mat) <- pa_top$pathway

cat("Pathways tested:", nrow(pa_mat), "  Samples:", ncol(pa_mat),
    " (", sum(treat == "CON"), "CON,", sum(treat == "SCFP"), "SCFP)\n\n")

# LEfSe-style scoring
lefse_results <- rbindlist(lapply(seq_len(nrow(pa_mat)), function(i) {
    v <- pa_mat[i, ]
    con  <- v[treat == "CON"]
    scfp <- v[treat == "SCFP"]
    if ((sum(con > 0) + sum(scfp > 0)) < 4) return(NULL)
    # Step 1: Kruskal-Wallis (equivalent to Wilcoxon for 2 groups)
    kw <- suppressWarnings(kruskal.test(v ~ factor(treat)))
    if (is.na(kw$p.value)) return(NULL)
    # Step 2: signed LDA effect size
    mean_con  <- mean(con)
    mean_scfp <- mean(scfp)
    mean_diff <- mean_scfp - mean_con
    # LEfSe-style log10 transformation of scaled mean difference
    lda_score <- sign(mean_diff) * log10(1 + SCALE_FACTOR * abs(mean_diff))
    data.table(
        pathway   = rownames(pa_mat)[i],
        CON_mean  = mean_con,
        SCFP_mean = mean_scfp,
        log2_fc   = log2((mean_scfp + 1e-12) / (mean_con + 1e-12)),
        kw_p      = kw$p.value,
        lda_score = lda_score,
        biomarker_class = ifelse(lda_score > 0, "SCFP", "CON")
    )
}))

# FDR adjusted KW p
lefse_results[, kw_q := p.adjust(kw_p, method = "BH")]
setorder(lefse_results, kw_p)

# LEfSe-style biomarker selection adapted for small n:
# Step 1: Kruskal-Wallis filter at p < 0.10
# Step 2: From remaining, take top N by absolute LDA score
candidates <- lefse_results[kw_p < KW_PVAL_THRESHOLD]
cat("Pathways passing KW p <", KW_PVAL_THRESHOLD, ":", nrow(candidates), "\n")
candidates[, abs_lda := abs(lda_score)]
setorder(candidates, -abs_lda)
biomarkers <- head(candidates, TOP_N_BIOMARKERS)
setorder(biomarkers, lda_score)
cat("\n=== Top", TOP_N_BIOMARKERS, "MetaCyc biomarkers (LEfSe-style) ===\n")
print(biomarkers[, .(pathway, CON_mean, SCFP_mean, log2_fc, kw_p,
                     lda_score, biomarker_class)])
cat("\nBiomarkers by class: CON =",
    sum(biomarkers$biomarker_class == "CON"),
    " SCFP =", sum(biomarkers$biomarker_class == "SCFP"), "\n")

fwrite(lefse_results, file.path(OUT_DIR, "metacyc_lda_all.tsv"), sep = "\t")
fwrite(biomarkers, file.path(OUT_DIR, "metacyc_lda_biomarkers.tsv"), sep = "\t")

# Permissive view: top 15 by absolute LDA score, regardless of KW p (gives CON
# biomarkers a chance to surface for the visualization).
permissive <- copy(lefse_results)
permissive[, abs_lda := abs(lda_score)]
setorder(permissive, -abs_lda)
permissive_top <- head(permissive, TOP_N_BIOMARKERS)
setorder(permissive_top, lda_score)
cat("\n=== Permissive top", TOP_N_BIOMARKERS, "by |LDA| (no KW filter) ===\n")
print(permissive_top[, .(pathway, CON_mean, SCFP_mean, kw_p,
                          lda_score, biomarker_class)])
cat("Permissive: CON =", sum(permissive_top$biomarker_class == "CON"),
    "  SCFP =", sum(permissive_top$biomarker_class == "SCFP"), "\n")
fwrite(permissive_top,
       file.path(OUT_DIR, "metacyc_lda_top15_permissive.tsv"), sep = "\t")

cat("\nDone. Outputs in:", OUT_DIR, "\n")
