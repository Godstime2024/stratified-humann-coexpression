# =============================================================================
# functional_permanova.R
# PERMANOVA on CLR-transformed unstratified KO matrix (CON vs SCFP) plus
# PCA ordination.
# =============================================================================

suppressPackageStartupMessages({
    library(data.table)
    library(vegan)
})

CLR_FILE <- "results_ko/preprocess/unstratified_clr.tsv"
META     <- data.table(
    sample_id = c("RUM_KSU_17","RUM_KSU_18","RUM_KSU_19","RUM_KSU_20","RUM_KSU_21",
                  "RUM_KSU_23","RUM_KSU_43","RUM_KSU_44","RUM_KSU_45","RUM_KSU_46",
                  "RUM_KSU_47","RUM_KSU_48"),
    treatment = c("CON","SCFP","SCFP","CON","CON",
                  "SCFP","CON","SCFP","SCFP","CON","CON","SCFP")
)
OUT_DIR <- "results_ko/permanova_func"
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

cat("Reading CLR matrix...\n")
clr <- fread(CLR_FILE)
# Detect first non-sample column (assume row IDs in first col)
id_col <- colnames(clr)[1]
sample_cols <- intersect(colnames(clr), META$sample_id)
cat("Samples found in matrix:", length(sample_cols), "\n")
clr_mat <- t(as.matrix(clr[, ..sample_cols]))  # samples x features
rownames(clr_mat) <- sample_cols
cat("KO features:", ncol(clr_mat), "  samples:", nrow(clr_mat), "\n")

meta_use <- META[match(rownames(clr_mat), sample_id)]

# Distance matrix (Euclidean on CLR = Aitchison distance, the standard for
# compositional functional data)
d <- vegdist(clr_mat, method = "euclidean")

# PERMANOVA
set.seed(1)
perm <- adonis2(d ~ treatment, data = meta_use, permutations = 9999,
                method = "euclidean")
cat("\n=== PERMANOVA on KO Aitchison distance ===\n")
print(perm)
fwrite(as.data.table(perm, keep.rownames = "term"),
       file.path(OUT_DIR, "permanova_KO.tsv"), sep = "\t")

# PERMDISP (homogeneity of dispersions) — to know whether group separation
# reflects location shift or dispersion difference
disp <- betadisper(d, meta_use$treatment)
disp_perm <- permutest(disp, permutations = 9999)
cat("\n=== PERMDISP (group dispersion homogeneity) ===\n")
print(disp_perm)

# PCA (eigendecomposition on the CLR matrix, the natural ordination for
# Aitchison-distanced compositional data)
pca <- prcomp(clr_mat, center = TRUE, scale. = FALSE)
ve <- pca$sdev^2 / sum(pca$sdev^2)
scores <- as.data.table(pca$x[, 1:4])
scores[, sample_id := rownames(pca$x)]
scores[, treatment := meta_use$treatment[match(sample_id, meta_use$sample_id)]]
fwrite(scores, file.path(OUT_DIR, "pca_scores.tsv"), sep = "\t")
cat("\n=== Variance explained (PC1-PC4) ===\n")
print(round(ve[1:4] * 100, 2))
fwrite(data.table(component = paste0("PC", seq_along(ve)),
                  variance_explained = ve),
       file.path(OUT_DIR, "pca_variance.tsv"), sep = "\t")

# Top contributing KOs to PC1 + PC2
load <- pca$rotation[, 1:2]
rn <- rownames(load)
if (is.null(rn)) rn <- paste0("feat_", seq_len(nrow(load)))
ko_load <- data.table(feature = rn, PC1 = load[,1], PC2 = load[,2])
ko_load[, abs1 := abs(PC1)]
setorder(ko_load, -abs1)
fwrite(head(ko_load[, .(feature, PC1, PC2)], 50),
       file.path(OUT_DIR, "top_PC1_loadings.tsv"), sep = "\t")

cat("\nDone. Outputs in:", OUT_DIR, "\n")
