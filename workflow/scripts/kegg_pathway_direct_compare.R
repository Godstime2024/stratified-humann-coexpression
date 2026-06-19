# =============================================================================
# kegg_pathway_direct_compare.R
# Direct CON vs SCFP test on KEGG pathways from the precomputed matrix.
# Independent of WGCNA and ANCOM-BC2; uses Wilcoxon on per-sample pathway counts.
#
# Special focus on:
#   ko00680  Methane metabolism (the methanogenesis pathway)
#   ko00720  Carbon fixation in prokaryotes (methanogen CO2 fixation)
#   ko00640  Propanoate metabolism (VFA - propionate)
#   ko00650  Butanoate metabolism (VFA - butyrate)
#   ko00010  Glycolysis
#   ko00190  Oxidative phosphorylation
# =============================================================================

suppressPackageStartupMessages({ library(data.table) })

KEGG_FILE <- "/mnt/d/Wisconsin_data1/merged_results/subset/kegg_pathway_matrix_subset.tsv"
OUT_DIR   <- "results_ko/kegg_pathway_direct"
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

meta <- data.table(
    sample_id = c("RUM_KSU_17","RUM_KSU_18","RUM_KSU_19","RUM_KSU_20","RUM_KSU_21",
                  "RUM_KSU_23","RUM_KSU_42","RUM_KSU_43","RUM_KSU_44","RUM_KSU_45",
                  "RUM_KSU_46","RUM_KSU_47","RUM_KSU_48"),
    treatment = c("CON","SCFP","SCFP","CON","CON",
                  "SCFP","CON","CON","SCFP","SCFP",
                  "CON","CON","SCFP")
)

pw <- fread(KEGG_FILE)
cat("Matrix:", nrow(pw), "samples x", ncol(pw) - 1, "columns\n")

# Dedupe ko* vs map* (same pathway under two prefixes)
pw_cols <- setdiff(colnames(pw), "Sample")
ko_cols <- pw_cols[grepl("^ko", pw_cols)]
cat("Using ko* prefix (dropping duplicated map*) -", length(ko_cols), "pathways\n")
pw <- pw[, c("Sample", ko_cols), with = FALSE]

pw <- merge(meta, pw, by.x = "sample_id", by.y = "Sample")

# Drop any sample with zero library size (would produce NaN on normalization)
pw_mat_pre <- as.matrix(pw[, ..ko_cols])
lib_pre <- rowSums(pw_mat_pre)
cat("\nLibrary sizes (pre-filter):\n"); print(setNames(lib_pre, pw$sample_id))
keep <- lib_pre > 0
if (sum(!keep) > 0) {
    cat("\nDropping zero-library samples:", paste(pw$sample_id[!keep], collapse = ", "), "\n")
    pw <- pw[keep, ]
}
samples <- pw$sample_id
treatment <- pw$treatment

# Library size normalization (TSS) to make per-sample comparable
pw_mat <- as.matrix(pw[, ..ko_cols])
lib_size <- rowSums(pw_mat)
cat("\nLibrary sizes (after filter):\n"); print(setNames(lib_size, samples))
pw_norm <- sweep(pw_mat, 1, lib_size, "/")  # relative abundance
cat("\nGroup sizes after filter:\n"); print(table(treatment))

cat("\n=== Direct test: methane and related pathways (relative abundance) ===\n")
focus <- c("ko00680" = "Methane metabolism",
           "ko00720" = "Carbon fixation in prokaryotes",
           "ko00640" = "Propanoate metabolism",
           "ko00650" = "Butanoate metabolism",
           "ko00010" = "Glycolysis",
           "ko00190" = "Oxidative phosphorylation",
           "ko00500" = "Starch and sucrose metabolism",
           "ko00520" = "Amino sugar and nucleotide sugar metabolism",
           "ko00770" = "Pantothenate and CoA biosynthesis",
           "ko02040" = "Flagellar assembly",
           "ko00550" = "Peptidoglycan biosynthesis")

focus_dt <- rbindlist(lapply(names(focus), function(p) {
    if (!(p %in% colnames(pw_norm))) return(NULL)
    vals <- pw_norm[, p]
    con  <- vals[treatment == "CON"]
    scfp <- vals[treatment == "SCFP"]
    wt   <- suppressWarnings(wilcox.test(con, scfp))
    raw_con  <- pw_mat[treatment == "CON",  p]
    raw_scfp <- pw_mat[treatment == "SCFP", p]
    data.table(
        pathway     = p,
        name        = focus[[p]],
        CON_rel_mean  = signif(mean(con), 3),
        SCFP_rel_mean = signif(mean(scfp), 3),
        log2_fc       = signif(log2((mean(scfp) + 1e-9) / (mean(con) + 1e-9)), 3),
        CON_raw_mean  = round(mean(raw_con), 1),
        SCFP_raw_mean = round(mean(raw_scfp), 1),
        W             = unname(wt$statistic),
        p             = wt$p.value
    )
}))
focus_dt[, q := p.adjust(p, method = "BH")]
setorder(focus_dt, p)
print(focus_dt)
fwrite(focus_dt, file.path(OUT_DIR, "focus_pathway_compare.tsv"), sep = "\t")

# Full pathway scan
cat("\n=== Full scan: top 25 pathways by raw p (relative abundance) ===\n")
full_dt <- rbindlist(lapply(ko_cols, function(p) {
    vals <- pw_norm[, p]
    if (any(is.na(vals))) return(NULL)
    con  <- vals[treatment == "CON"]
    scfp <- vals[treatment == "SCFP"]
    if ((sum(con) + sum(scfp)) == 0) return(NULL)
    wt <- suppressWarnings(wilcox.test(con, scfp))
    data.table(
        pathway = p,
        CON_rel_mean  = signif(mean(con), 3),
        SCFP_rel_mean = signif(mean(scfp), 3),
        log2_fc       = signif(log2((mean(scfp) + 1e-9) / (mean(con) + 1e-9)), 3),
        W = unname(wt$statistic),
        p = wt$p.value
    )
}))
full_dt[, q := p.adjust(p, method = "BH")]
setorder(full_dt, p)
print(head(full_dt, 25))
fwrite(full_dt, file.path(OUT_DIR, "all_pathway_compare.tsv"), sep = "\t")

# Methane metabolism per-sample values for plot/headline
cat("\n=== ko00680 Methane metabolism per sample (relative abundance) ===\n")
me <- data.table(
    sample_id = samples,
    treatment = treatment,
    methane_raw = pw_mat[, "ko00680"],
    methane_rel = signif(pw_norm[, "ko00680"], 4)
)
setorder(me, treatment, sample_id)
print(me)
cat("\nCON mean rel:", signif(mean(me[treatment == "CON",  methane_rel]), 4),
    "  SCFP mean rel:", signif(mean(me[treatment == "SCFP", methane_rel]), 4), "\n")
fwrite(me, file.path(OUT_DIR, "methane_per_sample.tsv"), sep = "\t")

cat("\nOutputs in:", OUT_DIR, "\n")
