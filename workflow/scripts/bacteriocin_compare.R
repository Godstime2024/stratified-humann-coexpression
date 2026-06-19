# =============================================================================
# bacteriocin_compare.R
# Compare bacteriocin gene burden between CON and SCFP.
# =============================================================================

suppressPackageStartupMessages({ library(data.table) })

BAC_FILE <- "/mnt/d/Wisconsin_data1/amr_bagel_results_hpc/merged/bacteriocin_screen.tsv"
OUT_DIR  <- "results_ko/amr_compare"
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

meta <- data.table(
    sample_id = c("RUM_KSU_17","RUM_KSU_18","RUM_KSU_19","RUM_KSU_20","RUM_KSU_21",
                  "RUM_KSU_23","RUM_KSU_42","RUM_KSU_43","RUM_KSU_44","RUM_KSU_45",
                  "RUM_KSU_46","RUM_KSU_47","RUM_KSU_48"),
    treatment = c("CON","SCFP","SCFP","CON","CON",
                  "SCFP","CON","CON","SCFP","SCFP",
                  "CON","CON","SCFP")
)

bac <- fread(BAC_FILE)
cat("Total bacteriocin hit rows:", nrow(bac), "\n")
cat("Bacteriocin families:\n"); print(table(bac$Bacteriocin_Family))

# ---- 1. Per sample family count -------------------------------------------
fam_count <- bac[, .(n_hits = .N), by = .(Sample, Bacteriocin_Family)]
fam_wide  <- dcast(fam_count, Sample ~ Bacteriocin_Family, value.var = "n_hits", fill = 0)
fam_wide  <- merge(meta, fam_wide, by.x = "sample_id", by.y = "Sample")
fam_cols  <- setdiff(colnames(fam_wide), c("sample_id", "treatment"))
fam_wide[, total_bac := rowSums(.SD), .SDcols = fam_cols]

cat("\n=== Per sample bacteriocin totals ===\n")
print(fam_wide[, .(sample_id, treatment, total_bac)])

cat("\n=== Mean by group ===\n")
print(fam_wide[, .(mean = mean(total_bac), sd = sd(total_bac), n = .N), by = treatment])

w <- wilcox.test(total_bac ~ treatment, data = fam_wide)
cat("\nWilcoxon total bacteriocin hits CON vs SCFP: W =", w$statistic,
    " p =", signif(w$p.value, 4), "\n")

# Per family test
cat("\n=== Per family Wilcoxon ===\n")
fam_stats <- rbindlist(lapply(fam_cols, function(f) {
    con  <- fam_wide[treatment == "CON",  get(f)]
    scfp <- fam_wide[treatment == "SCFP", get(f)]
    if (sum(con) + sum(scfp) == 0) return(NULL)
    w <- suppressWarnings(wilcox.test(con, scfp))
    data.table(
        family    = f,
        CON_mean  = mean(con),
        SCFP_mean = mean(scfp),
        diff      = mean(con) - mean(scfp),
        W         = unname(w$statistic),
        p         = w$p.value
    )
}))
fam_stats[, q := p.adjust(p, method = "BH")]
setorder(fam_stats, p)
print(fam_stats)
fwrite(fam_stats, file.path(OUT_DIR, "bacteriocin_family_compare.tsv"), sep = "\t")

# ---- 2. Per keyword (specific bacteriocin gene) test ----------------------
kw_count <- bac[, .(n_hits = .N), by = .(Sample, Keywords_Hit)]
kw_wide  <- dcast(kw_count, Sample ~ Keywords_Hit, value.var = "n_hits", fill = 0)
kw_wide  <- merge(meta, kw_wide, by.x = "sample_id", by.y = "Sample")
kw_cols  <- setdiff(colnames(kw_wide), c("sample_id", "treatment"))

cat("\n=== Per keyword Wilcoxon (top 25 by p, requiring >= 5 samples with hits) ===\n")
kw_stats <- rbindlist(lapply(kw_cols, function(k) {
    con  <- kw_wide[treatment == "CON",  get(k)]
    scfp <- kw_wide[treatment == "SCFP", get(k)]
    if ((sum(con > 0) + sum(scfp > 0)) < 5) return(NULL)
    w <- suppressWarnings(wilcox.test(con, scfp))
    data.table(
        keyword   = k,
        CON_mean  = round(mean(con), 2),
        SCFP_mean = round(mean(scfp), 2),
        diff      = round(mean(con) - mean(scfp), 2),
        n_with_hit = sum(con > 0) + sum(scfp > 0),
        W         = unname(w$statistic),
        p         = w$p.value
    )
}))
kw_stats[, q := p.adjust(p, method = "BH")]
setorder(kw_stats, p)
print(head(kw_stats, 25))
fwrite(kw_stats, file.path(OUT_DIR, "bacteriocin_gene_compare.tsv"), sep = "\t")

cat("\nOutputs in:", OUT_DIR, "\n")
