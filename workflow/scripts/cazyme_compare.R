# =============================================================================
# cazyme_compare.R
# Compare CAZyme family abundance between CON and SCFP.
# Special focus on fibrolytic GH families to validate the MEsalmon xynD finding.
# =============================================================================

suppressPackageStartupMessages({ library(data.table) })

CAZ_FILE <- "/mnt/d/Wisconsin_data1/merged_results/subset/cazyme_family_matrix_subset.tsv"
OUT_DIR  <- "results_ko/cazyme_compare"
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

meta <- data.table(
    sample_id = c("RUM_KSU_17","RUM_KSU_18","RUM_KSU_19","RUM_KSU_20","RUM_KSU_21",
                  "RUM_KSU_23","RUM_KSU_42","RUM_KSU_43","RUM_KSU_44","RUM_KSU_45",
                  "RUM_KSU_46","RUM_KSU_47","RUM_KSU_48"),
    treatment = c("CON","SCFP","SCFP","CON","CON",
                  "SCFP","CON","CON","SCFP","SCFP",
                  "CON","CON","SCFP")
)

# CAZyme matrix is families x samples; transpose
caz <- fread(CAZ_FILE)
fam_col <- "CAZyme_Family"
samples <- setdiff(colnames(caz), fam_col)
caz_long <- melt(caz, id.vars = fam_col, variable.name = "sample_id",
                 value.name = "count", variable.factor = FALSE)
caz_wide <- dcast(caz_long, sample_id ~ get(fam_col), value.var = "count", fill = 0)
caz_wide <- merge(meta, caz_wide, by = "sample_id")
fam_names <- setdiff(colnames(caz_wide), c("sample_id", "treatment"))

cat("Total CAZyme families:", length(fam_names), "\n")
cat("Samples:\n"); print(table(caz_wide$treatment))

# Per sample CAZyme total
caz_wide[, total_caz := rowSums(.SD), .SDcols = fam_names]
cat("\n=== Total CAZyme abundance per sample ===\n")
print(caz_wide[, .(sample_id, treatment, total_caz)])
cat("\n=== Mean by group ===\n")
print(caz_wide[, .(mean = mean(total_caz), sd = sd(total_caz), n = .N), by = treatment])
w_total <- wilcox.test(total_caz ~ treatment, data = caz_wide)
cat("\nWilcoxon total CAZymes CON vs SCFP: W =", w_total$statistic,
    " p =", signif(w_total$p.value, 4), "\n")

# Family-class summary (GH, GT, PL, CE, AA, CBM)
caz_long2 <- merge(meta, caz_long, by = "sample_id")
caz_long2[, class := sub("[0-9].*", "", get(fam_col))]
caz_long2[, class := ifelse(grepl("^GH", get(fam_col)), "GH",
                     ifelse(grepl("^GT", get(fam_col)), "GT",
                     ifelse(grepl("^PL", get(fam_col)), "PL",
                     ifelse(grepl("^CE", get(fam_col)), "CE",
                     ifelse(grepl("^AA", get(fam_col)), "AA",
                     ifelse(grepl("^CBM", get(fam_col)), "CBM", "other"))))))]
class_sum <- caz_long2[, .(count = sum(count)), by = .(sample_id, treatment, class)]
class_wide <- dcast(class_sum, sample_id + treatment ~ class, value.var = "count", fill = 0)
cat("\n=== Per CAZyme class totals ===\n")
print(class_wide)

cat("\n=== Per CAZyme class Wilcoxon ===\n")
class_cols <- setdiff(colnames(class_wide), c("sample_id", "treatment"))
class_stats <- rbindlist(lapply(class_cols, function(cl) {
    con  <- class_wide[treatment == "CON",  get(cl)]
    scfp <- class_wide[treatment == "SCFP", get(cl)]
    if (sum(con) + sum(scfp) == 0) return(NULL)
    wt <- suppressWarnings(wilcox.test(con, scfp))
    data.table(class = cl, CON_mean = mean(con), SCFP_mean = mean(scfp),
               diff = mean(con) - mean(scfp), W = unname(wt$statistic), p = wt$p.value)
}))
class_stats[, q := p.adjust(p, method = "BH")]
setorder(class_stats, p)
print(class_stats)
fwrite(class_stats, file.path(OUT_DIR, "cazyme_class_compare.tsv"), sep = "\t")

# Family level test
cat("\n=== Per family Wilcoxon (top 25 by p, requiring >= 5 samples with hits) ===\n")
fam_stats <- rbindlist(lapply(fam_names, function(f) {
    con  <- caz_wide[treatment == "CON",  get(f)]
    scfp <- caz_wide[treatment == "SCFP", get(f)]
    n_hit <- sum(con > 0) + sum(scfp > 0)
    if (n_hit < 5) return(NULL)
    wt <- suppressWarnings(wilcox.test(con, scfp))
    data.table(family = f, CON_mean = round(mean(con), 1), SCFP_mean = round(mean(scfp), 1),
               diff = round(mean(con) - mean(scfp), 1), n_with_hit = n_hit,
               W = unname(wt$statistic), p = wt$p.value)
}))
fam_stats[, q := p.adjust(p, method = "BH")]
setorder(fam_stats, p)
print(head(fam_stats, 25))
fwrite(fam_stats, file.path(OUT_DIR, "cazyme_family_compare.tsv"), sep = "\t")

# Focus on fibrolytic GH families (cellulose + hemicellulose)
fibrolytic <- c("GH5","GH8","GH9","GH10","GH11","GH26","GH28","GH43","GH44","GH45","GH48",
                "GH51","GH53","GH67","GH74","GH95","GH115","GH130","GH3","GH16","GH30")
cat("\n=== Fibrolytic GH families (cellulose + hemicellulose degradation) ===\n")
fib_stats <- fam_stats[family %in% fibrolytic]
print(fib_stats)
fwrite(fib_stats, file.path(OUT_DIR, "fibrolytic_GH_compare.tsv"), sep = "\t")

# Direction summary
cat("\n=== Direction of differential CAZymes (n=88 families with >= 5 samples) ===\n")
cat("Higher in CON (positive diff):", sum(fam_stats$diff > 0), "\n")
cat("Higher in SCFP (negative diff):", sum(fam_stats$diff < 0), "\n")
cat("Equal:", sum(fam_stats$diff == 0), "\n")

cat("\nOutputs in:", OUT_DIR, "\n")
