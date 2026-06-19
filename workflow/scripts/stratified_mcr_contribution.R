# =============================================================================
# stratified_mcr_contribution.R
# Decompose mcrA/B/C/G (and other methanogenesis KOs) stratified abundance by
# contributing taxon, CON vs SCFP. Demonstrates Nanopore species-level resolution
# of the methanogenesis suppression signal.
# =============================================================================
suppressPackageStartupMessages({ library(data.table) })

STRAT_FILE <- "results_ko/preprocess/stratified_filtered.tsv"
OUT_DIR    <- "results_ko/stratified_mcr"
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

meta <- data.table(
    sample_id = c("RUM_KSU_17","RUM_KSU_18","RUM_KSU_19","RUM_KSU_20","RUM_KSU_21",
                  "RUM_KSU_23","RUM_KSU_43","RUM_KSU_44","RUM_KSU_45","RUM_KSU_46",
                  "RUM_KSU_47","RUM_KSU_48"),
    treatment = c("CON","SCFP","SCFP","CON","CON",
                  "SCFP","CON","SCFP","SCFP","CON","CON","SCFP")
)

# Target methanogenesis KOs + Wood-Ljungdahl acsB for contrast
TARGET_KOS <- c(
    "K00399"="mcrA  methyl-CoM reductase alpha",
    "K00400"="mcrC  methyl-CoM reductase C",
    "K00401"="mcrB  methyl-CoM reductase beta",
    "K00402"="mcrG  methyl-CoM reductase gamma",
    "K14080"="mtrA  N5-methyl-tetrahydromethanopterin",
    "K14138"="acsB  acetyl-CoA synthase (Wood-Ljungdahl)"
)

strat <- fread(STRAT_FILE)
sample_cols <- intersect(colnames(strat), meta$sample_id)
cat("Samples in stratified matrix:", length(sample_cols), "\n")

# Focus on target KOs
sub <- strat[function_name %in% names(TARGET_KOS)]
cat("Stratified rows for target KOs:", nrow(sub), "\n")

# Long format
sub_long <- melt(sub, id.vars = c("feature", "function_name", "taxon"),
                 measure.vars = sample_cols,
                 variable.name = "sample_id", value.name = "abundance",
                 variable.factor = FALSE)
sub_long <- merge(sub_long, meta, by = "sample_id")

# Per-taxon, per-KO, per-treatment mean
agg <- sub_long[, .(mean_abund = mean(abundance)),
                by = .(function_name, taxon, treatment)]
agg[, gene := TARGET_KOS[function_name]]

# Total per-KO per-treatment for share calc
totals <- agg[, .(total = sum(mean_abund)), by = .(function_name, treatment)]
agg <- merge(agg, totals, by = c("function_name", "treatment"))
agg[, share := mean_abund / total]
agg[total == 0, share := 0]

fwrite(agg, file.path(OUT_DIR, "per_taxon_per_treatment.tsv"), sep = "\t")

# Top contributors overall (across both treatments)
overall_taxa <- agg[, .(total = sum(mean_abund)), by = .(function_name, taxon)]
setorder(overall_taxa, function_name, -total)
top_per_ko <- overall_taxa[, head(.SD, 8), by = function_name]
fwrite(top_per_ko, file.path(OUT_DIR, "top_contributors_per_KO.tsv"), sep = "\t")
cat("\n=== Top 8 contributing taxa per KO ===\n")
print(top_per_ko)

# CON vs SCFP comparison of Methanobrevibacter share
cat("\n=== Methanobrevibacter share per KO, CON vs SCFP ===\n")
mbrev <- agg[grepl("Methanobrevibacter", taxon)]
mbrev_share <- mbrev[, .(mbrev_share = sum(share),
                         mbrev_abund = sum(mean_abund)),
                     by = .(function_name, treatment)]
mbrev_share[, gene := TARGET_KOS[function_name]]
fwrite(mbrev_share, file.path(OUT_DIR, "methanobrevibacter_share.tsv"), sep = "\t")
print(mbrev_share[order(function_name, treatment)])

# Per-sample Methanobrevibacter mcr abundance for Wilcoxon
mcrkos <- c("K00399", "K00400", "K00401", "K00402")
mcr_sub <- sub_long[function_name %in% mcrkos & grepl("Methanobrevibacter", taxon)]
mcr_per_sample <- mcr_sub[, .(mbrev_mcr_total = sum(abundance)),
                          by = .(sample_id, treatment)]
cat("\n=== Per sample Methanobrevibacter total mcr abundance ===\n")
print(mcr_per_sample[order(treatment, sample_id)])
w <- wilcox.test(mbrev_mcr_total ~ treatment, data = mcr_per_sample)
cat("\nWilcoxon Methanobrevibacter total mcr CON vs SCFP: W =", w$statistic,
    " p =", signif(w$p.value, 4), "\n")
fwrite(mcr_per_sample, file.path(OUT_DIR, "mbrev_per_sample_mcr.tsv"), sep = "\t")

cat("\nDone. Outputs in:", OUT_DIR, "\n")
