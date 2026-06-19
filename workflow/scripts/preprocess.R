# =============================================================================
# preprocess.R
# Split stratified vs unstratified; prevalence filter; CLR transform;
# collapse stratified rows into taxa-level matrix
# =============================================================================

suppressPackageStartupMessages({
    library(tidyverse)
    library(data.table)
    library(compositions)
})

# ---- Snakemake I/O ---------------------------------------------------------
log_con <- file(snakemake@log[[1]], open = "wt")
sink(log_con, type = "output")
sink(log_con, type = "message")

input_csv     <- snakemake@input[["csv"]]
input_meta    <- snakemake@input[["meta"]]
out_strat     <- snakemake@output[["strat"]]
out_unstrat   <- snakemake@output[["unstrat"]]
out_strat_clr <- snakemake@output[["strat_clr"]]
out_unstrat_clr <- snakemake@output[["unstrat_clr"]]
out_taxa      <- snakemake@output[["taxa_collapsed"]]
out_taxa_clr  <- snakemake@output[["taxa_clr"]]
out_summary   <- snakemake@output[["summary"]]

p_sep         <- snakemake@params[["input_sep"]]
strat_sep     <- snakemake@params[["stratify_sep"]]
min_samples   <- as.integer(snakemake@params[["min_samples"]])
min_value     <- as.numeric(snakemake@params[["min_value"]])
pseudo_strat  <- snakemake@params[["pseudo_strategy"]]
pseudo_fixed  <- as.numeric(snakemake@params[["pseudo_fixed"]])

# Ensure output dirs
for (f in c(out_strat, out_unstrat, out_strat_clr, out_unstrat_clr,
            out_taxa, out_taxa_clr, out_summary)) {
    dir.create(dirname(f), recursive = TRUE, showWarnings = FALSE)
}

cat("[preprocess] Reading:", input_csv, "\n")
raw <- fread(input_csv, sep = p_sep, header = TRUE, check.names = FALSE)
colnames(raw)[1] <- "feature"
cat("[preprocess] Dimensions:", nrow(raw), "rows x", ncol(raw), "cols\n")

sample_cols <- setdiff(colnames(raw), "feature")
cat("[preprocess] Sample columns (", length(sample_cols), "):\n", sep = "")
cat(paste(" -", sample_cols), sep = "\n")

# ---- Split stratified vs unstratified --------------------------------------
is_strat <- grepl(strat_sep, raw$feature, fixed = TRUE)
strat <- raw[is_strat, ]
unstrat <- raw[!is_strat, ]
cat("[preprocess] Stratified rows:", nrow(strat), "\n")
cat("[preprocess] Unstratified rows:", nrow(unstrat), "\n")

# Add function + taxon columns to stratified
strat_split <- tstrsplit(strat$feature, strat_sep, fixed = TRUE)
strat$function_name <- strat_split[[1]]
strat$taxon <- strat_split[[2]]
setcolorder(strat, c("feature", "function_name", "taxon", sample_cols))

# ---- Prevalence filter -----------------------------------------------------
prev_filter <- function(df, samples, min_n, min_v) {
    mat <- as.matrix(df[, ..samples])
    prev <- rowSums(mat > min_v, na.rm = TRUE)
    df[prev >= min_n, ]
}

strat_f <- prev_filter(strat, sample_cols, min_samples, min_value)
unstrat_f <- prev_filter(unstrat, sample_cols, min_samples, min_value)
cat("[preprocess] After prevalence filter:\n")
cat("  - Stratified:", nrow(strat_f), "\n")
cat("  - Unstratified:", nrow(unstrat_f), "\n")

fwrite(strat_f, out_strat, sep = "\t")
fwrite(unstrat_f, out_unstrat, sep = "\t")

# ---- CLR transform ---------------------------------------------------------
clr_transform <- function(df, samples, strategy = "half_min", fixed = 1e-7) {
    mat <- as.matrix(df[, ..samples])
    if (strategy == "half_min") {
        nz <- mat[mat > 0]
        pseudo <- min(nz, na.rm = TRUE) / 2
    } else {
        pseudo <- fixed
    }
    mat_p <- mat + pseudo
    # CLR per sample (column): log(x) - mean(log(x))
    log_mat <- log(mat_p)
    clr_mat <- sweep(log_mat, 2, colMeans(log_mat), FUN = "-")
    clr_df <- as.data.table(clr_mat)
    cbind(df[, .(feature)], clr_df)
}

strat_clr <- clr_transform(strat_f, sample_cols, pseudo_strat, pseudo_fixed)
unstrat_clr <- clr_transform(unstrat_f, sample_cols, pseudo_strat, pseudo_fixed)
fwrite(strat_clr, out_strat_clr, sep = "\t")
fwrite(unstrat_clr, out_unstrat_clr, sep = "\t")
cat("[preprocess] CLR transform complete.\n")

# ---- Collapse stratified into taxa-level matrix ----------------------------
# Sum across functions for each taxon, then prevalence filter
taxa_mat <- strat_f[, lapply(.SD, sum, na.rm = TRUE),
                    by = taxon, .SDcols = sample_cols]
mat <- as.matrix(taxa_mat[, ..sample_cols])
prev <- rowSums(mat > min_value, na.rm = TRUE)
taxa_mat <- taxa_mat[prev >= min_samples, ]
cat("[preprocess] Collapsed taxa:", nrow(taxa_mat), "\n")
fwrite(taxa_mat, out_taxa, sep = "\t")

# CLR for taxa
nz <- mat[mat > 0]
pseudo <- if (pseudo_strat == "half_min") min(nz, na.rm = TRUE) / 2 else pseudo_fixed
log_mat <- log(mat + pseudo)
clr_mat <- sweep(log_mat, 2, colMeans(log_mat), FUN = "-")
taxa_clr <- cbind(taxa_mat[, .(taxon)], as.data.table(clr_mat))
fwrite(taxa_clr, out_taxa_clr, sep = "\t")

# ---- Summary ---------------------------------------------------------------
summary_dt <- data.table(
    metric = c("input_rows", "input_samples",
               "stratified_rows_raw", "stratified_rows_filtered",
               "unstratified_rows_raw", "unstratified_rows_filtered",
               "collapsed_taxa", "min_samples_threshold", "min_value_threshold",
               "pseudocount_strategy"),
    value  = c(nrow(raw), length(sample_cols),
               nrow(strat), nrow(strat_f),
               nrow(unstrat), nrow(unstrat_f),
               nrow(taxa_mat), min_samples, min_value, pseudo_strat)
)
fwrite(summary_dt, out_summary, sep = "\t")
cat("[preprocess] Summary written:", out_summary, "\n")

cat("[preprocess] Done.\n")
sink(type = "message"); sink(type = "output"); close(log_con)
