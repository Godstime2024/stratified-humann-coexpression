# =============================================================================
# ancombc.R
# ANCOM-BC2 compositional differential abundance analysis for:
#   - Unstratified function table (KO/COG features)
#   - Taxa-collapsed table
# Compares treatment groups defined in metadata.
# =============================================================================

suppressPackageStartupMessages({
    library(tidyverse)
    library(data.table)
    library(ANCOMBC)
    library(TreeSummarizedExperiment)
    library(ggplot2)
    library(ggrepel)
})

log_con <- file(snakemake@log[[1]], open = "wt")
sink(log_con, type = "output"); sink(log_con, type = "message")

# Inputs are FILTERED but pre-CLR tables (raw-ish abundance)
unstrat_path <- snakemake@input[["unstrat"]]
taxa_path    <- snakemake@input[["taxa"]]
meta_path    <- snakemake@input[["meta"]]
out_ko       <- snakemake@output[["ko_results"]]
out_taxa     <- snakemake@output[["taxa_results"]]
out_vko      <- snakemake@output[["volcano_ko"]]
out_vtaxa    <- snakemake@output[["volcano_taxa"]]
group_col    <- snakemake@params[["group_col"]]
fdr_cutoff   <- as.numeric(snakemake@params[["fdr_cutoff"]])
lfc_cutoff   <- as.numeric(snakemake@params[["lfc_cutoff"]])
prv_cut      <- as.numeric(snakemake@params[["prv_cut"]])

dir.create(dirname(out_ko), recursive = TRUE, showWarnings = FALSE)

meta <- fread(meta_path)
setnames(meta, 1, "sample_id")

# Helper: run ANCOM-BC2 on a feature matrix
run_ancombc <- function(feat_dt, feat_col, scale_factor = 1e6) {
    feat_id <- feat_dt[[feat_col]]
    sample_cols <- setdiff(colnames(feat_dt), feat_col)

    # Build count matrix (features in rows, samples in columns)
    mat <- as.matrix(feat_dt[, ..sample_cols])
    rownames(mat) <- feat_id

    # Scale to "counts" so ANCOM-BC2 is happy (it expects integer-ish counts)
    counts <- round(mat * scale_factor)
    counts[counts < 0] <- 0
    storage.mode(counts) <- "integer"

    # Align metadata to sample columns
    meta_ord <- meta[match(sample_cols, sample_id)]
    if (any(is.na(meta_ord$sample_id))) {
        stop("Metadata mismatch with sample columns")
    }
    coldata <- DataFrame(meta_ord)
    rownames(coldata) <- meta_ord$sample_id

    tse <- TreeSummarizedExperiment(
        assays = list(counts = counts),
        colData = coldata
    )

    cat("[ancombc] Running on", nrow(counts), "features,",
        ncol(counts), "samples\n")
    set.seed(42)
    res <- ancombc2(
        data = tse,
        assay_name = "counts",
        fix_formula = group_col,
        p_adj_method = "BH",
        prv_cut = prv_cut,
        lib_cut = 0,
        s0_perc = 0.05,
        group = group_col,
        struc_zero = TRUE,
        neg_lb = FALSE,
        alpha = fdr_cutoff,
        verbose = FALSE
    )
    return(res)
}

# Helper: extract result table for the group contrast
extract_results <- function(res, group_col) {
    out <- as.data.table(res$res)
    # Find the contrast column for group_col
    coef_col <- grep(paste0("^lfc_", group_col), colnames(out), value = TRUE)
    se_col   <- grep(paste0("^se_", group_col), colnames(out), value = TRUE)
    p_col    <- grep(paste0("^p_", group_col), colnames(out), value = TRUE)
    q_col    <- grep(paste0("^q_", group_col), colnames(out), value = TRUE)
    diff_col <- grep(paste0("^diff_", group_col), colnames(out), value = TRUE)

    if (length(coef_col) == 0) {
        warning("Could not find lfc column for ", group_col)
        return(NULL)
    }

    res_dt <- data.table(
        feature = out$taxon,
        lfc = out[[coef_col[1]]],
        se  = out[[se_col[1]]],
        p   = out[[p_col[1]]],
        q   = out[[q_col[1]]],
        diff_significant = out[[diff_col[1]]]
    )
    res_dt[order(q)]
}

# Helper: volcano plot
make_volcano <- function(res_dt, title_txt, fdr, lfc, out_path) {
    if (is.null(res_dt) || nrow(res_dt) == 0) {
        pdf(out_path, width = 6, height = 4)
        plot.new(); title(paste("No results:", title_txt))
        dev.off()
        return()
    }

    df <- copy(res_dt)
    df[, neg_log_q := -log10(pmax(q, 1e-10))]
    df[, sig := ifelse(q < fdr & abs(lfc) > lfc,
                       ifelse(lfc > 0, "Up in TRT", "Down in TRT"),
                       "Not significant")]
    df[, label := ifelse(sig != "Not significant" & rank(-neg_log_q) <= 15,
                         feature, NA_character_)]

    p <- ggplot(df, aes(x = lfc, y = neg_log_q, color = sig)) +
        geom_point(alpha = 0.6, size = 1.5) +
        geom_hline(yintercept = -log10(fdr), linetype = "dashed", color = "grey50") +
        geom_vline(xintercept = c(-lfc, lfc), linetype = "dashed", color = "grey50") +
        scale_color_manual(values = c("Up in TRT" = "#d73027",
                                       "Down in TRT" = "#4575b4",
                                       "Not significant" = "grey70")) +
        ggrepel::geom_text_repel(aes(label = label), size = 2.8,
                                  max.overlaps = 20, na.rm = TRUE) +
        labs(title = title_txt,
             x = "Log fold change (ANCOM-BC2)",
             y = "-log10(FDR q-value)",
             color = NULL) +
        theme_minimal(base_size = 11) +
        theme(legend.position = "bottom")

    ggsave(out_path, p, width = 8, height = 6, units = "in")
}

# ---- Run on unstratified function table -----------------------------------
cat("\n=== ANCOM-BC2 on functions (KO/COG) ===\n")
unstrat <- fread(unstrat_path)
# The preprocess output stores "feature" as the first column for unstrat
feat_col <- if ("feature" %in% colnames(unstrat)) "feature" else colnames(unstrat)[1]
res_func <- tryCatch(run_ancombc(unstrat, feat_col),
                     error = function(e) {
                         cat("[ancombc] Function run failed:",
                             conditionMessage(e), "\n"); NULL
                     })
if (!is.null(res_func)) {
    func_dt <- extract_results(res_func, group_col)
    fwrite(func_dt, out_ko, sep = "\t")
    cat("[ancombc] KO/function results:", nrow(func_dt),
        "| significant (q <", fdr_cutoff, "):",
        sum(func_dt$q < fdr_cutoff, na.rm = TRUE), "\n")
    make_volcano(func_dt, "ANCOM-BC2: functions (KO or COG)",
                 fdr_cutoff, lfc_cutoff, out_vko)
} else {
    fwrite(data.table(message = "ANCOM-BC2 function run failed"),
           out_ko, sep = "\t")
    pdf(out_vko); plot.new(); title("ANCOM-BC2 function run failed"); dev.off()
}

# ---- Run on taxa-collapsed table ------------------------------------------
cat("\n=== ANCOM-BC2 on taxa ===\n")
taxa <- fread(taxa_path)
feat_col_t <- if ("taxon" %in% colnames(taxa)) "taxon" else colnames(taxa)[1]
res_taxa <- tryCatch(run_ancombc(taxa, feat_col_t),
                     error = function(e) {
                         cat("[ancombc] Taxa run failed:",
                             conditionMessage(e), "\n"); NULL
                     })
if (!is.null(res_taxa)) {
    taxa_dt <- extract_results(res_taxa, group_col)
    fwrite(taxa_dt, out_taxa, sep = "\t")
    cat("[ancombc] Taxa results:", nrow(taxa_dt),
        "| significant (q <", fdr_cutoff, "):",
        sum(taxa_dt$q < fdr_cutoff, na.rm = TRUE), "\n")
    make_volcano(taxa_dt, "ANCOM-BC2: taxa", fdr_cutoff, lfc_cutoff, out_vtaxa)
} else {
    fwrite(data.table(message = "ANCOM-BC2 taxa run failed"),
           out_taxa, sep = "\t")
    pdf(out_vtaxa); plot.new(); title("ANCOM-BC2 taxa run failed"); dev.off()
}

cat("\n[ancombc] Done.\n")
sink(type = "message"); sink(type = "output"); close(log_con)
