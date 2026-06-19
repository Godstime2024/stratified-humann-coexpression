# =============================================================================
# mesalmon_deepdive.R
# Dissect the MEsalmon WGCNA module (top SCFP-suppressed module).
#
# Usage:
#   Rscript workflow/scripts/mesalmon_deepdive.R results_ko salmon
#
# Outputs to: results_ko/mesalmon_deepdive/
#   - module_kos_annotated.tsv  : 88 KOs in the module with KEGG names + ANCOM-BC2 stats
#   - module_taxa_contributors.tsv : taxa contributing stratified abundance to those KOs
#   - module_taxa_top.tsv : top contributing taxa ranked by total abundance
#   - methanogen_overlap.tsv : intersection with Methanobrevibacter / mcr genes
# =============================================================================

suppressPackageStartupMessages({
    library(data.table)
})

args <- commandArgs(trailingOnly = TRUE)
RES        <- if (length(args) >= 1) args[1] else "results_ko"
MODULE_TGT <- if (length(args) >= 2) args[2] else "salmon"

OUT <- file.path(RES, "mesalmon_deepdive")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

cat("==== MEsalmon deep dive ====\n")
cat("Result dir:", RES, "\n")
cat("Module:    ", MODULE_TGT, "\n\n")

# ---- 1. Pull module KOs ---------------------------------------------------
mod <- fread(file.path(RES, "wgcna/module_assignment.tsv"))
ko_ids <- mod[module == MODULE_TGT, feature]
cat("Module size:", length(ko_ids), "KOs\n\n")

# ---- 2. Annotate with KO names -------------------------------------------
ko_names <- NULL
ko_names_path <- ".kegg_cache/ko_names.tsv"
if (file.exists(ko_names_path)) {
    ko_names <- fread(ko_names_path)
}

ko_dt <- data.table(feature = ko_ids)
if (!is.null(ko_names)) {
    ko_dt <- merge(ko_dt, ko_names, by.x = "feature", by.y = "ko",
                   all.x = TRUE, sort = FALSE)
}

# ---- 3. Layer in ANCOM-BC2 stats -----------------------------------------
ancom <- fread(file.path(RES, "ancombc/ancombc_function_results.tsv"))
ko_dt <- merge(ko_dt, ancom[, .(feature, lfc, p, q)],
               by = "feature", all.x = TRUE, sort = FALSE)
setorder(ko_dt, p, na.last = TRUE)

fwrite(ko_dt, file.path(OUT, "module_kos_annotated.tsv"), sep = "\t")
cat("Top 20 KOs in MEsalmon by ANCOM-BC2 raw p-value:\n")
print(head(ko_dt, 20))

cat("\nDirection of LFC across module KOs (with ANCOM-BC2 results):\n")
print(table(sign(ko_dt$lfc), useNA = "always"))
cat("(negative LFC = lower in SCFP)\n\n")

# ---- 4. Taxa contributors via stratified file -----------------------------
strat <- fread(file.path(RES, "preprocess/stratified_filtered.tsv"))

# File already has function_name (KO id) and taxon columns parsed
# Restrict to MEsalmon KOs
strat_mod <- strat[function_name %in% ko_ids]
cat("Stratified rows mapped to MEsalmon KOs:", nrow(strat_mod), "\n")

if (nrow(strat_mod) > 0) {
    # Sample columns are the numeric ones (everything except feature/function_name/taxon)
    sample_cols <- setdiff(colnames(strat_mod), c("feature", "function_name", "taxon"))
    strat_mod[, total_abund := rowSums(.SD), .SDcols = sample_cols]

    taxa_summary <- strat_mod[, .(
        n_kos          = uniqueN(function_name),
        total_abund    = sum(total_abund),
        mean_per_ko    = sum(total_abund) / uniqueN(function_name)
    ), by = taxon]
    setorder(taxa_summary, -total_abund)

    fwrite(taxa_summary, file.path(OUT, "module_taxa_contributors.tsv"), sep = "\t")
    cat("\nTop 20 taxa contributing to MEsalmon KOs:\n")
    print(head(taxa_summary, 20))

    # Save the top-25 only (compact view)
    fwrite(head(taxa_summary, 25),
           file.path(OUT, "module_taxa_top.tsv"), sep = "\t")

    # ---- 5. Methanogen / mcr overlap -------------------------------------
    mcr_kos <- intersect(ko_ids, c("K00399", "K00400", "K00401", "K00402"))
    methanogen_taxa <- grep("Methanobrev|Methanoculleus|Methanosphaera|Methanomicro|Methano",
                            taxa_summary$taxon, value = TRUE, ignore.case = TRUE)

    cat("\n=== Methanogenesis overlap ===\n")
    cat("mcr genes in MEsalmon:", if (length(mcr_kos)) paste(mcr_kos, collapse = ", ") else "NONE", "\n")
    cat("Methanogen taxa contributing to MEsalmon KOs:",
        if (length(methanogen_taxa)) length(methanogen_taxa) else 0, "\n")
    if (length(methanogen_taxa)) {
        print(taxa_summary[taxon %in% methanogen_taxa])
    }

    fwrite(data.table(
        mcr_kos_in_module = paste(mcr_kos, collapse = ", "),
        n_methanogen_taxa = length(methanogen_taxa),
        methanogen_taxa = paste(methanogen_taxa, collapse = "; ")
    ), file.path(OUT, "methanogen_overlap.tsv"), sep = "\t")
}

cat("\n==== Done. Outputs in:", OUT, "====\n")
