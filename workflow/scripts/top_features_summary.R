# =============================================================================
# top_features_summary.R
# Standalone post-pipeline script that extracts the top N features from every
# analysis output into a single readable summary.
#
# Usage:
#   conda activate r_microbiome
#   Rscript workflow/scripts/top_features_summary.R results_ko 15
#
# Args:
#   1: result directory (e.g., results_ko or results)
#   2: top N to extract (default 15)
# =============================================================================

suppressPackageStartupMessages({
    library(data.table)
})

args <- commandArgs(trailingOnly = TRUE)
RES <- if (length(args) >= 1) args[1] else "results_ko"
TOP_N <- if (length(args) >= 2) as.integer(args[2]) else 15

cat("====================================================\n")
cat("Top features summary\n")
cat("Result dir:", RES, "\n")
cat("Top N:     ", TOP_N, "\n")
cat("====================================================\n\n")

OUT_DIR <- file.path(RES, "top_features_summary")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

# ---- Load KEGG annotation cache if available ------------------------------
ko_names_path <- ".kegg_cache/ko_names.tsv"
ko_names <- NULL
if (file.exists(ko_names_path)) {
    ko_names <- fread(ko_names_path)
    cat("[annotate] Loaded", nrow(ko_names), "KO names from cache\n\n")
} else if (file.exists(".kegg_cache/ko_names.rds")) {
    ko_names <- readRDS(".kegg_cache/ko_names.rds")
    cat("[annotate] Loaded", nrow(ko_names), "KO names from cache\n\n")
} else {
    cat("[annotate] No KO annotation cache found at .kegg_cache/ko_names.tsv\n")
    cat("[annotate] (Run KEGG enrichment first, or KO IDs will appear without names.)\n\n")
}

annotate_kos <- function(dt, ko_col = "feature") {
    # Add a "ko_name" column next to the KO ID column
    if (is.null(ko_names) || !ko_col %in% colnames(dt)) return(dt)
    out <- merge(dt, ko_names, by.x = ko_col, by.y = "ko",
                 all.x = TRUE, sort = FALSE)
    out
}

# ---- 1. TOP TREATMENT-CORRELATED MODULES ----------------------------------
cat("=== 1. Top WGCNA modules by treatment correlation ===\n\n")
trait_path <- file.path(RES, "wgcna/module_trait_correlation.tsv")
mod_assign_path <- file.path(RES, "wgcna/module_assignment.tsv")

if (file.exists(trait_path)) {
    tcor <- fread(trait_path)
    # First column = module, last column = treatment (TRT)
    cor_col <- colnames(tcor)[ncol(tcor)]
    tcor[, abs_r := abs(get(cor_col))]
    setorder(tcor, -abs_r)
    top_mods <- head(tcor, TOP_N)

    # Add module size
    if (file.exists(mod_assign_path)) {
        mod <- fread(mod_assign_path)
        sizes <- mod[, .N, by = module]
        sizes[, module_label := paste0("ME", module)]
        top_mods <- merge(top_mods, sizes[, .(module = module_label, N)],
                          by.x = "module", by.y = "module", all.x = TRUE)
        setorder(top_mods, -abs_r)
    }

    fwrite(top_mods, file.path(OUT_DIR, "top_modules.tsv"), sep = "\t")
    print(top_mods)
    cat("\n")
} else {
    cat("  Not found:", trait_path, "\n\n")
}

# ---- 2. TOP KEGG PATHWAYS ---------------------------------------------------
cat("=== 2. Top KEGG pathway enrichments ===\n\n")
kegg_path <- file.path(RES, "kegg_enrichment/module_pathway_table.tsv")

if (file.exists(kegg_path)) {
    kegg <- fread(kegg_path)
    if ("p.adjust" %in% colnames(kegg)) {
        setorder(kegg, p.adjust)
        top_kegg <- head(kegg, TOP_N)
        top_kegg_brief <- top_kegg[, .(module, ID, name = name, p.adjust, Count, n_module_kos)]
        fwrite(top_kegg_brief, file.path(OUT_DIR, "top_kegg_pathways.tsv"), sep = "\t")
        print(top_kegg_brief)
    } else {
        cat("  KEGG enrichment likely skipped (no p.adjust column)\n")
    }
    cat("\n")
} else {
    cat("  Not found:", kegg_path, "\n\n")
}

# ---- 3. TOP DIABLO LOADINGS (functions + taxa) ---------------------------
cat("=== 3. Top DIABLO loadings (taxa + functions) ===\n\n")
diablo_path <- file.path(RES, "diablo/diablo_loadings.tsv")

if (file.exists(diablo_path)) {
    diablo <- fread(diablo_path)
    if ("loading_comp1" %in% colnames(diablo)) {
        diablo[, abs_load := abs(loading_comp1)]
        setorder(diablo, -abs_load)

        # Split by block
        for (b in unique(diablo$block)) {
            cat("--- Top", TOP_N, "in block:", b, "---\n")
            top_b <- head(diablo[block == b], TOP_N)
            # Annotate KOs (function block); leave taxa block as-is
            if (b != "taxa") top_b <- annotate_kos(top_b, "feature")
            disp <- if ("ko_name" %in% colnames(top_b))
                top_b[, .(feature, ko_name, loading_comp1)]
            else
                top_b[, .(feature, loading_comp1)]
            print(disp)
            fwrite(top_b, file.path(OUT_DIR, paste0("top_diablo_", b, ".tsv")), sep = "\t")
            cat("\n")
        }
    }
} else {
    cat("  Not found:", diablo_path, "\n\n")
}

# ---- 4. TOP ANCOM-BC2 FEATURES BY RAW p-VALUE ----------------------------
cat("=== 4. Top ANCOM-BC2 features by raw p-value ===\n\n")
ancom_func_path <- file.path(RES, "ancombc/ancombc_function_results.tsv")
ancom_taxa_path <- file.path(RES, "ancombc/ancombc_taxa_results.tsv")

for (label_file in list(
    list("functions (KO/COG)", ancom_func_path),
    list("taxa", ancom_taxa_path)
)) {
    label <- label_file[[1]]
    p <- label_file[[2]]
    cat("---", label, "---\n")
    if (file.exists(p)) {
        a <- fread(p)
        if ("p" %in% colnames(a)) {
            setorder(a, p)
            top_a <- head(a, TOP_N)
            # Annotate KOs for the function block; skip for taxa
            if (label != "taxa") top_a <- annotate_kos(top_a, "feature")
            disp <- if ("ko_name" %in% colnames(top_a))
                top_a[, .(feature, ko_name, lfc, p, q)]
            else
                top_a[, .(feature, lfc, p, q)]
            print(disp)
            fwrite(top_a, file.path(OUT_DIR,
                                    paste0("top_ancombc_", gsub(" .*", "", label), ".tsv")),
                   sep = "\t")
        } else {
            cat("  No p-value column; analysis may have been skipped\n")
        }
    } else {
        cat("  Not found:", p, "\n")
    }
    cat("\n")
}

# ---- 5. TOP REDUNDANCY (most generalist + most specialist) ---------------
cat("=== 5. Top redundancy: most generalist vs most specialist functions ===\n\n")
red_path <- file.path(RES, "redundancy/redundancy_summary.tsv")

if (file.exists(red_path)) {
    red <- fread(red_path)
    if ("mean_n_taxa" %in% colnames(red)) {
        # Annotate function_name as KO if matches KO regex
        red_ann <- annotate_kos(red, "function_name")

        setorder(red_ann, -mean_n_taxa)
        cat("--- Top", TOP_N, "MOST generalist (high redundancy) ---\n")
        top_gen <- head(red_ann, TOP_N)
        print(top_gen)

        cat("\n--- Top", TOP_N, "MOST specialist (low redundancy) ---\n")
        setorder(red_ann, mean_n_taxa)
        top_spec <- head(red_ann[mean_n_taxa > 0], TOP_N)
        print(top_spec)

        fwrite(top_gen, file.path(OUT_DIR, "top_generalist_functions.tsv"), sep = "\t")
        fwrite(top_spec, file.path(OUT_DIR, "top_specialist_functions.tsv"), sep = "\t")
    }
    cat("\n")
} else {
    cat("  Not found:", red_path, "\n\n")
}

# ---- 6. WRITE COMBINED SUMMARY HEADLINE -----------------------------------
cat("=== 6. Headline summary ===\n\n")
headline <- data.table(
    section = c("Top treatment-correlated module",
                "Top KEGG pathway",
                "Top DIABLO function loading",
                "Top DIABLO taxa loading",
                "Top ANCOM-BC2 function",
                "Top ANCOM-BC2 taxa"),
    finding = c(
        if (exists("top_mods") && nrow(top_mods) > 0)
            sprintf("%s (r=%s, n=%s KOs)",
                    top_mods$module[1], signif(top_mods[[cor_col]][1], 3),
                    if ("N" %in% colnames(top_mods)) top_mods$N[1] else NA) else NA,
        if (exists("top_kegg") && nrow(top_kegg) > 0)
            sprintf("%s in %s (FDR=%s)",
                    top_kegg$name[1], top_kegg$module[1],
                    signif(top_kegg$p.adjust[1], 3)) else NA,
        if (exists("diablo") && nrow(diablo) > 0) {
            f_only <- diablo[block != "taxa"]
            if (nrow(f_only) > 0)
                sprintf("%s (loading=%s)", f_only$feature[1],
                        signif(f_only$loading_comp1[1], 3)) else NA
        } else NA,
        if (exists("diablo") && nrow(diablo) > 0) {
            t_only <- diablo[block == "taxa"]
            if (nrow(t_only) > 0)
                sprintf("%s (loading=%s)", t_only$feature[1],
                        signif(t_only$loading_comp1[1], 3)) else NA
        } else NA,
        # ANCOM-BC2 function: top by raw p-value (q may not pass FDR with small n)
        if (file.exists(ancom_func_path)) {
            af <- fread(ancom_func_path)
            if ("p" %in% colnames(af) && nrow(af) > 0) {
                setorder(af, p)
                af_ann <- annotate_kos(af, "feature")
                nm <- if ("ko_name" %in% colnames(af_ann) &&
                          !is.na(af_ann$ko_name[1]) && nzchar(af_ann$ko_name[1]))
                    sprintf("%s (%s)", af_ann$feature[1], af_ann$ko_name[1])
                else af_ann$feature[1]
                sprintf("%s  LFC=%s  raw p=%s  q=%s",
                        nm,
                        signif(af_ann$lfc[1], 3),
                        signif(af_ann$p[1], 3),
                        signif(af_ann$q[1], 3))
            } else NA
        } else NA,
        # ANCOM-BC2 taxa: top by raw p-value
        if (file.exists(ancom_taxa_path)) {
            at <- fread(ancom_taxa_path)
            if ("p" %in% colnames(at) && nrow(at) > 0) {
                setorder(at, p)
                sprintf("%s  LFC=%s  raw p=%s  q=%s",
                        at$feature[1],
                        signif(at$lfc[1], 3),
                        signif(at$p[1], 3),
                        signif(at$q[1], 3))
            } else NA
        } else NA
    )
)
fwrite(headline, file.path(OUT_DIR, "headline_summary.tsv"), sep = "\t")
print(headline)

cat("\n====================================================\n")
cat("All summary files saved in:", OUT_DIR, "\n")
cat("====================================================\n")
