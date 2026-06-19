# =============================================================================
# Install Bioconductor packages from within the activated env.
# Run AFTER `mamba env create -f workflow/envs/r_microbiome_lite.yaml`
# Usage:
#   conda activate r_microbiome
#   Rscript workflow/envs/install_bioc_packages.R
# =============================================================================

# Pin Bioconductor 3.18 (matches R 4.3)
if (!"BiocManager" %in% installed.packages()) install.packages("BiocManager")
BiocManager::install(version = "3.18", ask = FALSE, update = TRUE)

pkgs <- c(
    "mixOmics",
    "ComplexHeatmap",
    "ANCOMBC",
    "TreeSummarizedExperiment",
    "KEGGREST",
    "clusterProfiler",
    "mia"
)

for (p in pkgs) {
    cat("=== Installing", p, "===\n")
    tryCatch(
        BiocManager::install(p, ask = FALSE, update = FALSE),
        error = function(e) cat("FAILED:", p, "->", conditionMessage(e), "\n")
    )
}

# Verify
cat("\n=== Verifying ===\n")
for (p in pkgs) {
    ok <- requireNamespace(p, quietly = TRUE)
    cat(sprintf("  %-30s %s\n", p, ifelse(ok, "OK", "MISSING")))
}
