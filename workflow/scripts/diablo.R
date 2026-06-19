# =============================================================================
# diablo.R
# DIABLO (Data Integration Analysis for Biomarker discovery using a Latent
# cOmponents method) — cross-block sparse PLS-DA integrating taxa and functions
# against the treatment group.
# =============================================================================

suppressPackageStartupMessages({
    library(tidyverse)
    library(data.table)
    library(mixOmics)
})

log_con <- file(snakemake@log[[1]], open = "wt")
sink(log_con, type = "output"); sink(log_con, type = "message")

taxa_path    <- snakemake@input[["taxa_clr"]]
unstrat_path <- snakemake@input[["unstrat_clr"]]
meta_path    <- snakemake@input[["meta"]]
out_model    <- snakemake@output[["model"]]
out_loads    <- snakemake@output[["loadings"]]
out_scores   <- snakemake@output[["scores"]]
out_circos   <- snakemake@output[["circos"]]
ncomp        <- as.integer(snakemake@params[["ncomp"]])
design_val   <- as.numeric(snakemake@params[["design_value"]])
cutoff       <- as.numeric(snakemake@params[["cutoff"]])
group_col    <- snakemake@params[["group_col"]]

dir.create(dirname(out_model), recursive = TRUE, showWarnings = FALSE)

taxa <- fread(taxa_path)
func <- fread(unstrat_path)
meta <- fread(meta_path)

sample_cols <- intersect(colnames(taxa), colnames(func))
sample_cols <- setdiff(sample_cols, c("taxon", "feature"))
cat("[diablo] Shared samples:", length(sample_cols), "\n")

# Build matrices: samples in rows
taxa_mat <- t(as.matrix(taxa[, ..sample_cols]))
colnames(taxa_mat) <- taxa$taxon
func_mat <- t(as.matrix(func[, ..sample_cols]))
colnames(func_mat) <- func$feature

# Align to metadata + group
meta_align <- meta[match(rownames(taxa_mat), meta[[1]]), ]
Y <- factor(meta_align[[group_col]])
cat("[diablo] Group factor:", paste(levels(Y), collapse = ", "), "\n")

X <- list(taxa = taxa_mat, function_terms = func_mat)
design <- matrix(design_val, ncol = length(X), nrow = length(X),
                 dimnames = list(names(X), names(X)))
diag(design) <- 0

# Choose number of features per block (modest given small n)
keepX <- list(
    taxa = rep(min(30, ncol(taxa_mat) - 1), ncomp),
    function_terms = rep(min(30, ncol(func_mat) - 1), ncomp)
)

set.seed(123)
diablo <- block.splsda(X = X, Y = Y, ncomp = ncomp,
                        keepX = keepX, design = design)

saveRDS(diablo, out_model)

# Loadings (selected features)
load_dt <- list()
for (b in names(X)) {
    L <- selectVar(diablo, block = b, comp = 1)
    if (!is.null(L[[b]]$value)) {
        load_dt[[b]] <- data.table(
            block = b,
            feature = rownames(L[[b]]$value),
            loading_comp1 = L[[b]]$value$value.var,
            stringsAsFactors = FALSE
        )
    }
}
load_dt <- rbindlist(load_dt, fill = TRUE)
fwrite(load_dt, out_loads, sep = "\t")

# Score plot
pdf(out_scores, width = 9, height = 4)
plotIndiv(diablo, ind.names = TRUE, legend = TRUE,
          title = "DIABLO score plots (each block, comp 1 vs comp 2)")
dev.off()

# Circos plot
pdf(out_circos, width = 8, height = 8)
tryCatch(
    circosPlot(diablo, cutoff = cutoff, line = TRUE,
               color.blocks = c("steelblue", "darkorange"),
               color.cor = c("darkred", "darkgreen"),
               size.labels = 0.8),
    error = function(e) {
        plot.new()
        title("Circos plot failed; try lowering cutoff in config")
    }
)
dev.off()

cat("[diablo] Done.\n")
sink(type = "message"); sink(type = "output"); close(log_con)
