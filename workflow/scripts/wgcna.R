# =============================================================================
# wgcna.R
# Function-level co-expression network using WGCNA on CLR-transformed
# unstratified function table.
# =============================================================================

suppressPackageStartupMessages({
    library(tidyverse)
    library(data.table)
    library(WGCNA)
    library(pheatmap)
})
options(stringsAsFactors = FALSE)
allowWGCNAThreads()

log_con <- file(snakemake@log[[1]], open = "wt")
sink(log_con, type = "output"); sink(log_con, type = "message")

unstrat_path <- snakemake@input[["unstrat_clr"]]
meta_path    <- snakemake@input[["meta"]]
out_modules  <- snakemake@output[["modules"]]
out_eigen    <- snakemake@output[["eigengenes"]]
out_tcor     <- snakemake@output[["trait_cor"]]
out_dendro   <- snakemake@output[["dendro"]]
out_heatmap  <- snakemake@output[["heatmap"]]

power_param  <- snakemake@params[["power"]]
min_mod_size <- as.integer(snakemake@params[["min_mod_size"]])
merge_cut    <- as.numeric(snakemake@params[["merge_cut"]])
net_type     <- snakemake@params[["net_type"]]
group_col    <- snakemake@params[["group_col"]]
cont_traits  <- snakemake@params[["continuous"]]

dir.create(dirname(out_modules), recursive = TRUE, showWarnings = FALSE)

clr <- fread(unstrat_path)
meta <- fread(meta_path)

feats <- clr$feature
mat <- as.matrix(clr[, -1])
rownames(mat) <- feats
sample_cols <- colnames(mat)

# Samples in rows for WGCNA
expr <- t(mat)
cat("[wgcna] Input dims (samples x features):", dim(expr), "\n")

# Drop near-constant features (avoid SD = 0 errors)
sds <- apply(expr, 2, sd)
expr <- expr[, sds > 1e-6]
cat("[wgcna] After SD filter:", dim(expr), "\n")

# ---- Soft threshold --------------------------------------------------------
if (power_param == "auto") {
    powers <- 1:20
    sft <- pickSoftThreshold(expr, powerVector = powers, networkType = net_type,
                             verbose = 0)
    power_use <- sft$powerEstimate
    if (is.na(power_use)) power_use <- 6
    cat("[wgcna] Picked soft threshold:", power_use, "\n")
} else {
    power_use <- as.integer(power_param)
}

# ---- Block-wise module detection -------------------------------------------
net <- blockwiseModules(
    expr,
    power = power_use,
    TOMType = if (net_type == "signed") "signed" else "unsigned",
    networkType = net_type,
    minModuleSize = min_mod_size,
    mergeCutHeight = merge_cut,
    numericLabels = FALSE,
    saveTOMs = FALSE,
    verbose = 0
)

modules <- data.table(feature = colnames(expr),
                      module = net$colors)
fwrite(modules, out_modules, sep = "\t")
cat("[wgcna] Module assignments:", length(unique(net$colors)), "modules\n")

# Eigengenes
ME <- moduleEigengenes(expr, colors = net$colors)$eigengenes
ME_dt <- data.table(sample = rownames(ME), as.data.table(ME))
fwrite(ME_dt, out_eigen, sep = "\t")

# ---- Module-trait correlation ----------------------------------------------
# Build trait matrix: convert group_col to one-hot, append continuous traits
meta_align <- meta[match(sample_cols, meta[[1]]), ]
group_factor <- factor(meta_align[[group_col]])
group_dummy <- model.matrix(~ group_factor - 1)
colnames(group_dummy) <- levels(group_factor)

trait_mat <- group_dummy
for (ct in cont_traits) {
    if (ct %in% colnames(meta_align)) {
        trait_mat <- cbind(trait_mat, setNames(list(as.numeric(meta_align[[ct]])), ct))
    }
}
trait_mat <- as.matrix(as.data.frame(trait_mat))
rownames(trait_mat) <- meta_align[[1]]

cor_mat <- cor(ME, trait_mat, use = "pairwise.complete.obs")
pval_mat <- corPvalueStudent(cor_mat, nrow(expr))

tcor <- as.data.table(cor_mat, keep.rownames = "module")
fwrite(tcor, out_tcor, sep = "\t")

# ---- Plots -----------------------------------------------------------------
pdf(out_dendro, width = 12, height = 6)
plotDendroAndColors(net$dendrograms[[1]],
                    net$colors[net$blockGenes[[1]]],
                    "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05,
                    main = "Function clustering and module assignment")
dev.off()

text_mat <- paste0(signif(cor_mat, 2), "\n(p=", signif(pval_mat, 2), ")")
dim(text_mat) <- dim(cor_mat)

pdf(out_heatmap, width = 8, height = 8)
par(mar = c(6, 8, 3, 3))
labeledHeatmap(Matrix = cor_mat,
               xLabels = colnames(cor_mat),
               yLabels = rownames(cor_mat),
               ySymbols = rownames(cor_mat),
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = text_mat,
               setStdMargins = FALSE,
               cex.text = 0.7,
               zlim = c(-1, 1),
               main = "Module-trait relationships")
dev.off()

cat("[wgcna] Done.\n")
sink(type = "message"); sink(type = "output"); close(log_con)
