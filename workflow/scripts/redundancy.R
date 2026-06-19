# =============================================================================
# redundancy.R
# Compute functional redundancy: how many taxa contribute to each function
# per sample. Identify specialists (few contributors) vs generalists (many).
# =============================================================================

suppressPackageStartupMessages({
    library(tidyverse)
    library(data.table)
    library(ggplot2)
})

log_con <- file(snakemake@log[[1]], open = "wt")
sink(log_con, type = "output"); sink(log_con, type = "message")

strat_path <- snakemake@input[["strat"]]
meta_path  <- snakemake@input[["meta"]]
out_table  <- snakemake@output[["table"]]
out_sum    <- snakemake@output[["summary"]]
out_plot   <- snakemake@output[["plot"]]
group_col  <- snakemake@params[["group_col"]]

dir.create(dirname(out_table), recursive = TRUE, showWarnings = FALSE)

strat <- fread(strat_path)
meta  <- fread(meta_path)

sample_cols <- setdiff(colnames(strat), c("feature", "function_name", "taxon"))
cat("[redundancy] Samples:", length(sample_cols), "\n")

# For each function-sample combo, count taxa with abundance > 0
strat_long <- melt(strat, id.vars = c("feature", "function_name", "taxon"),
                   measure.vars = sample_cols,
                   variable.name = "sample", value.name = "abundance")
strat_long[, present := abundance > 0]

red_dt <- strat_long[, .(n_taxa = sum(present)),
                     by = .(function_name, sample)]
red_wide <- dcast(red_dt, function_name ~ sample,
                  value.var = "n_taxa", fill = 0)
fwrite(red_wide, out_table, sep = "\t")
cat("[redundancy] Redundancy table written.\n")

# Per-function summary
red_summary <- red_dt[, .(mean_n_taxa = mean(n_taxa),
                          sd_n_taxa = sd(n_taxa),
                          cv_n_taxa = sd(n_taxa) / mean(n_taxa)),
                      by = function_name][order(-mean_n_taxa)]
fwrite(red_summary, out_sum, sep = "\t")

# ---- Plot: top 25 most redundant + top 25 most specialist functions --------
top_gen <- head(red_summary, 25)[, type := "Generalist (high redundancy)"]
top_spec <- tail(red_summary, 25)[, type := "Specialist (low redundancy)"]
plot_dt <- rbind(top_gen, top_spec)
plot_dt[, function_name := factor(function_name, levels = function_name)]

p <- ggplot(plot_dt, aes(x = function_name, y = mean_n_taxa, fill = type)) +
    geom_col() +
    geom_errorbar(aes(ymin = mean_n_taxa - sd_n_taxa,
                      ymax = mean_n_taxa + sd_n_taxa), width = 0.3) +
    coord_flip() +
    scale_fill_manual(values = c("Generalist (high redundancy)" = "#2c7fb8",
                                  "Specialist (low redundancy)" = "#d95f0e")) +
    labs(title = "Functional redundancy across samples",
         x = NULL, y = "Mean number of contributing taxa", fill = NULL) +
    theme_minimal(base_size = 10) +
    theme(legend.position = "bottom",
          axis.text.y = element_text(size = 8))

ggsave(out_plot, p, width = 9, height = 11, units = "in")
cat("[redundancy] Plot saved:", out_plot, "\n")
cat("[redundancy] Done.\n")

sink(type = "message"); sink(type = "output"); close(log_con)
