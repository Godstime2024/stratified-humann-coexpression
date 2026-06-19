# =============================================================================
# Final HTML report
# =============================================================================

rule report:
    input:
        # Preprocess
        strat = f"{OUTDIR}/preprocess/stratified_filtered.tsv",
        unstrat = f"{OUTDIR}/preprocess/unstratified_filtered.tsv",
        summary = f"{OUTDIR}/preprocess/preprocess_summary.tsv",
        # Redundancy
        red_table = f"{OUTDIR}/redundancy/redundancy_table.tsv",
        red_plot = f"{OUTDIR}/redundancy/plot_specialist_vs_generalist.pdf",
        # WGCNA
        wgcna_modules = f"{OUTDIR}/wgcna/module_assignment.tsv",
        wgcna_dendro = f"{OUTDIR}/wgcna/plot_module_dendrogram.pdf",
        wgcna_heatmap = f"{OUTDIR}/wgcna/plot_module_trait_heatmap.pdf",
        # Bipartite
        bip_metrics = f"{OUTDIR}/bipartite/specialization_metrics.tsv",
        bip_plot = f"{OUTDIR}/bipartite/plot_bipartite_network.pdf",
        # DIABLO
        diablo_scores = f"{OUTDIR}/diablo/diablo_scores.pdf",
        diablo_circos = f"{OUTDIR}/diablo/plot_circos.pdf",
        # Procrustes
        proc_result = f"{OUTDIR}/procrustes/procrustes_test.tsv",
        proc_plot = f"{OUTDIR}/procrustes/plot_procrustes.pdf",
        # KEGG enrichment
        kegg_table = f"{OUTDIR}/kegg_enrichment/module_pathway_table.tsv",
        kegg_heatmap = f"{OUTDIR}/kegg_enrichment/plot_module_pathway_heatmap.pdf",
        # ANCOM-BC2
        ancom_func = f"{OUTDIR}/ancombc/ancombc_function_results.tsv",
        ancom_taxa = f"{OUTDIR}/ancombc/ancombc_taxa_results.tsv",
        ancom_vfunc = f"{OUTDIR}/ancombc/plot_volcano_function.pdf",
        ancom_vtaxa = f"{OUTDIR}/ancombc/plot_volcano_taxa.pdf",
        # Top features summary (annotated)
        top_modules = f"{OUTDIR}/top_features_summary/top_modules.tsv",
        top_kegg = f"{OUTDIR}/top_features_summary/top_kegg_pathways.tsv",
        top_headline = f"{OUTDIR}/top_features_summary/headline_summary.tsv",
        # Template
        rmd = "workflow/report/report.Rmd"
    output:
        html = f"{OUTDIR}/report.html"
    log:
        f"{OUTDIR}/logs/report.log"
    conda:
        "../envs/r_microbiome.yaml"
    threads: 1
    shell:
        """
        Rscript -e "proj <- getwd(); \
                    dir.create(file.path(proj, '.rmarkdown_intermediates'), showWarnings = FALSE); \
                    rmarkdown::render(input = file.path(proj, '{input.rmd}'), \
                                      output_file = 'report.html', \
                                      output_dir = file.path(proj, '{OUTDIR}'), \
                                      knit_root_dir = proj, \
                                      intermediates_dir = file.path(proj, '.rmarkdown_intermediates'), \
                                      params = list(outdir = file.path(proj, '{OUTDIR}')))" \
                                      > {log} 2>&1
        """
