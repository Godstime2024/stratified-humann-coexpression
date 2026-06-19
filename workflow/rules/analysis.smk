# =============================================================================
# Analysis rules: redundancy, WGCNA, bipartite, DIABLO, Procrustes
# =============================================================================

# ---- 1. Functional redundancy ----------------------------------------------
rule redundancy:
    input:
        strat = f"{OUTDIR}/preprocess/stratified_filtered.tsv",
        meta = config["metadata_tsv"]
    output:
        table = f"{OUTDIR}/redundancy/redundancy_table.tsv",
        summary = f"{OUTDIR}/redundancy/redundancy_summary.tsv",
        plot = f"{OUTDIR}/redundancy/plot_specialist_vs_generalist.pdf"
    params:
        group_col = config["group_column"]
    log:
        f"{OUTDIR}/logs/redundancy.log"
    conda:
        "../envs/r_microbiome.yaml"
    threads: 1
    script:
        "../scripts/redundancy.R"

# ---- 2. WGCNA co-expression modules ----------------------------------------
rule wgcna:
    input:
        unstrat_clr = f"{OUTDIR}/preprocess/unstratified_clr.tsv",
        meta = config["metadata_tsv"]
    output:
        modules = f"{OUTDIR}/wgcna/module_assignment.tsv",
        eigengenes = f"{OUTDIR}/wgcna/module_eigengenes.tsv",
        trait_cor = f"{OUTDIR}/wgcna/module_trait_correlation.tsv",
        dendro = f"{OUTDIR}/wgcna/plot_module_dendrogram.pdf",
        heatmap = f"{OUTDIR}/wgcna/plot_module_trait_heatmap.pdf"
    params:
        power = config["wgcna"]["power"],
        min_mod_size = config["wgcna"]["min_module_size"],
        merge_cut = config["wgcna"]["merge_cut_height"],
        net_type = config["wgcna"]["network_type"],
        group_col = config["group_column"],
        continuous = config["continuous_traits"]
    log:
        f"{OUTDIR}/logs/wgcna.log"
    conda:
        "../envs/r_microbiome.yaml"
    threads: config["threads_default"]
    script:
        "../scripts/wgcna.R"

# ---- 3. Bipartite taxon-function network -----------------------------------
rule bipartite:
    input:
        strat = f"{OUTDIR}/preprocess/stratified_filtered.tsv"
    output:
        matrix = f"{OUTDIR}/bipartite/bipartite_matrix.tsv",
        metrics = f"{OUTDIR}/bipartite/specialization_metrics.tsv",
        plot = f"{OUTDIR}/bipartite/plot_bipartite_network.pdf"
    params:
        top_taxa = config["bipartite"]["top_n_taxa"],
        top_func = config["bipartite"]["top_n_functions"]
    log:
        f"{OUTDIR}/logs/bipartite.log"
    conda:
        "../envs/r_microbiome.yaml"
    threads: 1
    script:
        "../scripts/bipartite.R"

# ---- 4. DIABLO integrative analysis (taxa + functions vs group) ------------
rule diablo:
    input:
        taxa_clr = f"{OUTDIR}/preprocess/taxa_clr.tsv",
        unstrat_clr = f"{OUTDIR}/preprocess/unstratified_clr.tsv",
        meta = config["metadata_tsv"]
    output:
        model = f"{OUTDIR}/diablo/diablo_model.rds",
        loadings = f"{OUTDIR}/diablo/diablo_loadings.tsv",
        scores = f"{OUTDIR}/diablo/diablo_scores.pdf",
        circos = f"{OUTDIR}/diablo/plot_circos.pdf"
    params:
        ncomp = config["diablo"]["ncomp"],
        design_value = config["diablo"]["design_value"],
        cutoff = config["diablo"]["cutoff_circos"],
        group_col = config["group_column"]
    log:
        f"{OUTDIR}/logs/diablo.log"
    conda:
        "../envs/r_microbiome.yaml"
    threads: config["threads_default"]
    script:
        "../scripts/diablo.R"

# ---- 5a. KEGG pathway enrichment per WGCNA module --------------------------
rule kegg_enrichment:
    input:
        modules = f"{OUTDIR}/wgcna/module_assignment.tsv"
    output:
        table = f"{OUTDIR}/kegg_enrichment/module_pathway_table.tsv",
        summary = f"{OUTDIR}/kegg_enrichment/pathway_enrichment_summary.tsv",
        heatmap = f"{OUTDIR}/kegg_enrichment/plot_module_pathway_heatmap.pdf"
    params:
        cache_dir = ".kegg_cache",
        top_n_paths = config.get("kegg", {}).get("top_n_paths_per_module", 10),
        pval_cutoff = config.get("kegg", {}).get("pval_cutoff", 0.05)
    log:
        f"{OUTDIR}/logs/kegg_enrichment.log"
    conda:
        "../envs/r_microbiome.yaml"
    threads: 1
    script:
        "../scripts/kegg_enrichment.R"

# ---- 5b. ANCOM-BC2 differential abundance ----------------------------------
rule ancombc:
    input:
        unstrat = f"{OUTDIR}/preprocess/unstratified_filtered.tsv",
        taxa    = f"{OUTDIR}/preprocess/taxa_collapsed.tsv",
        meta    = config["metadata_tsv"]
    output:
        ko_results   = f"{OUTDIR}/ancombc/ancombc_function_results.tsv",
        taxa_results = f"{OUTDIR}/ancombc/ancombc_taxa_results.tsv",
        volcano_ko   = f"{OUTDIR}/ancombc/plot_volcano_function.pdf",
        volcano_taxa = f"{OUTDIR}/ancombc/plot_volcano_taxa.pdf"
    params:
        group_col  = config["group_column"],
        fdr_cutoff = config.get("ancombc", {}).get("fdr_cutoff", 0.05),
        lfc_cutoff = config.get("ancombc", {}).get("lfc_cutoff", 1.0),
        prv_cut    = config.get("ancombc", {}).get("prv_cut", 0.1)
    log:
        f"{OUTDIR}/logs/ancombc.log"
    conda:
        "../envs/r_microbiome.yaml"
    threads: config["threads_default"]
    script:
        "../scripts/ancombc.R"

# ---- 5c. Top features summary (KO-annotated post-pipeline summary) ---------
rule top_features:
    input:
        modules = f"{OUTDIR}/wgcna/module_assignment.tsv",
        wgcna_trait = f"{OUTDIR}/wgcna/module_trait_correlation.tsv",
        kegg = f"{OUTDIR}/kegg_enrichment/module_pathway_table.tsv",
        diablo = f"{OUTDIR}/diablo/diablo_loadings.tsv",
        ancom_func = f"{OUTDIR}/ancombc/ancombc_function_results.tsv",
        ancom_taxa = f"{OUTDIR}/ancombc/ancombc_taxa_results.tsv",
        red = f"{OUTDIR}/redundancy/redundancy_summary.tsv"
    output:
        modules = f"{OUTDIR}/top_features_summary/top_modules.tsv",
        kegg = f"{OUTDIR}/top_features_summary/top_kegg_pathways.tsv",
        headline = f"{OUTDIR}/top_features_summary/headline_summary.tsv"
    params:
        top_n = config.get("top_n_features", 15),
        outdir = OUTDIR
    log:
        f"{OUTDIR}/logs/top_features.log"
    conda:
        "../envs/r_microbiome.yaml"
    threads: 1
    shell:
        """
        Rscript workflow/scripts/top_features_summary.R {params.outdir} {params.top_n} \
            > {log} 2>&1
        """

# ---- 6. Procrustes test (taxa structure vs function structure) -------------
rule procrustes:
    input:
        taxa_clr = f"{OUTDIR}/preprocess/taxa_clr.tsv",
        unstrat_clr = f"{OUTDIR}/preprocess/unstratified_clr.tsv",
        meta = config["metadata_tsv"]
    output:
        result = f"{OUTDIR}/procrustes/procrustes_test.tsv",
        plot = f"{OUTDIR}/procrustes/plot_procrustes.pdf"
    params:
        permutations = config["procrustes"]["permutations"],
        group_col = config["group_column"]
    log:
        f"{OUTDIR}/logs/procrustes.log"
    conda:
        "../envs/r_microbiome.yaml"
    threads: 1
    script:
        "../scripts/procrustes.R"
