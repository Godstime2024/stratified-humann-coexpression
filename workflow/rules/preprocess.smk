# =============================================================================
# Preprocessing rules: split stratified vs unstratified, filter, CLR transform
# =============================================================================

rule preprocess:
    input:
        csv = config["input_csv"],
        meta = config["metadata_tsv"]
    output:
        strat = f"{OUTDIR}/preprocess/stratified_filtered.tsv",
        unstrat = f"{OUTDIR}/preprocess/unstratified_filtered.tsv",
        strat_clr = f"{OUTDIR}/preprocess/stratified_clr.tsv",
        unstrat_clr = f"{OUTDIR}/preprocess/unstratified_clr.tsv",
        taxa_collapsed = f"{OUTDIR}/preprocess/taxa_collapsed.tsv",
        taxa_clr = f"{OUTDIR}/preprocess/taxa_clr.tsv",
        summary = f"{OUTDIR}/preprocess/preprocess_summary.tsv"
    params:
        input_sep = config["input_sep"],
        stratify_sep = config["stratify_sep"],
        min_samples = config["min_samples"],
        min_value = config["min_value"],
        pseudo_strategy = config["pseudocount_strategy"],
        pseudo_fixed = config["pseudocount_fixed"]
    log:
        f"{OUTDIR}/logs/preprocess.log"
    conda:
        "../envs/r_microbiome.yaml"
    threads: config["threads_default"]
    script:
        "../scripts/preprocess.R"
