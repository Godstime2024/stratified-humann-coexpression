# stratified-humann-coexpression

A Snakemake workflow for species-stratified functional co-expression analysis
of long-read Nanopore rumen metagenomes. The pipeline integrates HUMAnN3
species-stratified KEGG Orthology profiling, WGCNA coexpression module
discovery, ANCOM-BC2 bias-corrected compositional differential abundance
testing, PERMANOVA on Aitchison distances, KEGG pathway abundance and
enrichment analysis, MetaCyc pathway analysis with LEfSe-style LDA biomarker
discovery, stratified gene-to-taxon contribution analysis, and CAZyme
cross-validation.

This workflow was developed to characterize the mechanism by which
Saccharomyces cerevisiae fermentation product (SCFP) modulates the rumen
microbiome of cattle. It is generalizable to any two-group rumen functional
metagenomics comparison.

## Citation

If you use this pipeline, please cite both the companion paper and the
pipeline release:

> Taiwo G, et al. 2026. Saccharomyces cerevisiae fermentation product
> redirects rumen hydrogen disposal from methanogenesis to Wood-Ljungdahl
> acetogenesis in cattle. mSystems [in preparation].

> Taiwo G. 2026. stratified-humann-coexpression: a Snakemake workflow for
> species-stratified functional metagenomic analysis of Nanopore rumen
> metagenomes (v1.0.0). Zenodo. https://doi.org/10.5281/zenodo.XXXXXXX

## Quick start

```bash
# 1. Clone
git clone https://github.com/Godstime2024/stratified-humann-coexpression.git
cd stratified-humann-coexpression

# 2. Install Snakemake + mamba
conda install -n base -c conda-forge mamba
mamba create -n snakemake -c bioconda -c conda-forge snakemake
conda activate snakemake

# 3. Place inputs
#    data/Specie_stratified_ko_normalized.csv     (or COG equivalent)
#    metadata/sample_metadata.tsv

# 4. Run KO config
snakemake --configfile config/config_ko.yaml --cores 4 --use-conda

# 5. Or COG config
snakemake --configfile config/config_cog.yaml --cores 4 --use-conda
```

## Inputs

| File | Description |
|---|---|
| `data/Specie_stratified_ko_normalized.csv` | Sample by KO matrix, taxon-stratified (rows = `function|taxon`) |
| `metadata/sample_metadata.tsv` | Per-sample metadata: sample_id, treatment, animal_id |
| `config/config_ko.yaml` or `config/config_cog.yaml` | Pipeline configuration |

Metadata schema:

```
sample_id    treatment    animal_id    notes
SAMPLE_01    CON          A01
SAMPLE_02    SCFP         A02
```

## Pipeline outputs

| Directory | Contents |
|---|---|
| `results_ko/preprocess/` | Filtered + CLR-transformed matrices |
| `results_ko/ancombc/` | ANCOM-BC2 differential abundance results (taxa + KOs) |
| `results_ko/wgcna/` | Coexpression modules, eigengenes, module-trait correlations |
| `results_ko/kegg_enrichment/` | Per-module KEGG pathway enrichment |
| `results_ko/kegg_pathway_direct/` | Direct pathway abundance comparison |
| `results_ko/metacyc/` | MetaCyc pathway differential abundance |
| `results_ko/metacyc_lda/` | LEfSe-style LDA biomarker discovery |
| `results_ko/h2_sink/` | H2 sink and SCFA-producing KO panel results |
| `results_ko/stratified_mcr/` | Stratified mcr contribution by taxon |
| `results_ko/permanova_func/` | PERMANOVA + PCA ordination |
| `results_ko/cazyme_compare/` | CAZyme cross-check |
| `results_ko/top_features_summary/` | Annotated top features per analysis |
| `results_ko/amr_compare/` | AMR cross-check (optional) |
| `results_ko/bipartite/` | Bipartite taxon-function network |
| `results_ko/procrustes/` | Procrustes test (taxa vs function structure) |
| `results_ko/report.html` | Auto-generated R Markdown report |

## Analyses included

| Analysis | Implementation |
|---|---|
| Preprocessing (filter + CLR) | `workflow/rules/preprocess.smk` |
| ANCOM-BC2 differential abundance | `workflow/scripts/ancombc.R` |
| WGCNA coexpression modules | `workflow/scripts/wgcna.R` |
| KEGG pathway enrichment per module | `workflow/scripts/kegg_enrichment.R` |
| Direct KEGG pathway abundance test | `workflow/scripts/kegg_pathway_direct_compare.R` |
| MetaCyc pathway analysis | `workflow/scripts/metacyc_compare.R` |
| LEfSe-style LDA on MetaCyc | `workflow/scripts/metacyc_lda.R` |
| H2 sink and SCFA KO panel | `workflow/scripts/h2_sink_compare.R` |
| Stratified mcr taxonomic contribution | `workflow/scripts/stratified_mcr_contribution.R` |
| Functional PERMANOVA + PCA | `workflow/scripts/functional_permanova.R` |
| CAZyme cross-check | `workflow/scripts/cazyme_compare.R` |
| AMR cross-check (optional) | `workflow/scripts/amr_compare.R` |
| Bacteriocin cross-check (optional) | `workflow/scripts/bacteriocin_compare.R` |
| Top features summary with KEGG annotation | `workflow/scripts/top_features_summary.R` |
| MEsalmon module deep-dive | `workflow/scripts/mesalmon_deepdive.R` |
| Bipartite network | `workflow/scripts/bipartite.R` |
| Procrustes test | `workflow/scripts/procrustes.R` |
| HTML report (R Markdown) | `workflow/report/report.Rmd` |
| Paper figure generation | `workflow/scripts/make_paper_figures.py` |

## Run on HPC (SLURM)

```bash
snakemake --cluster "sbatch --time=2:00:00 --mem=16G --cpus-per-task={threads}" \
          --jobs 6 --use-conda \
          --configfile config/config_ko.yaml
```

## Software dependencies

Managed via conda environments in `workflow/envs/`. Snakemake auto-installs
per-rule when invoked with `--use-conda`. Key packages:

- Snakemake ≥ 7.0
- R ≥ 4.4 with: data.table, vegan, WGCNA, ANCOM-BC2, mixOmics, MASS,
  KEGGREST, clusterProfiler, SummarizedExperiment, mia
- Python ≥ 3.10 with: pandas, numpy, matplotlib, python-docx (for figure
  and document generation)
- HUMAnN3 ≥ 3.6 (if reprocessing reads)
- Kraken2 (if reprocessing reads)

## Workflow overview

```
                 ┌─────────────────────┐
                 │  Stratified KO TSV  │
                 │  + Metadata TSV     │
                 └──────────┬──────────┘
                            │
            ┌───────────────┼───────────────┐
            ▼               ▼               ▼
      Preprocess       Functional       Taxonomic
      (CLR + filter)   analysis         analysis
            │               │               │
            └───────┬───────┴───────┬───────┘
                    ▼               ▼
             ANCOM-BC2 +       MicrobiomeAnalyst
             WGCNA +           LEfSe biomarkers
             KEGG pathway      (taxonomic LDA)
             (KO level)
                    │               │
                    └───────┬───────┘
                            ▼
                    Integrated report
                    (HTML + figures + tables)
```

## Companion paper

This pipeline implements the analytical workflow described in:

> Taiwo G, et al. 2026. Saccharomyces cerevisiae fermentation product
> redirects rumen hydrogen disposal from methanogenesis to Wood-Ljungdahl
> acetogenesis in cattle. mSystems [in preparation].

A second manuscript using the same data with an extended AMR-BAGEL resistome
analysis is in preparation.

## License

MIT License. See `LICENSE`.

## Contact

**Godstime Taiwo**
Department of Animal and Range Sciences
New Mexico State University, Las Cruces, New Mexico, USA
gtaiwo@nmsu.edu

## Acknowledgments

This work uses the bioBakery toolkit (HUMAnN3), MicrobiomeAnalyst 2.0,
ANCOM-BC2, WGCNA, vegan, mixOmics, BAGEL, and clusterProfiler frameworks.
Please cite these tools when using the pipeline; full references are in the
companion paper.
