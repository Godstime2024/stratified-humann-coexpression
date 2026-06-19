"""
Generate the three new figures for Paper 1 additions:
P1_fig6: PCA + PERMANOVA (functional ordination)
P1_fig7: MetaCyc top differential pathways
P1_fig8: Stratified mcr taxonomic contribution
"""
import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl

mpl.rcParams.update({
    "font.family": "DejaVu Sans",
    "font.size": 11,
    "axes.titlesize": 12,
    "axes.labelsize": 11,
    "xtick.labelsize": 10,
    "ytick.labelsize": 10,
    "figure.dpi": 200,
    "savefig.dpi": 300,
    "savefig.bbox": "tight",
})

CON_COLOR = "#1f77b4"
SCFP_COLOR = "#d62728"
FIG_DIR = "workflow/figures"
os.makedirs(FIG_DIR, exist_ok=True)

# =============================================================================
# Fig P1.6: PCA scatter, CON vs SCFP, with PERMANOVA stats overlaid
# =============================================================================
pca = pd.read_csv("results_ko/permanova_func/pca_scores.tsv", sep="\t")
ve = pd.read_csv("results_ko/permanova_func/pca_variance.tsv", sep="\t")
pc1_pct = ve.loc[0, "variance_explained"] * 100
pc2_pct = ve.loc[1, "variance_explained"] * 100

fig, ax = plt.subplots(figsize=(6, 5))
for grp, color in [("CON", CON_COLOR), ("SCFP", SCFP_COLOR)]:
    sub = pca[pca.treatment == grp]
    ax.scatter(sub.PC1, sub.PC2, s=130, color=color, edgecolor="black",
               linewidth=0.8, label=grp, alpha=0.85)
    for _, r in sub.iterrows():
        ax.annotate(r.sample_id.replace("RUM_KSU_", ""), (r.PC1, r.PC2),
                    fontsize=8, xytext=(6, 4), textcoords="offset points")

ax.axhline(0, color="gray", linewidth=0.5, linestyle="--")
ax.axvline(0, color="gray", linewidth=0.5, linestyle="--")
ax.set_xlabel(f"PC1 ({pc1_pct:.1f}%)")
ax.set_ylabel(f"PC2 ({pc2_pct:.1f}%)")
ax.legend(loc="upper right")
# PERMANOVA inset
ax.text(0.02, 0.02,
        "PERMANOVA (Aitchison)\nF = 0.86\n$R^2$ = 0.080\np = 0.71",
        transform=ax.transAxes, fontsize=10,
        bbox=dict(boxstyle="round,pad=0.5", facecolor="white",
                  edgecolor="gray", alpha=0.9),
        verticalalignment="bottom")
ax.spines[["top", "right"]].set_visible(False)
fig.tight_layout()
fig.savefig(f"{FIG_DIR}/P1_fig6_pca_permanova.png")
plt.close(fig)
print("Saved P1_fig6_pca_permanova.png")

# =============================================================================
# Fig P1.7: MetaCyc top differential pathways
# =============================================================================
mc = pd.read_csv("results_ko/metacyc/metacyc_all_pathways.tsv", sep="\t")
mc_top = mc.head(12).copy()
# Shorten pathway labels
mc_top["short"] = mc_top.pathway.apply(
    lambda s: (s.split(": ", 1)[1] if ": " in s else s)[:55])
mc_top["id"] = mc_top.pathway.apply(lambda s: s.split(": ", 1)[0] if ": " in s else "")
mc_top["label"] = mc_top["short"] + " (" + mc_top["id"] + ")"

fig, ax = plt.subplots(figsize=(11, 5.5))
y = np.arange(len(mc_top))
colors = [SCFP_COLOR if x > 0 else CON_COLOR for x in mc_top.log2_fc]
# Positive log2_fc = UP in SCFP, color SCFP; negative = UP in CON
colors = [SCFP_COLOR if x > 0 else CON_COLOR for x in mc_top.log2_fc]
ax.barh(y, mc_top.log2_fc, color=colors, edgecolor="black", linewidth=0.6)
ax.axvline(0, color="black", lw=0.7)
ax.set_yticks(y)
ax.set_yticklabels(mc_top.label, fontsize=9)
ax.invert_yaxis()
ax.set_xlabel("MetaCyc pathway log2 fold change\n(SCFP relative to CON; positive = UP in SCFP)")
vmin, vmax = mc_top.log2_fc.min(), mc_top.log2_fc.max()
rng = max(vmax - vmin, 0.1)
x_p = vmax + 0.18 * rng
for i, p in enumerate(mc_top.p):
    ax.text(x_p, i, f"p={p:.3g}", va="center", ha="left", fontsize=9)
ax.set_xlim(vmin - 0.15 * rng, x_p + 0.6 * rng)
ax.spines[["top", "right"]].set_visible(False)
fig.tight_layout()
fig.savefig(f"{FIG_DIR}/P1_fig7_metacyc.png")
plt.close(fig)
print("Saved P1_fig7_metacyc.png")

# =============================================================================
# Fig P1.8: Stratified mcr taxonomic contribution
# =============================================================================
strat = pd.read_csv("results_ko/stratified_mcr/per_taxon_per_treatment.tsv", sep="\t")
strat = strat[strat.function_name.isin(["K00399", "K00400", "K00401", "K00402"])].copy()
strat["gene"] = strat.function_name.map({
    "K00399": "mcrA (K00399)",
    "K00400": "mcrC (K00400)",
    "K00401": "mcrB (K00401)",
    "K00402": "mcrG (K00402)"
})

# Stacked bars: x = gene * treatment, stack = taxon contribution
pivot = strat.pivot_table(
    index=["gene", "treatment"], columns="taxon",
    values="mean_abund", aggfunc="sum", fill_value=0
)
pivot = pivot * 1e5  # scale to ppm for readability

# Order taxa by total contribution
taxa_order = pivot.sum(axis=0).sort_values(ascending=False).index.tolist()
pivot = pivot[taxa_order]

# Build category labels: 4 genes x 2 treatments = 8 bars
fig, ax = plt.subplots(figsize=(11, 5.5))
n_groups = 4
n_treat = 2
bar_width = 0.38
gene_names = ["mcrA (K00399)", "mcrB (K00401)", "mcrC (K00400)", "mcrG (K00402)"]
x_pos = np.arange(n_groups)

# Color palette per taxon (use Methanobrevibacter-distinguishing palette)
taxon_colors = {}
palette = plt.cm.tab10.colors
for i, t in enumerate(taxa_order):
    taxon_colors[t] = palette[i % len(palette)]

labeled = set()
for j, treat in enumerate(["CON", "SCFP"]):
    offsets = x_pos + (j - 0.5) * bar_width
    bottoms = np.zeros(n_groups)
    for taxon in taxa_order:
        vals = []
        for gene in gene_names:
            try:
                v = pivot.loc[(gene, treat), taxon]
            except KeyError:
                v = 0
            vals.append(v)
        vals = np.array(vals)
        if np.all(vals == 0):
            bottoms += vals
            continue
        label = (taxon.replace("_", " ")
                 if taxon not in labeled else None)
        labeled.add(taxon)
        edge = "black" if treat == "CON" else "gray"
        ax.bar(offsets, vals, bottom=bottoms, width=bar_width,
               color=taxon_colors[taxon], edgecolor=edge, linewidth=0.6,
               label=label)
        bottoms += vals
    # Track bar tops for treatment labels above the bars
    bar_top_y = {}
    for k, xp in enumerate(offsets):
        bar_top_y[xp] = bottoms[k]
    # Add CON/SCFP labels above each bar
    for xp, top_v in bar_top_y.items():
        ax.text(xp, top_v + 0.07, treat, ha="center", va="bottom",
                fontsize=9, fontweight="bold",
                color=CON_COLOR if treat == "CON" else SCFP_COLOR)

# Gene labels on the x-axis
ax.set_xticks(x_pos)
ax.set_xticklabels(gene_names, fontsize=10)
ax.set_xlabel("")
ax.set_ylabel("Stratified abundance (per 100,000)")

# Single deduped legend
handles, labels = ax.get_legend_handles_labels()
ax.legend(handles, labels, title="Contributing taxon",
          loc="upper left", bbox_to_anchor=(1.01, 1.0),
          fontsize=9, title_fontsize=10, frameon=False)
ax.spines[["top", "right"]].set_visible(False)
fig.tight_layout()
fig.savefig(f"{FIG_DIR}/P1_fig8_stratified_mcr.png")
plt.close(fig)
print("Saved P1_fig8_stratified_mcr.png")

print("\nAll three new figures saved to:", FIG_DIR)
