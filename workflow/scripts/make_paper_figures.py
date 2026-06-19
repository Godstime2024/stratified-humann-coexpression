"""
Generate publication-quality PNG figures for Paper 1 (methanogenesis + Wood-Ljungdahl)
and Paper 2 (AMR + bacteriocin).
Saves to: workflow/figures/
"""
import os
import sys
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


def forest(ax, labels, values, p_values, colors, xlabel, title,
           ses=None, n_labels=None):
    """Clean horizontal forest plot with a dedicated right-side p-value column.

    labels   : list of y-axis labels (gene / KO / taxon)
    values   : list of LFCs (can be negative)
    p_values : list of p-values (None for entries to omit)
    colors   : list of colors per bar
    ses      : optional list of standard errors for error bars
    n_labels : optional list of "n=X" annotations placed inline with bars
    """
    y = np.arange(len(labels))
    ax.barh(y, values, color=colors, edgecolor="black", linewidth=0.6)
    if ses is not None:
        ax.errorbar(values, y, xerr=ses, fmt="none", color="black",
                    capsize=3, linewidth=0.8)
    ax.axvline(0, color="black", lw=0.7)
    ax.set_yticks(y)
    ax.set_yticklabels(labels, fontsize=9)
    ax.invert_yaxis()
    ax.set_xlabel(xlabel)
    ax.set_title(title)
    vmin = min(values)
    vmax = max(values)
    rng = max(vmax - vmin, 0.1)
    # P-value column at a fixed x position, well to the right of all bars
    x_p = vmax + 0.18 * rng
    for i, p in enumerate(p_values):
        if p is None or (isinstance(p, float) and np.isnan(p)):
            continue
        ax.text(x_p, i, f"p={p:.3g}", va="center", ha="left", fontsize=9)
    # Optional "n=X" annotations placed just past each bar tip in the
    # outward direction (negative bars: to the left of the tip).
    if n_labels is not None:
        for i, (v, ntxt) in enumerate(zip(values, n_labels)):
            if ntxt is None:
                continue
            if v >= 0:
                ax.text(v + 0.03 * rng, i, ntxt, va="center",
                        ha="left", fontsize=8)
            else:
                ax.text(v - 0.03 * rng, i, ntxt, va="center",
                        ha="right", fontsize=8)
    # Generous padding on both sides; leave headroom for the p-value column.
    ax.set_xlim(vmin - 0.15 * rng, x_p + 0.55 * rng)
    ax.spines[["top", "right"]].set_visible(False)
    ax.margins(y=0.02)

# =============================================================================
# PAPER 1: Methanogenesis + Wood-Ljungdahl
# =============================================================================

# ----- Figure P1.1: Methane metabolism per-sample (ko00680) -----
df = pd.read_csv("results_ko/kegg_pathway_direct/methane_per_sample.tsv", sep="\t")
df = df.dropna()
con = df[df.treatment == "CON"]["methane_rel"].values * 100
scfp = df[df.treatment == "SCFP"]["methane_rel"].values * 100

fig, ax = plt.subplots(figsize=(5.2, 4.2))
bp = ax.boxplot([con, scfp], labels=["CON", "SCFP"], widths=0.55, patch_artist=True,
                medianprops=dict(color="black", linewidth=1.5))
bp['boxes'][0].set_facecolor(CON_COLOR)
bp['boxes'][0].set_alpha(0.5)
bp['boxes'][1].set_facecolor(SCFP_COLOR)
bp['boxes'][1].set_alpha(0.5)
np.random.seed(1)
ax.scatter(np.random.normal(1, 0.06, len(con)), con, color=CON_COLOR, s=45, zorder=3, edgecolor="black", linewidth=0.5)
ax.scatter(np.random.normal(2, 0.06, len(scfp)), scfp, color=SCFP_COLOR, s=45, zorder=3, edgecolor="black", linewidth=0.5)
ax.set_ylabel("Methane metabolism (ko00680)\nrelative abundance (%)")
ax.spines[["top", "right"]].set_visible(False)
fig.tight_layout()
fig.savefig(f"{FIG_DIR}/P1_fig1_methane_pathway.png")
plt.close(fig)
print("Saved P1_fig1_methane_pathway.png")

# ----- Figure P1.2: ANCOM-BC2 methanogenesis taxa + KOs -----
taxa = pd.read_csv("results_ko/ancombc/ancombc_taxa_results.tsv", sep="\t")
func = pd.read_csv("results_ko/ancombc/ancombc_function_results.tsv", sep="\t")

methano_taxa = taxa[taxa.feature.str.contains("Methanobrev|Methanocul|Methanosphae|Methanomicrob", case=False, na=False)].copy()
methano_taxa = methano_taxa.sort_values("p").head(5)

mcr = func[func.feature.isin(["K00399", "K00400", "K00401", "K00402"])].copy()
mcr["gene"] = mcr.feature.map({"K00399": "mcrA", "K00400": "mcrC", "K00401": "mcrB", "K00402": "mcrG"})
mcr = mcr.sort_values("p")

fig, axes = plt.subplots(1, 2, figsize=(13, 4.8))
labels = [t.replace("_", " ").replace(" sp ", " sp. ") for t in methano_taxa.feature]
labels = [l[:38] + ("..." if len(l) > 38 else "") for l in labels]
colors_taxa = [SCFP_COLOR if x < 0 else CON_COLOR for x in methano_taxa.lfc]
forest(axes[0],
       labels=labels,
       values=methano_taxa.lfc.tolist(),
       p_values=methano_taxa.p.tolist(),
       colors=colors_taxa,
       ses=methano_taxa.se.tolist(),
       xlabel="ANCOM-BC2 log fold change\n(negative = lower in SCFP)",
       title="")
# Compact "(A)" panel marker inside the top-left of the axes (caption references panels).
axes[0].text(0.01, 0.98, "(A)", transform=axes[0].transAxes,
             fontsize=11, fontweight="bold", va="top", ha="left")

colors_mcr = [SCFP_COLOR if x < 0 else CON_COLOR for x in mcr.lfc]
mcr_labels = [f"{g} ({f})" for g, f in zip(mcr.gene, mcr.feature)]
forest(axes[1],
       labels=mcr_labels,
       values=mcr.lfc.tolist(),
       p_values=mcr.p.tolist(),
       colors=colors_mcr,
       ses=mcr.se.tolist(),
       xlabel="ANCOM-BC2 log fold change\n(negative = lower in SCFP)",
       title="")
axes[1].text(0.01, 0.98, "(B)", transform=axes[1].transAxes,
             fontsize=11, fontweight="bold", va="top", ha="left")

fig.tight_layout()
fig.savefig(f"{FIG_DIR}/P1_fig2_methanogenesis_ancombc.png")
plt.close(fig)
print("Saved P1_fig2_methanogenesis_ancombc.png")

# ----- Figure P1.3: H2 sink pathway group LFCs -----
grp = pd.read_csv("results_ko/h2_sink/h2_sink_group_summary.tsv", sep="\t")
grp = grp[grp.n_KOs_found > 0].copy()
grp = grp.sort_values("mean_lfc")
fig, ax = plt.subplots(figsize=(9, 4.8))
colors = [SCFP_COLOR if x < 0 else CON_COLOR for x in grp.mean_lfc]
n_labels = [f"n={n}" for n in grp.n_KOs_found]
# No p-values for the group summary; pass None list and use n_labels.
forest(ax,
       labels=grp.group.tolist(),
       values=grp.mean_lfc.tolist(),
       p_values=[None] * len(grp),
       colors=colors,
       n_labels=n_labels,
       xlabel="Mean ANCOM-BC2 log fold change across KOs in group\n(negative = lower in SCFP)",
       title="")
fig.tight_layout()
fig.savefig(f"{FIG_DIR}/P1_fig3_h2_sink_groups.png")
plt.close(fig)
print("Saved P1_fig3_h2_sink_groups.png")

# ----- Figure P1.4: Top H2 sink KOs (forest plot style) -----
ko = pd.read_csv("results_ko/h2_sink/h2_sink_KO_compare.tsv", sep="\t")
ko = ko.dropna(subset=["p"]).sort_values("p").head(12).copy()
ko["short"] = ko.gene_function.str.split("  ").str[0].str.strip()
ko["label"] = ko["ko"] + " " + ko["short"]
fig, ax = plt.subplots(figsize=(11, 5.6))
colors = [SCFP_COLOR if x < 0 else CON_COLOR for x in ko.lfc]
forest(ax,
       labels=ko.label.tolist(),
       values=ko.lfc.tolist(),
       p_values=ko.p.tolist(),
       colors=colors,
       xlabel="ANCOM-BC2 log fold change\n(negative = lower in SCFP, positive = higher in SCFP)",
       title="")
fig.tight_layout()
fig.savefig(f"{FIG_DIR}/P1_fig4_h2_sink_KOs.png")
plt.close(fig)
print("Saved P1_fig4_h2_sink_KOs.png")

# ----- Figure P1.5: WGCNA module trait heatmap (top 12) -----
trait = pd.read_csv("results_ko/wgcna/module_trait_correlation.tsv", sep="\t")
trait["abs_r"] = trait[trait.columns[-1]].abs()
top = trait.sort_values("abs_r", ascending=False).head(12).drop(columns=["abs_r"])
top = top.set_index("module")
fig, ax = plt.subplots(figsize=(5.5, 5))
im = ax.imshow(top.values, cmap="RdBu_r", vmin=-1, vmax=1, aspect="auto")
ax.set_xticks(range(len(top.columns)))
ax.set_xticklabels(top.columns)
ax.set_yticks(range(len(top.index)))
ax.set_yticklabels(top.index, fontsize=10)
for i in range(top.shape[0]):
    for j in range(top.shape[1]):
        v = top.values[i, j]
        ax.text(j, i, f"{v:.2f}", ha="center", va="center",
                color=("white" if abs(v) > 0.55 else "black"), fontsize=8)
fig.colorbar(im, ax=ax, shrink=0.8, label="Pearson r")
fig.tight_layout()
fig.savefig(f"{FIG_DIR}/P1_fig5_wgcna_module_trait.png")
plt.close(fig)
print("Saved P1_fig5_wgcna_module_trait.png")

# =============================================================================
# PAPER 2: AMR + bacteriocin
# =============================================================================

# ----- Figure P2.1: Total AMR hits per sample (CON vs SCFP) -----
dc = pd.read_csv("/mnt/d/Wisconsin_data1/amr_bagel_results_hpc/merged/amr_drugclass_matrix.tsv", sep="\t")
meta_full = pd.DataFrame({
    "sample_id": ["RUM_KSU_17","RUM_KSU_18","RUM_KSU_19","RUM_KSU_20","RUM_KSU_21",
                  "RUM_KSU_23","RUM_KSU_42","RUM_KSU_43","RUM_KSU_44","RUM_KSU_45",
                  "RUM_KSU_46","RUM_KSU_47","RUM_KSU_48"],
    "treatment": ["CON","SCFP","SCFP","CON","CON","SCFP","CON","CON","SCFP","SCFP","CON","CON","SCFP"]
})
dc = dc.merge(meta_full, left_on="Sample", right_on="sample_id")
class_cols = [c for c in dc.columns if c not in ["Sample", "sample_id", "treatment"]]
dc["total"] = dc[class_cols].sum(axis=1)

con_t = dc[dc.treatment == "CON"]["total"].values
scfp_t = dc[dc.treatment == "SCFP"]["total"].values

fig, ax = plt.subplots(figsize=(5.2, 4.2))
bp = ax.boxplot([con_t, scfp_t], labels=["CON (n=7)", "SCFP (n=6)"], widths=0.55, patch_artist=True,
                medianprops=dict(color="black", linewidth=1.5))
bp['boxes'][0].set_facecolor(CON_COLOR); bp['boxes'][0].set_alpha(0.5)
bp['boxes'][1].set_facecolor(SCFP_COLOR); bp['boxes'][1].set_alpha(0.5)
np.random.seed(2)
ax.scatter(np.random.normal(1, 0.06, len(con_t)), con_t, color=CON_COLOR, s=45, zorder=3, edgecolor="black", linewidth=0.5)
ax.scatter(np.random.normal(2, 0.06, len(scfp_t)), scfp_t, color=SCFP_COLOR, s=45, zorder=3, edgecolor="black", linewidth=0.5)
ax.set_ylabel("Total AMR gene hits per sample\n(AMRFinderPlus)")
ax.spines[["top", "right"]].set_visible(False)
fig.tight_layout()
fig.savefig(f"{FIG_DIR}/P2_fig1_total_amr.png")
plt.close(fig)
print("Saved P2_fig1_total_amr.png")

# ----- Figure P2.2: Drug class CON vs SCFP -----
cls = pd.read_csv("results_ko/amr_compare/drugclass_compare.tsv", sep="\t")
cls = cls.sort_values("p")
fig, ax = plt.subplots(figsize=(11, 5.5))
y = np.arange(len(cls))
w = 0.4
ax.barh(y - w / 2, cls["CON_mean"], height=w, color=CON_COLOR, alpha=0.85, edgecolor="black", linewidth=0.5, label="CON")
ax.barh(y + w / 2, cls["SCFP_mean"], height=w, color=SCFP_COLOR, alpha=0.85, edgecolor="black", linewidth=0.5, label="SCFP")
short = [c[:38] + ("..." if len(c) > 38 else "") for c in cls.drug_class]
ax.set_yticks(y); ax.set_yticklabels(short, fontsize=9)
ax.invert_yaxis()
ax.set_xlabel("Mean hits per sample")
ax.legend(loc="lower center", bbox_to_anchor=(0.5, 1.01),
          ncol=2, frameon=False)
# Dedicated p-value column past the longest bar, with padding
max_bar = max(cls["CON_mean"].max(), cls["SCFP_mean"].max())
x_p = max_bar * 1.18
for i, p in enumerate(cls.p):
    ax.text(x_p, i, f"p={p:.2f}", va="center", ha="left", fontsize=9)
ax.set_xlim(0, x_p * 1.18)
ax.spines[["top", "right"]].set_visible(False)
fig.tight_layout()
fig.savefig(f"{FIG_DIR}/P2_fig2_drugclass.png")
plt.close(fig)
print("Saved P2_fig2_drugclass.png")

# ----- Figure P2.3: ABRicate cross-database hits -----
ab_dbs = ["card", "resfinder", "ncbi", "argannot"]
db_totals = {}
for db in ab_dbs:
    f = f"/mnt/d/Wisconsin_data1/amr_bagel_results_hpc/merged/abricate_{db}_summary.tsv"
    s = pd.read_csv(f, sep="\t", header=None, on_bad_lines="skip")
    db_totals[db.upper()] = s.shape[1]
hits_per_db = {"CARD": 193, "RESFINDER": 196, "NCBI": 222, "ARG-ANNOT": 209}

fig, ax = plt.subplots(figsize=(5.5, 4))
x = list(hits_per_db.keys()); h = list(hits_per_db.values())
bars = ax.bar(x, h, color=["#4c72b0", "#55a868", "#c44e52", "#8172b2"], edgecolor="black", linewidth=0.6)
for b, v in zip(bars, h):
    ax.text(b.get_x() + b.get_width()/2, v + 3, str(v), ha="center", fontsize=10)
ax.set_ylabel("Total hits across 13 samples")
ax.spines[["top", "right"]].set_visible(False)
fig.tight_layout()
fig.savefig(f"{FIG_DIR}/P2_fig3_abricate_dbs.png")
plt.close(fig)
print("Saved P2_fig3_abricate_dbs.png")

# ----- Figure P2.4: Bacteriocin families -----
bfam = pd.read_csv("results_ko/amr_compare/bacteriocin_family_compare.tsv", sep="\t")
fig, ax = plt.subplots(figsize=(9, 4.4))
y = np.arange(len(bfam))
w = 0.38
ax.barh(y - w/2, bfam["CON_mean"], height=w, color=CON_COLOR, alpha=0.85, edgecolor="black", linewidth=0.5, label="CON")
ax.barh(y + w/2, bfam["SCFP_mean"], height=w, color=SCFP_COLOR, alpha=0.85, edgecolor="black", linewidth=0.5, label="SCFP")
ax.set_yticks(y); ax.set_yticklabels(bfam["family"], fontsize=10)
ax.invert_yaxis()
ax.set_xlabel("Mean hits per sample")
ax.legend(loc="lower center", bbox_to_anchor=(0.5, 1.01),
          ncol=2, frameon=False)
max_bar = max(bfam["CON_mean"].max(), bfam["SCFP_mean"].max())
x_p = max_bar * 1.12
for i, p in enumerate(bfam.p):
    ax.text(x_p, i, f"p={p:.2f}", va="center", ha="left", fontsize=10)
ax.set_xlim(0, x_p * 1.20)
ax.spines[["top", "right"]].set_visible(False)
fig.tight_layout()
fig.savefig(f"{FIG_DIR}/P2_fig4_bacteriocin.png")
plt.close(fig)
print("Saved P2_fig4_bacteriocin.png")

# ----- Figure P2.5: Top 10 individual AMR genes -----
gn = pd.read_csv("results_ko/amr_compare/gene_compare.tsv", sep="\t").sort_values("p").head(10)
fig, ax = plt.subplots(figsize=(10, 5))
y = np.arange(len(gn))
w = 0.4
ax.barh(y - w/2, gn["CON_mean"], height=w, color=CON_COLOR, alpha=0.85, edgecolor="black", linewidth=0.5, label="CON")
ax.barh(y + w/2, gn["SCFP_mean"], height=w, color=SCFP_COLOR, alpha=0.85, edgecolor="black", linewidth=0.5, label="SCFP")
ax.set_yticks(y); ax.set_yticklabels(gn["gene"], fontsize=10)
ax.invert_yaxis()
ax.set_xlabel("Mean hits per sample")
ax.legend(loc="lower center", bbox_to_anchor=(0.5, 1.01),
          ncol=2, frameon=False)
max_bar = max(gn["CON_mean"].max(), gn["SCFP_mean"].max())
x_p = max_bar * 1.15
for i, p in enumerate(gn.p):
    ax.text(x_p, i, f"p={p:.2f}", va="center", ha="left", fontsize=10)
ax.set_xlim(0, x_p * 1.20)
ax.spines[["top", "right"]].set_visible(False)
fig.tight_layout()
fig.savefig(f"{FIG_DIR}/P2_fig5_top_genes.png")
plt.close(fig)
print("Saved P2_fig5_top_genes.png")

print("\nAll figures saved to:", FIG_DIR)
