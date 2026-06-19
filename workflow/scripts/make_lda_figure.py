"""Generate LEfSe-style LDA effect-size biomarker figure for MetaCyc pathways."""
import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl

mpl.rcParams.update({
    "font.family": "DejaVu Sans",
    "font.size": 11,
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

df = pd.read_csv("results_ko/metacyc_lda/metacyc_lda_top15_permissive.tsv",
                 sep="\t")
df["short"] = df.pathway.apply(
    lambda s: (s.split(": ", 1)[1] if ": " in s else s)[:55])
df["id"] = df.pathway.apply(
    lambda s: s.split(": ", 1)[0] if ": " in s else "")
df["label"] = df["short"] + " (" + df["id"] + ")"
# Sort so SCFP positive at top, CON negative at bottom (classical LEfSe layout)
df = df.sort_values("lda_score").reset_index(drop=True)

fig, ax = plt.subplots(figsize=(9.5, 6))
y = np.arange(len(df))
colors = [SCFP_COLOR if x > 0 else CON_COLOR for x in df.lda_score]
ax.barh(y, df.lda_score, color=colors, edgecolor="black", linewidth=0.6)
ax.axvline(0, color="black", lw=0.7)
ax.set_yticks(y)
ax.set_yticklabels(df.label, fontsize=9)
ax.set_xlabel("LDA effect size (p $\\leq$ 0.10)")

vmin, vmax = df.lda_score.min(), df.lda_score.max()
rng = vmax - vmin
ax.set_xlim(vmin - 0.10 * rng, vmax + 0.10 * rng)

# Legend
con_patch = plt.Rectangle((0, 0), 1, 1, color=CON_COLOR, label="CON biomarker")
scfp_patch = plt.Rectangle((0, 0), 1, 1, color=SCFP_COLOR, label="SCFP biomarker")
ax.legend(handles=[con_patch, scfp_patch], loc="lower center",
          bbox_to_anchor=(0.5, 1.01), ncol=2, frameon=False)
ax.spines[["top", "right"]].set_visible(False)
fig.tight_layout()
fig.savefig(f"{FIG_DIR}/P1_fig9_metacyc_lda.png")
plt.close(fig)
print("Saved P1_fig9_metacyc_lda.png")
