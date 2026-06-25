# /home/tmhaxs421/brannanlab/tmhaxs421/CLIP/TERT_DELTA/analysis/ACLY/presentation.ipynb
# you can copy and paste this code into a .py file and run it, but it's really meant to be run as a Jupyter notebook for easier tweaking of the config and visualization

# ── CONFIG ────────────────────────────────────────────────────────────────────
list_of_files = [
    "/home/tmhaxs421/brannanlab/tmhaxs421/CLIP/TERT_DELTA/TERT_ACLY/09a_annotated_normalized_clipper_peaks/ACLY_IP1_S11_R1_001.peakClusters.normed.compressed.bed.compressed.sorted.annotated.bed",
    "/home/tmhaxs421/brannanlab/tmhaxs421/CLIP/TERT_DELTA/TERT_ACLY/09a_annotated_normalized_clipper_peaks/ACLY_IP2_S12_R1_001.peakClusters.normed.compressed.bed.compressed.sorted.annotated.bed",
    "/home/tmhaxs421/brannanlab/tmhaxs421/CLIP/TERT_DELTA/TERT_ACLY/10_idr/ACLY_IP1_S11_R1_001_ACLY_IP2_S12_R1_001_reproducible_peaks.sorted.annotated.bed",
]

mapping = {
    "ACLY_IP1_S11_R1_001.peakClusters.normed.compressed.bed.compressed.sorted.annotated.bed": "ACLY_IP1",
    "ACLY_IP2_S12_R1_001.peakClusters.normed.compressed.bed.compressed.sorted.annotated.bed": "ACLY_IP2",
    "ACLY_IP1_S11_R1_001_ACLY_IP2_S12_R1_001_reproducible_peaks.sorted.annotated.bed":        "ACLY_IDR",
}

order = ["ACLY_IP1", "ACLY_IP2", "ACLY_IDR"]

save_dir  = "/home/tmhaxs421/brannanlab/tmhaxs421/CLIP/TERT_DELTA/analysis/ACLY"
save_name = "ACLY_feature_type_distribution.png"

pval_threshold  = 3
log2fc_threshold = 3
label_pct_min   = 5.0   # min % width to print a label inside a segment

feature_type_colors = {
    "5utr":                      "#E5F5E0",
    "CDS":                       "#31A354",
    "3utr":                      "#006D2C",
    "proxintron500":              "#9ECAE1",
    "distintron500":              "#6BAED6",
    "proxnoncoding_intron500":   "#2171B5",
    "distnoncoding_intron500":   "#08519C",
    "noncoding_exon":            "#D55E00",
    "miRNA":                     "#CC79A7",
    "intergenic":                "#999999",
    "stop_codon":                "#000000",
}
# ─────────────────────────────────────────────────────────────────────────────

import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker

plt.rcParams.update({
    "font.size":         8,
    "axes.linewidth":    0.8,
    "xtick.major.width": 0.8,
    "ytick.major.width": 0.8,
    "xtick.major.size":  3,
    "ytick.major.size":  3,
    "pdf.fonttype":      42,
    "svg.fonttype":      "none",
})

# ── Load & filter ─────────────────────────────────────────────────────────────
df_list = []
for f in list_of_files:
    print(f"Reading: {f}")
    df = pd.read_csv(f, sep="\t", header=None,
                     names=['chrom', 'start', 'end', '-log10(pval)', 'mean(log2FC)',
                            'strand', 'gene_id', 'gene_name', 'feature_type', 'transcript_id'])
    df['sample'] = os.path.basename(f)
    df = df[(df['-log10(pval)'] > pval_threshold) & (df['mean(log2FC)'] > log2fc_threshold)]
    df_list.append(df)

idr_main_df = pd.concat(df_list, ignore_index=True)
idr_main_df['sample'] = idr_main_df['sample'].map(mapping)
idr_main_df['sample'] = pd.Categorical(idr_main_df['sample'], categories=order, ordered=True)

# ── Build count / percent tables ──────────────────────────────────────────────
ct = (
    pd.crosstab(idr_main_df["sample"], idr_main_df["feature_type"])
      .reindex(order)
      .fillna(0)
      .astype(int)
)

stack_order = [c for c in feature_type_colors if c in ct.columns] + \
              [c for c in ct.columns if c not in feature_type_colors]
ct     = ct.reindex(columns=stack_order)
ct_pct = ct.div(ct.sum(axis=1), axis=0).replace([np.inf, -np.inf], np.nan).fillna(0) * 100

gene_counts = (
    idr_main_df.groupby('sample', observed=True)['gene_id']
    .apply(lambda s: len({
        gene.strip().split('.')[0]
        for genes in s.dropna()
        for gene in genes.split(',')
    }))
    .reindex(order)
)

# Reverse so first sample plots at top
ct_pct_plot       = ct_pct.iloc[::-1].copy()
totals            = ct.iloc[::-1].sum(axis=1).values
gene_counts_plot  = gene_counts.iloc[::-1].values
n_samples         = len(ct_pct_plot)
y                 = np.arange(n_samples)

# ── Plot ──────────────────────────────────────────────────────────────────────
fig, (ax_pct, ax_n, ax_g) = plt.subplots(
    1, 3,
    figsize=(9.5, 0.65 * n_samples + 1.4),
    sharey=True,
    gridspec_kw={"width_ratios": [3.5, 1, 1], "wspace": 0.1},
)

bar_height = 0.55
lefts = np.zeros(n_samples)

for col in ct_pct_plot.columns:
    vals  = ct_pct_plot[col].values
    color = feature_type_colors.get(col, "#BDBDBD")
    ax_pct.barh(y, vals, left=lefts, height=bar_height, color=color, linewidth=0, label=col)
    for i, (v, l) in enumerate(zip(vals, lefts)):
        if v >= label_pct_min:
            ax_pct.text(l + v / 2, i, f"{v:.1f}%",
                        ha="center", va="center", fontsize=6.5, color="black")
    lefts += vals

ax_pct.set_xlim(0, 100)
ax_pct.set_xlabel("Percent of Peaks", fontsize=8, labelpad=4)
ax_pct.spines["top"].set_visible(False)
ax_pct.spines["right"].set_visible(False)
ax_pct.set_yticks(y)
ax_pct.set_yticklabels(ct_pct_plot.index, fontsize=8)
ax_pct.tick_params(axis="y", length=0)

for ax, vals, label in [(ax_n, totals, "# Peaks"), (ax_g, gene_counts_plot, "# Genes")]:
    ax.barh(y, vals, height=bar_height, color="#888888", linewidth=0)
    for i, v in enumerate(vals):
        ax.text(v + vals.max() * 0.02, i, f"{int(v):,}",
                ha="left", va="center", fontsize=7, color="black")
    ax.set_xlabel(label, fontsize=8, labelpad=4)
    ax.set_xlim(0, vals.max() * 1.45)
    ax.xaxis.set_major_locator(mticker.MaxNLocator(3, integer=True))
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["left"].set_visible(False)
    ax.tick_params(axis="y", length=0)

ax_g.legend(
    *ax_pct.get_legend_handles_labels(),
    title="Feature type", title_fontsize=8, fontsize=7,
    bbox_to_anchor=(1.02, 1), loc="upper left",
    frameon=False, handlelength=1.2, handleheight=1.0,
)

fig.savefig(f"{save_dir}/{save_name}", dpi=600, bbox_inches="tight", facecolor="white")
print(f"Saved to {save_dir}/{save_name}")
plt.show()
