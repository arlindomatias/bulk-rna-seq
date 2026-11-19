"""Figures and overlap tables for the cross-dataset DEG table (meta_analysis).

Run through Snakemake (`script:` directive, rule `plots`). Reads results/all_degs.csv and each
project's top_table.csv (to know which genes were tested in each dataset) and writes:
  plots/deg_counts.png         DEGs up/down per dataset
  plots/shared_degs_lfc.png    log2FC of the genes that are DEGs in more than one dataset
  results/pairwise_overlap.csv DEG overlap for each pair of datasets, with a hypergeometric test
Genes are matched across datasets (and species) by upper-case symbol.
"""
import itertools
import math
import os

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

inp, out, cfg = snakemake.input, snakemake.output, snakemake.config  # noqa: F821
datasets = list(cfg["datasets"])

UP, DOWN = "#eb6834", "#2a78d6"
INK2 = "#52514e"
plt.rcParams.update({
    "figure.dpi": 150, "savefig.dpi": 150, "font.size": 10, "axes.titlesize": 11,
    "axes.edgecolor": "#cccccc", "axes.labelcolor": INK2, "xtick.color": INK2, "ytick.color": INK2,
    "axes.spines.top": False, "axes.spines.right": False, "axes.axisbelow": True,
    "legend.frameon": False,
})

degs = pd.read_csv(inp.degs)
degs["key"] = degs["SYMBOL"].str.upper()


def save(fig, path):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    fig.savefig(path, bbox_inches="tight", facecolor="white")
    plt.close(fig)


def hypergeom_sf(x, N, K, n):
    """P(X >= x) for X ~ Hypergeometric(population N, K successes, n draws)."""
    def logpmf(k):
        return (math.lgamma(K + 1) - math.lgamma(k + 1) - math.lgamma(K - k + 1)
                + math.lgamma(N - K + 1) - math.lgamma(n - k + 1) - math.lgamma(N - K - n + k + 1)
                - (math.lgamma(N + 1) - math.lgamma(n + 1) - math.lgamma(N - n + 1)))
    return min(1.0, sum(math.exp(logpmf(k)) for k in range(x, min(K, n) + 1)))


# ---------------------------------------------------------------- pairwise overlap
tested = {}
for name, path in zip(datasets, inp.top_tables):
    tt = pd.read_csv(path, usecols=["SYMBOL"])
    tested[name] = set(tt["SYMBOL"].dropna().str.upper())

rows = []
for a, b in itertools.combinations(datasets, 2):
    universe = tested[a] & tested[b]
    da = degs[(degs.dataset == a) & degs.key.isin(universe)].drop_duplicates("key").set_index("key")
    db = degs[(degs.dataset == b) & degs.key.isin(universe)].drop_duplicates("key").set_index("key")
    shared = da.index.intersection(db.index)
    same = int((np.sign(da.loc[shared, "log2FoldChange"]) == np.sign(db.loc[shared, "log2FoldChange"])).sum())
    N, K, n, x = len(universe), len(da), len(db), len(shared)
    rows.append({
        "dataset_a": a, "dataset_b": b, "genes_tested_in_both": N,
        "degs_a_in_universe": K, "degs_b_in_universe": n, "shared_degs": x,
        "expected_by_chance": K * n / N if N else np.nan,
        "hypergeometric_p": hypergeom_sf(x, N, K, n) if x else 1.0,
        "same_direction": same, "opposite_direction": x - same,
        "shared_genes": ";".join(sorted(da.loc[shared, "SYMBOL"])),
    })
pd.DataFrame(rows).to_csv(out.overlap, index=False)

# ---------------------------------------------------------------- DEG counts
# <figure name="deg_counts">
counts = (degs.groupby(["dataset", "direction"]).size().unstack(fill_value=0)
          .reindex(index=datasets, columns=["up", "down"], fill_value=0))
fig, ax = plt.subplots(figsize=(6, 3.2))
y = np.arange(len(datasets))
ax.barh(y - 0.2, counts["up"], 0.38, color=UP, label="up")
ax.barh(y + 0.2, counts["down"], 0.38, color=DOWN, label="down")
for i, d in enumerate(datasets):
    ax.text(counts.loc[d, "up"], i - 0.2, f" {counts.loc[d, 'up']}", va="center", fontsize=9, color=INK2)
    ax.text(counts.loc[d, "down"], i + 0.2, f" {counts.loc[d, 'down']}", va="center", fontsize=9, color=INK2)
ax.set_yticks(y, datasets)
ax.invert_yaxis()
ax.set_xlabel("DEGs")
ax.grid(axis="x", color="#ebebeb")
ax.legend(loc="lower right")
save(fig, out.deg_counts)
# </figure>

# ---------------------------------------------------------------- shared DEGs
# <figure name="shared_degs_lfc">
shared = degs[degs.n_datasets > 1]
mat = shared.pivot_table(index="key", columns="dataset", values="log2FoldChange", aggfunc="first")
mat = mat.reindex(columns=datasets)
labels = shared.drop_duplicates("key").set_index("key")["SYMBOL"]
order = mat.notna().astype(int).mul([2 ** i for i in range(len(datasets))]).sum(axis=1)
mat = mat.loc[sorted(mat.index, key=lambda k: (order[k], mat.loc[k].mean()))]
fig, ax = plt.subplots(figsize=(4.2, 0.28 * len(mat) + 1.2))
lim = np.nanmax(np.abs(mat.values)) if len(mat) else 1
im = ax.imshow(mat.values, cmap="RdBu_r", vmin=-lim, vmax=lim, aspect="auto")
for (i, j), v in np.ndenumerate(mat.values):
    if not np.isnan(v):
        ax.text(j, i, f"{v:.2f}", ha="center", va="center", fontsize=7,
                color="white" if abs(v) > 0.6 * lim else "black")
ax.set_xticks(range(len(datasets)), datasets, rotation=30, ha="right")
ax.set_yticks(range(len(mat)), [labels[k] for k in mat.index], fontstyle="italic")
ax.tick_params(length=0)
for s in ax.spines.values():
    s.set_visible(False)
fig.colorbar(im, ax=ax, label="log2FC (treated vs control)", shrink=0.6)
save(fig, out.shared_lfc)
# </figure>
