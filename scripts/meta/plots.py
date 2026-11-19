"""Figures and overlap tables for the cross-dataset DEG table (meta_analysis).

Run through Snakemake (`script:` directive, rule `plots`). Reads results/all_degs.csv,
results/gene_map.csv (the genes tested in each dataset) and results/meta_random_effects.csv, and writes:
  plots/deg_counts.png         DEGs up/down per dataset
  plots/shared_degs_lfc.png    log2FC of the genes that are DEGs in more than one dataset
  results/pairwise_overlap.csv DEG overlap for each pair of datasets, with a hypergeometric test
  plots/meta_volcano.png       pooled log2FC vs p of the random-effects meta-analysis
  plots/meta_forest.png        per-dataset and pooled log2FC of the top meta DEGs
Genes are matched across datasets (and species) by `key` from gene_map.csv: the mouse Ensembl ID, through
the 1:1 ortholog for rat genes. Genes without a key are left out of the overlap.
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
degs["label"] = degs["key_symbol"].fillna(degs["SYMBOL"])
gene_map = pd.read_csv(inp.gene_map, dtype=str)


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
tested = {d: set(gene_map.loc[gene_map.dataset == d, "key"].dropna()) for d in datasets}

rows = []
for a, b in itertools.combinations(datasets, 2):
    universe = tested[a] & tested[b]
    da = degs[(degs.dataset == a) & degs.key.isin(universe)].set_index("key")
    db = degs[(degs.dataset == b) & degs.key.isin(universe)].set_index("key")
    assert da.index.is_unique and db.index.is_unique
    shared = da.index.intersection(db.index)
    same = int((np.sign(da.loc[shared, "log2FoldChange"]) == np.sign(db.loc[shared, "log2FoldChange"])).sum())
    N, K, n, x = len(universe), len(da), len(db), len(shared)
    rows.append({
        "dataset_a": a, "dataset_b": b, "genes_tested_in_both": N,
        "degs_a_in_universe": K, "degs_b_in_universe": n, "shared_degs": x,
        "expected_by_chance": K * n / N if N else np.nan,
        "hypergeometric_p": hypergeom_sf(x, N, K, n) if x else 1.0,
        "same_direction": same, "opposite_direction": x - same,
        "shared_genes": ";".join(sorted(da.loc[shared, "label"])),
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
labels = shared.drop_duplicates("key").set_index("key")["label"]
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

# ---------------------------------------------------------------- random-effects meta-analysis
meta = pd.read_csv(inp.meta)
mcfg = cfg["meta"]
DS_COLORS = dict(zip(datasets, ["#2a78d6", "#1baf7a", "#8a5cc2", "#d6a32a", "#d6336c", "#495057"]))

# <figure name="meta_volcano">
fig, ax = plt.subplots(figsize=(6, 4.6))
y = -np.log10(meta["pvalue"].clip(lower=1e-300))
sig = meta["meta_deg"].astype(bool)
ax.scatter(meta.loc[~sig, "estimate"], y[~sig], s=5, color="#c4c4c4", alpha=0.6, linewidths=0, label="not significant")
for d, col in (("up", UP), ("down", DOWN)):
    m = sig & (meta["direction"] == d)
    ax.scatter(meta.loc[m, "estimate"], y[m], s=10, color=col, linewidths=0, label=f"meta DEG {d} ({int(m.sum())})")
padj_line = meta.loc[meta["padj"] < mcfg["padj"], "pvalue"].max()
if pd.notna(padj_line):
    ax.axhline(-np.log10(padj_line), color="#9a9a9a", lw=0.8, ls="--")
for v in (-mcfg["lfc"], mcfg["lfc"]):
    ax.axvline(v, color="#9a9a9a", lw=0.8, ls="--")
# Label the most significant meta DEGs, skipping labels that would overlap one already placed.
placed = []
xr = ax.get_xlim()[1] - ax.get_xlim()[0]
yr = ax.get_ylim()[1] - ax.get_ylim()[0]
for _, r in meta[sig].nsmallest(25, "pvalue").iterrows():
    x0, y0 = r["estimate"], -np.log10(r["pvalue"])
    if any(abs(x0 - x1) < 0.14 * xr and abs(y0 - y1) < 0.035 * yr for x1, y1 in placed):
        continue
    placed.append((x0, y0))
    ax.annotate(r["key_symbol"], (x0, y0), fontsize=7, fontstyle="italic",
                xytext=(3, 2), textcoords="offset points", color=INK2)
ax.set_xlabel("pooled log2FC (random effects)")
ax.set_ylabel("-log10 p")
ax.grid(color="#ebebeb")
ax.legend(loc="upper left", fontsize=8, markerscale=1.5)
save(fig, out.meta_volcano)
# </figure>

# <figure name="meta_forest">
n_top = int(mcfg["forest_top_n"])
top = meta[sig].nsmallest(n_top, "pvalue")
if top.empty:
    top = meta.nsmallest(n_top, "pvalue")
top = top.iloc[::-1]
fig, ax = plt.subplots(figsize=(6.4, 0.34 * len(top) + 1.4))
offsets = np.linspace(-0.27, 0.27, len(datasets) + 1)
in_legend = set()
for i, (_, r) in enumerate(top.iterrows()):
    for j, d in enumerate(datasets):
        lfc, se = r.get(f"log2FC_{d}"), r.get(f"lfcSE_{d}")
        if pd.notna(lfc):
            ax.errorbar(lfc, i + offsets[j], xerr=1.96 * se, fmt="o", ms=3, color=DS_COLORS[d],
                        elinewidth=0.8, capsize=0, label=None if d in in_legend else d)
            in_legend.add(d)
    yy = i + offsets[-1]
    ax.fill([r["ci_lb"], r["estimate"], r["ci_ub"], r["estimate"]], [yy, yy + 0.1, yy, yy - 0.1],
            color="black", label="pooled (95% CI)" if i == 0 else None)
ax.axvline(0, color="#9a9a9a", lw=0.8)
ax.set_yticks(range(len(top)), [f"{s}  (I² {i2:.0f}%)" for s, i2 in zip(top["key_symbol"], top["I2"])],
              fontstyle="italic", fontsize=8)
ax.set_xlabel("log2FC ± 95% CI (treated vs control)")
ax.grid(axis="x", color="#ebebeb")
handles, labels = ax.get_legend_handles_labels()
order = [labels.index("pooled (95% CI)")] + [labels.index(d) for d in datasets if d in labels]
ax.legend([handles[k] for k in order], [labels[k] for k in order], loc="upper center",
          bbox_to_anchor=(0.5, -0.6 / (0.34 * len(top) + 1.4)), ncol=min(len(order), 4), fontsize=8)
save(fig, out.meta_forest)
# </figure>
