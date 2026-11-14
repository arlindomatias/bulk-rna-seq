"""QC figures for the raw-read processing step (FastQC, Salmon and the tximport count table).

Usage: python3 qc_plots.py <samples.tsv> <qc_summary.tsv> <counts_gene.tsv> <tx2gene.tsv> \
           <fastqc_dir> <plots_dir> <results_dir>

Writes PNG figures to <plots_dir> and the tables behind them to <results_dir>.
"""
import os
import re
import sys
import zipfile

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

samples_tsv, qc_tsv, counts_tsv, tx2gene_tsv, fastqc_dir, plots_dir, results_dir = sys.argv[1:8]
os.makedirs(plots_dir, exist_ok=True)
os.makedirs(results_dir, exist_ok=True)

# Categorical slots 1-2 of the validated palette (pass the CVD checks); neutral grey for "rest".
GROUP_ORDER = ["Saline", "AC-T"]
GROUP_COL = {"Saline": "#2a78d6", "AC-T": "#eb6834"}
CAT_COLS = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4"]
GREY = "#b5b4ae"
INK, INK2 = "#0b0b0b", "#52514e"
plt.rcParams.update({
    "figure.dpi": 150, "savefig.dpi": 150, "font.size": 10, "axes.titlesize": 11,
    "axes.edgecolor": "#cccccc", "axes.labelcolor": INK2, "xtick.color": INK2, "ytick.color": INK2,
    "axes.spines.top": False, "axes.spines.right": False, "axes.grid": True,
    "grid.color": "#ebebeb", "grid.linewidth": 0.8, "axes.axisbelow": True,
    "legend.frameon": False, "lines.linewidth": 1.6,
})


def save(fig, name):
    fig.savefig(os.path.join(plots_dir, name), bbox_inches="tight", facecolor="white")
    plt.close(fig)


def group_legend(ax, **kw):
    handles = [plt.Line2D([], [], marker="o", ls="", color=GROUP_COL[g], ms=7, label=g) for g in GROUP_ORDER]
    ax.legend(handles=handles, **kw)


samples = pd.read_csv(samples_tsv, sep="\t", dtype={"mouse_id": str})
samples["group"] = pd.Categorical(samples["condition"], GROUP_ORDER, ordered=True)
samples = samples.sort_values(["group", "mouse_id"]).reset_index(drop=True)
samples["label"] = samples["mouse_id"] + " · " + samples["condition"]
order = samples["run"].tolist()
label = dict(zip(samples["run"], samples["label"]))
group = dict(zip(samples["run"], samples["condition"]))

# ---------------------------------------------------------------- Salmon: reads and mapping
qc = pd.read_csv(qc_tsv, sep="\t").set_index("run").loc[order].reset_index()
qc["unmapped"] = qc["num_processed"] - qc["num_mapped"]

# <figure name="qc_reads_mapping">
fig, ax = plt.subplots(figsize=(8, 4.6))
y = np.arange(len(qc))[::-1]
ax.barh(y, qc["num_mapped"] / 1e6, color="#2a78d6", height=0.7, label="Mapped (Salmon)",
        edgecolor="white", linewidth=1)
ax.barh(y, qc["unmapped"] / 1e6, left=qc["num_mapped"] / 1e6, color=GREY, height=0.7,
        label="Not mapped", edgecolor="white", linewidth=1)
ax.set_yticks(y, [label[r] for r in qc["run"]])
ax.set_xlabel("Read pairs (millions)")
ax.grid(axis="y", visible=False)
for yi, (tot, pct) in zip(y, zip(qc["num_processed"] / 1e6, qc["percent_mapped"])):
    ax.text(tot + 0.4, yi, f"{pct:.1f}%", va="center", fontsize=8.5, color=INK2)
ax.set_xlim(0, qc["num_processed"].max() / 1e6 * 1.12)
ax.legend(loc="upper center", bbox_to_anchor=(0.5, 1.1), ncol=2)
save(fig, "qc_reads_mapping.png")
# </figure>

# ---------------------------------------------------------------- FastQC
STATUS_COL = {"PASS": "#0ca30c", "WARN": "#fab219", "FAIL": "#d03b3b"}


def read_fastqc(zpath):
    """Returns ({module: rows}, {module: status}) from a FastQC zip."""
    with zipfile.ZipFile(zpath) as z:
        root = z.namelist()[0].split("/")[0]
        data = z.read(f"{root}/fastqc_data.txt").decode()
        summary = z.read(f"{root}/summary.txt").decode()
    modules, cur, header = {}, None, None
    for line in data.splitlines():
        if line.startswith(">>END_MODULE"):
            cur = None
        elif line.startswith(">>"):
            cur = line[2:].split("\t")[0]
            modules[cur] = []
            header = None
        elif cur is not None:
            if line.startswith("#"):
                if header is None:
                    header = line[1:].split("\t")
                    modules[cur + "#header"] = header
            else:
                modules[cur].append(line.split("\t"))
    status = {l.split("\t")[1]: l.split("\t")[0] for l in summary.splitlines() if l.strip()}
    return modules, status


def pos_mid(b):
    nums = [int(x) for x in re.findall(r"\d+", b)]
    return sum(nums) / len(nums)


fq = {}
for run in order:
    for read in (1, 2):
        z = os.path.join(fastqc_dir, f"{run}_{read}_fastqc.zip")
        if os.path.exists(z):
            fq[(run, read)] = read_fastqc(z)

if fq:
    # <figure name="fastqc_per_base_quality">
    # Mean quality per position
    fig, axes = plt.subplots(1, 2, figsize=(11, 3.8), sharey=True)
    for read, ax in zip((1, 2), axes):
        for run in order:
            if (run, read) not in fq:
                continue
            rows = fq[(run, read)][0]["Per base sequence quality"]
            x = [pos_mid(r[0]) for r in rows]
            ax.plot(x, [float(r[1]) for r in rows], color=GROUP_COL[group[run]], alpha=0.75, lw=1.2)
        ax.axhspan(28, 42, color="#0ca30c", alpha=0.06, lw=0)
        ax.set_title(f"Read {read}", loc="left", color=INK)
        ax.set_xlabel("Position in read (bp)")
    axes[0].set_ylabel("Mean Phred quality")
    axes[0].set_ylim(0, 41)
    group_legend(axes[1], loc="lower left")
    save(fig, "fastqc_per_base_quality.png")
    # </figure>

    # <figure name="fastqc_gc_content">
    # GC content distribution
    fig, axes = plt.subplots(1, 2, figsize=(11, 3.8), sharey=True)
    for read, ax in zip((1, 2), axes):
        for run in order:
            if (run, read) not in fq:
                continue
            rows = fq[(run, read)][0]["Per sequence GC content"]
            gc = np.array([float(r[0]) for r in rows])
            cnt = np.array([float(r[1]) for r in rows])
            ax.plot(gc, 100 * cnt / cnt.sum(), color=GROUP_COL[group[run]], alpha=0.75, lw=1.2)
        ax.set_title(f"Read {read}", loc="left", color=INK)
        ax.set_xlabel("Mean GC content of the read (%)")
    axes[0].set_ylabel("Reads (%)")
    group_legend(axes[1], loc="upper right")
    save(fig, "fastqc_gc_content.png")
    # </figure>

    # <figure name="fastqc_adapter_content">
    # Adapter content (max over adapter types per position)
    fig, axes = plt.subplots(1, 2, figsize=(11, 3.4), sharey=True)
    for read, ax in zip((1, 2), axes):
        for run in order:
            if (run, read) not in fq:
                continue
            rows = fq[(run, read)][0]["Adapter Content"]
            x = [pos_mid(r[0]) for r in rows]
            ax.plot(x, [max(float(v) for v in r[1:]) for r in rows], color=GROUP_COL[group[run]], alpha=0.75, lw=1.2)
        ax.set_title(f"Read {read}", loc="left", color=INK)
        ax.set_xlabel("Position in read (bp)")
    axes[0].set_ylabel("Reads with adapter (%)")
    axes[0].set_ylim(bottom=0)
    group_legend(axes[1], loc="upper left")
    save(fig, "fastqc_adapter_content.png")
    # </figure>

    # <figure name="fastqc_module_status">
    # Module status grid (PASS / WARN / FAIL), with letters so it never relies on colour alone
    keys = [k for k in fq]
    mods = list(fq[keys[0]][1].keys())
    fig, ax = plt.subplots(figsize=(11, 0.42 * len(mods) + 1.6))
    for j, k in enumerate(keys):
        for i, m in enumerate(mods):
            st = fq[k][1].get(m, "")
            ax.add_patch(plt.Rectangle((j + 0.04, i + 0.04), 0.92, 0.92, color=STATUS_COL.get(st, "#eeeeee"), lw=0))
            ax.text(j + 0.5, i + 0.5, st[:1], ha="center", va="center", fontsize=7.5,
                    color="white" if st in ("PASS", "FAIL") else INK)
    ax.set_xlim(0, len(keys))
    ax.set_ylim(len(mods), 0)
    ax.set_xticks(np.arange(len(keys)) + 0.5, [f"{label[r].split(' ·')[0]} R{rd}" for r, rd in keys],
                  rotation=90, fontsize=8)
    ax.set_yticks(np.arange(len(mods)) + 0.5, mods, fontsize=8.5)
    ax.grid(False)
    for s in ax.spines.values():
        s.set_visible(False)
    ax.tick_params(length=0)
    handles = [plt.Rectangle((0, 0), 1, 1, color=c, label=f"{s[0]} = {s}") for s, c in STATUS_COL.items()]
    ax.legend(handles=handles, loc="upper center", bbox_to_anchor=(0.5, 1.0 + 1.0 / (len(mods) + 3)), ncol=3)
    save(fig, "fastqc_module_status.png")
    # </figure>

    stat_rows = [{"run": r, "read": rd, "module": m, "status": s}
                 for (r, rd), (_, st) in fq.items() for m, s in st.items()]
    pd.DataFrame(stat_rows).to_csv(os.path.join(results_dir, "fastqc_module_status.tsv"), sep="\t", index=False)

# ---------------------------------------------------------------- count table
counts = pd.read_csv(counts_tsv, sep="\t").set_index("gene_id")
names = counts.pop("gene_name")
counts = counts[order]

# Gene category from transcript biotypes (the GENCODE FASTA header has no gene biotype).
tx = pd.read_csv(tx2gene_tsv, sep="\t")
bt = tx.groupby("gene_id")["transcript_biotype"].agg(set)


def category(gene_id):
    if str(names.get(gene_id, "")).startswith("mt-"):
        return "Mitochondrial"
    b = bt.get(gene_id, set())
    if "protein_coding" in b:
        return "Protein coding"
    if "lncRNA" in b:
        return "lncRNA"
    if b & {"rRNA", "rRNA_pseudogene", "Mt_rRNA"}:
        return "rRNA"
    return "Other"


CATS = ["Protein coding", "lncRNA", "Mitochondrial", "rRNA", "Other"]
cat = pd.Series([category(g) for g in counts.index], index=counts.index)
comp = counts.groupby(cat).sum().reindex(CATS).fillna(0)
frac = 100 * comp / comp.sum()
frac.T.rename(index=label).round(3).to_csv(os.path.join(results_dir, "count_composition.tsv"), sep="\t")

# <figure name="count_composition">
fig, ax = plt.subplots(figsize=(8, 4.6))
left = np.zeros(len(order))
y = np.arange(len(order))[::-1]
for c, col in zip(CATS, CAT_COLS[:4] + [GREY]):
    v = frac.loc[c, order].values
    ax.barh(y, v, left=left, color=col, height=0.7, label=c, edgecolor="white", linewidth=1)
    left += v
ax.set_yticks(y, [label[r] for r in order])
ax.set_xlim(0, 100)
ax.set_xlabel("Share of assigned reads (%)")
ax.grid(axis="y", visible=False)
ax.legend(loc="upper center", bbox_to_anchor=(0.5, 1.12), ncol=5, fontsize=9)
save(fig, "count_composition.png")
# </figure>

# <figure name="detected_genes">
# Detected genes
detected = (counts >= 10).sum()
fig, ax = plt.subplots(figsize=(5, 3.8))
for i, g in enumerate(GROUP_ORDER):
    runs = [r for r in order if group[r] == g]
    jitter = np.linspace(-0.12, 0.12, len(runs))
    ax.scatter(i + jitter, detected[runs], s=46, color=GROUP_COL[g], edgecolor="white", linewidth=1.5, zorder=3)
ax.set_xticks(range(len(GROUP_ORDER)), GROUP_ORDER)
ax.set_xlim(-0.6, len(GROUP_ORDER) - 0.4)
ax.set_ylabel("Genes with ≥ 10 counts")
ax.grid(axis="x", visible=False)
save(fig, "detected_genes.png")
# </figure>

# <figure name="top_expressed_genes">
# Most abundant genes
share = counts.div(counts.sum()) * 100
top = share.mean(axis=1).sort_values(ascending=False).head(20)
fig, ax = plt.subplots(figsize=(7, 5.2))
yy = np.arange(len(top))[::-1]
ax.barh(yy, top.values, color="#2a78d6", height=0.65)
ax.set_yticks(yy, [names[g] for g in top.index], fontsize=9)
ax.set_xlabel("Mean share of assigned reads (%)")
ax.grid(axis="y", visible=False)
save(fig, "top_expressed_genes.png")
# </figure>
pd.DataFrame({"gene_id": top.index, "gene_name": [names[g] for g in top.index],
              "mean_percent_of_reads": top.values.round(3)}).to_csv(
    os.path.join(results_dir, "top_expressed_genes.tsv"), sep="\t", index=False)

# log-CPM of genes with >= 10 counts in every sample
keep = (counts >= 10).all(axis=1)
lcpm = np.log2(counts[keep].div(counts.sum()) * 1e6 + 1)
corr = np.corrcoef(lcpm.T.values)
pd.DataFrame(corr, index=[label[r] for r in order], columns=[label[r] for r in order]).round(4).to_csv(
    os.path.join(results_dir, "sample_correlation_logcpm.tsv"), sep="\t")

# Sequential blue ramp (single hue, light -> dark) for a magnitude
from matplotlib.colors import LinearSegmentedColormap

# <figure name="sample_correlation">
blues = LinearSegmentedColormap.from_list("blues", ["#cde2fb", "#86b6ef", "#3987e5", "#1c5cab", "#0d366b"])
fig, ax = plt.subplots(figsize=(7.6, 6.4))
im = ax.imshow(corr, cmap=blues, vmin=np.floor(corr.min() * 100) / 100, vmax=1)
ax.set_xticks(range(len(order)), [label[r] for r in order], rotation=90, fontsize=8.5)
ax.set_yticks(range(len(order)), [label[r] for r in order], fontsize=8.5)
ax.grid(False)
for i in range(len(order)):
    for j in range(len(order)):
        v = corr[i, j]
        ax.text(j, i, f"{v:.3f}", ha="center", va="center", fontsize=6.5,
                color="white" if (v - im.norm.vmin) / (1 - im.norm.vmin) > 0.4 else INK)
cb = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.03)
cb.set_label("Pearson r (log2 CPM)")
cb.outline.set_visible(False)
save(fig, "sample_correlation.png")
# </figure>

# <figure name="pca_logcpm">
# PCA on the 1,000 most variable genes
var = lcpm.var(axis=1).sort_values(ascending=False)
X = lcpm.loc[var.index[:1000]].T.values
X = X - X.mean(axis=0)
U, S, Vt = np.linalg.svd(X, full_matrices=False)
pcs = U * S
pve = 100 * S ** 2 / (S ** 2).sum()
fig, ax = plt.subplots(figsize=(6.2, 4.8))
for i, r in enumerate(order):
    ax.scatter(pcs[i, 0], pcs[i, 1], s=60, color=GROUP_COL[group[r]], edgecolor="white", linewidth=2, zorder=3)
    ax.annotate(label[r].split(" ·")[0], (pcs[i, 0], pcs[i, 1]), xytext=(5, 4), textcoords="offset points",
                fontsize=8, color=INK2)
ax.set_xlabel(f"PC1 ({pve[0]:.1f}%)")
ax.set_ylabel(f"PC2 ({pve[1]:.1f}%)")
group_legend(ax, loc="best")
save(fig, "pca_logcpm.png")
# </figure>
pd.DataFrame({"run": order, "label": [label[r] for r in order], "PC1": pcs[:, 0].round(3),
              "PC2": pcs[:, 1].round(3)}).to_csv(os.path.join(results_dir, "pca_logcpm.tsv"), sep="\t", index=False)

print(f"QC figures written to {plots_dir}")
