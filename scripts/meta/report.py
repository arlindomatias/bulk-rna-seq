"""HTML report for meta_analysis: DEGs of every DE project of the repository in one table.

Run through Snakemake (`script:` directive, rule `report`). Uses the figures in plots/ and the
tables in results/; every number is computed from them.
"""
import os
import sys

import pandas as pd

# Shared report module: <repository>/scripts/html_report.py (this script lives in scripts/meta/).
sys.path.insert(0, os.path.join(snakemake.scriptdir, ".."))  # noqa: F821
from html_report import Report  # noqa: E402

inp, out, cfg = snakemake.input, snakemake.output, snakemake.config  # noqa: F821
ds = cfg["datasets"]

degs = pd.read_csv(inp.degs)
overlap = pd.read_csv(inp.overlap)
summaries = {name: pd.read_csv(path, sep="\t").iloc[0] for name, path in zip(ds, inp.summaries)}

gene_map = pd.read_csv(inp.gene_map, dtype=str)
ortho = cfg["orthologs"]
degs["label"] = degs["key_symbol"].fillna(degs["SYMBOL"])
n_genes = degs["key"].nunique() + int(degs["key"].isna().sum())
shared = degs[degs.n_datasets > 1]
n_shared = shared["key"].nunique()
max_n = int(degs["n_datasets"].max())

r = Report("Chemotherapy DEGs across datasets",
           subtitle="Meta-analysis · DEG table of every DE project in the repository",
           code_files=list(inp.code))

r.kpis([
    (len(ds), "datasets"),
    (f"{len(degs):,}", "rows (gene × dataset)"),
    (f"{n_genes:,}", "distinct genes"),
    (n_shared, "genes that are DEGs in > 1 dataset"),
    (max_n, "most datasets for one gene"),
])

# ---------------------------------------------------------------- datasets
r.h2("Datasets")
rows = []
links = []
for name, d in ds.items():
    s = summaries[name]
    crit = f"{s['deg_criterion']} < {s['deg_cutoff']}, |log2FC| > {s['lfc_cutoff']}"
    rows.append({"dataset": name, "design": d["description"], "species": d["species"],
                 "method": d["method"], "contrast": s["contrast"], "samples": int(s["samples"]),
                 "genes tested": int(s["genes_tested"]), "DEG criterion": crit,
                 "up": int(s["deg_up"]), "down": int(s["deg_down"])})
    links.append(f"* **{name}**: [project report]({os.path.join('..', d['report'])})")
fallback = [n for n in ds if summaries[n]["deg_criterion"] == "pvalue"]
r.md(f"""
This project does not fit any model. It collects the differentially expressed genes (DEGs) that each DE
project of the repository already reports (`results/top_table.csv` and `results/de_summary.tsv`), using the
same DEG rule as that project, and puts them in one table: `results/all_degs.csv`.

{chr(10).join(links)}
""")
r.table(pd.DataFrame(rows), "Datasets and DEG criteria (from each project's `de_summary.tsv`).")
if fallback:
    r.md(f"""
> **Exploratory DEGs:** {", ".join(fallback)} has no gene at padj < {summaries[fallback[0]]['padj_cutoff']}, so its DEGs
> use the unadjusted p-value cutoff. Many of them are expected by chance; overlaps with these datasets are weak evidence.
""")

# ---------------------------------------------------------------- table
r.h2("DEG table")
n_no_entrez = int(degs["ENTREZID"].isna().sum())
n_no_ens = int(degs["ENSEMBL"].isna().sum())
r.md(f"""
`results/all_degs.csv` has one row per gene and dataset where the gene is a DEG: a gene that is a DEG in
two datasets has two rows. Columns:

* `SYMBOL`, `ENTREZID`, `ENSEMBL`: the project's own identifier (Ensembl ID or symbol) plus the others from
  `org.Mm.eg.db` / `org.Rn.eg.db`. Several IDs for one gene are joined with `;`.
  {n_no_entrez} rows have no Entrez ID and {n_no_ens} have no Ensembl ID.
* `dataset`, `species`, `contrast`, `method`, `deg_criterion`.
* `direction`, `log2FoldChange` (positive = higher in the treated group), `stat`, `stat_type`
  (DESeq2 Wald statistic or limma moderated t), `pvalue`, `padj`.
* `key`, `key_symbol`: the cross-dataset gene key (mouse Ensembl ID and its symbol; see below).
* `n_datasets`, `datasets`: in how many and which datasets the gene is a DEG, matched by `key`.

Effect sizes from DESeq2 and limma are not on the same scale, and the treatments differ
(AC-T in GSE295281, doxorubicin alone in the others).
""")
r.figure(inp.deg_counts, "DEGs per dataset, by direction.", width="60%")

# ---------------------------------------------------------------- gene matching
r.h2("Matching genes across datasets and species")
status_names = {
    "ensembl": "own Ensembl ID", "symbol": "symbol → one Ensembl ID",
    "alias": "outdated symbol (alias) → one Ensembl ID",
    "ortholog_1to1": "1:1 mouse ortholog", "ambiguous_symbol": "symbol → several Ensembl IDs",
    "no_ensembl": "symbol without Ensembl ID", "not_1to1": "only 1:n or n:m orthologs",
    "no_ortholog": "no mouse ortholog", "duplicate_key": "key shared by several genes",
}
st = (gene_map.groupby(["dataset", "map_status"]).size().unstack(fill_value=0)
      .reindex(index=list(ds), fill_value=0))
st = st[[c for c in status_names if c in st.columns]]
st.insert(0, "genes tested", st.sum(axis=1))
st.insert(1, "with key", gene_map.dropna(subset=["key"]).groupby("dataset").size().reindex(list(ds)))
st = st.rename(columns=status_names).reset_index()
keys = {d: set(gene_map.loc[(gene_map.dataset == d) & gene_map.key.notna(), "key"]) for d in ds}
n_all = len(set.intersection(*keys.values()))
deg_nokey = int(degs["key"].isna().sum())
r.md(f"""
Every gene tested in each dataset (`results/gene_map.csv`) gets a `key`: the Ensembl gene ID of the
mouse gene. Mouse genes use their own Ensembl ID. GSE178812 has symbols, some of them outdated (*Ctgf*
for *Ccn2*), so they are mapped through `org.Mm.eg.db`: first as current symbols, then as aliases, and only
when they give exactly one Ensembl ID (an alias whose gene is already in the dataset under its current
symbol is dropped). Rat genes use their **1:1 mouse ortholog** from Ensembl Compara
release {ortho['release']} (`ortholog_one2one`), the release of the GSE236404 rat IDs (Ensembl 104, rn6).
Genes without a unique key are not matched across datasets. This replaces matching by upper-case
symbol, which pairs genes that merely share a name and misses orthologs named differently.

{n_all:,} keys are tested in all {len(ds)} datasets. {deg_nokey} of the {len(degs):,} DEG rows have no key.
""")
r.table(st, "Genes tested per dataset and how they got their key.", digits=0)

# ---------------------------------------------------------------- overlap
r.h2("Overlap between datasets")
r.md("""
For each pair of datasets, the universe is the genes (keys) tested in both. The p-value is the
hypergeometric probability of at least that many shared DEGs if the two DEG lists were independent.
The DEG counts in this table only include DEGs whose key was tested in both datasets, so they are lower
than the totals above.
""")
show = overlap.drop(columns="shared_genes").rename(columns=lambda c: c.replace("_", " "))
r.table(show, "Pairwise DEG overlap (`results/pairwise_overlap.csv`).", digits=3)

conc = shared.groupby("key")["direction"].nunique()
same, opp = int((conc == 1).sum()), int((conc > 1).sum())
in_all = sorted(shared.loc[shared.n_datasets == len(ds), "label"].unique())
in_all_text = "none" if not in_all else f"{len(in_all)} ({', '.join(f'*{g}*' for g in in_all)})"
r.md(f"""
{n_shared} genes are DEGs in more than one dataset; {in_all_text} in all {len(ds)}.
{same} change in the same direction in every dataset where they are DEGs and {opp} change in opposite directions.
""")
r.figure(inp.shared_lfc, "log2FC of the genes that are DEGs in more than one dataset (blank = not a DEG there).",
         width="55%")
cols = ["label", "SYMBOL", "dataset", "direction", "log2FoldChange", "stat", "pvalue", "padj"]
r.table(shared[cols], "Genes that are DEGs in more than one dataset.", digits=3, collapsed=True)

# ---------------------------------------------------------------- appendix
r.appendix(
    outputs=[
        ("results/gene_map.csv", "Every tested gene of every dataset with its cross-dataset key"),
        ("results/all_degs.csv", "DEGs of every dataset, one row per gene × dataset"),
        ("results/pairwise_overlap.csv", "Pairwise DEG overlap and hypergeometric test"),
        ("plots/deg_counts.png", "DEGs per dataset"),
        ("plots/shared_degs_lfc.png", "log2FC of the shared DEGs"),
        ("reports/sessionInfo.txt", "R session of the table step"),
    ],
    env_files=list(inp.envs),
    config_files=[inp.config],
    session_file=inp.session,
)
r.write(out.html)
