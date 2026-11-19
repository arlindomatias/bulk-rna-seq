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

key = degs["SYMBOL"].str.upper()
n_genes = key.nunique()
shared = degs[degs.n_datasets > 1]
n_shared = shared["SYMBOL"].str.upper().nunique()
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
* `n_datasets`, `datasets`: in how many and which datasets the gene is a DEG. Genes are matched by
  upper-case symbol, so mouse and rat orthologs with the same symbol match.

Effect sizes from DESeq2 and limma are not on the same scale, and the treatments differ
(AC-T in GSE295281, doxorubicin alone in the others).
""")
r.figure(inp.deg_counts, "DEGs per dataset, by direction.", width="60%")

# ---------------------------------------------------------------- overlap
r.h2("Overlap between datasets")
r.md("""
For each pair of datasets, the universe is the genes tested in both (matched by symbol). The p-value is the
hypergeometric probability of at least that many shared DEGs if the two DEG lists were independent.
The DEG counts in this table only include DEGs whose symbol was tested in both datasets, so they are lower
than the totals above.
""")
show = overlap.drop(columns="shared_genes").rename(columns=lambda c: c.replace("_", " "))
r.table(show, "Pairwise DEG overlap (`results/pairwise_overlap.csv`).", digits=3)

conc = shared.groupby(shared["SYMBOL"].str.upper())["direction"].nunique()
same, opp = int((conc == 1).sum()), int((conc > 1).sum())
r.md(f"""
{n_shared} genes are DEGs in more than one dataset; {"none" if max_n < len(ds) else "some"} in all {len(ds)}.
{same} change in the same direction in every dataset where they are DEGs and {opp} change in opposite directions.
""")
r.figure(inp.shared_lfc, "log2FC of the genes that are DEGs in more than one dataset (blank = not a DEG there).",
         width="55%")
cols = ["SYMBOL", "dataset", "direction", "log2FoldChange", "stat", "pvalue", "padj"]
r.table(shared[cols], "Genes that are DEGs in more than one dataset.", digits=3, collapsed=True)

# ---------------------------------------------------------------- appendix
r.appendix(
    outputs=[
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
