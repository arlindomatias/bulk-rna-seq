# Chemotherapy DEGs across datasets (meta-analysis)

One table with the differentially expressed genes (DEGs) of every DE project in this repository, and
their overlap between datasets.

**Report:** [arlindomatias.github.io/bulk-rna-seq/meta_analysis/reports/report.html](https://arlindomatias.github.io/bulk-rna-seq/meta_analysis/reports/report.html)

| Dataset | Project | Species | Method |
|---|---|---|---|
| GSE295281 | [`../bulk_rna_seq_tximport`](../bulk_rna_seq_tximport) | mouse | DESeq2 |
| GSE178812 | [`../analysis_from_geo_data/GSE178812`](../analysis_from_geo_data/GSE178812) | mouse | DESeq2 |
| GSE236404 | [`../analysis_from_geo_data/GSE236404`](../analysis_from_geo_data/GSE236404) | rat | limma-trend |

Nothing is refitted here: the workflow reads each project's `results/top_table.csv` and
`results/de_summary.tsv` and applies the same DEG rule as that project (`deg_criterion`, `deg_cutoff`,
`lfc_cutoff`; GSE236404 uses the exploratory unadjusted `p < 0.01` because no gene passes the FDR). The run
stops if the DEG counts differ from `de_summary.tsv`. Rerun the DE projects first when their results change.

## Outputs

| File | Contents |
|---|---|
| `results/all_degs.csv` | One row per gene and dataset where it is a DEG: `SYMBOL`, `ENTREZID`, `ENSEMBL`, `dataset`, `species`, `contrast`, `method`, `direction`, `log2FoldChange`, `stat`, `stat_type`, `pvalue`, `padj`, `deg_criterion`, `n_datasets`, `datasets` |
| `results/pairwise_overlap.csv` | Shared DEGs of each pair of datasets (universe = genes tested in both), hypergeometric p, same/opposite direction |
| `plots/deg_counts.png`, `plots/shared_degs_lfc.png` | DEGs per dataset; log2FC of the genes that are DEGs in more than one dataset |
| `reports/report.html`, `reports/sessionInfo.txt` | Report and R session |

Identifiers missing from a project's table come from `org.Mm.eg.db` / `org.Rn.eg.db`; several IDs for one
gene are joined with `;`. Genes are matched across datasets (and species) by upper-case symbol.
`stat` is the DESeq2 Wald statistic or the limma moderated t (`stat_type`), so it is not on the same scale
across methods.

## Structure

```
.
├── Snakefile          # collect_degs (R) -> plots (Python) -> report (Python)
├── config.yaml        # datasets: project folder, species, OrgDb, gene ID type, method
├── envs/              # limma (R + org.Mm.eg.db + org.Rn.eg.db, pinned), report (Python)
├── results/, plots/, reports/
```

Scripts: [`../scripts/meta/`](../scripts/meta) (`collect_degs.R`, `plots.py`, `report.py`).

## Running

```bash
conda activate genomic-variants   # Snakemake 9
snakemake -n --use-conda
snakemake --use-conda --cores 1
```
