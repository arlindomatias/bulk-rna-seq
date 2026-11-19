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
| `results/gene_map.csv` | Every gene tested in every dataset with its cross-dataset `key` (mouse Ensembl ID), `key_symbol` and `map_status` |
| `results/all_degs.csv` | One row per gene and dataset where it is a DEG: `key`, `key_symbol`, `SYMBOL`, `ENTREZID`, `ENSEMBL`, `dataset`, `species`, `contrast`, `method`, `direction`, `log2FoldChange`, `stat`, `stat_type`, `pvalue`, `padj`, `deg_criterion`, `n_datasets`, `datasets` |
| `results/meta_random_effects.csv` | Random-effects meta-analysis of every gene tested in ≥ 2 datasets: pooled log2FC, 95% CI, z, p, BH padj, `meta_deg`, Knapp-Hartung p, τ², I², Q, leave-one-out, per-dataset log2FC/SE, DEG status in each project |
| `results/pairwise_overlap.csv` | Shared DEGs of each pair of datasets (universe = keys tested in both), hypergeometric p, same/opposite direction |
| `plots/deg_counts.png`, `plots/shared_degs_lfc.png` | DEGs per dataset; log2FC of the genes that are DEGs in more than one dataset |
| `plots/meta_volcano.png`, `plots/meta_forest.png` | Meta-analysis volcano; forest plot of the top meta DEGs |
| `reports/report.html`, `reports/sessionInfo*.txt` | Report and R sessions |

Identifiers missing from a project's table come from `org.Mm.eg.db` / `org.Rn.eg.db`; several IDs for one
gene are joined with `;`.

**Meta-analysis.** For every gene (key) tested in at least `meta.min_datasets` datasets, the
`log2FoldChange` and `lfcSE` of each dataset are combined by `metafor::rma` (random effects, τ² by REML, z
test, BH over genes). A meta DEG has padj < `meta.padj` and |pooled log2FC| > `meta.lfc`. Robustness
columns: Knapp-Hartung p (very conservative with 2–3 studies) and leave-one-out (largest p and sign flip when
each dataset is dropped; with 2 datasets, the p of the weaker dataset alone). The log2FC scales differ
slightly (DESeq2 on counts, limma on log2(normalized + 1)) and the treatments differ (AC-T vs doxorubicin).

**Gene matching.** Genes are matched across datasets and species by `key`, the Ensembl ID of the mouse
gene. Mouse genes use their own ID (GSE178812 symbols are mapped through `org.Mm.eg.db` when a symbol has
exactly one Ensembl ID); rat genes use their 1:1 mouse ortholog (`ortholog_one2one`) from Ensembl Compara
release 104, the release of the GSE236404 rat IDs (rn6). Genes without a unique key (no 1:1 ortholog,
ambiguous symbol, or a key shared by two genes of one dataset) are not matched. The `orthologs` rule
downloads the Compara tables once into `cache/compara/` (internet needed on the first run; delete the folder
to refresh). Ensembl BioMart is not used: it did not answer in November 2026.
`stat` is the DESeq2 Wald statistic or the limma moderated t (`stat_type`), so it is not on the same scale
across methods.

## Structure

```
.
├── Snakefile          # orthologs (download) -> gene_map (R) -> collect_degs (R) -> meta_random_effects (R) -> plots, report (Python)
├── config.yaml        # datasets (project folder, species, OrgDb, gene ID type, method); orthologs (Compara); meta
├── cache/compara/     # Ensembl Compara homology tables (downloaded, git-ignored)
├── envs/              # limma (R + org.Mm.eg.db + org.Rn.eg.db), metafor (R + metafor), both pinned; report (Python)
├── results/, plots/, reports/
```

Scripts: [`../scripts/meta/`](../scripts/meta) (`gene_map.R`, `collect_degs.R`, `meta_random_effects.R`, `plots.py`, `report.py`).

## Running

```bash
conda activate genomic-variants   # Snakemake 9
snakemake -n --use-conda
snakemake --use-conda --cores 1
```
