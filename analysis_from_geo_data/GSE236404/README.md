# Doxorubicin-induced transcriptional changes in the rat dentate gyrus (GSE236404)

Differential expression and functional enrichment of the public dataset **GSE236404** (NCBI GEO,
PMID 38225896): dentate gyrus of the hippocampus of young (6-month) male F344 rats treated with
doxorubicin (**DOX**, n = 8) or vehicle (**VEH**, n = 10). The series also has a doxorubicin + ABT-263
(senolytic) group, which is not used here.

**Report:** [arlindomatias.github.io/bulk-rna-seq/analysis_from_geo_data/GSE236404/reports/report.html](https://arlindomatias.github.io/bulk-rna-seq/analysis_from_geo_data/GSE236404/reports/report.html)

## Why limma instead of DESeq2

GEO only has a **normalized** (non-integer) count table for this series; the raw reads are only in SRA
(PRJNA990762, Ion Torrent). DESeq2 and voom need raw counts and library sizes, so the DE step is
**limma-trend** on log2(normalized counts + 1) (`lmFit` + `eBayes(trend = TRUE, robust = TRUE)`). Enrichment,
STRING and the report are the same shared scripts used by [`../GSE178812`](../GSE178812) and
[`../../bulk_rna_seq_tximport`](../../bulk_rna_seq_tximport). Fold changes and p-values are not numerically
comparable to the DESeq2 projects, although the cutoffs are the same.

## Repository structure

```
.
├── Snakefile                         # Snakemake workflow that runs the whole pipeline
├── config.yaml                       # Dataset settings (input, groups, cutoffs, filter, enrichment, de.method)
├── envs/
│   ├── limma.yaml                    # Conda environment `limma` (rna-seq stack + org.Rn.eg.db)
│   ├── limma.linux-64.pin.txt        # Exact package versions (linux-64)
│   └── report.yaml                   # Python environment of the report rule
├── data/
│   ├── GSE236404_DOX_Normalized_Count-6-23-23.csv.gz      # Normalized table from GEO
│   ├── GSE236404_DG-DOX-Sample-Attributes-12-16-22.csv.gz # Sample sheet from GEO
│   ├── expression.csv                # Gene-level normalized counts (rows sharing a gene ID summed)
│   ├── colData.csv                   # Sample metadata used in the model
│   └── se_raw.rds, fit.rds           # SummarizedExperiment and limma fit (generated, not versioned)
├── results/
│   ├── top_table.csv                 # limma results for all tested genes
│   ├── de_summary.tsv                # Number of DEGs and cutoffs
│   ├── log2_normalized_counts.csv    # log2(normalized + 1) of the tested genes
│   ├── GSEA_GO_all_full.csv          # GSEA (GO: BP, CC, MF)
│   ├── GSEA_KEGG_full.csv            # GSEA (KEGG, rno)
│   ├── enrichment_STRING_full.csv    # STRING functional enrichment (empty when there are < 2 DEGs)
│   ├── edges_PPI.csv, nodes_PPI.csv  # PPI tables (empty when there are < 2 DEGs)
├── plots/                            # QC, DE and enrichment figures (pathview maps in plots/pathview/)
├── reports/report.html               # HTML report of the whole analysis (rule `report`)
└── cache/                            # STRING and KEGG downloads (generated, not versioned)
```

## Pipeline

| Rule | Script (`../../scripts/de/`) | Main outputs |
|---|---|---|
| `prepare` | `01_prepare_geo_normalized.R` | Downloads the GEO table and sample sheet; `data/expression.csv`, `data/se_raw.rds`, `data/colData.csv` |
| `differential` | `02_differential_expression_limma.R` | `data/fit.rds`, `results/top_table.csv`, `results/de_summary.tsv`, `results/log2_normalized_counts.csv`, QC/DE plots |
| `kegg_annotation` | (curl) | KEGG gene–pathway table (`rno`) in `cache/kegg/` (downloaded once) |
| `enrichment` | `03_enrichment.R` | `results/GSEA_*.csv`, `results/*_PPI.csv`, `results/enrichment_STRING_full.csv`, `reports/sessionInfo.txt`, enrichment plots |
| `report` | `report.py` | `reports/report.html` |

`00_common.R` is sourced by every step (logging, config, sample metadata). Logs go to `logs/{rule}.log`.

## Methods summary

- **Input:** the deposited table has 24,012 rows for 17,910 Ensembl gene IDs (several features for some
  genes); rows that share a gene ID are summed. Values are on a common normalized scale (14.8–20.3 million
  per sample).
- **Filtering:** genes with at least 10 normalized counts in at least 8 samples (the smallest group).
- **Differential expression:** limma-trend, design `~ group` (`VEH` as reference). Contrast **DOX vs VEH**:
  a positive log2 fold change means higher expression in DOX. DEGs: `padj < 0.05` and `|log2FC| > 0.5`;
  when no gene has `padj < 0.05`, the unadjusted `p < 0.01` (`de.pvalue_fallback`) is used instead (exploratory).
- **GSEA:** genes ranked by the limma moderated t statistic. GO (BP, CC, MF; Ensembl IDs, `org.Rn.eg.db`)
  and KEGG (`rno`, cached annotation) with `clusterProfiler`, Benjamini-Hochberg adjustment and a fixed seed.
- **STRING:** v12.0, *Rattus norvegicus* (taxon 10116), combined score ≥ 700, all tested genes as background.
  Only run when at least 2 DEGs map to STRING (here: the exploratory p < 0.01 DEGs).
- **Reproducibility:** fixed seed (`seed` in `config.yaml`); R/package versions pinned in
  `envs/limma.linux-64.pin.txt`.

## Results

DOX vs VEH, 14,569 genes tested:

- **No gene reaches `padj < 0.05`** (smallest padj 0.16): within-group variability is high and the PCA does
  not separate the groups. Following the repository rule, the DEGs therefore use the **unadjusted
  `p < 0.01`** (and `|log2FC| > 0.5`): **169 exploratory DEGs** (49 up, 120 down in DOX). About 146 genes are
  expected below p = 0.01 by chance (294 observed without the fold-change cutoff), so many are false
  positives. Top genes: *Slc35d3*, *Hhatl*, *Fcrl2*, *Ttc23*, *Dcst1*, *Dbp*, *Gfap*, *Fn1*, *Col3a1*.
- **STRING:** 163 of the 169 DEGs mapped, 42 with at least one interaction (84 edges); enriched terms are
  mostly secreted / extracellular-region proteins and response to lipid.
- **GSEA** (independent of the DEG cutoff): 293 GO terms (105 enriched in DOX, 188 in VEH) and 27 KEGG
  pathways (4 / 23) with `padj < 0.05`. Ribosome, oxidative phosphorylation and fatty-acid metabolism are lower
  in DOX; mRNA processing / splicing (spliceosome) are higher. Pathview maps: `rno03010` (ribosome),
  `rno05171` (coronavirus disease), `rno00190` (oxidative phosphorylation).
- **Sample weights:** limma `arrayWeights` (`de.array_weights`) was tested and is **off**. The weights
  (0.38–3.47) followed PC1, a hidden factor that does not match the groups (D5, D7 cluster with VEH; V10, V11
  with DOX), and gave 2,974 genes with padj < 0.05 with half of all genes at p < 0.1, i.e. inflated p-values.

## Running

```bash
conda activate genomic-variants      # Snakemake 9
cd bulk-rna-seq/analysis_from_geo_data/GSE236404
snakemake -n --use-conda             # dry-run
snakemake --use-conda --cores 2
```

The `enrichment` step needs internet access (STRING, pathview); the first run also downloads the GEO files and
the KEGG annotation.

## Data and citation

The analysis uses the public dataset [GSE236404](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE236404)
from NCBI GEO (PMID 38225896). Please cite the original authors when using these data.
