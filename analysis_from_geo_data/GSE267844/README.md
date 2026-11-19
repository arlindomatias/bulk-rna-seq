# Docetaxel-induced transcriptional changes in the rat hippocampus (GSE267844)

Differential expression and functional enrichment of the public dataset **GSE267844** (NCBI GEO,
BioProject PRJNA1113364): hippocampus of rats treated with docetaxel (**DTX**, n = 3) or saline (**NC**,
n = 3). The series also has a docetaxel + edaravone dexborneol group (DTX-ED), which is not used here.

**Report:** [arlindomatias.github.io/bulk-rna-seq/analysis_from_geo_data/GSE267844/reports/report.html](https://arlindomatias.github.io/bulk-rna-seq/analysis_from_geo_data/GSE267844/reports/report.html)

## Why limma instead of DESeq2

GEO only has an **FPKM** table for this series (`GSE267844_2024-5-17_genes_fpkm_expression.xlsx`); the raw
reads are only in SRA. DESeq2 and voom need raw counts, so the DE step is **limma-trend** on log2(FPKM + 1)
(`lmFit` + `eBayes(trend = TRUE, robust = TRUE)`), as in [`../GSE236404`](../GSE236404). Enrichment, STRING
and the report are the shared scripts of the other DE projects. Fold changes and p-values are not numerically
comparable to the DESeq2 projects, although the cutoffs are the same.

## Repository structure

```
.
├── Snakefile                         # Snakemake workflow that runs the whole pipeline
├── config.yaml                       # Dataset settings (input, sample patterns, groups, cutoffs, filter, enrichment)
├── envs/
│   ├── limma.yaml                    # Conda environment `limma` (same as ../GSE236404)
│   ├── limma.linux-64.pin.txt        # Exact package versions (linux-64)
│   └── report.yaml                   # Python environment of the report rule
├── data/
│   ├── GSE267844_2024-5-17_genes_fpkm_expression.xlsx  # FPKM table from GEO
│   ├── expression.csv                # Gene-level FPKM of the analysed samples
│   ├── colData.csv                   # Sample metadata used in the model
│   └── se_raw.rds, fit.rds           # SummarizedExperiment and limma fit (generated, not versioned)
├── results/                          # top_table.csv, de_summary.tsv, log2 values, GSEA, STRING, PPI tables
├── plots/                            # QC, DE and enrichment figures (pathview maps in plots/pathview/)
├── reports/report.html               # HTML report of the whole analysis (rule `report`)
└── cache/                            # STRING and KEGG downloads (generated, not versioned)
```

## Pipeline

| Rule | Script (`../../scripts/de/`) | Main outputs |
|---|---|---|
| `prepare` | `01_prepare_geo_normalized.R` | Downloads the GEO table; samples from its `FPKM.*` columns (`geo.sample_patterns`); `data/expression.csv`, `data/se_raw.rds`, `data/colData.csv` |
| `differential` | `02_differential_expression_limma.R` | `data/fit.rds`, `results/top_table.csv`, `results/de_summary.tsv`, `results/log2_normalized_counts.csv`, QC/DE plots |
| `kegg_annotation` | (curl) | KEGG gene–pathway table (`rno`) in `cache/kegg/` (downloaded once) |
| `enrichment` | `03_enrichment.R` | `results/GSEA_*.csv`, `results/*_PPI.csv`, `results/enrichment_STRING_full.csv`, `reports/sessionInfo.txt`, enrichment plots |
| `report` | `report.py` | `reports/report.html` |

## Methods summary

- **Input:** 24,858 Ensembl rat genes (rn6-era IDs) × 9 samples; the 3 DTX-ED samples are dropped.
- **Filtering:** genes with FPKM ≥ 1 in at least 3 samples (the group size).
- **Differential expression:** limma-trend, design `~ group` (`NC` as reference). Contrast **DTX vs NC**:
  a positive log2 fold change means higher expression in DTX. DEGs: `padj < 0.05` and `|log2FC| > 0.5`;
  when no gene has `padj < 0.05`, the unadjusted `p < 0.01` is used instead (exploratory).
- **GSEA:** genes ranked by the limma moderated t statistic; GO (`org.Rn.eg.db`) and KEGG (`rno`, cached).
- **STRING:** v12.0, *Rattus norvegicus* (taxon 10116), combined score ≥ 700, all tested genes as background.

## Results

DTX vs NC, 13,869 genes tested:

- **No gene-level signal.** The PCA does not separate the groups and the smallest padj is 0.94. The DEGs
  therefore use the exploratory **unadjusted `p < 0.01`**: 65 genes (31 up, 34 down). Fewer genes fall below
  p = 0.01 (86) than expected by chance (~139), so these DEGs are essentially noise; the genes highlighted
  by the original study (*Dbh*, *Fgb*, *Muc6*, *Tnfrsf9*) are not reproduced.
- **GSEA** is more consistent: 689 GO terms (71 up / 618 down in DTX) and 46 KEGG pathways (6 / 40) at
  `padj < 0.05`. Lower in DTX: integrin signaling, focal adhesion, ECM–receptor interaction, cell
  motility/migration; higher: spliceosome / RNA splicing, proteasome, cell cycle.

## Running

```bash
conda activate genomic-variants      # Snakemake 9
cd bulk-rna-seq/analysis_from_geo_data/GSE267844
snakemake -n --use-conda             # dry-run
snakemake --use-conda --cores 2
```

The `enrichment` step needs internet access (STRING, pathview).

## Data and citation

The analysis uses the public dataset [GSE267844](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE267844)
from NCBI GEO. Please cite the original authors when using these data.
