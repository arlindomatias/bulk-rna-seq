# Doxorubicin-induced transcriptional changes: bulk RNA-seq analysis (GSE178812)

A reproducible R pipeline for bulk RNA-seq differential expression and functional enrichment analysis of the public dataset **GSE178812** (NCBI GEO), investigating transcriptional alterations induced by doxorubicin, with a focus on mechanisms relevant to chemotherapy-induced cognitive impairment ("chemobrain").

## Overview

The pipeline downloads and prepares the data, performs quality control and differential expression analysis with DESeq2, and runs functional enrichment analyses (GO, KEGG and GSEA), including pathway-level visualization.

## Repository structure

```
.
├── data/
│   ├── colData.csv                  # Sample metadata
│   ├── counts.csv                   # Raw count matrix
│   ├── keytable.csv                 # Gene identifier mapping
│   ├── top_table.csv                # Differential expression results
│   └── GSE178812_Doxo-counts.xlsx   # Original count table from GEO
├── plots/                           # QC, DE and enrichment figures
├── scripts/
│   ├── 00_packages.R
│   ├── 01_download_and_prepare.R
│   ├── 02_differential_expression_and_qc.R
│   ├── 03_enrichment.R
│   └── 04_run_pipeline.R
└── GSE178812_series_matrix.txt.gz   # Series matrix from GEO
```

## Pipeline

| Script | Description |
|---|---|
| `00_packages.R` | Installs and loads the required packages |
| `01_download_and_prepare.R` | Retrieves data from GEO and prepares the count matrix and sample metadata |
| `02_differential_expression_and_qc.R` | Quality control (PCA, density, boxplot, dispersion, Cook's distance) and differential expression with DESeq2 |
| `03_enrichment.R` | GO and KEGG over-representation analysis, GSEA, and pathway visualization with pathview |
| `04_run_pipeline.R` | Runs all steps in sequence |

## Outputs

- **Quality control:** PCA (VST), count density, boxplot, dispersion plot, Cook's distance
<img width="500" alt="boxplot" src="https://github.com/user-attachments/assets/bb0a4e49-9cbf-48d2-8e07-2a539bd6287c" />


- **Differential expression:** volcano plot, MA plot, heatmap, FDR histogram, QQ plot
<img width="500" alt="volcano" src="https://github.com/user-attachments/assets/52692b2d-18ab-4610-b139-f3e4234ab2da" />

  
- **Functional enrichment:** GO and KEGG dotplots, enrichment maps, gene-concept networks, ridgeplots, GSEA plots
<img width="500" alt="dotplot_GO" src="https://github.com/user-attachments/assets/bfb27a9e-41e2-497a-aab0-d221e8aad49b" />

- **Pathway visualization:** pathview maps (e.g. mmu04668, TNF signaling pathway)
<img width="500" alt="gsea_kegg" src="https://github.com/user-attachments/assets/c6e4201f-6ff2-4797-b672-87173a0fd0cc" />


## Requirements

- R (>= 4.0)
- Main packages: `DESeq2`, `clusterProfiler`, `enrichplot`, `pathview`, `ggplot2` (the full list is in [`sessionInfo.txt`](sessionInfo.txt)

## Usage

```bash
git clone git@github.com:arlindomatias/bulk-rna-seq.git
cd bulk-rna-seq
Rscript scripts/04_run_pipeline.R
```

Figures are saved to `plots/` and result tables to `data/`.

## Data

The analysis uses the public dataset [GSE178812](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE178812) from NCBI GEO. Please cite the original authors when using these data.

## Author

Arlindo César Matias Pereira
