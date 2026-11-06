# Doxorubicin-induced transcriptional changes: bulk RNA-seq analysis (GSE178812)

A reproducible R pipeline for bulk RNA-seq differential expression and functional enrichment analysis of the public dataset **GSE178812** (NCBI GEO), investigating transcriptional alterations induced by doxorubicin, with a focus on mechanisms relevant to chemotherapy-induced cognitive impairment ("chemobrain").

## Overview

The pipeline downloads and prepares the data, performs quality control and differential expression analysis with DESeq2, and runs functional enrichment analyses (GO, KEGG and GSEA), including pathway-level visualization.

## Repository structure

```
.
├── Makefile                          # Runs the whole pipeline
├── data/
│   ├── colData.csv                   # Sample metadata
│   ├── counts.csv                    # Raw count matrix (gene symbols)
│   ├── dds.rds                       # DESeq2 object (generated, not versioned)
│   └── GSE178812_Doxo-counts.xlsx    # Original count table from GEO
├── results/
│   ├── top_table.csv                 # DESeq2 results for all tested genes
│   ├── GSEA_GO_all_full.csv          # GSEA (GO: BP, CC, MF)
│   ├── GSEA_KEGG_full.csv            # GSEA (KEGG)
│   ├── enrichment_STRING_full.csv    # STRING functional enrichment
│   ├── edges_PPI.csv                 # PPI edges (Cytoscape-ready)
│   └── nodes_PPI.csv                 # PPI node attributes (log2FC, padj)
├── plots/                            # QC, DE, enrichment and network figures
├── scripts/
│   ├── 00_packages.R
│   ├── 01_download_and_prepare.R
│   ├── 02_differential_expression_and_qc.R
│   └── 03_enrichment.R
├── cache/                            # STRING downloads (generated, not versioned)
├── sessionInfo.txt                   # R and package versions used
└── GSE178812_series_matrix.txt.gz    # Series matrix from GEO
```
## Pipeline

Each script reads its inputs from disk and writes its outputs to disk, so steps can be run independently or chained with `make`.

| Script | Description | Main outputs |
|---|---|---|
| `00_packages.R` | Installs, loads and configures packages and project folders (sourced by all steps) | `sessionInfo.txt` |
| `01_download_and_prepare.R` | Retrieves data from GEO and prepares the count matrix and sample metadata | `data/counts.csv`, `data/colData.csv` |
| `02_differential_expression_and_qc.R` | Quality control and differential expression with DESeq2 | `data/dds.rds`, `results/top_table.csv`, QC and DE plots |
| `03_enrichment.R` | GSEA (GO, KEGG), pathview, STRING PPI network and STRING enrichment | `results/GSEA_*.csv`, `results/*_PPI.csv`, `results/enrichment_STRING_full.csv`, enrichment plots |

## Methods summary

- **Filtering:** genes with at least 10 counts in at least 3 samples are kept.
- **Differential expression:** DESeq2 with design `~ group` (`VEH` as reference) and `fitType = "local"`. Contrast: **CHEM vs VEH**, so a positive log2 fold change means higher expression in CHEM. Differentially expressed genes (DEGs) are defined as `padj < 0.05` and `|log2FC| > 0.5`.
- **GSEA:** genes are ranked by the DESeq2 Wald statistic. GO (BP, CC and MF) and KEGG (`mmu`) are tested with `clusterProfiler`, with Benjamini-Hochberg adjustment. A positive NES indicates enrichment in CHEM and a negative NES in VEH.
- **STRING:** DEGs are mapped to STRING v12.0 (mouse, taxon 10090, combined score >= 700). Functional enrichment uses all tested genes as the background. Literature (`PMID`) and network-neighbor categories are excluded from the significant-term tables.
- **Reproducibility:** a fixed random seed is set for GSEA permutations, and the R/package versions are recorded in `sessionInfo.txt`.

## Outputs

- **Quality control:** PCA (VST), normalized count density and boxplot, dispersion plot, Cook's distance
<img width="500" alt="boxplot" src="https://github.com/user-attachments/assets/bb0a4e49-9cbf-48d2-8e07-2a539bd6287c" />

- **Differential expression:** volcano plot, MA plot, heatmap, FDR histogram, QQ plot
<img width="500" alt="volcano" src="https://github.com/user-attachments/assets/52692b2d-18ab-4610-b139-f3e4234ab2da" />

- **Functional enrichment:** GO and KEGG dotplots, enrichment maps, gene-concept networks, ridgeplots, GSEA plots
<img width="500" alt="dotplot_GO" src="https://github.com/user-attachments/assets/bfb27a9e-41e2-497a-aab0-d221e8aad49b" />

- **Network analysis:** STRING PPI network of DEGs, STRING enrichment (Reactome, WikiPathways, mammalian phenotype), and edge/node tables for Cytoscape

- **Pathway visualization:** pathview maps (e.g. mmu04668, TNF signaling pathway) and GSEA running-score plots
<img width="500" alt="gsea_kegg" src="https://github.com/user-attachments/assets/c6e4201f-6ff2-4797-b672-87173a0fd0cc" />

## Requirements

- R (>= 4.0)
- GNU Make (>= 4.3)
- Internet connection (GEO, STRING and KEGG are queried at run time)
- Main packages: `DESeq2`, `clusterProfiler`, `enrichplot`, `pathview`, `STRINGdb`, `ReactomePA`, `ggplot2` (missing packages are installed automatically by `00_packages.R`; the full list of versions is in [`sessionInfo.txt`](sessionInfo.txt))

## Usage

```bash
git clone git@github.com:arlindomatias/bulk-rna-seq.git
cd bulk-rna-seq
make
```

Run all commands from the repository root.
Figures are saved to `plots/` and result tables to `results/`.

## Data

The analysis uses the public dataset [GSE178812](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE178812) from NCBI GEO. Please cite the original authors when using these data.

## Author

Arlindo César Matias Pereira
