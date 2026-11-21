# Cisplatin and the mouse hippocampus: single-cell pseudobulk (GSE216146)

Differential expression and functional enrichment of the public single-cell RNA-seq dataset **GSE216146**
(NCBI GEO; Park et al. 2024, PMID 38446899): hippocampi of 9-week-old female C57BL/6J mice injected with
cisplatin (**CIS**, chemobrain model) or PBS (**PBS**), collected on day 22 (10x Genomics, Cell Ranger 5.0.0,
mm10), n = 3 per group. The study tested 40 Hz sensory stimulation (GENUS, 1 h a day for 21 days) in both
groups; only the **unstimulated** groups are used here (cisplatin vs PBS).

**Report:** [arlindomatias.github.io/bulk-rna-seq/analysis_from_geo_data/GSE216146/reports/report.html](https://arlindomatias.github.io/bulk-rna-seq/analysis_from_geo_data/GSE216146/reports/report.html)

## Pseudobulk of the hippocampus

The counts of the cells of each sample are summed and analysed like bulk RNA-seq, so the result can be
compared with the bulk hippocampus datasets of this repository (`../../scripts/de/01_prepare_scrna.R`,
settings in `scrna:` of `config.yaml`):

- Input: the raw UMI counts deposited in GEO (cells × genes, gene symbols) and the authors' cell table
  (QC-filtered cells, no predicted doublets, sample, batch and cell type `MajorCellType`). The matrix has no
  barcodes; its row order is checked against the cell table (detected genes per cell vs `n_genes_by_counts`,
  r = 0.99999).
- **Choroid plexus and ependymal cells are left out.** Their share follows the dissection (choroid plexus
  2–38% of the cells of a sample) and a few choroid plexus genes (e.g. *Ttr*) would dominate the pseudobulk.
  `plots/cell_composition.png` shows the composition of every sample.
- Samples were sequenced in two batches (replicate 1 = D20, one mouse per sample; replicates 2–3 = D21, two
  mice pooled per sample), each with every group, so the model is `~ batch + group`.

## Repository structure

```
.
├── Snakefile                         # Snakemake workflow that runs the whole pipeline
├── config.yaml                       # Dataset settings (input files, scrna columns, groups, batch, cutoffs, enrichment)
├── envs/
│   ├── rna-seq.yaml                  # Conda environment `rna-seq` (same as ../GSE178812)
│   ├── rna-seq.linux-64.pin.txt      # Exact package versions (linux-64)
│   └── report.yaml                   # Python environment of the report rule
├── data/
│   ├── GSE216146_series_matrix.txt.gz  # GEO sample sheet (accessions, titles)
│   ├── pseudobulk_counts.csv         # Pseudobulk counts (gene symbol × sample)
│   ├── colData.csv                   # Sample metadata used in the model
│   └── dds_raw.rds, dds.rds          # DESeq2 objects (generated, not versioned)
├── results/                          # pseudobulk_summary.csv, cell_composition.csv, top_table.csv, de_summary.tsv, GSEA, STRING, PPI tables
├── plots/                            # cell_composition.png, QC, DE and enrichment figures (pathview maps in plots/pathview/)
├── reports/report.html               # HTML report of the whole analysis (rule `report`)
└── cache/                            # STRING and KEGG downloads (generated, not versioned)
```

The GEO supplementary files (count matrix 127 MB, gene list, cell table, and the authors' `.h5ad`, which is not
used) are not in the repository: they live in `$DATA/GSE216146/`.

## Pipeline

| Rule | Script (`../../scripts/de/`) | Main outputs |
|---|---|---|
| `prepare` | `01_prepare_scrna.R` | Pseudobulk counts; `data/dds_raw.rds`, `data/colData.csv`, `results/pseudobulk_summary.csv`, `results/cell_composition.csv`, `plots/cell_composition.png` |
| `differential` | `02_differential_expression_and_qc.R` | `data/dds.rds`, `results/top_table.csv`, `results/de_summary.tsv`, QC/DE plots |
| `kegg_annotation` | (curl) | KEGG gene–pathway table (`mmu`) in `cache/kegg/` (downloaded once) |
| `enrichment` | `03_enrichment.R` | `results/GSEA_*.csv`, `results/*_PPI.csv`, `results/enrichment_STRING_full.csv`, `reports/sessionInfo.txt`, enrichment plots |
| `report` | `report.py` | `reports/report.html` |

## Results

CIS vs PBS, 14,789 genes tested (≥ 10 pseudobulk counts in ≥ 3 samples), 965–2,612 cells and 7.0–15.9 million
UMIs per sample:

- The groups separate on the PCA; PC2 follows the batch (D20 vs D21).
- **723 DEGs** at `padj < 0.05` and `|log2FC| > 0.5` (413 up, 310 down in CIS); 1,329 genes at `padj < 0.05`.
  Top genes: *Vwf*, *Depp1*, *Sult1a1*, *Adm*, *Cxcl12*, *Itgad* up; *Hist1h1e*, *Cryab*, *Shtn1*, *Map1b*
  down.
- **GSEA:** 1,528 GO terms (1,009 up / 519 down in CIS) and 133 KEGG pathways (115 / 18) at `padj < 0.05`.
  Higher in CIS: blood vessel development and angiogenesis, extracellular space, cytosolic ribosome, HIF-1
  signaling; lower: synaptic membrane, synapse organization, axon, glutamatergic synapse.
- STRING: 667 of the DEGs mapped, 986 edges (score ≥ 700).

Caveats: the pseudobulk mixes cell types in the proportions captured. Endothelial cells are a larger share of the
cells in CIS (9.6–14.7% vs 5.4–9.4% in PBS), so part of the vascular signal is composition. *Vwf* (4.8-fold),
though, rises more than the endothelial share. Choroid plexus genes (*Ttr*, *Folr1*) are still higher in CIS
even without choroid plexus cells, probably ambient RNA from CN2 (38% choroid plexus cells). Neurons are few
(4.6–7.0% of the cells used), as usual in single-cell suspensions of adult brain.

## Running

```bash
conda activate genomic-variants      # Snakemake 9
cd bulk-rna-seq/analysis_from_geo_data/GSE216146
snakemake -n --use-conda             # dry-run
snakemake --use-conda --cores 2
```

`$DATA` must point to the folder with `GSE216146/` (`GSE216146_expression_matrix.mtx.gz`,
`GSE216146_genes.tsv.gz`, `GSE216146_cell_metadata.tsv.gz`). The `prepare` step downloads the series matrix
and the `enrichment` step needs internet access (STRING, pathview).

## Data and citation

The analysis uses the public dataset [GSE216146](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE216146)
from NCBI GEO (PMID 38446899). Please cite the original authors when using these data.
