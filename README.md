# Bulk RNA-seq: chemotherapy-induced transcriptional changes in the rodent brain

Reproducible Snakemake + R workflows for bulk RNA-seq analysis of public GEO datasets on
chemotherapy-induced cognitive impairment ("chemobrain"). They cover everything from raw reads to
differential expression, functional enrichment and protein-protein interaction (PPI) networks.

The repository contains eight projects. Each has its own Snakefile, conda environments and README:

| Folder | Dataset | Input | What it does |
|---|---|---|---|
| [`raw_reads_processing/`](raw_reads_processing) | [GSE295281](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE295281) | SRA runs | Download → FASTQ → FastQC/MultiQC → Salmon → tximport (gene-level counts) |
| [`bulk_rna_seq_tximport/`](bulk_rna_seq_tximport) | GSE295281 | tximport object from `raw_reads_processing` | DESeq2, QC, GSEA (GO/KEGG), pathview, STRING; configured through `config.yaml` |
| [`analysis_from_geo_data/GSE178812/`](analysis_from_geo_data/GSE178812) | [GSE178812](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE178812) | Count matrix deposited in GEO | DESeq2, QC, GSEA (GO/KEGG), pathview, STRING |
| [`analysis_from_geo_data/GSE236404/`](analysis_from_geo_data/GSE236404) | [GSE236404](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE236404) | Normalized count table deposited in GEO (rat) | limma-trend, QC, GSEA (GO/KEGG), pathview, STRING |
| [`analysis_from_geo_data/GSE267844/`](analysis_from_geo_data/GSE267844) | [GSE267844](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE267844) | FPKM table deposited in GEO (rat) | limma-trend, QC, GSEA (GO/KEGG), pathview, STRING |
| [`analysis_from_geo_data/GSE325554/`](analysis_from_geo_data/GSE325554) | [GSE325554](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE325554) | Visium spatial counts (mouse) | Hippocampus pseudobulk, then DESeq2, QC, GSEA (GO/KEGG), pathview, STRING |
| [`analysis_from_geo_data/GSE216146/`](analysis_from_geo_data/GSE216146) | [GSE216146](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE216146) | 10x single-cell counts (mouse) | Hippocampus pseudobulk (choroid plexus left out), then DESeq2 with batch, QC, GSEA (GO/KEGG), pathview, STRING |
| [`meta_analysis/`](meta_analysis) | all five | DE results of the five DE projects | DEG table and overlap across datasets, genes matched by 1:1 mouse–rat orthologs; random-effects meta-analysis per gene |

```
SRA (GSE295281) ──► raw_reads_processing ──► data/txi_gene.rds ──► bulk_rna_seq_tximport ──► DE + enrichment
GEO counts (GSE178812) ──────────────────────────────────────────► GSE178812 ─► DE + enrichment
GEO normalized counts (GSE236404) ───────────────────────────────► GSE236404 ─► DE (limma) + enrichment
GEO FPKM (GSE267844) ────────────────────────────────────────────► GSE267844 ─► DE (limma) + enrichment
GEO Visium (GSE325554) ──► hippocampus pseudobulk ───────────────► GSE325554 ─► DE + enrichment
GEO scRNA-seq (GSE216146) ──► sample pseudobulk ─────────────────► GSE216146 ─► DE + enrichment
DE results of the five projects ─────────────────────────────────► meta_analysis ─► DEG table, overlap, meta-analysis
```

## Reports

Each project has an HTML report that walks through the analysis from the input data to the last step,
with figures and tables:

| Project | Report |
|---|---|
| `raw_reads_processing` | [SRA → FastQC → Salmon → count table](https://arlindomatias.github.io/bulk-rna-seq/raw_reads_processing/reports/report.html) |
| `bulk_rna_seq_tximport` | [GSE295281: differential expression and enrichment](https://arlindomatias.github.io/bulk-rna-seq/bulk_rna_seq_tximport/reports/report.html) |
| `GSE178812` | [GSE178812: differential expression and enrichment](https://arlindomatias.github.io/bulk-rna-seq/analysis_from_geo_data/GSE178812/reports/report.html) |
| `GSE236404` | [GSE236404: differential expression and enrichment](https://arlindomatias.github.io/bulk-rna-seq/analysis_from_geo_data/GSE236404/reports/report.html) |
| `GSE267844` | [GSE267844: differential expression and enrichment](https://arlindomatias.github.io/bulk-rna-seq/analysis_from_geo_data/GSE267844/reports/report.html) |
| `GSE325554` | [GSE325554: differential expression and enrichment](https://arlindomatias.github.io/bulk-rna-seq/analysis_from_geo_data/GSE325554/reports/report.html) |
| `GSE216146` | [GSE216146: differential expression and enrichment](https://arlindomatias.github.io/bulk-rna-seq/analysis_from_geo_data/GSE216146/reports/report.html) |
| `meta_analysis` | [Chemotherapy DEGs across datasets](https://arlindomatias.github.io/bulk-rna-seq/meta_analysis/reports/report.html) |

The reports are built by the `report` rule of each Snakefile (`reports/report.html`) and published with
GitHub Pages.

## Project layout

Every project follows the same folder layout (other folders, such as `logs/` or `cache/`, are added when needed).
The scripts are not inside the projects: they live in the repository-level [`scripts/`](scripts) folder,
shared by the projects (see below).

| Folder | Contents |
|---|---|
| `data/` | Inputs (count tables, sample metadata), intermediate objects (`dds.rds`, `txi_gene.rds`) and small reference files |
| `results/` | Result tables only (`.csv`/`.tsv`) |
| `plots/` | All figures |
| `envs/` | Conda environments; each `envs/<name>.yaml` creates the environment `<name>` and has a `<name>.linux-64.pin.txt` lock when pinned |
| `reports/` | `report.html`: self-contained HTML report explaining the results with figures, from the first to the last step (built in Python, so it works for projects without R) |

Large files (SRA, FASTQ, genomes, indices, Salmon outputs) are never stored in the repository: they live under `$DATA`.

### Shared scripts

```
scripts/
├── html_report.py                         # HTML report builder used by every project's report rule
├── de/                                    # DE projects (bulk_rna_seq_tximport and the GEO projects)
│   ├── 00_common.R                        # config, logging, sample metadata (sourced by every step)
│   ├── 01_prepare_tximport.R              # tximport object -> DESeq2 object
│   ├── 01_prepare_geo.R                   # GEO count table + series matrix -> DESeq2 object
│   ├── 01_prepare_geo_normalized.R        # GEO normalized/FPKM table (+ sample sheet) -> SummarizedExperiment
│   ├── 01_prepare_visium.R                # Visium Space Ranger outputs -> region pseudobulk -> DESeq2 object
│   ├── 01_prepare_scrna.R                 # single-cell counts + cell table -> sample pseudobulk -> DESeq2 object
│   ├── 02_differential_expression_and_qc.R  # DESeq2
│   ├── 02_differential_expression_limma.R   # limma-trend (normalized values only)
│   ├── 03_enrichment.R                    # GSEA GO/KEGG, pathview, STRING
│   └── report.py
├── raw/                                   # raw_reads_processing
│   ├── tx2gene.py, qc_summary.py, tximport.R
│   ├── qc_plots.py
│   └── report.py
└── meta/                                  # meta_analysis
    ├── gene_map.R                         # cross-dataset gene key (mouse Ensembl; rat via 1:1 Compara orthologs)
    ├── collect_degs.R                     # DEGs of every DE project -> results/all_degs.csv
    ├── meta_random_effects.R              # random-effects meta-analysis per gene (metafor)
    ├── plots.py                           # DEG counts, pairwise overlap, shared-DEG heatmap
    └── report.py
```

The DE projects run the same scripts (GSE236404 and GSE267844 swap the prepare and DE steps for their limma
versions, `de.method: limma`; GSE325554 and GSE216146 have their own Visium and single-cell prepare steps); everything specific to a dataset (input, groups, gene identifier type,
cutoffs, filter, enrichment settings, seed) is in its `config.yaml`. Each project keeps its own complete
`Snakefile`.

## Datasets

- **GSE295281**: hippocampus of 12-week-old female C57BL/6J mice, one month after AC-T chemotherapy
  (doxorubicin + cyclophosphamide, then paclitaxel) vs saline. n = 6 saline vs 5 AC-T, paired-end,
  NextSeq 2000, about 31–41 M read pairs per sample. PMID 41519984.
- **GSE178812**: hippocampus of young C57BL/6J mice, doxorubicin (CHEM) vs sham (VEH). The series also has a
  doxorubicin + MitoQ group. PMID 34147619. The analysis starts from the count table the authors deposited in GEO.
- **GSE236404**: dentate gyrus of young (6-month) male F344 rats, doxorubicin (DOX, n = 8) vs vehicle (VEH, n = 10).
  The series also has a doxorubicin + ABT-263 (senolytic) group. PMID 38225896. GEO only has normalized counts,
  so the analysis uses limma-trend.
- **GSE267844**: hippocampus of rats, docetaxel (DTX, n = 3) vs saline (NC, n = 3). The series also has a
  docetaxel + edaravone dexborneol group. GEO only has FPKM, so the analysis uses limma-trend on log2(FPKM + 1).
- **GSE325554**: brains of p16-3MR mice 6 months after paclitaxel (PTX, n = 2) vs vehicle (Ctrl, n = 2),
  10x Visium spatial transcriptomics. The series also has PTX + ABT263 and PTX + ganciclovir groups. The
  hippocampal spots of each section are summed (pseudobulk) and analysed with DESeq2.
  Its spatial, spot-level analysis is in the separate [spatial-transcriptomics](https://github.com/arlindomatias/spatial-transcriptomics) repository.
- **GSE216146**: hippocampus of 9-week-old female C57BL/6J mice, cisplatin (CIS, n = 3) vs PBS (n = 3), day 22,
  10x single-cell RNA-seq. The series also tests 40 Hz sensory stimulation (GENUS) in both groups; the
  stimulated groups are not used. PMID 38446899. The cells of each sample (except choroid plexus and ependymal
  cells) are summed (pseudobulk) and analysed with DESeq2, with batch as a covariate.

## Projects

### 1. `raw_reads_processing`: from SRA to a gene count table

- Downloads `.sra` files from the NCBI AWS mirror and checks them with `vdb-validate`, then runs `fasterq-dump`.
  FASTQ files are temporary: they go to a local SSD scratch folder and are deleted after FastQC and Salmon.
- Builds a Salmon index on the GENCODE mouse transcriptome (release M35, `--gencode`, k = 31).
- Runs `salmon quant` with automatic library type detection, `--gcBias` and `--seqBias`. A sample fails the
  run if its mapping rate is below `min_mapping_rate` (30%). All samples mapped at about 92% (library type ISR);
  see `results/qc_summary.tsv`.
- Summarizes to gene level with tximport and writes counts, length-scaled TPM counts, TPM and lengths
  (`results/*.tsv`) plus `data/txi_gene.rds`. It also writes a MultiQC report to `reports/`.

Large files (SRA, reference, index, Salmon output) are written under `$DATA/<data_subdir>`, so the `DATA`
environment variable must be set:

```bash
cd raw_reads_processing
export DATA=/path/to/large/storage
snakemake -n --use-conda
snakemake --use-conda --cores 8 --resources download=3 scratch_dump=1
```

### 2. `bulk_rna_seq_tximport`: DE and enrichment from tximport

Runs the shared DE scripts (`scripts/de/`) on the tximport object of `raw_reads_processing`. Everything specific
to a dataset lives in `config.yaml`: groups, covariates, excluded samples, filters, cutoffs and enrichment settings.

- `DESeqDataSetFromTximport` (with transcript-length offsets), Ensembl gene IDs, apeglm-shrunken log2FC.
- KEGG annotation is downloaded once into `cache/kegg/`, then used for GSEA with `TERM2GENE` (both DE projects).
- QC: sample-distance heatmap, labelled PCA and a marker-gene plot (*Xist*, *Ddx3y*, glial markers).
- Current result (AC-T vs saline, 18,767 genes tested): **262 up- and 408 down-regulated DEGs**
  (`padj < 0.05`, `|log2FC| > 0.5`).

```bash
cd bulk_rna_seq_tximport
snakemake --use-conda --cores 2
# another dataset:
snakemake --use-conda --cores 2 --configfile other.yaml
```

### 3. `GSE178812`: DE and enrichment from GEO counts

Same scripts as `bulk_rna_seq_tximport`, with a different first step (`01_prepare_geo.R`): it downloads the
GSE178812 count table and series matrix, takes the groups from the GEO sample titles and uses gene symbols as
identifiers. Then DESeq2 (`~ group`, VEH as reference), GSEA (GO and KEGG), pathview, a STRING v12 PPI
network and STRING enrichment. Current result: 7 up- and 51 down-regulated DEGs (CHEM vs VEH).

```bash
cd analysis_from_geo_data/GSE178812
snakemake --use-conda --cores 2
```

### 4. `GSE236404`: DE (limma) and enrichment from GEO normalized counts

GEO only has a normalized (non-integer) count table for this series, so DESeq2 cannot be used.
`01_prepare_geo_normalized.R` downloads it with the sample sheet, drops the DOX + ABT-263 group and sums rows
that share an Ensembl gene ID; `02_differential_expression_limma.R` fits limma-trend on log2(normalized + 1).
Enrichment and report are the shared scripts (rat: `org.Rn.eg.db`, KEGG `rno`, STRING taxon 10116).
Current result (DOX vs VEH, 14,569 genes tested): no gene at `padj < 0.05`, so the DEGs use the unadjusted
`p < 0.01` (169 exploratory DEGs: 49 up, 120 down); GSEA finds 293 GO terms and 27 KEGG pathways.

```bash
cd analysis_from_geo_data/GSE236404
snakemake --use-conda --cores 2
```

### 5. `GSE267844`: DE (limma) and enrichment from a GEO FPKM table

Same pipeline as `GSE236404`: `01_prepare_geo_normalized.R` reads the GEO `.xlsx` table and takes the samples
from its `FPKM.*` columns (`geo.sample_patterns`), and `02_differential_expression_limma.R` fits limma-trend
on log2(FPKM + 1). No gene reaches `padj < 0.05`, so the DEGs use the exploratory unadjusted `p < 0.01`.

```bash
cd analysis_from_geo_data/GSE267844
snakemake --use-conda --cores 2
```

### 6. `GSE325554`: Visium hippocampus pseudobulk

`01_prepare_visium.R` reads the Space Ranger outputs from `$DATA/GSE325554/GSE325554_RAW`, selects the
hippocampal spots of each section by marker genes (spatially connected high-score spots, see
`plots/region_spots.png`) and sums their counts; the pseudobulk then goes through the shared DESeq2,
enrichment and report scripts. With n = 2 per group the results are exploratory.

This is the **bulk (pseudobulk)** view of the dataset, comparable with the other hippocampus projects. The
**spatial** analysis of the same sections (spot level, all brain regions, niches and hotspots) lives in the
[spatial-transcriptomics](https://github.com/arlindomatias/spatial-transcriptomics) repository ([report](https://arlindomatias.github.io/spatial-transcriptomics/GSE325554/reports/report.html)).

```bash
cd analysis_from_geo_data/GSE325554
snakemake --use-conda --cores 2
```

### 7. `GSE216146`: single-cell hippocampus pseudobulk

`01_prepare_scrna.R` reads the raw 10x counts and the authors' cell table from `$DATA/GSE216146`, leaves out
choroid plexus and ependymal cells (their share follows the dissection: `plots/cell_composition.png`) and sums
the counts of the remaining cells per sample; DESeq2 uses `~ batch + group` (the samples were sequenced in two
batches, each with every group). The meta-analysis does not include this dataset yet.

```bash
cd analysis_from_geo_data/GSE216146
snakemake --use-conda --cores 2
```

### 8. `meta_analysis`: DEGs across datasets and meta-analysis

Collects the DEGs of the five DE projects (same DEG rule as each project, read from its `de_summary.tsv`)
into `results/all_degs.csv`: one row per gene and dataset where it is a DEG, with symbol, Entrez ID, Ensembl ID,
log2FC, statistic (Wald or moderated t), p-value and padj. Genes are matched across datasets and species by
the mouse Ensembl ID (rat genes through their 1:1 ortholog from Ensembl Compara 104; `results/gene_map.csv`).
`results/pairwise_overlap.csv` tests the DEG overlap of each pair of datasets (hypergeometric), and
`results/meta_random_effects.csv` combines the log2FC of every gene tested in at least two datasets in a
random-effects model (metafor). Run it after the DE projects:

```bash
cd meta_analysis
snakemake --use-conda --cores 1
```

## Shared methods

The DE pipelines use the same core settings (GSE236404: limma-trend instead of DESeq2, moderated t for GSEA,
rat annotation):

- **Filtering:** genes with ≥ 10 counts in at least *k* samples (3 for GSE178812; the smallest group size for GSE295281).
- **DESeq2:** `fitType = "local"`. DEGs are genes with `padj < 0.05` and `|log2FC| > 0.5`. When no gene has
  `padj < 0.05`, the DEGs use the unadjusted `p < 0.01` instead (`de.pvalue_fallback`, exploratory); the report
  and `results/de_summary.tsv` (`deg_criterion`, `deg_cutoff`) say which rule was used.
  A positive log2FC means higher expression in the treated group.
- **GSEA:** genes are ranked by the Wald statistic. GO (BP, CC, MF) and KEGG `mmu` are tested with
  `clusterProfiler`, with Benjamini-Hochberg adjustment and a fixed seed. Positive NES = enriched in the treated group.
- **STRING:** v12.0, *Mus musculus* (taxon 10090), combined score ≥ 700, all tested genes as background.
  Edge and node tables are written for Cytoscape.

Outputs go to `results/` (tables) and `plots/` (PCA, density, boxplot, dispersion, Cook's distance, volcano,
MA plot, heatmap, FDR histogram, QQ plot, GO/KEGG dotplots, emaplots, cnetplots, ridgeplots, GSEA plots,
pathview maps, PPI network, Reactome and WikiPathways enrichment).

## Requirements

- Linux (x86-64), Conda/Mamba and Snakemake ≥ 9. Each rule creates its own environment with `--use-conda`.
- Pinned versions: Salmon 2.11.0, sra-tools 3.4.1, FastQC 0.13.0, MultiQC 1.35, R 4.5 / tximport 1.38.2.
  The R/Bioconductor stack for DE is pinned in `analysis_from_geo_data/GSE178812/envs/rna-seq.linux-64.pin.txt` and
  `bulk_rna_seq_tximport/envs/deseq2.linux-64.pin.txt` (`analysis_from_geo_data/GSE236404/envs/limma.linux-64.pin.txt`
  adds `org.Rn.eg.db`).
- Internet access (SRA, GENCODE, GEO, KEGG REST, STRING).
- For `raw_reads_processing`: about 25 GB of fast scratch space per sample running at the same time, plus
  room under `$DATA` for SRA files, the Salmon index and quantifications.

## Reproducing the full GSE295281 analysis

```bash
git clone git@github.com:arlindomatias/bulk-rna-seq.git
cd bulk-rna-seq

export DATA=/path/to/large/storage
(cd raw_reads_processing  && snakemake --use-conda --cores 8 --resources download=3 scratch_dump=1)
(cd bulk_rna_seq_tximport && snakemake --use-conda --cores 2)
```

`*.rds` files are not versioned (see `.gitignore`), so `txi_gene.rds` has to be regenerated by
`raw_reads_processing` before `bulk_rna_seq_tximport` can run on a fresh clone.

## Data and citation

All datasets are public on NCBI GEO. Please cite the original authors of
[GSE295281](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE295281),
[GSE178812](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE178812),
[GSE236404](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE236404),
[GSE267844](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE267844),
[GSE325554](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE325554) and
[GSE216146](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE216146) when using these data.

## License

[MIT](LICENSE) © 2026 Arlindo César Matias Pereira
