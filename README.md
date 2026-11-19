# Bulk RNA-seq: chemotherapy-induced transcriptional changes in the rodent brain

Reproducible Snakemake + R workflows for bulk RNA-seq analysis of public GEO datasets on
chemotherapy-induced cognitive impairment ("chemobrain"). They cover everything from raw reads to
differential expression, functional enrichment and protein-protein interaction (PPI) networks.

The repository contains five projects. Each has its own Snakefile, conda environments and README:

| Folder | Dataset | Input | What it does |
|---|---|---|---|
| [`raw_reads_processing/`](raw_reads_processing) | [GSE295281](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE295281) | SRA runs | Download → FASTQ → FastQC/MultiQC → Salmon → tximport (gene-level counts) |
| [`bulk_rna_seq_tximport/`](bulk_rna_seq_tximport) | GSE295281 | tximport object from `raw_reads_processing` | DESeq2, QC, GSEA (GO/KEGG), pathview, STRING; configured through `config.yaml` |
| [`analysis_from_geo_data/GSE178812/`](analysis_from_geo_data/GSE178812) | [GSE178812](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE178812) | Count matrix deposited in GEO | DESeq2, QC, GSEA (GO/KEGG), pathview, STRING |
| [`analysis_from_geo_data/GSE236404/`](analysis_from_geo_data/GSE236404) | [GSE236404](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE236404) | Normalized count table deposited in GEO (rat) | limma-trend, QC, GSEA (GO/KEGG), pathview, STRING |
| [`meta_analysis/`](meta_analysis) | all three | DE results of the three DE projects | One table with every DEG (symbol, Entrez, Ensembl, log2FC, stat, p, padj per dataset) and their overlap, genes matched by 1:1 mouse–rat orthologs |

```
SRA (GSE295281) ──► raw_reads_processing ──► data/txi_gene.rds ──► bulk_rna_seq_tximport ──► DE + enrichment
GEO counts (GSE178812) ──────────────────────────────────────────► GSE178812 ─► DE + enrichment
GEO normalized counts (GSE236404) ───────────────────────────────► GSE236404 ─► DE (limma) + enrichment
DE results of the three projects ────────────────────────────────► meta_analysis ─► DEG table + overlap
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
├── de/                                    # DE projects (bulk_rna_seq_tximport, GSE178812, GSE236404)
│   ├── 00_common.R                        # config, logging, sample metadata (sourced by every step)
│   ├── 01_prepare_tximport.R              # tximport object -> DESeq2 object
│   ├── 01_prepare_geo.R                   # GEO count table + series matrix -> DESeq2 object
│   ├── 01_prepare_geo_normalized.R        # GEO normalized table + sample sheet -> SummarizedExperiment
│   ├── 02_differential_expression_and_qc.R  # DESeq2
│   ├── 02_differential_expression_limma.R   # limma-trend (normalized counts only)
│   ├── 03_enrichment.R                    # GSEA GO/KEGG, pathview, STRING
│   └── report.py
├── raw/                                   # raw_reads_processing
│   ├── tx2gene.py, qc_summary.py, tximport.R
│   ├── qc_plots.py
│   └── report.py
└── meta/                                  # meta_analysis
    ├── gene_map.R                         # cross-dataset gene key (mouse Ensembl; rat via 1:1 Compara orthologs)
    ├── collect_degs.R                     # DEGs of every DE project -> results/all_degs.csv
    ├── plots.py                           # DEG counts, pairwise overlap, shared-DEG heatmap
    └── report.py
```

The DE projects run the same scripts (GSE236404 swaps the prepare and DE steps for their limma versions,
`de.method: limma`); everything specific to a dataset (input, groups, gene identifier type,
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

### 5. `meta_analysis`: DEGs across datasets

Collects the DEGs of the three DE projects (same DEG rule as each project, read from its `de_summary.tsv`)
into `results/all_degs.csv`: one row per gene and dataset where it is a DEG, with symbol, Entrez ID, Ensembl ID,
log2FC, statistic (Wald or moderated t), p-value and padj. Missing IDs come from `org.Mm.eg.db` /
`org.Rn.eg.db`. `results/pairwise_overlap.csv` tests the DEG overlap of each pair of datasets
(hypergeometric, genes matched by symbol). Run it after the DE projects:

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
[GSE178812](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE178812) and
[GSE236404](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE236404) when using these data.

## License

[MIT](LICENSE) © 2026 Arlindo César Matias Pereira
