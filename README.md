# Bulk RNA-seq: chemotherapy-induced transcriptional changes in the mouse brain

Reproducible Snakemake + R workflows for bulk RNA-seq analysis of public GEO datasets on
chemotherapy-induced cognitive impairment ("chemobrain"). They cover everything from raw reads to
differential expression, functional enrichment and protein-protein interaction (PPI) networks.

The repository contains three projects. Each has its own Snakefile, conda environments and README:

| Folder | Dataset | Input | What it does |
|---|---|---|---|
| [`raw_reads_processing/`](raw_reads_processing) | [GSE295281](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE295281) | SRA runs | Download → FASTQ → FastQC/MultiQC → Salmon → tximport (gene-level counts) |
| [`bulk_rna_seq_tximport/`](bulk_rna_seq_tximport) | GSE295281 | tximport object from `raw_reads_processing` | DESeq2, QC, GSEA (GO/KEGG), pathview, STRING; configured through `config.yaml` |
| [`analysis_from_geo_data/`](analysis_from_geo_data) | [GSE178812](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE178812) | Count matrix deposited in GEO | DESeq2, QC, GSEA (GO/KEGG), pathview, STRING |

```
SRA (GSE295281) ──► raw_reads_processing ──► txi_gene.rds ──► bulk_rna_seq_tximport ──► DE + enrichment
GEO counts (GSE178812) ─────────────────────────────────────► analysis_from_geo_data ─► DE + enrichment
```

## Datasets

- **GSE295281**: hippocampus of 12-week-old female C57BL/6J mice, one month after AC-T chemotherapy
  (doxorubicin + cyclophosphamide, then paclitaxel) vs saline. n = 6 saline vs 5 AC-T, paired-end,
  NextSeq 2000, about 31–41 M read pairs per sample. PMID 41519984.
- **GSE178812**: hippocampus of young C57BL/6J mice, doxorubicin (CHEM) vs sham (VEH). The series also has a
  doxorubicin + MitoQ group. PMID 34147619. The analysis starts from the count table the authors deposited in GEO.

## Projects

### 1. `raw_reads_processing`: from SRA to a gene count table

- Downloads `.sra` files from the NCBI AWS mirror and checks them with `vdb-validate`, then runs `fasterq-dump`.
  FASTQ files are temporary: they go to a local SSD scratch folder and are deleted after FastQC and Salmon.
- Builds a Salmon index on the GENCODE mouse transcriptome (release M35, `--gencode`, k = 31).
- Runs `salmon quant` with automatic library type detection, `--gcBias` and `--seqBias`. A sample fails the
  run if its mapping rate is below `min_mapping_rate` (30%). All samples mapped at about 92% (library type ISR);
  see `results/qc_summary.tsv`.
- Summarizes to gene level with tximport and writes counts, length-scaled TPM counts, TPM and lengths
  (`results/*.tsv`) plus `results/txi_gene.rds`. It also writes a MultiQC report.

Large files (SRA, reference, index, Salmon output) are written under `$DATA/<data_subdir>`, so the `DATA`
environment variable must be set:

```bash
cd raw_reads_processing
export DATA=/path/to/large/storage
snakemake -n --use-conda
snakemake --use-conda --cores 8 --resources download=3 scratch_dump=1
```

### 2. `bulk_rna_seq_tximport`: DE and enrichment from tximport

A configurable version of the GEO pipeline. Everything specific to a dataset lives in `config.yaml`: groups,
covariates, excluded samples, filters, cutoffs and enrichment settings.

- `DESeqDataSetFromTximport` (with transcript-length offsets), Ensembl gene IDs, apeglm-shrunken log2FC.
- KEGG annotation is downloaded once into `cache/kegg/`, then used for GSEA with `TERM2GENE`.
- Extra QC: sample-distance heatmap, labelled PCA, marker-gene plot (*Xist*, *Ddx3y*, glial markers).
- Current result (AC-T vs saline, 18,767 genes tested): **262 up- and 408 down-regulated DEGs**
  (`padj < 0.05`, `|log2FC| > 0.5`).

```bash
cd bulk_rna_seq_tximport
snakemake --use-conda --cores 2
# another dataset:
snakemake --use-conda --cores 2 --configfile other.yaml
```

### 3. `analysis_from_geo_data`: DE and enrichment from GEO counts

The original pipeline. It downloads the GSE178812 count table and series matrix, then runs DESeq2
(`~ group`, VEH as reference), GSEA (GO and KEGG), pathview, a STRING v12 PPI network and STRING enrichment.

```bash
cd analysis_from_geo_data
snakemake --use-conda --cores 1
```

## Shared methods

Both DE pipelines use the same core settings:

- **Filtering:** genes with ≥ 10 counts in at least *k* samples (3 for GSE178812; the smallest group size for GSE295281).
- **DESeq2:** `fitType = "local"`. DEGs are genes with `padj < 0.05` and `|log2FC| > 0.5`.
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
  The R/Bioconductor stack for DE is pinned in `environment.linux-64.pin.txt` and
  `bulk_rna_seq_tximport/envs/deseq2.linux-64.pin.txt`.
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

Both datasets are public on NCBI GEO. Please cite the original authors of
[GSE295281](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE295281) and
[GSE178812](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE178812) when using these data.

## License

[MIT](LICENSE) © 2026 Arlindo César Matias Pereira
