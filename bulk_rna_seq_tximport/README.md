# Bulk RNA-seq from Salmon/tximport: differential expression and enrichment

Differential expression and enrichment of **Salmon quantifications summarized with tximport**, configured per
dataset in `config.yaml`. It runs the same scripts (`../scripts/de/`) as
[`analysis_from_geo_data`](../analysis_from_geo_data), which starts from a GEO count table instead.

Current dataset: **GSE295281**: hippocampus of 12-week-old female C57BL/6J mice, one month after AC-T
chemotherapy (doxorubicin + cyclophosphamide, then paclitaxel) vs saline (n = 5 vs 6). The input comes from
[`../raw_reads_processing`](../raw_reads_processing) (Salmon 2.11.0, GENCODE M35).

**Report:** [arlindomatias.github.io/bulk-rna-seq/bulk_rna_seq_tximport/reports/report.html](https://arlindomatias.github.io/bulk-rna-seq/bulk_rna_seq_tximport/reports/report.html)

## Differences from `analysis_from_geo_data`

Both projects run the same scripts. Only the input differs, and that sets a few options in `config.yaml`:

| | `analysis_from_geo_data` | this project |
|---|---|---|
| Input | GEO xlsx of counts by gene symbol (`01_prepare_geo.R`) | `txi_gene.rds` (tximport) + `samples.tsv` + `tx2gene.tsv` (`01_prepare_tximport.R`) |
| DESeq2 object | `DESeqDataSetFromMatrix` | `DESeqDataSetFromTximport` (uses transcript-length offsets) |
| Gene IDs (`gene_id_type`) | `SYMBOL` (duplicates dropped) | `ENSEMBL` without version; symbol kept as annotation |
| GSEA keys | `SYMBOL` for GO, Entrez for KEGG | `ENSEMBL` for GO, Entrez for KEGG (ties resolved by largest \|stat\|) |
| Filter | ≥ 10 counts in 3 samples | ≥ 10 counts in the smallest group size |
| Marker-gene plot | not configured | *Xist*, *Ddx3y*, glial markers |

## Usage

```bash
conda activate genomic-variants        # Snakemake 9
snakemake -n --use-conda               # dry-run
snakemake --use-conda --cores 2
```

The R environment is built from `envs/deseq2.linux-64.pin.txt` (exact package versions, copied from
`../analysis_from_geo_data/envs/rna-seq.linux-64.pin.txt`); `envs/deseq2.yaml` (env `deseq2`) is the unpinned fallback. Step 03 needs internet access (STRING, pathview).
The KEGG REST API can be very slow (about 10 min for the 1 MB mouse gene-pathway list on 2026-11-12), so the
annotation is cached in `cache/kegg/`; delete that folder to refresh it.

### Another dataset

Copy `config.yaml`, point `input:` to the tximport object and sample table, set the condition column and
group labels (syntactic R names, e.g. `ACT` instead of `AC-T`) and run with
`snakemake --use-conda --cores 2 --configfile other.yaml`. Outputs go to `results/` and `plots/`, so use
a separate checkout or working directory per dataset.

## Pipeline

| Rule | Script | Main outputs |
|---|---|---|
| `prepare` | `../scripts/de/01_prepare_tximport.R` | `data/dds_raw.rds`, `data/colData.csv` |
| `differential` | `../scripts/de/02_differential_expression_and_qc.R` | `data/dds.rds`, `results/top_table.csv`, `results/de_summary.tsv`, `results/normalized_counts.csv`, QC/DE plots |
| `kegg_annotation` | (curl) | `cache/kegg/{org}_link_pathway.tsv`, `cache/kegg/{org}_list_pathway.tsv` |
| `enrichment` | `../scripts/de/03_enrichment.R` | `results/GSEA_GO_all_full.csv`, `results/GSEA_KEGG_full.csv`, `results/enrichment_STRING_full.csv`, `results/{edges,nodes}_PPI.csv`, `reports/sessionInfo.txt`, enrichment plots |

`../scripts/de/00_common.R` is sourced by every step (logging, config, colors). Logs are written to `logs/{rule}.log`.

## Methods

- **Filtering:** genes with ≥ 10 counts in at least as many samples as the smallest group (5).
- **Differential expression:** DESeq2, design `~ group` with the reference group from `config.yaml`,
  `fitType = "local"`. DEGs: `padj < 0.05` and `|log2FC| > 0.5`. `top_table.csv` also has the apeglm-shrunken
  log2FC (used in the MA plot).
- **GSEA:** genes ranked by the Wald statistic; GO (BP, CC, MF) and KEGG with `clusterProfiler`,
  BH adjustment, fixed seed. Positive NES = enriched in the treatment group.
- **STRING:** v12.0, mouse, combined score ≥ 700; one gene per symbol (lowest padj); background = all tested genes.

## Notes for GSE295281

- The first digit of the mouse ID (2/3 = saline, 5/6 = AC-T) matches the `ACTGolgi_N` suffix in the GEO count
  matrix: probably cage or processing batch. It is **nested within the group**, so it cannot be added as a
  covariate.
- SRR33247678 (saline, mouse 300) has almost no *Xist*, although all mice are female and no Y-linked gene is
  expressed. Check `plots/marker_genes.png`, `plots/pca_vst.png` and `plots/sample_distances.png`; it can be
  removed with `exclude_samples` in a sensitivity analysis.

## Data

GEO [GSE295281](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE295281) (PMID 41519984). Please cite the
original authors.
