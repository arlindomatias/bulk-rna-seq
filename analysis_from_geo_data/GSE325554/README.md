# Paclitaxel and the mouse hippocampus: Visium pseudobulk (GSE325554)

Differential expression and functional enrichment of the public spatial transcriptomics dataset
**GSE325554** (NCBI GEO, BioProject PRJNA1440463): coronal brain hemisections of p16-3MR mice treated with
paclitaxel (**PTX**, 5 mg/kg/day, 2 cycles) or vehicle (**Ctrl**), collected 6 months after treatment
(10x Visium, Space Ranger 1.2.1), n = 2 per group. The series also has PTX + ABT263 and PTX + ganciclovir
(senolytic) groups, which are not used here.

**Report:** [arlindomatias.github.io/bulk-rna-seq/analysis_from_geo_data/GSE325554/reports/report.html](https://arlindomatias.github.io/bulk-rna-seq/analysis_from_geo_data/GSE325554/reports/report.html)

> **Bulk-style analysis here; spatial analysis in another repository.** This project sums the hippocampal spots
> of each section (**pseudobulk**) and analyses them like bulk RNA-seq (DESeq2), so it can be compared with the
> other hippocampus datasets of this repository; the spot positions are only used to select the region. The
> **spatial** analysis of the same sections (PTX vs Ctrl, spot level, every brain region: anatomical domains,
> spatially variable genes, niches, senescence hotspots, DE inside each domain) is in
> [spatial-transcriptomics](https://github.com/arlindomatias/spatial-transcriptomics/tree/main/GSE325554)
> ([report](https://arlindomatias.github.io/spatial-transcriptomics/GSE325554/reports/report.html)).

## Hippocampus pseudobulk

Each Visium section covers several brain regions and is cut at a different level and orientation, so the
counts of a whole section are not comparable between mice. The other datasets of this repository are
hippocampus, so only the **hippocampal spots** are used (`../../scripts/de/01_prepare_visium.R`, settings in
`visium:` of `config.yaml`):

1. per spot, the log-normalized expression of hippocampal markers (*Prox1*, *C1ql2*, *Dsp*, *Spink8*,
   *Fibcd1*, *Zbtb20*, *Hpca*) is z-scored across the section; score = mean z;
2. spots with score > 0.75 are grouped into connected components on the Visium grid; components with ≥ 10
   spots are kept (a coronal section can cut the hippocampus in two pieces);
3. the selection grows by one ring of neighbouring spots (layers with fewer marker transcripts).

The UMI counts of the selected spots are summed per section and analysed with the shared DESeq2 pipeline.
`plots/region_spots.png` shows the selection on each section; `results/region_spots.csv` has every spot.

## Repository structure

```
.
├── Snakefile                         # Snakemake workflow that runs the whole pipeline
├── config.yaml                       # Dataset settings (raw folder, sample patterns, visium, groups, cutoffs, enrichment)
├── envs/
│   ├── rna-seq.yaml                  # Conda environment `rna-seq` (same as ../GSE178812)
│   ├── rna-seq.linux-64.pin.txt      # Exact package versions (linux-64)
│   └── report.yaml                   # Python environment of the report rule
├── data/
│   ├── pseudobulk_counts.csv         # Hippocampus pseudobulk counts (Ensembl ID × section)
│   ├── colData.csv                   # Sample metadata used in the model
│   └── dds_raw.rds, dds.rds          # DESeq2 objects (generated, not versioned)
├── results/                          # region_spots.csv, pseudobulk_summary.csv, top_table.csv, de_summary.tsv, GSEA, STRING, PPI tables
├── plots/                            # region_spots.png, QC, DE and enrichment figures (pathview maps in plots/pathview/)
├── reports/report.html               # HTML report of the whole analysis (rule `report`)
└── cache/                            # STRING and KEGG downloads (generated, not versioned)
```

Raw Space Ranger files (count matrices, spot positions, images; 371 MB) are not in the repository: they live
in `$DATA/GSE325554/GSE325554_RAW` (GEO `GSE325554_RAW.tar`, extracted).

## Pipeline

| Rule | Script (`../../scripts/de/`) | Main outputs |
|---|---|---|
| `prepare` | `01_prepare_visium.R` | Hippocampal spots, pseudobulk counts; `data/dds_raw.rds`, `data/colData.csv`, `results/region_spots.csv`, `results/pseudobulk_summary.csv`, `plots/region_spots.png` |
| `differential` | `02_differential_expression_and_qc.R` | `data/dds.rds`, `results/top_table.csv`, `results/de_summary.tsv`, QC/DE plots |
| `kegg_annotation` | (curl) | KEGG gene–pathway table (`mmu`) in `cache/kegg/` (downloaded once) |
| `enrichment` | `03_enrichment.R` | `results/GSEA_*.csv`, `results/*_PPI.csv`, `results/enrichment_STRING_full.csv`, `reports/sessionInfo.txt`, enrichment plots |
| `report` | `report.py` | `reports/report.html` |

## Results

PTX vs Ctrl, 11,632 genes tested (≥ 10 pseudobulk counts in ≥ 2 sections):

- **Spot selection:** 123–205 hippocampal spots per section (0.7–4.2 million UMIs). **Ctrl1 is a weak
  section**: median 2.6 k UMIs per selected spot vs 11–17 k in the others, and it separates from the other
  three on PC1 (71% of the variance). With n = 2 per group, one weak control drives much of the contrast.
- **4 DEGs** at `padj < 0.05` (3 up, 1 down): *Tac2*, *Hbb-bs*, *Hbb-bt* up and *Phyhip* down in PTX.
  Haemoglobin genes point to blood content, and *Tac2* to neighbouring tissue in the selection, rather than
  to a hippocampal response to paclitaxel. Fewer genes fall below p = 0.01 (52) than expected by chance
  (~116): power is very low.
- **GSEA:** 362 GO terms (182 up / 180 down in PTX) and 44 KEGG pathways (18 / 26) at `padj < 0.05`.
  Higher in PTX: ribosome, oxidative phosphorylation, electron transport chain; lower: synapse and dendrite
  organization, FoxO signaling, dopaminergic synapse. Ribosome/OXPHOS genes also track sequencing depth, so
  the Ctrl1 difference may contribute.
- Treat this dataset as exploratory.

## Running

```bash
conda activate genomic-variants      # Snakemake 9
cd bulk-rna-seq/analysis_from_geo_data/GSE325554
snakemake -n --use-conda             # dry-run
snakemake --use-conda --cores 2
```

`$DATA` must point to the folder with `GSE325554/GSE325554_RAW`. The `enrichment` step needs internet access
(STRING, pathview).

## Data and citation

The analysis uses the public dataset [GSE325554](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE325554)
from NCBI GEO. Please cite the original authors when using these data.
