# ============================================================
# SCRIPT 01 (single-cell projects): cell × gene count matrix -> sample pseudobulk -> DESeqDataSet
# Writes the same outputs as the other 01 scripts: data/dds_raw.rds (mcols gene_id, SYMBOL) and
# data/colData.csv, so 02 and 03 do not depend on the data source.
#
# Input (GEO supplementary files, read from $DATA/<input.raw_dir>):
#   * input.matrix_file: MatrixMarket cells × genes raw counts (rows in the order of the cell table);
#   * input.genes_file: one gene per row (first column = gene name, order of the matrix columns);
#   * input.cells_file: cell metadata (first column = barcode) with the sample, condition, batch and cell type
#     columns named in config `scrna:`.
# Cells of the types in scrna.exclude_cell_types are dropped (e.g. choroid plexus, whose share depends on the
# dissection); the counts of the remaining cells are summed per sample (pseudobulk) and analysed as bulk
# RNA-seq. GEO accessions and titles come from the series matrix (downloaded into data/ if missing).
# Also writes results/pseudobulk_summary.csv, results/cell_composition.csv and plots/cell_composition.png.
# ============================================================
source(file.path(snakemake@scriptdir, "00_common.R"))
suppressPackageStartupMessages({
  library(DESeq2)
  library(Matrix)
})
icfg <- cfg$input
scfg <- cfg$scrna

data_root <- Sys.getenv("DATA")
if (!nzchar(data_root)) stop("Set the DATA environment variable (heavy files live in $DATA/<dataset>/)")
raw_file <- function(f) {
  p <- file.path(data_root, icfg$raw_dir, f)
  if (!file.exists(p)) stop("Missing input file: ", p)
  p
}

### 1. Samples (series matrix) ---------------------------------------------------
series <- snakemake@output$series
if (!file.exists(series) || file.size(series) == 0) {
  message("Downloading ", icfg$series_matrix_url)
  download.file(icfg$series_matrix_url, series, mode = "wb", quiet = TRUE)
}
sm <- readLines(gzfile(series))
sm_field <- function(tag) {
  lines <- grep(paste0("^!", tag, "\t"), sm, value = TRUE)
  lapply(lines, function(l) gsub('^"|"$', "", strsplit(l, "\t")[[1]][-1]))
}
geo <- tibble(GEO_ID = sm_field("Sample_geo_accession")[[1]], title = sm_field("Sample_title")[[1]])
# The sample id of the cell table is one of the Sample_description fields (e.g. D21-2746).
desc <- sm_field("Sample_description")

### 2. Cells and counts -----------------------------------------------------------
cells <- read_tsv(raw_file(icfg$cells_file), show_col_types = FALSE, progress = FALSE,
                  name_repair = "unique_quiet")
names(cells)[1] <- "barcode"
need <- unlist(scfg[c("sample_column", "condition_column", "batch_column", "cell_type_column")])
stopifnot(all(need %in% names(cells)))
genes <- read_tsv(raw_file(icfg$genes_file), show_col_types = FALSE, progress = FALSE,
                  name_repair = "unique_quiet")[[1]]
stopifnot(!anyDuplicated(genes))

# Matrix::readMM does not accept the "unsigned-integer" field type: read the triplets directly.
hdr <- readLines(gzfile(raw_file(icfg$matrix_file)), n = 30)
dims_line <- hdr[!startsWith(hdr, "%")][1]
dims <- as.integer(strsplit(trimws(dims_line), "\\s+")[[1]])
tr <- read_delim(raw_file(icfg$matrix_file), delim = " ", skip = which(hdr == dims_line),
                 col_names = c("i", "j", "x"), col_types = "iid", progress = FALSE)
stopifnot(nrow(tr) == dims[3], all(tr$x == round(tr$x)))
if (dims[1] == nrow(cells) && dims[2] == length(genes)) {
  m <- sparseMatrix(i = tr$j, j = tr$i, x = tr$x, dims = dims[2:1])        # cells × genes -> genes × cells
} else if (dims[1] == length(genes) && dims[2] == nrow(cells)) {
  m <- sparseMatrix(i = tr$i, j = tr$j, x = tr$x, dims = dims)
} else {
  stop("Matrix dimensions ", paste(dims[1:2], collapse = " × "), " do not match ", nrow(cells),
       " cells and ", length(genes), " genes")
}
rm(tr)
dimnames(m) <- list(genes, cells$barcode)
# Row order check: genes detected per cell must follow the QC column of the cell table, when present.
if ("n_genes_by_counts" %in% names(cells)) {
  r <- cor(colSums(m > 0), cells$n_genes_by_counts)
  message("Correlation of detected genes per cell with n_genes_by_counts: ", signif(r, 6))
  if (r < 0.99) stop("Cell order of the matrix does not match the cell table")
}
message("Count matrix: ", nrow(m), " genes × ", ncol(m), " cells")

sample_of <- cells[[scfg$sample_column]]
cond_of   <- cells[[scfg$condition_column]]
type_of   <- cells[[scfg$cell_type_column]]

sample_tab <- cells |>
  distinct(sample = .data[[scfg$sample_column]], condition = .data[[scfg$condition_column]],
           batch = .data[[scfg$batch_column]])
stopifnot(!anyDuplicated(sample_tab$sample))
desc_hit <- vapply(desc, function(d) mean(sample_tab$sample %in% d), numeric(1))
if (!any(desc_hit == 1)) stop("No Sample_description field of the series matrix matches the sample ids")
geo$sample <- desc[[which(desc_hit == 1)[1]]]
sample_tab <- left_join(sample_tab, geo, by = "sample") |>
  mutate(replicate = sub(".*replicate\\s*", "", title),
         label = paste0(condition, replicate)) |>
  arrange(condition, replicate)

excluded <- unlist(cfg$geo$exclude_conditions)
if (length(excluded)) {
  message("Samples excluded by condition (", paste(excluded, collapse = ", "), "): ",
          paste(sample_tab$sample[sample_tab$condition %in% excluded], collapse = ", "))
}
samples <- sample_tab[!sample_tab$condition %in% excluded, ]

### 3. Pseudobulk per sample ----------------------------------------------------------
drop_types <- unlist(scfg$exclude_cell_types)
stopifnot(all(drop_types %in% type_of))
use <- sample_of %in% samples$sample & !type_of %in% drop_types
message("Cell types excluded: ", paste(drop_types, collapse = ", "), " (",
        sum(sample_of %in% samples$sample & type_of %in% drop_types), " cells)")

grp <- factor(sample_of[use], levels = samples$sample)
agg <- sparseMatrix(i = seq_along(grp), j = as.integer(grp), x = 1,
                    dims = c(length(grp), nlevels(grp)), dimnames = list(NULL, levels(grp)))
counts <- as.matrix(m[, use] %*% agg)
storage.mode(counts) <- "integer"

umi <- colSums(m)
summary_tab <- samples |>
  rowwise() |>
  mutate(cells_total = sum(sample_of == sample),
         cells_excluded = sum(sample_of == sample & type_of %in% drop_types),
         cells_used = sum(sample_of == sample & !type_of %in% drop_types),
         umi_used = sum(umi[sample_of == sample & !type_of %in% drop_types]),
         median_umi_per_cell = median(umi[sample_of == sample & !type_of %in% drop_types])) |>
  ungroup() |>
  dplyr::select(sample, GEO_ID, title, condition, batch, cells_total, cells_excluded, cells_used,
                umi_used, median_umi_per_cell)
print(as.data.frame(summary_tab))

composition <- tibble(sample = sample_of, cell_type = type_of) |>
  filter(sample %in% samples$sample) |>
  dplyr::count(sample, cell_type, name = "cells") |>
  group_by(sample) |>
  mutate(fraction = cells / sum(cells)) |>
  ungroup() |>
  mutate(excluded = cell_type %in% drop_types)

### 4. DESeqDataSet ----------------------------------------------------------------
coldata <- build_coldata(samples, colnames(counts))
counts <- counts[, coldata$sample_id, drop = FALSE]

dds <- DESeqDataSetFromMatrix(countData = counts, colData = coldata, design = build_design())
mcols(dds) <- cbind(mcols(dds), data.frame(gene_id = rownames(dds), SYMBOL = rownames(dds)))

write_csv(data.frame(gene_id = rownames(counts), SYMBOL = rownames(counts), counts, check.names = FALSE),
          snakemake@output$counts)
write_csv(summary_tab, snakemake@output$summary)
write_csv(composition, snakemake@output$composition)
saveRDS(dds, snakemake@output$dds)
write_csv(coldata, snakemake@output$coldata)
message("DESeqDataSet: ", nrow(dds), " genes x ", ncol(dds), " samples")

### 5. Cell-type composition per sample --------------------------------------------------
# <figure name="cell_composition">
lab <- setNames(paste0(samples$label, "\n", samples$sample), samples$sample)
comp_plot <- composition |>
  mutate(sample = factor(sample, levels = coldata$sample_id),
         cell_type = ifelse(excluded, paste0(cell_type, " (excluded)"), cell_type))
p <- ggplot(comp_plot, aes(x = sample, y = fraction, fill = cell_type)) +
  geom_col(colour = "white", linewidth = 0.2) +
  scale_x_discrete(labels = lab) +
  scale_y_continuous(labels = scales::percent) +
  labs(x = NULL, y = "Share of cells", fill = "Cell type (authors' annotation)",
       title = "Cell-type composition per sample") +
  theme_bw(base_size = 11) +
  theme(panel.grid.major.x = element_blank())
ggsave(snakemake@output$composition_plot, p, width = 9, height = 5, dpi = 300)
# </figure>
