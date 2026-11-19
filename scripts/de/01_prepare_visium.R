# ============================================================
# SCRIPT 01 (Visium projects): Space Ranger outputs -> region pseudobulk -> DESeqDataSet
# Writes the same outputs as the other 01 scripts: data/dds_raw.rds (mcols gene_id, SYMBOL) and
# data/colData.csv, so 02 and 03 do not depend on the data source.
#
# Each Visium section covers several brain regions, at a different level in each sample, so the counts of
# a whole section are not comparable between samples. Only the spots of one region (config visium.region,
# e.g. hippocampus) are used:
#   1. per spot, log-normalized expression (log1p of counts per 10,000) of the region markers
#      (visium.markers), z-scored across the spots of the section; score = mean z of the markers;
#   2. spots with score > visium.score_threshold, grouped into spatially connected components on the
#      Visium hexagonal grid; components with >= visium.min_component spots are kept (a coronal section can
#      cut the region in two pieces), isolated spots are dropped;
#   3. the selection grows by visium.dilate rings of neighbouring spots (the region's layers with fewer
#      marker transcripts, e.g. stratum radiatum).
# The UMI counts of the selected spots are summed per sample (pseudobulk) and analysed as bulk RNA-seq.
# Raw files are read from $DATA/<input.raw_dir> (Space Ranger: *_matrix.mtx.gz, *_features.tsv.gz,
# *_barcodes.tsv.gz, *_tissue_positions_list.csv.gz, *_scalefactors_json.json.gz,
# *_tissue_lowres_image.png.gz). Genes are Ensembl IDs (gene_id_type: ENSEMBL).
# Also writes results/region_spots.csv (every spot with its score and selection), results/
# pseudobulk_summary.csv and plots/region_spots.png (selected spots on each section).
# ============================================================
source(file.path(snakemake@scriptdir, "00_common.R"))
suppressPackageStartupMessages({
  library(DESeq2)
  library(Matrix)
})
icfg <- cfg$input
vcfg <- cfg$visium

data_root <- Sys.getenv("DATA")
if (!nzchar(data_root)) stop("Set the DATA environment variable (heavy files live in $DATA/<dataset>/)")
raw_dir <- file.path(data_root, icfg$raw_dir)
stopifnot(dir.exists(raw_dir))
raw_file <- function(s, suffix) file.path(raw_dir, paste0(s, "_", suffix))

### 1. Samples -------------------------------------------------------------------
# File prefixes look like GSM9608398_Ctrl1: GEO accession + sample name.
prefixes <- sort(unique(sub("_matrix\\.mtx\\.gz$", "", list.files(raw_dir, "_matrix\\.mtx\\.gz$"))))
names_ <- sub("^GSM[0-9]+_", "", prefixes)
patterns <- unlist(cfg$geo$sample_patterns)
condition <- vapply(names_, function(n) {
  hit <- names(patterns)[vapply(patterns, grepl, logical(1), x = n)]
  if (length(hit)) hit[1] else NA_character_
}, character(1), USE.NAMES = FALSE)
if (anyNA(condition)) stop("Samples without a match in geo.sample_patterns: ",
                           paste(names_[is.na(condition)], collapse = ", "))
samples <- tibble(GEO_ID = sub("_.*", "", prefixes), sample = names_, condition = condition,
                  prefix = prefixes)
excluded <- unlist(cfg$geo$exclude_conditions)
if (length(excluded)) {
  message("Samples excluded by condition (", paste(excluded, collapse = ", "), "): ",
          paste(samples$sample[samples$condition %in% excluded], collapse = ", "))
  samples <- samples[!samples$condition %in% excluded, ]
}

### 2. Region spots and pseudobulk per sample ------------------------------------
# Neighbours on the Visium grid (array_row, array_col): same row +-2 columns, adjacent rows +-1 column.
NB <- rbind(c(0, 2), c(0, -2), c(1, 1), c(1, -1), c(-1, 1), c(-1, -1))
spot_key <- function(r, c) paste(r, c)
neighbours <- function(keys) {
  rc <- do.call(rbind, lapply(strsplit(keys, " "), as.integer))
  unique(unlist(lapply(seq_len(nrow(NB)), function(i) spot_key(rc[, 1] + NB[i, 1], rc[, 2] + NB[i, 2]))))
}
components <- function(keys) {
  left <- keys
  comps <- list()
  while (length(left)) {
    comp <- left[1]
    frontier <- comp
    left <- left[-1]
    while (length(frontier)) {
      nb <- intersect(neighbours(frontier), left)
      comp <- c(comp, nb)
      left <- setdiff(left, nb)
      frontier <- nb
    }
    comps[[length(comps) + 1]] <- comp
  }
  comps
}

features <- NULL
pseudo <- list()
spots <- list()
summary_rows <- list()
for (i in seq_len(nrow(samples))) {
  s <- samples$prefix[i]
  m <- readMM(raw_file(s, "matrix.mtx.gz"))
  feat <- read.delim(raw_file(s, "features.tsv.gz"), header = FALSE,
                     col.names = c("gene_id", "SYMBOL", "type"))
  if (is.null(features)) features <- feat
  stopifnot(identical(feat$gene_id, features$gene_id))
  bc <- readLines(raw_file(s, "barcodes.tsv.gz"))
  pos <- read.csv(raw_file(s, "tissue_positions_list.csv.gz"), header = FALSE,
                  col.names = c("barcode", "in_tissue", "array_row", "array_col", "pxl_row", "pxl_col"))
  pos <- pos[match(bc, pos$barcode), ]
  stopifnot(!anyNA(pos$barcode), all(pos$in_tissue == 1))
  sf <- jsonlite::fromJSON(raw_file(s, "scalefactors_json.json.gz"))

  umi <- colSums(m)
  mk <- match(unlist(vcfg$markers), feat$SYMBOL)
  if (anyNA(mk)) stop("Markers not in the features: ", paste(unlist(vcfg$markers)[is.na(mk)], collapse = ", "))
  logn <- log1p(t(t(as.matrix(m[mk, , drop = FALSE])) / umi * 1e4))
  z <- t(scale(t(logn)))
  score <- colMeans(z)

  keys <- spot_key(pos$array_row, pos$array_col)
  comps <- components(keys[score > vcfg$score_threshold])
  sizes <- vapply(comps, length, integer(1))
  sel <- unlist(comps[sizes >= vcfg$min_component])
  for (k in seq_len(vcfg$dilate)) sel <- union(sel, intersect(neighbours(sel), keys))
  selected <- keys %in% sel

  pseudo[[samples$sample[i]]] <- rowSums(m[, selected, drop = FALSE])
  spots[[i]] <- tibble(sample = samples$sample[i], barcode = bc, array_row = pos$array_row,
                       array_col = pos$array_col,
                       x_lowres = pos$pxl_col * sf$tissue_lowres_scalef,
                       y_lowres = pos$pxl_row * sf$tissue_lowres_scalef,
                       umi = umi, score = score, selected = selected)
  summary_rows[[i]] <- tibble(sample = samples$sample[i], GEO_ID = samples$GEO_ID[i],
                              condition = samples$condition[i], spots_in_tissue = length(bc),
                              components_kept = sum(sizes >= vcfg$min_component),
                              spots_selected = sum(selected), umi_selected = sum(umi[selected]),
                              median_umi_per_spot = median(umi[selected]))
  message(samples$sample[i], ": ", length(bc), " spots, components ", paste(sort(sizes, TRUE)[1:min(3, length(sizes))], collapse = "/"),
          " -> ", sum(selected), " selected spots, ", sum(umi[selected]), " UMIs")
}

counts <- do.call(cbind, pseudo)
rownames(counts) <- features$gene_id
counts <- counts[features$type == "Gene Expression", , drop = FALSE]
storage.mode(counts) <- "integer"
symbols <- features$SYMBOL[match(rownames(counts), features$gene_id)]

### 3. DESeqDataSet ----------------------------------------------------------------
samples$label <- samples$sample
coldata <- build_coldata(samples, colnames(counts))
counts <- counts[, coldata$sample_id, drop = FALSE]

dds <- DESeqDataSetFromMatrix(countData = counts, colData = coldata, design = build_design())
mcols(dds) <- cbind(mcols(dds), data.frame(gene_id = rownames(dds), SYMBOL = symbols))

spots_tab <- bind_rows(spots)
write_csv(data.frame(gene_id = rownames(counts), SYMBOL = symbols, counts, check.names = FALSE),
          snakemake@output$counts)
write_csv(spots_tab, snakemake@output$spots)
write_csv(bind_rows(summary_rows), snakemake@output$summary)
saveRDS(dds, snakemake@output$dds)
write_csv(coldata, snakemake@output$coldata)
message("DESeqDataSet: ", nrow(dds), " genes x ", ncol(dds), " samples")

### 4. Selected spots on each section ------------------------------------------------
read_gz_png <- function(path) {
  con <- gzfile(path, "rb")
  on.exit(close(con))
  png::readPNG(readBin(con, "raw", 5e7))
}

# <figure name="region_spots">
n <- nrow(samples)
png(snakemake@output$spots_plot, width = 4 * n, height = 4.6, units = "in", res = 150)
op <- par(mfrow = c(1, n), mar = c(0.5, 0.5, 2.5, 0.5))
for (i in seq_len(n)) {
  img <- read_gz_png(raw_file(samples$prefix[i], "tissue_lowres_image.png.gz"))
  st <- spots_tab[spots_tab$sample == samples$sample[i], ]
  plot(NA, xlim = c(0, dim(img)[2]), ylim = c(dim(img)[1], 0), asp = 1, axes = FALSE, xlab = "", ylab = "",
       main = sprintf("%s (%s): %d spots", samples$sample[i], samples$condition[i], sum(st$selected)))
  rasterImage(img, 0, dim(img)[1], dim(img)[2], 0)
  points(st$x_lowres[!st$selected], st$y_lowres[!st$selected], pch = 16, cex = 0.45, col = "#4dabf799")
  points(st$x_lowres[st$selected], st$y_lowres[st$selected], pch = 16, cex = 0.7, col = "#e8590c")
}
par(op)
dev.off()
# </figure>
