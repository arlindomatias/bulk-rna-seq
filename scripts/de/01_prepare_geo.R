# ============================================================
# SCRIPT 01 (GEO count-table projects): GEO count table + series matrix -> DESeqDataSet
# Writes the same outputs as 01_prepare_tximport.R: data/dds_raw.rds (mcols gene_id, SYMBOL) and
# data/colData.csv, so 02 and 03 do not depend on the data source. Genes are identified by the
# symbols of the deposited table (gene_id = SYMBOL, config gene_id_type: SYMBOL).
# ============================================================
source(file.path(snakemake@scriptdir, "00_common.R"))
suppressPackageStartupMessages({
  library(DESeq2)
  library(GEOquery)
  library(readxl)
})
icfg <- cfg$input
gcfg <- cfg$geo

### 1. Count table deposited in GEO ------------------------------------------
counts_file <- icfg$counts_file
if (!file.exists(counts_file)) {
  message("Downloading the GEO count table: ", icfg$counts_url)
  dir.create(dirname(counts_file), recursive = TRUE, showWarnings = FALSE)
  download.file(icfg$counts_url, destfile = counts_file, mode = "wb")
}
counts_raw <- as.data.frame(read_excel(counts_file))
colnames(counts_raw)[1] <- "SYMBOL"

### 2. Sample metadata from the series matrix -------------------------------
series_file <- icfg$series_matrix
if (file.exists(series_file)) {
  gse <- getGEO(filename = series_file, GSEMatrix = TRUE)
} else {
  message("Downloading the series matrix of ", icfg$geo_accession)
  gse <- getGEO(icfg$geo_accession, GSEMatrix = TRUE, destdir = dirname(series_file))
  if (is.list(gse)) gse <- gse[[1]]   # online retrieval returns a list
}
info <- if (inherits(gse, "SummarizedExperiment")) as.data.frame(colData(gse)) else pData(gse)

# Group of each sample from its GEO title (first matching pattern); excluded groups dropped.
if (!is.null(gcfg$exclude_title_pattern)) {
  dropped <- grepl(gcfg$exclude_title_pattern, info$title)
  message("Samples excluded by title pattern '", gcfg$exclude_title_pattern, "': ", sum(dropped))
  info <- info[!dropped, ]
}
patterns <- unlist(gcfg$title_patterns)
group_raw <- vapply(info$title, function(t) {
  hit <- names(patterns)[vapply(patterns, grepl, logical(1), x = t)]
  if (length(hit)) hit[1] else NA_character_
}, character(1), USE.NAMES = FALSE)

samples <- data.frame(
  GEO_ID    = rownames(info),
  condition = sub(":.*", "", info$title),   # column name of the sample in the count table
  group     = group_raw,
  title     = info$title
)
samples <- samples[!is.na(samples$group), ]

### 3. Gene identifiers --------------------------------------------------------
missing_symbols <- sum(is.na(counts_raw$SYMBOL) | counts_raw$SYMBOL == "")
duplicated_symbols <- length(unique(counts_raw$SYMBOL[duplicated(counts_raw$SYMBOL)]))
message("Missing symbols: ", missing_symbols, " | duplicated symbols: ", duplicated_symbols,
        " (first occurrence kept)")
counts_raw <- counts_raw[!duplicated(counts_raw$SYMBOL) & !is.na(counts_raw$SYMBOL) &
                           counts_raw$SYMBOL != "", ]

stopifnot(all(samples$condition %in% colnames(counts_raw)))
counts <- as.matrix(counts_raw[, samples$condition])
rownames(counts) <- counts_raw$SYMBOL
colnames(counts) <- samples$GEO_ID
stopifnot(all(counts == round(counts)), !any(counts < 0), !anyNA(counts))
storage.mode(counts) <- "integer"

### 4. DESeqDataSet ------------------------------------------------------------
coldata <- build_coldata(samples, colnames(counts))
counts <- counts[, coldata$sample_id, drop = FALSE]

dds <- DESeqDataSetFromMatrix(countData = counts, colData = coldata, design = build_design())
mcols(dds) <- cbind(mcols(dds), data.frame(gene_id = rownames(dds), SYMBOL = rownames(dds)))

write_csv(data.frame(SYMBOL = rownames(counts), counts, check.names = FALSE), snakemake@output$counts)
saveRDS(dds, snakemake@output$dds)
write_csv(coldata, snakemake@output$coldata)
message("DESeqDataSet: ", nrow(dds), " genes x ", ncol(dds), " samples")
