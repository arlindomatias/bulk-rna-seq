# ============================================================
# SCRIPT 01 (GEO normalized-table projects): normalized expression table + sample sheet from GEO
# -> SummarizedExperiment for limma (02_differential_expression_limma.R).
#
# For series that only deposit normalized (non-integer) counts, so DESeq2 cannot be used.
# Writes data/se_raw.rds (assay "normalized", rowData gene_id and SYMBOL), data/expression.csv
# (gene-level table) and data/colData.csv (built with build_coldata()), so 03 does not depend on the
# data source. Rows that share a gene ID (several features per gene in the deposited table) are summed:
# the values are on the same per-sample scale, so the sum is the gene-level normalized count.
# ============================================================
source(file.path(snakemake@scriptdir, "00_common.R"))
suppressPackageStartupMessages(library(SummarizedExperiment))
icfg <- cfg$input

download_once <- function(url, path) {
  if (!file.exists(path)) {
    message("Downloading ", url)
    dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
    download.file(url, destfile = paste0(path, ".part"), mode = "wb")
    file.rename(paste0(path, ".part"), path)
  }
  path
}

### 1. Normalized table and sample sheet -------------------------------------
expr_raw <- read_csv(download_once(icfg$expression_url, icfg$expression_file), show_col_types = FALSE)
samples  <- read_csv(download_once(icfg$samples_url, icfg$samples_file), show_col_types = FALSE)

# Conditions outside the analysis (e.g. a third treatment group) are dropped before build_coldata().
cond_col <- cfg$samples$condition_column
excluded <- unlist(cfg$geo$exclude_conditions)
if (length(excluded)) {
  drop <- samples[[cond_col]] %in% excluded
  message("Samples excluded by condition (", paste(excluded, collapse = ", "), "): ", sum(drop))
  samples <- samples[!drop, ]
}

### 2. Gene-level matrix -------------------------------------------------------
id_col  <- icfg$gene_id_column
sym_col <- icfg$symbol_column
sample_ids <- samples[[cfg$samples$id_column]]
stopifnot(all(sample_ids %in% names(expr_raw)))

vals <- as.matrix(expr_raw[, sample_ids])
stopifnot(!anyNA(vals), !any(vals < 0))
n_rows <- nrow(expr_raw)
gene_ids <- expr_raw[[id_col]]
expr <- rowsum(vals, group = gene_ids, reorder = TRUE)
symbols <- expr_raw[[sym_col]][match(rownames(expr), gene_ids)]
message("Rows in the deposited table: ", n_rows, " -> genes: ", nrow(expr),
        " (", sum(duplicated(gene_ids)), " rows shared a gene ID and were summed)")
message("Library sizes (sum of normalized values): ",
        paste(range(round(colSums(expr))), collapse = " - "))

### 3. SummarizedExperiment -------------------------------------------------------
coldata <- build_coldata(samples, colnames(expr))
expr <- expr[, coldata$sample_id, drop = FALSE]

se <- SummarizedExperiment(
  assays  = list(normalized = expr),
  rowData = DataFrame(gene_id = rownames(expr), SYMBOL = symbols),
  colData = DataFrame(coldata, check.names = FALSE)
)

write_csv(data.frame(gene_id = rownames(expr), SYMBOL = symbols, expr, check.names = FALSE),
          snakemake@output$expression)
saveRDS(se, snakemake@output$se)
write_csv(coldata, snakemake@output$coldata)
message("SummarizedExperiment: ", nrow(se), " genes x ", ncol(se), " samples")
