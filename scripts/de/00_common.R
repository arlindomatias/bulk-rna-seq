# ============================================================
# SCRIPT 00: shared setup (sourced by every step)
# Shared by every DE project (bulk_rna_seq_tximport, GSE178812); all dataset-specific
# values come from the project's config.yaml (`snakemake@config`).
# Scripts run through the Snakemake `script:` directive, which creates the `snakemake` object.
# ============================================================
if (!exists("snakemake")) {
  stop("Run this script through Snakemake (`snakemake --use-conda --cores 1`).")
}

# R messages and output go to the rule log.
log_con <- file(snakemake@log[[1]], open = "wt")
sink(log_con)
sink(log_con, type = "message")

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(ggplot2)
})

cfg <- snakemake@config

for (d in c("results", "plots")) dir.create(d, showWarnings = FALSE)

options(timeout = 600, lifecycle_verbosity = "quiet")
set.seed(cfg$seed %||% 123)
if (!interactive()) pdf(NULL)

ref_level <- cfg$groups$reference
trt_level <- cfg$groups$treatment
group_colors <- unlist(cfg$groups$colors)

# Gene identifier type of the count table: ENSEMBL (tximport) or SYMBOL (GEO symbol tables).
gene_id_type <- cfg$gene_id_type %||% "ENSEMBL"
stopifnot(gene_id_type %in% c("ENSEMBL", "SYMBOL"))

# DEG criterion: padj < de$padj; when no gene reaches it, the unadjusted p-value with the stricter
# cutoff de$pvalue_fallback (default 0.01) is used instead (exploratory; report.py applies the same rule).
# Returns the column and its cutoff.
deg_criterion <- function(tt) {
  if (any(tt$padj < cfg$de$padj, na.rm = TRUE)) {
    list(col = "padj", cut = cfg$de$padj)
  } else {
    list(col = "pvalue", cut = cfg$de$pvalue_fallback %||% 0.01)
  }
}

# Readable gene label: symbol, or the Ensembl ID when there is no symbol.
gene_label <- function(symbol, gene_id) {
  ifelse(is.na(symbol) | symbol == "", gene_id, symbol)
}

# Sample metadata -> colData used by DESeq2 (shared by the prepare scripts).
# samples: one row per sample; ids: sample ids present in the count data.
build_coldata <- function(samples, ids) {
  id_col    <- cfg$samples$id_column
  cond_col  <- cfg$samples$condition_column
  label_col <- cfg$samples$label_column %||% id_col
  labels    <- unlist(cfg$groups$labels)

  missing_in_counts <- setdiff(samples[[id_col]], ids)
  missing_in_table  <- setdiff(ids, samples[[id_col]])
  if (length(missing_in_counts) || length(missing_in_table)) {
    stop("Samples differ between the sample table and the counts. Only in table: ",
         paste(missing_in_counts, collapse = ", "), " | only in counts: ",
         paste(missing_in_table, collapse = ", "))
  }

  unknown <- setdiff(unique(samples[[cond_col]]), names(labels))
  if (length(unknown)) stop("Conditions without a label in config groups$labels: ",
                            paste(unknown, collapse = ", "))
  bad_labels <- labels[make.names(labels) != labels]
  if (length(bad_labels)) stop("Group labels must be syntactic R names: ",
                               paste(bad_labels, collapse = ", "))
  stopifnot(ref_level %in% labels, trt_level %in% labels)

  exclude <- unlist(cfg$exclude_samples)
  if (length(exclude)) {
    stopifnot(all(exclude %in% samples[[id_col]]))
    message("Excluding samples: ", paste(exclude, collapse = ", "))
  }

  coldata <- samples |>
    filter(!.data[[id_col]] %in% exclude) |>
    dplyr::select(!matches("fastq|md5|ftp|url", ignore.case = TRUE)) |>
    mutate(
      sample_id = .data[[id_col]],
      label     = as.character(.data[[label_col]]),
      group     = factor(labels[.data[[cond_col]]], levels = c(ref_level, trt_level))
    ) |>
    filter(!is.na(group)) |>          # conditions outside the contrast
    arrange(group, sample_id) |>
    as.data.frame()
  rownames(coldata) <- coldata$sample_id

  covariates <- unlist(cfg$covariates)
  stopifnot(all(covariates %in% names(coldata)))
  for (cv in covariates) if (is.character(coldata[[cv]])) coldata[[cv]] <- factor(coldata[[cv]])

  message("Samples per group:")
  print(table(coldata$group))
  coldata
}

# Design formula: covariates from the config, then the group.
build_design <- function() {
  design <- as.formula(paste("~", paste(c(unlist(cfg$covariates), "group"), collapse = " + ")))
  message("Design: ", deparse(design))
  design
}
