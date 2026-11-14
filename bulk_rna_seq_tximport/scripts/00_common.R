# ============================================================
# SCRIPT 00: shared setup (sourced by every step)
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

# Readable gene label: symbol, or the Ensembl ID when there is no symbol.
gene_label <- function(symbol, gene_id) {
  ifelse(is.na(symbol) | symbol == "", gene_id, symbol)
}
