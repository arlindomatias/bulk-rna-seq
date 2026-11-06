# ============================================================
# SCRIPT 00: Project Configuration and Package Loading
# ============================================================

# CRAN packages
cran_packages <- c(
  "readr", "dplyr", "tidyr", "ggplot2", "gridExtra",
  "pheatmap", "ggrepel", "readxl", "tibble", "ggforce",
  "shadowtext", "stringr", "qqman", "patchwork", "ggtangle",
  "igraph"
)

# Bioconductor packages
bioc_packages <- c(
  "GEOquery", "DESeq2", "clusterProfiler", "org.Mm.eg.db",
  "EnhancedVolcano", "enrichplot", "pathview", "apeglm",
  "SummarizedExperiment", "STRINGdb", "ReactomePA"
)

# Load packages
invisible(lapply(
  c(cran_packages, bioc_packages),
  library,
  character.only = TRUE
))

# Project configuration
if (!file.exists("scripts/00_packages.R")) {
  setwd(path.expand("~/github/bulk-rna-seq"))
}
project_dir <- getwd()

for (d in c("data", "plots", "results", "cache")) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

geo_accession <- "GSE178812"

options(timeout = 600)
options(lifecycle_verbosity = "quiet")

set.seed(123)

# Graphics device for non-interactive execution
if (!interactive()) pdf(NULL)