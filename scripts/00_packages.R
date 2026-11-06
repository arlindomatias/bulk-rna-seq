# ============================================================
# SCRIPT 00: PACKAGE CHECK AND PROJECT CONFIGURATION
# ============================================================
### 1. Required packages
required_packages <- c(
  # CRAN
  "readr",
  "dplyr",
  "tidyr",
  "ggplot2",
  "gridExtra",
  "pheatmap",
  "ggrepel",
  "readxl",
  "tibble",
  "ggforce",
  "shadowtext",
  "stringr",
  "qqman",
  "patchwork",
  "ggtangle",
  "igraph",
  "ggridges",
  
  # Bioconductor
  "GEOquery",
  "DESeq2",
  "clusterProfiler",
  "org.Mm.eg.db",
  "EnhancedVolcano",
  "enrichplot",
  "pathview",
  "apeglm",
  "SummarizedExperiment",
  "STRINGdb",
  "ReactomePA"
)

# Check missing packages
missing_packages <- required_packages[
  !vapply(
    required_packages,
    requireNamespace,
    logical(1),
    quietly = TRUE
  )
]

if (length(missing_packages) > 0) {
  
  stop(
    "Missing packages: ",
    paste(missing_packages, collapse = ", "),
    "\nPlease recreate the Conda environment using 'conda env create -f environment.yml'",
  )
  
}

# Load packages
invisible(
  lapply(
    required_packages,
    library,
    character.only = TRUE
  )
)

### 2. Project configuration
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

if (!interactive()) pdf(NULL)