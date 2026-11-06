# Install BiocManager if necessary
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

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

# Install only what is missing
installed <- rownames(installed.packages())

missing_cran <- setdiff(cran_packages, installed)
  if (length(missing_cran) > 0) install.packages(missing_cran)

missing_bioc <- setdiff(bioc_packages, installed)
  if (length(missing_bioc) > 0) BiocManager::install(missing_bioc, update = FALSE, ask = FALSE)



# Load packages (each one only once)
invisible(lapply(c(cran_packages, bioc_packages),
                 library, character.only = TRUE))

rm(list = ls())

# Project configuration
# 1. Project configuration -----------------------------------
project_dir <- path.expand("~/github/bulk-rna-seq")
setwd(project_dir)

dir.create("data", showWarnings = FALSE)
dir.create("plots", showWarnings = FALSE)
dir.create("results", showWarnings = FALSE)

geo_accession <- "GSE178812"
options(timeout = 600)

options(device = "RStudioGD")

# Rscript has no default graphics device; open a null one so print(<ggplot>)
# works under make. ggsave()/png()/pdf() calls are unaffected.
if (!interactive()) pdf(NULL)

writeLines(capture.output(sessionInfo()), "sessionInfo.txt")

