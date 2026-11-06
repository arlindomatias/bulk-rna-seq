# Main RNA-seq pipeline
project_dir <- path.expand("~/github/bulk-rna-seq")
setwd(project_dir)

scripts <- c(
  "scripts/00_packages.R",
  "scripts/01_download_and_prepare.R",
  "scripts/02_differential_expression_and_qc.R",
  "scripts/03_enrichment.R"
)

for (script in scripts) {
  message("\nRunning: ", script)
  
  source(script)
  
  message("Completed: ", script)
}

message("\nRNA-seq pipeline completed successfully.")