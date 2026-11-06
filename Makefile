# Bulk RNA-seq analysis pipeline
# GSE178812

RSCRIPT = Rscript

SCRIPTS = scripts/00_packages.R \
          scripts/01_download_and_prepare.R \
          scripts/02_differential_expression_and_qc.R \
          scripts/03_enrichment.R

.PHONY: all packages prepare differential enrichment clean

# Execute the complete pipeline
all: packages prepare differential enrichment

# Step 00: Install/load required packages
packages:
	$(RSCRIPT) scripts/00_packages.R

# Step 01: Download and prepare data
prepare: packages
	$(RSCRIPT) scripts/01_download_and_prepare.R

# Step 02: Quality control and differential expression
differential: prepare
	$(RSCRIPT) scripts/02_differential_expression_and_qc.R

# Step 03: Functional enrichment analysis
enrichment: differential
	$(RSCRIPT) scripts/03_enrichment.R

# Remove generated results and plots
clean:
	rm -rf data/colData.csv data/counts.csv data/dds.rds
