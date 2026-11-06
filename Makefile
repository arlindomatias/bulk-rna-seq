# Bulk RNA-seq pipeline (GSE178812)
# Run from the repository root: make

RSCRIPT := Rscript
PKGS    := scripts/00_packages.R

.PHONY: all clean cleanall

all: results/GSEA_KEGG_full.csv results/enrichment_STRING_full.csv

# 01 - download and prepare counts / sample metadata
data/counts.csv data/colData.csv &: scripts/01_download_and_prepare.R $(PKGS)
	@mkdir -p data
	$(RSCRIPT) scripts/01_download_and_prepare.R

# 02 - differential expression and QC
data/dds.rds results/top_table.csv plots/volcano.png &: \
		scripts/02_differential_expression_and_qc.R \
		data/counts.csv data/colData.csv $(PKGS)
	@mkdir -p results plots
	$(RSCRIPT) scripts/02_differential_expression_and_qc.R

# 03 - network and enrichment (STRING, GSEA GO/KEGG)
results/GSEA_KEGG_full.csv results/enrichment_STRING_full.csv &: \
		scripts/03_enrichment.R results/top_table.csv data/dds.rds $(PKGS)
	@mkdir -p results plots
	$(RSCRIPT) scripts/03_enrichment.R

# Remove only intermediate objects (results/ and plots/ are kept)
clean:
	rm -f data/dds.rds

# Remove everything generated, INCLUDING results and plots
cleanall: clean
	rm -rf results plots