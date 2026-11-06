RSCRIPT := Rscript

# One representative output file per pipeline stage
COUNTS    := data/counts.csv
TOP_TABLE := data/top_table.csv
ENRICH    := plots/gsea_kegg.png

.PHONY: all clean help

all: $(ENRICH)

# Step 1: download and prepare the data
$(COUNTS): scripts/00_packages.R scripts/01_download_and_prepare.R
	$(RSCRIPT) scripts/01_download_and_prepare.R

# Step 2: quality control and differential expression
$(TOP_TABLE): $(COUNTS) data/colData.csv scripts/00_packages.R scripts/02_differential_expression_and_qc.R
	$(RSCRIPT) scripts/02_differential_expression_and_qc.R

# Step 3: functional enrichment
$(ENRICH): $(TOP_TABLE) scripts/00_packages.R scripts/03_enrichment.R
	$(RSCRIPT) scripts/03_enrichment.R

# Remove generated results (scripts are kept)
clean:
	rm -f $(TOP_TABLE) $(ENRICH)
	rm -f plots/*.png plots/*.pdf
	rm -rf plots/pathview

help:
	@echo "make        run the full pipeline (only outdated steps)"
	@echo "make -n     show what would run, without running it"
	@echo "make clean  delete generated results"