# ============================================================
# Bulk RNA-seq analysis pipeline
# Dataset: GSE178812
# ============================================================

RSCRIPT = "Rscript"


# ------------------------------------------------------------
# Complete pipeline
# ------------------------------------------------------------

rule all:
    input:
        # Prepared data
        "data/counts.csv",
        "data/colData.csv",

        # Differential expression
        "results/top_table.csv",

        # QC plots
        "plots/pca_vst.png",
        "plots/volcano.png",
        "plots/heatmap.png",
        "plots/fdr_histogram.png",
        "plots/boxplot.png",
        "plots/QQPlot.png",
        "plots/MAplot.png",
        "plots/dispersion_plot.png",
        "plots/density_counts.png",
        "plots/cooks_distance.png",

        # Functional enrichment
        "results/GSEA_GO_all_full.csv",
        "results/GSEA_KEGG_full.csv",
        "results/enrichment_STRING_full.csv",
        "results/edges_PPI.csv",
        "results/nodes_PPI.csv",
        "sessionInfo.txt"


# ------------------------------------------------------------
# Step 00: Packages and configuration
# ------------------------------------------------------------

rule packages:
    output:
        touch(".pipeline/packages.done")

    shell:
        """
        mkdir -p .pipeline
        {RSCRIPT} scripts/00_packages.R
        """


# ------------------------------------------------------------
# Step 01: Data acquisition and preparation
# ------------------------------------------------------------

rule prepare:
    input:
        packages = ".pipeline/packages.done"

    output:
        counts = "data/counts.csv",
        coldata = "data/colData.csv"

    shell:
        """
        {RSCRIPT} scripts/01_download_and_prepare.R
        """


# ------------------------------------------------------------
# Step 02: Quality control and differential expression
# ------------------------------------------------------------

rule differential:
    input:
        counts = "data/counts.csv",
        coldata = "data/colData.csv"

    output:
        top_table = "results/top_table.csv",

        # QC plots
        pca = "plots/pca_vst.png",
        volcano = "plots/volcano.png",
        heatmap = "plots/heatmap.png",
        fdr = "plots/fdr_histogram.png",
        boxplot = "plots/boxplot.png",
        qqplot = "plots/QQPlot.png",
        ma = "plots/MAplot.png",
        dispersion = "plots/dispersion_plot.png",
        density = "plots/density_counts.png",
        cooks = "plots/cooks_distance.png"

    shell:
        """
        {RSCRIPT} scripts/02_differential_expression_and_qc.R
        """


# ------------------------------------------------------------
# Step 03: Functional enrichment analysis
# ------------------------------------------------------------

rule enrichment:
    input:
        top_table = "results/top_table.csv",
        counts = "data/counts.csv",
        coldata = "data/colData.csv"

    output:
        go = "results/GSEA_GO_all_full.csv",
        kegg = "results/GSEA_KEGG_full.csv",
        string = "results/enrichment_STRING_full.csv",
        edges = "results/edges_PPI.csv",
        nodes = "results/nodes_PPI.csv",
        session = "sessionInfo.txt"

    shell:
        """
        {RSCRIPT} scripts/03_enrichment.R
        """