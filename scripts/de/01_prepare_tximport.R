# ============================================================
# SCRIPT 01 (tximport projects): tximport + sample table -> DESeqDataSet
# Writes the same outputs as 01_prepare_geo.R: data/dds_raw.rds (mcols gene_id, SYMBOL) and
# data/colData.csv, so 02 and 03 do not depend on the data source.
# ============================================================
source(file.path(snakemake@scriptdir, "00_common.R"))
suppressPackageStartupMessages(library(DESeq2))

### 1. Import --------------------------------------------------------------
txi <- readRDS(snakemake@input$txi)
stopifnot(all(c("counts", "abundance", "length") %in% names(txi)))
if (!identical(txi$countsFromAbundance, "no")) {
  stop("txi must be imported with countsFromAbundance = 'no' for DESeqDataSetFromTximport ",
       "(found: ", txi$countsFromAbundance, ").")
}

samples <- as.data.frame(read_tsv(snakemake@input$samples, show_col_types = FALSE))
tx2gene <- read_tsv(snakemake@input$tx2gene, show_col_types = FALSE)

### 2. Sample metadata ---------------------------------------------------
coldata <- build_coldata(samples, colnames(txi$counts))

### 3. Subset txi to the samples used, in colData order -------------
for (m in c("counts", "abundance", "length")) {
  txi[[m]] <- txi[[m]][, coldata$sample_id, drop = FALSE]
}

### 4. Gene annotation --------------------------------------------------
genes <- unique(tx2gene[, c("gene_id", "gene_name")])
stopifnot(!anyDuplicated(genes$gene_id))

annot <- data.frame(gene_id_version = rownames(txi$counts))
annot$gene_id <- sub("\\.[0-9]+(_PAR_Y)?$", "\\1", annot$gene_id_version)
annot$SYMBOL  <- genes$gene_name[match(annot$gene_id_version, genes$gene_id)]
stopifnot(!anyDuplicated(annot$gene_id))
message("Genes: ", nrow(annot), " | without symbol: ", sum(is.na(annot$SYMBOL)),
        " | symbols shared by >1 gene: ",
        sum(duplicated(annot$SYMBOL) | duplicated(annot$SYMBOL, fromLast = TRUE)))

### 5. DESeqDataSet --------------------------------------------------------
design <- build_design()

dds <- DESeqDataSetFromTximport(txi, colData = coldata, design = design)
rownames(dds) <- annot$gene_id
mcols(dds) <- cbind(mcols(dds), annot)

saveRDS(dds, snakemake@output$dds)
write_csv(coldata, snakemake@output$coldata)
message("DESeqDataSet: ", nrow(dds), " genes x ", ncol(dds), " samples")
