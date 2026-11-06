# ============================================================
# SCRIPT 02: RNA-SEQ QUALITY CONTROL AND DIFFERENTIAL ANALYSIS
# ============================================================
source("scripts/00_packages.R")

message("Importing count data and colData")
colData <- as.data.frame(read_csv("data/colData.csv", show_col_types = FALSE))
colData$group <- factor(colData$group, levels = c("VEH", "CHEM"))

counts <- as.data.frame(read_csv("data/counts.csv", show_col_types = FALSE))
rownames(counts) <- counts$SYMBOL
counts$SYMBOL <- NULL

colnames(counts) <- colData$GEO_ID

stopifnot(all(counts == round(counts)), !any(counts < 0), !anyNA(counts))

### 1. Create DDS object and Top Table --------------------------
# DESeq2 object: load from disk if it already exists, otherwise create and save it
message("Creating or loading DESeq2 object \n")
dds_file <- file.path("data", "dds.rds")

if (file.exists(dds_file)) {
  message("Loading existing DESeq2 object: ", dds_file)
  dds <- readRDS(dds_file)
} else {
  message("Creating DESeq2 object (this may take a while)...")
  dds <- DESeqDataSetFromMatrix(countData = counts, colData = colData, design = ~ group)
  keep <- rowSums(counts(dds) >= 10) >= 3
  dds  <- dds[keep, ]
  dds  <- DESeq(dds, fitType = "local")
}
saveRDS(dds, "data/dds.rds")
resultsNames(dds)

tt <- data.frame(results(dds, contrast = c("group", "CHEM", "VEH"))) # Top table
tt$SYMBOL <- rownames(tt)
tt_filtered <- subset(tt, padj < 0.05 & abs(log2FoldChange) > 0.5)

message("Saving top table")
write.csv(tt, "results/top_table.csv", row.names = F)  # Export top table


### 1. Create DEG and QC plots ---------------------------------------
message("Creating QC and DEG plots")

# PCA
vsd <- vst(dds, blind = TRUE)
pca_data   <- plotPCA(vsd, intgroup = "group", returnData = TRUE)
percentVar <- round(100 * attr(pca_data, "percentVar"))

pca <- ggplot(pca_data, aes(x = PC1, y = PC2, color = group)) +
  geom_point(size = 4) +
  scale_color_manual(values = c("#3399FF", "#CC0000", "#66A61E")) +
  labs(x = paste0("PC1: ", percentVar[1], "% variance"),
       y = paste0("PC2: ", percentVar[2], "% variance"),
       color = "Group") +
  coord_fixed() +
  theme_bw(base_size = 12)

ggsave("plots/pca_vst.png", plot = pca,
       width = 6, height = 5, dpi = 300)

# Volcano plot
keyvals <- ifelse(tt$padj < 0.05 & tt$log2FoldChange >  0.5, "#CC0000",
                  ifelse(tt$padj < 0.05 & tt$log2FoldChange < -0.5, "#3399FF", "gray"))
keyvals[is.na(keyvals)] <- "gray"
names(keyvals) <- ifelse(keyvals == "#CC0000", "Up",
                         ifelse(keyvals == "#3399FF", "Down", "NS"))

vp <- EnhancedVolcano(tt,
                     lab = tt$SYMBOL,
                     x = 'log2FoldChange',
                     y = 'padj',
                     xlab = bquote(~Log[2]~ 'fold change'),
                     ylab = bquote(~-Log[10]~ 'FDR'),
                     title = "Chemotherapy-induced gene expression change",
                     caption = paste("Total genes:", nrow(tt)),
                     pCutoff = 0.05,
                     FCcutoff = 0.5,
                     pointSize = 1.0,
                     labSize = 3,
                     xlim = c(-5, 5),
                     #ylim = c(0, 10),
                     colCustom = keyvals,    
                     colAlpha = 0.5,
                     legendLabSize = 7,
                     legendIconSize = 7,
                     drawConnectors = TRUE,
                     widthConnectors = 0.3,
                     max.overlaps = 15)

ggsave("plots/volcano.png", plot = vp, width = 8, height = 7,
       units = "in", dpi = 300)

# Heatmap of DEGs (VST-transformed counts, row-scaled)
deg_genes <- tt_filtered$SYMBOL[tt_filtered$SYMBOL %in% rownames(vsd)]
deg_mat <- assay(vsd)[deg_genes, , drop = FALSE]
heat_legend <- data.frame(Group = colData(dds)$group,
                          row.names = colnames(dds))
heat_colors <- list(Group = c(VEH = "#3399FF", CHEM = "#CC0000"))

hm <- pheatmap(deg_mat,
               annotation_col    = heat_legend,
               annotation_colors = heat_colors,
               show_rownames = TRUE,
               cluster_rows  = TRUE,
               cluster_cols  = TRUE,
               scale         = "row",
               fontsize      = 8,
               cellwidth     = 10,
               cellheight    = 12)

ggsave("plots/heatmap.png", plot = hm, 
       width = 4, height = 12, units = "in", dpi = 300, bg = "white")

# FDR histogram 
histo <- ggplot(tt, aes(x = padj)) +
  geom_histogram(binwidth = 0.05, fill = "skyblue", color = "black",
                 boundary = 0, na.rm = TRUE) +
  labs(title = "FDR distribution", x = "padj", y = "Number of genes") +
  theme_classic()

ggsave("plots/fdr_histogram.png", plot = histo,
       width = 8, height = 5, units = "in", dpi = 300)

# Box plot of normalized counts, colored by group
norm_counts <- counts(dds, normalized = TRUE)

meta <- data.frame(sample = colnames(dds),
                   group  = colData(dds)$group)

df_long <- as.data.frame(norm_counts) %>%
  rownames_to_column("SYMBOL") %>%
  pivot_longer(-SYMBOL, names_to = "sample", values_to = "expression") %>%
  left_join(meta, by = "sample") %>%
  mutate(sample = factor(sample, levels = meta$sample[order(meta$group)]))

bp <- ggplot(df_long, aes(x = sample, y = log2(expression + 1), fill = group)) +
  geom_boxplot(width = 0.6, outlier.shape = NA) +
  scale_fill_manual(values = c("VEH" = "#3399FF", "CHEM" = "#CC0000")) +
  labs(title = "Normalized expression",
       x = NULL, y = "log2(normalized counts + 1)", fill = "Group") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

ggsave("plots/boxplot.png", plot = bp,
       width = 8, height = 5, units = "in", dpi = 300)

# QQ-Plot
png("plots/QQPlot.png", 
    width = 6, height = 6, units = "in", res = 300)
qqman::qq(tt$pvalue)
dev.off()

# MA plot
resultsNames(dds)
res_shr <- lfcShrink(dds, coef = "group_CHEM_vs_VEH", type = "apeglm")

ma_df <- as.data.frame(res_shr)
ma_df <- ma_df[ma_df$baseMean > 0, ]

ma_df$sig <- ifelse(is.na(ma_df$padj) | ma_df$padj >= 0.05, "ns",
                    ifelse(ma_df$log2FoldChange > 0, "Up", "Down"))
ma_df$sig <- factor(ma_df$sig, levels = c("Up", "Down", "ns"))

p_ma <- ggplot(ma_df, aes(x = baseMean, y = log2FoldChange, color = sig)) +
  geom_point(size = 0.8, alpha = 0.6, na.rm = TRUE) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  scale_x_log10() +
  scale_color_manual(values = c("Up" = "#CC0000", "Down" = "#0066CC", "ns" = "grey70")) +
  labs(x = "Mean of normalized counts (log10)",
       y = "log2 fold change",
       color = NULL,
       title = "MA plot") +
  coord_cartesian(ylim = c(-3, 3)) +
  theme_classic(base_size = 12)

ggsave("plots/MAplot.png", plot = p_ma,
       width = 6, height = 5, dpi = 300)

# Dispersion plot
png("plots/dispersion_plot.png",
    width = 6, height = 5, units = "in", res = 300)
plotDispEsts(dds, main = "Dispersion estimates")
dev.off()

# Normalized counts density
norm_counts <- counts(dds, normalized = TRUE)

df_dens <- as.data.frame(log2(norm_counts + 1)) %>%
  pivot_longer(everything(), names_to = "amostra", values_to = "log2_expr")

df_dens <- as.data.frame(log2(norm_counts + 1)) %>%
  pivot_longer(everything(), names_to = "amostra", values_to = "log2_expr")
meta <- data.frame(amostra = colnames(dds), group = colData(dds)$group)
df_dens <- left_join(df_dens, meta, by = "amostra")

p_dens <- ggplot(df_dens, aes(x = log2_expr, color = amostra, linetype = group)) +
  geom_density(linewidth = 0.7) +
  labs(x = "log2(normalized counts + 1)", y = "Density",
       title = "Density of normalized counts") +
  theme_classic(base_size = 12)

ggsave("plots/density_counts.png", plot = p_dens,
       width = 7, height = 5, dpi = 300)

# Cooks Distance
png("plots/cooks_distance.png",
    width = 7, height = 5, units = "in", res = 300)
par(mar = c(8, 4.5, 3, 1))
boxplot(log10(assays(dds)[["cooks"]]), range = 0, las = 2,
        ylab = "log10(Cook's distance)", main = "Cook's distance")
dev.off()

# ============================================================
# SCRIPT 01: DATA ACQUISITION AND ORGANIZATION
# Dataset: GSE178812
# ============================================================
source("scripts/00_packages.R")

### 1. Configure used files -----------------------------------
matrix_file <- file.path("data/GSE178812_Doxo-counts.xlsx")

matrix_url <- "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE178nnn/GSE178812/suppl/GSE178812%5FDoxo%2Dcounts%2Exlsx"

if (!file.exists(matrix_file)) {
  
  message("Downloading GEO supplementary expression matrix...")
  
  download.file(
    url = matrix_url,
    destfile = matrix_file,
    mode = "wb"
  )
  
} else {
  
  message("Expression matrix already exists. Skipping download.")
  
}

if (!file.exists(matrix_file)) {
  stop("Expression matrix file was not downloaded successfully.")
}

counts_raw <- read_excel( # Import count matrix
  file.path("data/GSE178812_Doxo-counts.xlsx"))
colnames(counts_raw)[1] <- c("SYMBOL")

### 2. Retrieve GEO metadata -------------------------------------

local_file <- file.path(project_dir, paste0(geo_accession, "_series_matrix.txt.gz"))

# Use the local file if it exists; otherwise download from GEO
if (file.exists(local_file)) {
  message("Local file found: ", local_file)
  gse <- getGEO(filename = local_file, GSEMatrix = TRUE)
} else {
  message("Local file not found. Downloading from GEO: ", geo_accession)
  gse <- getGEO(geo_accession, GSEMatrix = TRUE, destdir = project_dir)
  if (is.list(gse)) gse <- gse[[1]]   # online retrieval returns a list
}

# Extract sample metadata according to the object class
if (inherits(gse, "SummarizedExperiment")) {
  colnames(rowData(gse)) <- make.names(colnames(rowData(gse)))
  info <- as.data.frame(colData(gse))
} else {
  fvarLabels(gse) <- make.names(fvarLabels(gse))
  info <- pData(gse)
}

# Remove non-used group
info <- info[!grepl("DMQ", info$title), ]

colData <- data.frame(
  GEO_ID = rownames(info),
  condition = sub(":.*", "", info$title))


### 3. Define experimental groups -------------------------------
group <- ifelse(
  grepl("Doxo-treated", info$title),
  "CHEM",
  ifelse(
    grepl("Sham", info$title),
    "VEH",
    NA_character_
  )) 

colData$group <- group
colData$group <- factor(colData$group, levels = c("VEH", "CHEM"))

counts_raw <- counts_raw %>% dplyr::select(SYMBOL, all_of(colData$condition))

### 4. Inspect gene identifiers ---------------------------------
# Missing identifiers.
missing_symbols <- sum(is.na(counts_raw$SYMBOL) | counts_raw$SYMBOL == "")

# Duplicated gene symbols.
duplicated_symbol_values <- length(unique(
  counts_raw$SYMBOL[
    duplicated(counts_raw$SYMBOL)
  ]
))

message("Missing symbols:", missing_symbols, "|" ,"Duplicated symbols:", duplicated_symbol_values, "\n")

counts <- counts_raw[!duplicated(counts_raw$SYMBOL),]
counts <- counts[!is.na(counts$SYMBOL),]
counts <- as.data.frame(counts)

message("Missing and duplicated rows removed \n")


# 5. Export organized data ------------------------------------
message("Exporting counts and colData \n")
write.csv(counts, "data/counts.csv", row.names = F)
write.csv(colData, "data/colData.csv", row.names = F)

gc()

# ============================================================

