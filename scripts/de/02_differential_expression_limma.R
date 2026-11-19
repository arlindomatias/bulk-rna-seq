# ============================================================
# SCRIPT 02 (limma): QC and differential expression from normalized counts
# For GEO series that only deposit normalized (non-integer) counts. limma-trend on
# log2(normalized + 1): lmFit + eBayes(trend = TRUE), so the prior variance follows the mean-variance
# trend of the log values (voom needs raw counts and library sizes, which are not available).
# Writes results/top_table.csv with the same columns 03_enrichment.R and report.py read from the
# DESeq2 version (gene_id, SYMBOL, baseMean, log2FoldChange, stat, pvalue, padj), plus AveExpr and B.
# ============================================================
source(file.path(snakemake@scriptdir, "00_common.R"))
suppressPackageStartupMessages({
  library(SummarizedExperiment)
  library(limma)
  library(pheatmap)
  library(EnhancedVolcano)
  library(qqman)
})
out <- snakemake@output

se <- readRDS(snakemake@input$se)
padj_cut <- cfg$de$padj
lfc_cut  <- cfg$de$lfc
coldata  <- as.data.frame(colData(se))

### 1. Filter and fit --------------------------------------------------------
norm <- assay(se, "normalized")
value_name <- cfg$input$value_name %||% "normalized counts"  # e.g. "FPKM": used in plot labels
min_samples <- cfg$filter$min_samples %||% min(table(coldata$group))
keep <- rowSums(norm >= cfg$filter$min_count) >= min_samples
message("Filter: >= ", cfg$filter$min_count, " ", value_name, " in >= ", min_samples,
        " samples -> ", sum(keep), " of ", nrow(se), " genes")
se   <- se[keep, ]
norm <- norm[keep, , drop = FALSE]
logexpr <- log2(norm + 1)

design_formula <- build_design()
design <- model.matrix(design_formula, data = coldata)
coef_name <- paste0("group", trt_level)
stopifnot(coef_name %in% colnames(design))

# Sample quality weights (arrayWeights): down-weights samples that are less consistent with their group.
if (isTRUE(cfg$de$array_weights)) {
  sample_weights <- arrayWeights(logexpr, design)
  message("Sample weights (arrayWeights): ",
          paste(sprintf("%s=%.2f", colnames(logexpr), sample_weights), collapse = ", "))
} else {
  sample_weights <- setNames(rep(1, ncol(logexpr)), colnames(logexpr))
}
fit <- lmFit(logexpr, design, weights = if (isTRUE(cfg$de$array_weights)) sample_weights else NULL)
fit <- eBayes(fit, trend = isTRUE(cfg$de$trend %||% TRUE), robust = isTRUE(cfg$de$robust %||% TRUE))
saveRDS(fit, out$fit)
message("Coefficients: ", paste(colnames(design), collapse = ", "),
        " | prior df (median): ", signif(median(fit$df.prior), 3))

res <- topTable(fit, coef = coef_name, number = Inf, sort.by = "none")
ann <- as.data.frame(rowData(se)[, c("gene_id", "SYMBOL")])
tt <- data.frame(ann,
                 baseMean       = rowMeans(norm),
                 AveExpr        = res$AveExpr,
                 log2FoldChange = res$logFC,
                 lfcSE          = fit$stdev.unscaled[, coef_name] * sqrt(fit$s2.post),
                 stat           = res$t,
                 pvalue         = res$P.Value,
                 padj           = res$adj.P.Val,
                 B              = res$B,
                 row.names = NULL, check.names = FALSE)
tt$label <- gene_label(tt$SYMBOL, tt$gene_id)
tt <- tt[order(tt$padj, tt$pvalue), ]
deg_rule <- deg_criterion(tt)
sig_col <- deg_rule$col
sig_cut <- deg_rule$cut
if (sig_col == "pvalue") message("No gene with padj < ", padj_cut,
                                 ": DEGs use the unadjusted p-value < ", sig_cut, " (exploratory)")
tt_filtered <- tt[!is.na(tt[[sig_col]]) & tt[[sig_col]] < sig_cut & abs(tt$log2FoldChange) > lfc_cut, ]

write_csv(tt, out$top_table)

de_summary <- tibble(
  contrast   = paste(trt_level, "vs", ref_level),
  design     = deparse(design_formula),
  samples    = ncol(se),
  genes_tested = nrow(tt),
  padj_cutoff  = padj_cut,
  deg_criterion = sig_col,
  deg_cutoff    = sig_cut,
  lfc_cutoff   = lfc_cut,
  deg_up   = sum(tt_filtered$log2FoldChange > 0),
  deg_down = sum(tt_filtered$log2FoldChange < 0),
  padj_only = sum(tt$padj < padj_cut, na.rm = TRUE)
)
write_tsv(de_summary, out$summary)
print(as.data.frame(de_summary))

write_csv(data.frame(gene_id = rownames(logexpr), SYMBOL = rowData(se)$SYMBOL,
                     round(logexpr, 4), check.names = FALSE),
          out$log_expression)

### 2. QC plots ---------------------------------------------------------------
message("Creating QC and DEG plots")
sample_labels <- setNames(coldata$label, coldata$sample_id)

# <figure name="pca">
# PCA of log2(normalized + 1), 500 most variable genes (as DESeq2::plotPCA)
top_var <- head(order(apply(logexpr, 1, var), decreasing = TRUE), 500)
pc <- prcomp(t(logexpr[top_var, ]))
percentVar <- round(100 * pc$sdev^2 / sum(pc$sdev^2))
pca_data <- data.frame(pc$x[, 1:2], group = coldata$group, label = sample_labels[rownames(pc$x)])

pca <- ggplot(pca_data, aes(x = PC1, y = PC2, color = group, label = label)) +
  geom_point(size = 4) +
  ggrepel::geom_text_repel(size = 3, show.legend = FALSE) +
  scale_color_manual(values = group_colors) +
  labs(x = paste0("PC1: ", percentVar[1], "% variance"),
       y = paste0("PC2: ", percentVar[2], "% variance"),
       color = "Group") +
  coord_fixed() +
  theme_bw(base_size = 12)
ggsave(out$pca, plot = pca, width = 6, height = 5, dpi = 300)
# </figure>

# <figure name="sample_distances">
# Sample-to-sample distances (log2 normalized counts)
annot_col  <- data.frame(Group = coldata$group, row.names = coldata$sample_id)
annot_cols <- list(Group = group_colors)
dist_mat <- as.matrix(dist(t(logexpr)))
dimnames(dist_mat) <- list(sample_labels[rownames(dist_mat)], sample_labels[colnames(dist_mat)])
rownames(annot_col) <- sample_labels[rownames(annot_col)]
pheatmap(dist_mat,
         clustering_distance_rows = as.dist(dist_mat),
         clustering_distance_cols = as.dist(dist_mat),
         annotation_col = annot_col, annotation_colors = annot_cols,
         color = colorRampPalette(c("#08306B", "#F7FBFF"))(100),
         main = paste0("Euclidean distance between samples (log2 ", value_name, ")"),
         filename = out$sample_dist, width = 7, height = 6)
# </figure>

# <figure name="marker_genes">
# Marker genes (e.g. sex and cell type markers) per sample
markers <- unlist(cfg$plots$marker_genes)
if (!length(markers)) {
  # The output is declared in every project; without configured markers it only says so.
  p_mk <- ggplot() + annotate("text", x = 0, y = 0, label = "No marker genes configured") + theme_void()
  ggsave(out$markers, plot = p_mk, width = 4, height = 1, dpi = 100, bg = "white")
} else {
  mk_idx <- which(rowData(se)$SYMBOL %in% markers)
  mk_df <- as.data.frame(logexpr[mk_idx, , drop = FALSE])
  mk_df$SYMBOL <- rowData(se)$SYMBOL[mk_idx]
  mk_df <- mk_df |>
    pivot_longer(-SYMBOL, names_to = "sample", values_to = "logexpr") |>
    mutate(label = factor(sample_labels[sample], levels = sample_labels),
           group = coldata$group[match(sample, coldata$sample_id)])
  not_found <- setdiff(markers, mk_df$SYMBOL)
  if (length(not_found)) message("Marker genes absent after filtering: ",
                                 paste(not_found, collapse = ", "))
  p_mk <- ggplot(mk_df, aes(x = label, y = logexpr, fill = group)) +
    geom_col() +
    facet_wrap(~ SYMBOL, scales = "free_y") +
    scale_fill_manual(values = group_colors) +
    labs(x = NULL, y = paste0("log2(", value_name, " + 1)"), fill = "Group",
         title = "Marker genes",
         caption = if (length(not_found)) paste("Filtered out / absent:",
                                                paste(not_found, collapse = ", "))) +
    theme_bw(base_size = 10) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
  ggsave(out$markers, plot = p_mk, width = 9, height = 6, dpi = 300)
}
# </figure>

# <figure name="volcano">
# Volcano plot
keyvals <- ifelse(tt[[sig_col]] < sig_cut & tt$log2FoldChange >  lfc_cut, "#CC0000",
                  ifelse(tt[[sig_col]] < sig_cut & tt$log2FoldChange < -lfc_cut, "#3399FF", "gray"))
keyvals[is.na(keyvals)] <- "gray"
names(keyvals) <- ifelse(keyvals == "#CC0000", "Up",
                         ifelse(keyvals == "#3399FF", "Down", "NS"))

vp <- EnhancedVolcano(tt,
                      lab = tt$label,
                      x = "log2FoldChange",
                      y = sig_col,
                      xlab = bquote(~Log[2]~ "fold change"),
                      ylab = if (sig_col == "padj") bquote(~-Log[10]~ "FDR") else
                        bquote(~-Log[10]~ "p-value (unadjusted)"),
                      title = cfg$title,
                      subtitle = paste(trt_level, "vs", ref_level, "(limma-trend)"),
                      caption = paste("Total genes:", nrow(tt)),
                      pCutoff = sig_cut,
                      FCcutoff = lfc_cut,
                      pointSize = 1.0,
                      labSize = 3,
                      colCustom = keyvals,
                      colAlpha = 0.5,
                      legendLabSize = 7,
                      legendIconSize = 7,
                      drawConnectors = TRUE,
                      widthConnectors = 0.3,
                      max.overlaps = 15)
ggsave(out$volcano, plot = vp, width = 8, height = 7, units = "in", dpi = 300)
# </figure>

# <figure name="heatmap">
# Heatmap: top genes by FDR (DEGs first; capped so the plot stays legible)
n_max <- cfg$plots$heatmap_max_genes
hm_genes <- head(tt_filtered$gene_id, n_max)
hm_title <- paste0("Top ", length(hm_genes), " of ", nrow(tt_filtered), " DEGs",
                   if (sig_col == "pvalue") paste0(" (unadjusted p < ", sig_cut, ")") else "")
if (length(hm_genes) < 2) {
  hm_genes <- head(tt$gene_id[!is.na(tt$padj)], n_max)
  hm_title <- paste0("No DEGs: top ", length(hm_genes), " genes by FDR")
}
hm_mat <- logexpr[hm_genes, , drop = FALSE]
rownames(hm_mat) <- tt$label[match(hm_genes, tt$gene_id)]
colnames(hm_mat) <- sample_labels[colnames(hm_mat)]
pheatmap(hm_mat,
         annotation_col    = annot_col,
         annotation_colors = annot_cols,
         show_rownames = TRUE,
         cluster_rows  = TRUE,
         cluster_cols  = TRUE,
         scale         = "row",
         fontsize      = 8,
         main          = hm_title,
         filename = out$heatmap, width = 7,
         height = min(4 + 0.13 * length(hm_genes), 14))
# </figure>

# <figure name="fdr_histogram">
# FDR histogram
histo <- ggplot(tt, aes(x = padj)) +
  geom_histogram(binwidth = 0.05, fill = "skyblue", color = "black",
                 boundary = 0, na.rm = TRUE) +
  labs(title = "FDR distribution", x = "padj", y = "Number of genes") +
  theme_classic()
ggsave(out$fdr, plot = histo, width = 8, height = 5, units = "in", dpi = 300)
# </figure>

# <figure name="boxplot">
# Box plot of log2 normalized counts, colored by group
meta <- data.frame(sample = coldata$sample_id, group = coldata$group,
                   label = factor(sample_labels, levels = sample_labels))
df_long <- as.data.frame(logexpr) |>
  rownames_to_column("gene_id") |>
  pivot_longer(-gene_id, names_to = "sample", values_to = "log2_expr") |>
  left_join(meta, by = "sample")

bp <- ggplot(df_long, aes(x = label, y = log2_expr, fill = group)) +
  geom_boxplot(width = 0.6, outlier.shape = NA) +
  scale_fill_manual(values = group_colors) +
  labs(title = "Normalized expression",
       x = NULL, y = paste0("log2(", value_name, " + 1)"), fill = "Group") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
ggsave(out$boxplot, plot = bp, width = 8, height = 5, units = "in", dpi = 300)
# </figure>

# <figure name="QQPlot">
# QQ-Plot
png(out$qqplot, width = 6, height = 6, units = "in", res = 300)
qqman::qq(tt$pvalue[!is.na(tt$pvalue)])
dev.off()
# </figure>

# <figure name="MAplot">
# MA plot: average log2 expression vs log2FC
ma_df <- tt
ma_df$sig <- ifelse(is.na(ma_df[[sig_col]]) | ma_df[[sig_col]] >= sig_cut, "ns",
                    ifelse(ma_df$log2FoldChange > 0, "Up", "Down"))
ma_df$sig <- factor(ma_df$sig, levels = c("Up", "Down", "ns"))

p_ma <- ggplot(ma_df, aes(x = AveExpr, y = log2FoldChange, color = sig)) +
  geom_point(size = 0.8, alpha = 0.6, na.rm = TRUE) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  scale_color_manual(values = c("Up" = "#CC0000", "Down" = "#0066CC", "ns" = "grey70")) +
  labs(x = paste0("Average log2(", value_name, " + 1)"),
       y = "log2 fold change",
       color = NULL,
       title = "MA plot") +
  coord_cartesian(ylim = c(-3, 3)) +
  theme_classic(base_size = 12)
ggsave(out$ma, plot = p_ma, width = 6, height = 5, dpi = 300)
# </figure>

# <figure name="mean_variance">
# limma-trend: residual standard deviation vs average expression, with the fitted prior trend
png(out$mean_variance, width = 6, height = 5, units = "in", res = 300)
plotSA(fit, main = "Mean-variance trend (limma-trend)")
dev.off()
# </figure>

# <figure name="sample_weights">
# Sample weights used in the fit (all 1 when de.array_weights is false)
w_df <- data.frame(label = factor(sample_labels[colnames(logexpr)], levels = sample_labels),
                   weight = as.numeric(sample_weights), group = coldata$group)
p_w <- ggplot(w_df, aes(x = label, y = weight, fill = group)) +
  geom_col() +
  geom_hline(yintercept = 1, linetype = "dashed") +
  scale_fill_manual(values = group_colors) +
  labs(x = NULL, y = "Sample weight", fill = "Group",
       title = if (isTRUE(cfg$de$array_weights)) "Sample quality weights (arrayWeights)" else
         "Equal sample weights (arrayWeights off)") +
  theme_classic(base_size = 12) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
ggsave(out$sample_weights, plot = p_w, width = 7, height = 4.5, dpi = 300)
# </figure>

# <figure name="density_counts">
# Normalized counts density
p_dens <- ggplot(df_long, aes(x = log2_expr, color = label, linetype = group)) +
  geom_density(linewidth = 0.7) +
  labs(x = paste0("log2(", value_name, " + 1)"), y = "Density", color = "Sample",
       title = paste0("Density of ", value_name)) +
  theme_classic(base_size = 12)
ggsave(out$density, plot = p_dens, width = 7, height = 5, dpi = 300)
# </figure>
