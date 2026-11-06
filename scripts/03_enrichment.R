# ============================================================
# SCRIPT 03: ENRICHMENT ANALYSIS
# Load organized data from Script 02 if necessary ---------
# ============================================================
source("scripts/00_packages.R")

cat("Importing count data and colData")
colData <- as.data.frame(read_csv("data/colData.csv", show_col_types = FALSE))
colData$group <- factor(colData$group, levels = c("VEH", "CHEM"))

counts <- as.data.frame(read_csv("data/counts.csv", show_col_types = FALSE))
rownames(counts) <- counts$SYMBOL
counts$SYMBOL <- NULL
colnames(counts) <- colData$GEO_ID

tt <- as.data.frame(read_csv("results/top_table.csv", show_col_types = FALSE))
tt_filtered <- subset(tt, padj < 0.05 & abs(log2FoldChange) > 0.5)
tt_filtered_up <- subset(tt, padj < 0.05 & log2FoldChange > 0.5)
tt_filtered_down <- subset(tt, padj < 0.05 & log2FoldChange < - 0.5)

### 1. Prepare count matrix with SYMBOLS -------------------------------------------
counts_raw <- counts
counts_raw$"SYMBOL" <- rownames(counts)
counts_raw <- merge(tt, counts_raw, by = "SYMBOL")

counts_sig <- counts_raw[counts_raw$SYMBOL %in% tt_filtered$SYMBOL, ]

# Prepare ranking list for GSEA
gene_list <- counts_raw$log2FoldChange
names(gene_list) <- counts_raw$SYMBOL

gene_list <- gene_list[!is.na(gene_list) & !is.na(names(gene_list))]
gene_list <- gene_list[!duplicated(names(gene_list))]
gene_list <- sort(gene_list, decreasing = TRUE)

### 2. GSEA Analysis ---------------------------------------------------------------
# GSEA with GO
cat("Performing GSEA")
gse_GO_all <- gseGO(geneList= gene_list, 
                    ont ="ALL", # BP, CC, MF
                    keyType = "SYMBOL", 
                    minGSSize = 3, 
                    maxGSSize = 800, 
                    pvalueCutoff = 0.05, 
                    verbose = TRUE, 
                    OrgDb = "org.Mm.eg.db",
                    pAdjustMethod = "fdr")
gse_GO_all <- pairwise_termsim(gse_GO_all) 

# Dotplot
dot <- dotplot(gse_GO_all,
                        showCategory = 10,
                        split = "ONTOLOGY",
                        font.size = 12) +
  facet_grid(ONTOLOGY ~ .sign, scales = "free_y", space = "free_y")

dot
ggsave("plots/dotplot_GO.png", plot = dot, width = 7, height = 12, dpi = 300)

# Ridgeplot
make_ridge <- function(obj, ont, n = 10) {
  sub <- obj
  sub@result <- obj@result[obj@result$ONTOLOGY == ont, ]
  ridgeplot(sub,
            showCategory = n,
            fill = "p.adjust",
            core_enrichment = TRUE,
            label_format = 50,
            decreasing = TRUE) +
    ggtitle(ont)
}

onts <- c("BP", "CC", "MF")
onts <- onts[onts %in% gse_GO_all@result$ONTOLOGY]
ridge_list <- lapply(onts, make_ridge, obj = gse_GO_all)
names(ridge_list) <- onts

ridge_all <- wrap_plots(ridge_list, ncol = 1)
ridge_all
ggsave("plots/ridge.png", plot = ridge_all, width = 7, height = 12, dpi = 300)

# Cnetplot
make_cnet <- function(obj, ont, n = 20) {
  sub <- obj
  sub@result <- obj@result[obj@result$ONTOLOGY == ont, ]
  cnetplot(sub,
           layout = igraph::layout_nicely,
           showCategory = n,
           color_category = "#E5C494",
           size_category = 1,
           color_item = "#B3B3B3",
           size_item = 1,
           color_edge = "grey",
           size_edge = 0.5,
           node_label = "category",
           foldChange = gene_list,
           hilight = "none",
           hilight_alpha = 0.3) +
    ggtitle(ont)
}
cnet_list <- lapply(onts, make_cnet, obj = gse_GO_all)
names(cnet_list) <- onts
cnet_all <- wrap_plots(cnet_list, ncol = 1)

cnet_all
ggsave("plots/cnet.png", plot = cnet_all, width = 10, height = 10, dpi = 300)

# Emaplot
make_ema <- function(obj, ont, n = 20) {
  sub <- obj
  sub@result <- obj@result[obj@result$ONTOLOGY == ont, ]
  emapplot(sub,
           layout = igraph::layout_with_kk,
           showCategory = n,
           color = "p.adjust",
           size_category = 1,
           min_edge = 0.2,
           color_edge = "grey",
           size_edge = 0.5,
           node_label = "category_grouped",
           label_format = 30,
           clusterFunction = stats::kmeans,
           nWords = 4,
           nCluster = NULL) +
    ggtitle(ont)
}
ema_list <- lapply(onts, make_ema, obj = gse_GO_all)
names(ema_list) <- onts
ema_all <- wrap_plots(ema_list, ncol = 1)

ema_all
ggsave("plots/emaplot.png", plot = ema_all, width = 9, height = 40, dpi = 300)

# Single pathway GSEA plot
go_sets <- data.frame(gse_GO_all@result[["Description"]])
id_go = 4 # Choose pathway of interest

gsea <- gseaplot2(gse_GO_all, 
          geneSetID = id_go,
          pvalue_table = TRUE,
          title = gse_GO_all$Description[id_go],
          color = c("#E495A5", "#86B875", "#7DB0DD"), 
          ES_geom = "line")
gsea
ggsave("plots/gsea.png", plot = gsea, width = 9, height = 7, dpi = 300)

### 3. Convert gene_list to Entrez ID for further analyzes ---------------------
map <- bitr(names(gene_list),
            fromType = "SYMBOL",
            toType   = "ENTREZID",
            OrgDb    = org.Mm.eg.db)

gene_list_entrez <- gene_list[map$SYMBOL]
names(gene_list_entrez) <- map$ENTREZID

gene_list_entrez <- gene_list_entrez[!is.na(names(gene_list_entrez))]
gene_list_entrez <- gene_list_entrez[!duplicated(names(gene_list_entrez))]
gene_list_entrez <- sort(gene_list_entrez, decreasing = TRUE)

gse_kegg <- gseKEGG(geneList = gene_list_entrez,
                    keyType = "kegg",
                    nPerm = 10000,
                    organism = "mmu",
                    minGSSize = 3,
                    maxGSSize = 800,
                    pvalueCutoff = 0.05,
                    pAdjustMethod = "none")
gse_kegg <- pairwise_termsim(gse_kegg)

# Dotplot KEGG
dot_k <- dotplot(gse_kegg,
                 showCategory = 10,
                 split = ".sign",
                 font.size = 12) +
  facet_grid(. ~ .sign)

dot_k
ggsave("plots/dotplot_kegg.png", plot = dot, width = 7, height = 12, dpi = 300)

# Emaplot KEGG
emakegg <- emapplot(gse_kegg)

emakegg
ggsave("plots/ema_kegg.png", plot = emakegg, width = 8, height = 9, 
       dpi = 300, bg = "white")

# Cnetplot KEGG
cnetkegg <- cnetplot(gse_kegg,
                     layout = igraph::layout_nicely,
                     showCategory = 10,
                     color_category = "#E5C494",
                     size_category = 1,
                     size_item = 1,
                     color_edge = "grey",
                     size_edge = 0.5,
                     node_label = "category",
                     foldChange = gene_list_entrez,
                     hilight = "none",
                     hilight_alpha = 0.3)
cnetkegg
ggsave("plots/cnet_kegg.png", plot = cnetkegg, width = 7, height = 5, 
       dpi = 300, bg = "white")

# Ridge KEGG
ridgekegg <- ridgeplot(gse_kegg,
          showCategory = 10,
          fill = "p.adjust",
          core_enrichment = TRUE,
          label_format = 30
) + labs(x = "enrichment distribution")

ridgekegg
ggsave("plots/ridge_kegg.png", plot = ridgekegg, width = 7, height = 5, dpi = 300)

# Pathview
kegg_sets <- data.frame(Description = gse_kegg@result[["Description"]],
                        ID = gse_kegg@result[["ID"]])

run_pathview <- function(pathway_id, out_dir = file.path(project_dir, "plots", "pathview"), ...) {
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  old_wd <- setwd(out_dir)
  on.exit(setwd(old_wd))  
  pathview(gene.data = gene_list_entrez, pathway.id = pathway_id,
           species = "mmu", gene.idtype = "ENTREZ", kegg.dir = out_dir, ...)
}
run_pathview("mmu04668", limit = list(gene = 1, cpd = 1),
             out.suffix = "log2FC")

id_kegg = 8
gsea_kegg <- gseaplot2(gse_kegg, 
          geneSetID = id_kegg,
          pvalue_table = TRUE,
          title = gse_kegg$Description[id_kegg],
          color = c("#E495A5", "#86B875", "#7DB0DD"), 
          ES_geom = "line")

gsea_kegg
ggsave("plots/gsea_kegg.png", plot = gsea_kegg, width = 9, height = 7, dpi = 300)

### 4. Save enrichment tables ---------------------------------------------------
cat("Saving enrichment tables")
fdr_cut <- 0.05

# Full result table (core_enrichment = leading-edge genes, separated by "/")
gse_df <- as.data.frame(gse_GO_all) |>
  mutate(direction = ifelse(NES > 0, "activated", "suppressed")) |>
  arrange(p.adjust)

write.csv(gse_df, "results/GSEA_GO_all_full.csv", row.names = FALSE)

# Significant terms only
gse_sig <- gse_df |> filter(p.adjust < fdr_cut)
write.csv(gse_sig, "results/GSEA_GO_all_significant.csv", row.names = FALSE)

# Convert leading-edge genes from Entrez IDs to gene symbols
gse_kegg_sym <- setReadable(gse_kegg, OrgDb = org.Mm.eg.db, keyType = "ENTREZID")

# Full result table (core_enrichment = leading-edge genes, separated by "/")
kegg_df <- as.data.frame(gse_kegg_sym) |>
  mutate(
    Description = sub(" - Mus musculus.*$", "", Description),  # drop organism suffix
    direction   = ifelse(NES > 0, "activated", "suppressed")
  ) |>
  arrange(p.adjust)

write.csv(kegg_df, "results/GSEA_KEGG_full.csv", row.names = FALSE)

### 4. String Analysis----------------------------------------------------------
cat("Performing STRING analysis")

string_db <- STRINGdb$new(version = "12.0", species = 10090, # 10090 = mouse
                          score_threshold = 700, # 400 = medium, 700 = high
                          network_type = "full", input_directory = "")

mapped <- string_db$map(tt_filtered, "SYMBOL", removeUnmappedRows = TRUE)
mapped <- mapped[!duplicated(mapped$SYMBOL), ]    # one STRING_id per gene
mapped <- mapped[!duplicated(mapped$STRING_id), ] # and one gene per STRING_id

message("Mapped to STRING: ", nrow(mapped), " of ", nrow(tt_filtered))
hits <- mapped$STRING_id

# ---- PPI network ------------------------------------------
colored <- string_db$add_diff_exp_color(mapped, logFcColStr = "log2FoldChange")
payload <- string_db$post_payload(colored$STRING_id, colors = colored$color)

png("plots/rede_PPI.png", width = 3000, height = 3000, res = 300)
string_db$plot_network(hits, payload_id = payload)
dev.off()

# Export for Cytoscape (edges and node attributes with gene symbols)
edges <- string_db$get_interactions(hits)
edges$from <- mapped$SYMBOL[match(edges$from, mapped$STRING_id)]
edges$to   <- mapped$SYMBOL[match(edges$to,   mapped$STRING_id)]
write.csv(edges, "results/edges_PPI.csv", row.names = FALSE)
write.csv(mapped[, c("SYMBOL", "log2FoldChange", "padj")],
          "results/nodes_PPI.csv", row.names = FALSE)

# Genes with no interaction at this score threshold (isolated nodes)
connected <- unique(c(edges$from, edges$to))
message("Genes with at least one interaction: ", length(connected),
        " of ", nrow(mapped))

# ---- Enrichment (background = all tested genes) -----------
universe <- string_db$map(tt, "SYMBOL", removeUnmappedRows = TRUE)
universe <- universe[!duplicated(universe$STRING_id), ]
string_db$set_background(universe$STRING_id)   # must be set before get_enrichment()

enrich <- string_db$get_enrichment(hits)
write.csv(enrich, "results/enrichment_STRING_full.csv", row.names = FALSE)

print(table(enrich$category))

# Drop literature/network categories, which are mostly noise here
keep_cats <- c("Process", "Function", "Component", "KEGG", "RCTM",
               "WikiPathways", "MPO")
enrich_sig <- enrich |>
  filter(category %in% keep_cats, fdr < fdr_cut)

write.csv(enrich_sig, "results/enrichment_STRING_significant.csv",
          row.names = FALSE)

# ---- Plot helper ------------------------------------------
plot_enrichment <- function(df, cat, title, n = 15) {
  d <- df |> filter(category == cat)
  if (nrow(d) == 0) {
    message("No significant terms for ", cat, "; skipping plot.")
    return(invisible(NULL))
  }
  d <- d |>
    arrange(fdr) |>
    slice_head(n = n) |>
    mutate(
      gene_ratio  = number_of_genes / length(hits),
      description = str_wrap(description, 50),
      description = factor(description, levels = rev(unique(description)))
    )
  
  ggplot(d, aes(x = gene_ratio, y = description,
                size = number_of_genes, color = fdr)) +
    geom_point() +
    scale_color_gradient(low = "red", high = "blue") +
    labs(x = "Gene ratio", y = NULL, color = "FDR",
         size = "Gene count", title = title) +
    theme_bw(base_size = 12)
}

# ---- Reactome, WikiPathways ------------------------------
reac <- plot_enrichment(enrich_sig, "RCTM", "Reactome enrichment")
wp   <- plot_enrichment(enrich_sig, "WikiPathways", "WikiPathways enrichment")

if (!is.null(reac)) {
  print(reac)
  ggsave("plots/reactome.png", reac, width = 8, height = 4, dpi = 300)
}
if (!is.null(wp)) {
  ggsave("plots/wikipathways.png", wp, width = 8, height = 4, dpi = 300)
}

# ---- Mammalian Phenotype (MPO) ---------------------------
# Few terms expected; report as a table instead of a plot when n < 3
mpo <- enrich_sig |> filter(category == "MPO")
print(mpo[, c("term", "description", "number_of_genes",
              "number_of_genes_in_background", "fdr")])

# ---- Combined plot (all significant categories) ----------
if (nrow(enrich_sig) > 0) {
  multi_df <- enrich_sig |>
    group_by(category) |>
    slice_min(fdr, n = 8, with_ties = FALSE) |>
    ungroup() |>
    mutate(description = str_wrap(description, 45))
  
  comb <- ggplot(multi_df, aes(x = -log10(fdr),
                               y = reorder(description, -fdr),
                               fill = category)) +
    geom_col() +
    facet_grid(category ~ ., scales = "free_y", space = "free_y") +
    labs(x = expression(-log[10](FDR)), y = NULL) +
    theme_bw(base_size = 11) +
    theme(legend.position = "none")
  
  ggsave("plots/enrichment_combined.png", comb,
         width = 9, height = 8, dpi = 300)
}
