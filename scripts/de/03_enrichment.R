# ============================================================
# SCRIPT 03: ENRICHMENT ANALYSIS
# GSEA uses the gene identifiers of the count table (config gene_id_type): Ensembl IDs without
# version for tximport projects (in GENCODE some symbols map to more than one gene), symbols for
# GEO symbol tables. KEGG uses Entrez IDs mapped from them. STRING always uses symbols.
# ============================================================
source(file.path(snakemake@scriptdir, "00_common.R"))
suppressPackageStartupMessages({
  library(clusterProfiler)
  library(enrichplot)
  library(pathview)
  library(STRINGdb)
  library(patchwork)
  library(stringr)
  library(igraph)
})
out <- snakemake@output
ecfg <- cfg$enrichment
suppressPackageStartupMessages(library(ecfg$orgdb, character.only = TRUE))
orgdb <- get(ecfg$orgdb)

dir_up   <- paste0("Enriched_", trt_level)
dir_down <- paste0("Enriched_", ref_level)

tt <- as.data.frame(read_csv(snakemake@input$top_table, show_col_types = FALSE))
deg_rule <- deg_criterion(tt)
sig_col <- deg_rule$col
tt_filtered <- tt[!is.na(tt[[sig_col]]) & tt[[sig_col]] < deg_rule$cut &
                    abs(tt$log2FoldChange) > cfg$de$lfc, ]
message("DEGs: ", nrow(tt_filtered), " (", sig_col, " < ", deg_rule$cut, ")")

write_empty <- function(path, what) {
  write_csv(data.frame(), path)
  message("No significant ", what, "; wrote an empty table: ", path)
}
save_plot <- function(path, expr, ...) {
  tryCatch(ggsave(path, plot = expr, dpi = 300, bg = "white", ...),
           error = function(e) message("Plot ", path, " failed: ", conditionMessage(e)))
}

### 1. Ranking for GSEA: Wald statistic -------------------------------------
gene_list <- setNames(tt$stat, tt$gene_id)
gene_list <- sort(gene_list[!is.na(gene_list)], decreasing = TRUE)
fc_list   <- setNames(tt$log2FoldChange, tt$gene_id)[names(gene_list)]  # colors only

n_mapped <- sum(names(gene_list) %in% keys(orgdb, gene_id_type))
message("Ranked genes: ", length(gene_list), " | in ", ecfg$orgdb, " (", gene_id_type, "): ", n_mapped)

n_gsea_plot <- ecfg$gsea_plot_top_n %||% 3   # pathways/terms in the single-set GSEA plots
gsea_args <- list(minGSSize = ecfg$min_gs_size, maxGSSize = ecfg$max_gs_size,
                  pvalueCutoff = ecfg$pvalue_cutoff, pAdjustMethod = "BH", eps = 0,
                  seed = TRUE)

### 2. GSEA GO -----------------------------------------------------------------
message("Performing GSEA (GO)")
gse_GO_all <- do.call(gseGO, c(list(geneList = gene_list, ont = "ALL", keyType = gene_id_type,
                                    OrgDb = orgdb, verbose = FALSE), gsea_args))
has_go <- !is.null(gse_GO_all) && nrow(gse_GO_all@result) > 0

if (has_go) {
  if (gene_id_type != "SYMBOL") gse_GO_all <- setReadable(gse_GO_all, OrgDb = orgdb, keyType = gene_id_type)
  gse_df <- as.data.frame(gse_GO_all) |>
    mutate(direction = ifelse(NES > 0, dir_up, dir_down)) |>
    arrange(p.adjust)
  write_csv(gse_df, out$go)
  message("Significant GO terms: ", nrow(gse_df))

  if (nrow(gse_GO_all@result) >= 2) gse_GO_all <- pairwise_termsim(gse_GO_all)

  # <figure name="dotplot_GO">
  # Top 10 terms per ontology *and* direction. showCategory = 10 with split = "ONTOLOGY" ranks by
  # p.adjust within each ontology only, so the much more significant activated terms hid all the
  # suppressed ones.
  top_go <- gse_GO_all@result |>
    mutate(sign = ifelse(NES > 0, "activated", "suppressed")) |>
    group_by(ONTOLOGY, sign) |>
    slice_min(p.adjust, n = 10, with_ties = FALSE) |>
    ungroup()
  go_dot <- gse_GO_all
  go_dot@result <- gse_GO_all@result[gse_GO_all@result$ID %in% top_go$ID, ]
  save_plot("plots/dotplot_GO.png",
            dotplot(go_dot, showCategory = nrow(go_dot@result), split = ".sign", font.size = 11,
                    label_format = 60) +
              facet_grid(ONTOLOGY ~ .sign, scales = "free_y", space = "free_y"),
            width = 13, height = 18)
  # </figure>

  onts <- names(which(table(gse_GO_all@result$ONTOLOGY) >= 2))
  sub_ont <- function(obj, ont) {
    obj@result <- obj@result[obj@result$ONTOLOGY == ont, ]
    obj
  }

  if (length(onts) > 0) {
    # <figure name="ridge">
    ridge_list <- lapply(onts, function(o)
      ridgeplot(sub_ont(gse_GO_all, o), showCategory = 10, fill = "p.adjust",
                core_enrichment = TRUE, label_format = 50, decreasing = TRUE) + ggtitle(o))
    save_plot("plots/ridge.png", wrap_plots(ridge_list, ncol = 1),
              width = 7, height = 4 * length(onts))
    # </figure>

    # <figure name="cnet">
    cnet_list <- lapply(onts, function(o)
      cnetplot(sub_ont(gse_GO_all, o), layout = igraph::layout_nicely, showCategory = 20,
               color_category = "#E5C494", size_category = 1, color_item = "#B3B3B3",
               size_item = 1, color_edge = "grey", size_edge = 0.5,
               node_label = "category", foldChange = fc_list,
               hilight = "none", hilight_alpha = 0.3) + ggtitle(o))
    save_plot("plots/cnet.png", wrap_plots(cnet_list, ncol = 1),
              width = 10, height = 6 * length(onts))
    # </figure>

    # <figure name="emaplot">
    ema_list <- lapply(onts, function(o)
      emapplot(sub_ont(gse_GO_all, o), layout = igraph::layout_with_kk, showCategory = 20,
               color = "p.adjust", size_category = 1, min_edge = 0.2,
               color_edge = "grey", size_edge = 0.5,
               node_label = "category",  # enrichplot >= 1.30: category, group, all or none
               node_label_size = 3, label_format = 30, clusterFunction = stats::kmeans, nWords = 4,
               nCluster = NULL) + ggtitle(o))
    save_plot("plots/emaplot.png", wrap_plots(ema_list, ncol = 1),
              width = 9, height = 8 * length(onts), limitsize = FALSE)
    # </figure>
  } else {
    message("No GO ontology has >= 2 significant terms; skipping ridge/cnet/ema plots.")
  }

  # <figure name="gsea">
  # Running score of the n most significant terms (result is sorted by p.adjust).
  # pvalue_table = FALSE: with several sets, enrichplot 1.30 mislabels the rows of that table
  # (names and colors do not match the curves); the legend is correct and the report has the p-values.
  top_ids <- seq_len(min(n_gsea_plot, nrow(gse_GO_all@result)))
  save_plot("plots/gsea.png",
            gseaplot2(gse_GO_all, geneSetID = top_ids, pvalue_table = FALSE,
                      title = paste("Top", length(top_ids), "GO terms"), ES_geom = "line"),
            width = 11, height = 8)
  # </figure>
} else {
  write_empty(out$go, "GO terms")
}

### 3. GSEA KEGG (Entrez IDs) ---------------------------------------------------
map <- bitr(names(gene_list), fromType = gene_id_type, toType = "ENTREZID", OrgDb = orgdb)
map$stat <- gene_list[map[[gene_id_type]]]
map <- map[order(-abs(map$stat)), ]                       # on conflicts, keep the largest |stat|
map <- map[!duplicated(map[[gene_id_type]]) & !duplicated(map$ENTREZID), ]
message("Genes with Entrez ID for KEGG: ", nrow(map), " of ", length(gene_list))

gene_list_entrez <- sort(setNames(map$stat, map$ENTREZID), decreasing = TRUE)
fc_entrez <- setNames(fc_list[map[[gene_id_type]]], map$ENTREZID)

# Local annotation (kegg_annotation rule) instead of gseKEGG downloading it on every run.
# Formats: "path:mmu00010<TAB>mmu:100042025" and "mmu00010<TAB>Glycolysis ... - Mus musculus (house mouse)".
kegg_link <- read_tsv(snakemake@input$kegg_link, col_names = c("term", "gene"),
                      col_types = "cc") |>
  mutate(term = sub("^path:", "", term), gene = sub("^[a-z]+:", "", gene))
kegg_names <- read_tsv(snakemake@input$kegg_names, col_names = c("term", "name"),
                       col_types = "cc") |>
  mutate(term = sub("^path:", "", term))
message("KEGG annotation: ", n_distinct(kegg_link$term), " pathways, ",
        n_distinct(kegg_link$gene), " genes (downloaded ",
        format(file.mtime(snakemake@input$kegg_link), "%Y-%m-%d"), ")")

message("Performing GSEA (KEGG)")
gse_kegg <- do.call(GSEA, c(list(geneList = gene_list_entrez, TERM2GENE = kegg_link,
                                 TERM2NAME = kegg_names, verbose = FALSE), gsea_args))
has_kegg <- !is.null(gse_kegg) && nrow(gse_kegg@result) > 0

if (has_kegg) {
  gse_kegg <- setReadable(gse_kegg, OrgDb = orgdb, keyType = "ENTREZID")
  # KEGG names end in " - <species> (<common name>)", e.g. " - Mus musculus (house mouse)".
  gse_kegg@result$Description <- sub(" - [^-]+ \\([^)]*\\)$", "", gse_kegg@result$Description)
  kegg_df <- as.data.frame(gse_kegg) |>
    mutate(direction = ifelse(NES > 0, dir_up, dir_down)) |>
    arrange(p.adjust)
  write_csv(kegg_df, out$kegg)
  message("Significant KEGG pathways: ", nrow(kegg_df))

  # <figure name="ema_kegg">
  if (nrow(gse_kegg@result) >= 2) {
    gse_kegg <- pairwise_termsim(gse_kegg)
    save_plot("plots/ema_kegg.png", emapplot(gse_kegg), width = 8, height = 9)
  }
  # </figure>
  # <figure name="dotplot_kegg">
  save_plot("plots/dotplot_kegg.png",
            dotplot(gse_kegg, showCategory = 10, split = ".sign", font.size = 12) +
              facet_grid(. ~ .sign),
            width = 9, height = 8)
  # </figure>
  # <figure name="cnet_kegg">
  save_plot("plots/cnet_kegg.png",
            cnetplot(gse_kegg, layout = igraph::layout_nicely, showCategory = 10,
                     color_category = "#E5C494", size_category = 1, size_item = 1,
                     color_edge = "grey", size_edge = 0.5, node_label = "category",
                     foldChange = fc_entrez, hilight = "none", hilight_alpha = 0.3),
            width = 9, height = 7)
  # </figure>
  # <figure name="ridge_kegg">
  save_plot("plots/ridge_kegg.png",
            ridgeplot(gse_kegg, showCategory = 10, fill = "p.adjust",
                      core_enrichment = TRUE, label_format = 30) +
              labs(x = "Enrichment distribution"),
            width = 7, height = 5)
  # </figure>
  # <figure name="gsea_kegg">
  top_ids <- seq_len(min(n_gsea_plot, nrow(gse_kegg@result)))
  save_plot("plots/gsea_kegg.png",
            gseaplot2(gse_kegg, geneSetID = top_ids, pvalue_table = FALSE,
                      title = paste("Top", length(top_ids), "KEGG pathways"), ES_geom = "line"),
            width = 11, height = 8)
  # </figure>
} else {
  write_empty(out$kegg, "KEGG pathways")
}

### 4. Pathview -----------------------------------------------------------------
pathways <- unlist(ecfg$pathview_pathways)
if (!length(pathways) && has_kegg) pathways <- head(kegg_df$ID, ecfg$pathview_top_n)

if (length(pathways)) {
  # pathview writes the colored map to the working directory and the downloaded KEGG files
  # (KGML + base PNG) to kegg.dir: figures go to plots/pathview/, downloads to cache/kegg/pathview/.
  # Create before normalizing: normalizePath() of a missing directory returns the relative path,
  # which stops being valid after setwd().
  dir.create(file.path("plots", "pathview"), recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path("cache", "kegg", "pathview"), recursive = TRUE, showWarnings = FALSE)
  pv_dir <- normalizePath(file.path("plots", "pathview"), mustWork = TRUE)
  kegg_dir <- normalizePath(file.path("cache", "kegg", "pathview"), mustWork = TRUE)
  old_wd <- setwd(pv_dir)
  old_opt <- options(timeout = 180)   # KGML/PNG files are small; do not wait 10 min per pathway
  on.exit(options(old_opt), add = TRUE)
  for (pw in pathways) {
    message("Pathview: ", pw)
    tryCatch(pathview(gene.data = fc_entrez, pathway.id = pw, species = ecfg$organism_kegg,
                      gene.idtype = "ENTREZ", kegg.dir = kegg_dir),
             error = function(e) message("pathview ", pw, " failed: ", conditionMessage(e)))
  }
  setwd(old_wd)
} else {
  message("No KEGG pathways selected for Pathview.")
}

### 5. STRING ------------------------------------------------------------------
message("Performing STRING analysis")
dir.create(ecfg$string_cache, recursive = TRUE, showWarnings = FALSE)
string_db <- STRINGdb$new(version = ecfg$string_version, species = ecfg$string_species,
                          score_threshold = ecfg$string_score, network_type = "full",
                          input_directory = ecfg$string_cache)

# One gene per symbol (lowest padj) before mapping.
degs_sym <- tt_filtered |>
  filter(!is.na(SYMBOL)) |>
  arrange(.data[[sig_col]]) |>
  distinct(SYMBOL, .keep_all = TRUE) |>
  as.data.frame()

hits <- character(0)
if (nrow(degs_sym) > 0) {
  mapped <- string_db$map(degs_sym, "SYMBOL", removeUnmappedRows = TRUE)
  mapped <- mapped[!duplicated(mapped$SYMBOL) & !duplicated(mapped$STRING_id), ]
  hits <- mapped$STRING_id
  message("Mapped to STRING: ", nrow(mapped), " of ", nrow(degs_sym))
}

if (length(hits) >= 2) {
  # PPI network (the image comes from the STRING server)
  tryCatch({
    # <figure name="rede_PPI">
    colored <- string_db$add_diff_exp_color(mapped, logFcColStr = "log2FoldChange")
    payload <- string_db$post_payload(colored$STRING_id, colors = colored$color)
    png("plots/rede_PPI.png", width = 3000, height = 3000, res = 300)
    string_db$plot_network(hits, payload_id = payload)
    dev.off()
    # </figure>
  }, error = function(e) message("STRING network plot failed: ", conditionMessage(e)))

  # Cytoscape tables
  edges <- string_db$get_interactions(hits)
  edges$from <- mapped$SYMBOL[match(edges$from, mapped$STRING_id)]
  edges$to   <- mapped$SYMBOL[match(edges$to,   mapped$STRING_id)]
  write_csv(edges, out$edges)
  write_csv(mapped[, unique(c("SYMBOL", "gene_id", "log2FoldChange", "padj", sig_col))], out$nodes)
  message("Genes with at least one interaction: ",
          length(unique(c(edges$from, edges$to))), " of ", nrow(mapped))

  # Enrichment, background = all tested genes
  universe <- tt |> filter(!is.na(SYMBOL)) |> distinct(SYMBOL) |> as.data.frame()
  universe <- string_db$map(universe, "SYMBOL", removeUnmappedRows = TRUE)
  string_db$set_background(unique(universe$STRING_id))
  enrich <- string_db$get_enrichment(hits)
  write_csv(enrich, out$string)

  keep_cats <- c("Process", "Function", "Component", "KEGG", "RCTM", "WikiPathways", "MPO")
  enrich_sig <- enrich |> filter(category %in% keep_cats, fdr < 0.05)

  # <figure name="plot_enrichment">
  plot_enrichment <- function(df, cat, title, n = 15) {
    d <- df |> filter(category == cat)
    if (nrow(d) == 0) {
      message("No significant terms for ", cat, "; skipping plot.")
      return(NULL)
    }
    d |>
      arrange(fdr) |>
      slice_head(n = n) |>
      mutate(gene_ratio  = number_of_genes / length(hits),
             description = str_wrap(description, 50),
             description = factor(description, levels = rev(unique(description)))) |>
      ggplot(aes(x = gene_ratio, y = description, size = number_of_genes, color = fdr)) +
      geom_point() +
      scale_color_gradient(low = "red", high = "blue") +
      labs(x = "Gene ratio", y = NULL, color = "FDR", size = "Gene count", title = title) +
      theme_bw(base_size = 12)
  }
  # </figure>
  # <figure name="reactome">
  reac <- plot_enrichment(enrich_sig, "RCTM", "Reactome enrichment")
  wp   <- plot_enrichment(enrich_sig, "WikiPathways", "WikiPathways enrichment")
  if (!is.null(reac)) save_plot("plots/reactome.png", reac, width = 8, height = 5)
  if (!is.null(wp))   save_plot("plots/wikipathways.png", wp, width = 8, height = 5)
  # </figure>

  mpo <- enrich_sig |> filter(category == "MPO")
  if (nrow(mpo)) print(mpo[, c("term", "description", "number_of_genes",
                               "number_of_genes_in_background", "fdr")])

  if (nrow(enrich_sig) > 0) {
    # <figure name="enrichment_combined">
    multi_df <- enrich_sig |>
      group_by(category) |>
      slice_min(fdr, n = 8, with_ties = FALSE) |>
      ungroup() |>
      mutate(description = str_wrap(description, 45))
    comb <- ggplot(multi_df, aes(x = -log10(fdr), y = reorder(description, -fdr),
                                 fill = category)) +
      geom_col() +
      facet_grid(category ~ ., scales = "free_y", space = "free_y") +
      labs(x = expression(-log[10](FDR)), y = NULL) +
      theme_bw(base_size = 11) +
      theme(legend.position = "none")
    save_plot("plots/enrichment_combined.png", comb, width = 9, height = 10)
    # </figure>
  }
} else {
  message("Fewer than 2 DEGs mapped to STRING; skipping network and enrichment.")
  for (p in c(out$edges, out$nodes, out$string)) write_empty(p, "STRING results")
}

writeLines(capture.output(sessionInfo()), out$session)
