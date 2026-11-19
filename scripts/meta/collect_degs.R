# ============================================================
# Collect the DEGs of every DE project into one table (meta_analysis/results/all_degs.csv).
# One row per gene per dataset where it is a DEG. The DEG criterion of each dataset is read from its
# de_summary.tsv (deg_criterion = padj or pvalue, deg_cutoff, lfc_cutoff), so it matches the project's
# own report. Missing identifiers are filled from the species OrgDb; several IDs are joined with ";".
# Genes are matched across datasets by `key` from results/gene_map.csv (mouse Ensembl ID; rat genes through
# their 1:1 ortholog); genes without a key are never matched.
# ============================================================
if (!exists("snakemake")) {
  stop("Run this script through Snakemake (`snakemake --use-conda --cores 1`).")
}

log_con <- file(snakemake@log[[1]], open = "wt")
sink(log_con)
sink(log_con, type = "message")

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(AnnotationDbi)
})

cfg <- snakemake@config
gene_map <- read_csv(snakemake@input$gene_map, show_col_types = FALSE,
                     col_types = cols(.default = col_character()))

# All values of `column` for each key, joined with ";" (NA when there is none).
map_ids <- function(orgdb, keys, keytype, column) {
  keys <- unique(keys[!is.na(keys) & keys != ""])
  valid <- intersect(keys, AnnotationDbi::keys(orgdb, keytype = keytype))
  if (!length(valid)) return(setNames(character(0), character(0)))
  tab <- suppressMessages(AnnotationDbi::select(orgdb, keys = valid, keytype = keytype, columns = column))
  tab <- tab[!is.na(tab[[column]]), ]
  vapply(split(tab[[column]], tab[[keytype]]),
         function(v) paste(sort(unique(v)), collapse = ";"), character(1))
}

lookup <- function(map, keys) unname(map[keys])

collect_one <- function(gse, d) {
  tt  <- read_csv(file.path(d$project, "results", "top_table.csv"), show_col_types = FALSE)
  sm  <- read_tsv(file.path(d$project, "results", "de_summary.tsv"), show_col_types = FALSE)
  crit <- sm$deg_criterion[1]
  cut  <- sm$deg_cutoff[1]
  lfc  <- sm$lfc_cutoff[1]
  stopifnot(crit %in% c("padj", "pvalue"))

  deg <- tt %>%
    filter(!is.na(.data[[crit]]), .data[[crit]] < cut, abs(log2FoldChange) > lfc)

  n_up <- sum(deg$log2FoldChange > 0)
  n_down <- sum(deg$log2FoldChange < 0)
  if (n_up != sm$deg_up[1] || n_down != sm$deg_down[1]) {
    stop(gse, ": DEG count (", n_up, "/", n_down, ") differs from de_summary.tsv (",
         sm$deg_up[1], "/", sm$deg_down[1], ")")
  }

  suppressPackageStartupMessages(library(d$orgdb, character.only = TRUE))
  orgdb <- get(d$orgdb)

  if (d$gene_id_type == "ENSEMBL") {
    ensembl <- deg$gene_id
    entrez  <- lookup(map_ids(orgdb, ensembl, "ENSEMBL", "ENTREZID"), ensembl)
    symbol  <- ifelse(is.na(deg$SYMBOL) | deg$SYMBOL == "",
                      lookup(map_ids(orgdb, ensembl, "ENSEMBL", "SYMBOL"), ensembl),
                      deg$SYMBOL)
  } else {
    symbol  <- deg$SYMBOL
    ensembl <- lookup(map_ids(orgdb, symbol, "SYMBOL", "ENSEMBL"), symbol)
    entrez  <- lookup(map_ids(orgdb, symbol, "SYMBOL", "ENTREZID"), symbol)
  }

  message(gse, ": ", nrow(deg), " DEGs (", n_up, " up / ", n_down, " down; ", crit, " < ", cut,
          ", |log2FC| > ", lfc, "); without Entrez: ", sum(is.na(entrez)),
          "; without Ensembl: ", sum(is.na(ensembl)))

  gm <- gene_map[gene_map$dataset == gse, ]
  gm_row <- match(deg$gene_id, gm$gene_id)
  stopifnot(!anyNA(gm_row))

  tibble(
    key            = gm$key[gm_row],
    key_symbol     = gm$key_symbol[gm_row],
    SYMBOL         = symbol,
    ENTREZID       = entrez,
    ENSEMBL        = ensembl,
    dataset        = gse,
    species        = d$species,
    contrast       = sm$contrast[1],
    method         = d$method,
    direction      = ifelse(deg$log2FoldChange > 0, "up", "down"),
    log2FoldChange = deg$log2FoldChange,
    stat           = deg$stat,
    stat_type      = d$stat_type,
    pvalue         = deg$pvalue,
    padj           = deg$padj,
    deg_criterion  = paste0(crit, " < ", cut, ", |log2FC| > ", lfc,
                            ifelse(crit == "pvalue", " (exploratory: no gene at FDR)", ""))
  )
}

all_degs <- bind_rows(Map(collect_one, names(cfg$datasets), cfg$datasets))

# Datasets in which each gene (key) is a DEG. Genes without a key only count their own dataset.
all_degs <- all_degs %>%
  mutate(.key = ifelse(is.na(key), paste(dataset, row_number()), key)) %>%
  group_by(.key) %>%
  mutate(n_datasets = n_distinct(dataset),
         datasets   = paste(sort(unique(dataset)), collapse = ";")) %>%
  ungroup() %>%
  dplyr::select(-.key) %>%
  arrange(desc(n_datasets), toupper(coalesce(key_symbol, SYMBOL)), dataset)

write_csv(all_degs, snakemake@output$table, na = "")
writeLines(capture.output(sessionInfo()), snakemake@output$session)
message("Rows: ", nrow(all_degs), "; without key: ", sum(is.na(all_degs$key)),
        "; genes DEG in more than one dataset: ", n_distinct(all_degs$key[all_degs$n_datasets > 1]))
