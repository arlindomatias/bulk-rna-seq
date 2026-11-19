# ============================================================
# Cross-dataset gene key (meta_analysis/results/gene_map.csv).
# One row per gene tested in each dataset (every row of its top_table.csv), with `key`: the Ensembl gene ID
# of the key species (config orthologs.key_species, mouse). Genes of the key species use their own
# Ensembl ID (symbols are mapped through the species OrgDb); genes of another species use their
# 1:1 ortholog (Ensembl Compara `ortholog_one2one`). Genes without a unique key get key = NA and are not
# matched across datasets. `map_status` says why:
#   ensembl          own Ensembl ID (key species)
#   symbol           symbol -> one Ensembl ID through the OrgDb (key species)
#   alias            outdated symbol -> one Ensembl ID through the OrgDb aliases (key species)
#   ortholog_1to1    1:1 ortholog of another species
#   ambiguous_symbol symbol -> several Ensembl IDs
#   no_ensembl       symbol without an Ensembl ID
#   not_1to1         only one-to-many / many-to-many orthologs
#   no_ortholog      no ortholog in the key species
#   duplicate_key    the key is shared by several genes of the same dataset (all of them get NA)
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
key_species <- cfg$orthologs$key_species
compara_name <- function(species) gsub(" ", "_", tolower(species))

# Compara homologies between `species` and the key species, from every downloaded file. A pair is stored
# once, in the file of either species, so both orientations are read.
read_homologies <- function(species) {
  bind_rows(lapply(unlist(snakemake@input$homologies), function(f) {
    h <- read_tsv(f, col_types = cols(.default = col_character()), na = c("", "NULL"),
                  col_select = c("gene_stable_id", "species", "homology_type",
                                 "homology_gene_stable_id", "homology_species", "is_high_confidence"))
    bind_rows(
      h %>% filter(species == !!species, homology_species == key_species) %>%
        transmute(gene_id = gene_stable_id, key = homology_gene_stable_id, homology_type, is_high_confidence),
      h %>% filter(species == key_species, homology_species == !!species) %>%
        transmute(gene_id = homology_gene_stable_id, key = gene_stable_id, homology_type, is_high_confidence)
    )
  })) %>% distinct(gene_id, key, homology_type, .keep_all = TRUE)
}

map_one <- function(gse, d) {
  tt <- read_csv(file.path(d$project, "results", "top_table.csv"), show_col_types = FALSE,
                 col_select = c("gene_id", "SYMBOL"))
  species <- compara_name(d$species)
  ids <- sub("\\.[0-9]+$", "", tt$gene_id)

  if (species == key_species) {
    if (d$gene_id_type == "ENSEMBL") {
      key <- ids
      status <- rep("ensembl", length(ids))
    } else {
      suppressPackageStartupMessages(library(d$orgdb, character.only = TRUE))
      orgdb <- get(d$orgdb)
      valid <- intersect(unique(ids), AnnotationDbi::keys(orgdb, keytype = "SYMBOL"))
      # Unique Ensembl ID of each key (NA when none or several); `found` = keys with any Ensembl ID.
      unique_ensembl <- function(keys, keytype) {
        valid <- intersect(unique(keys), AnnotationDbi::keys(orgdb, keytype = keytype))
        if (!length(valid)) return(list(one = character(0), found = character(0)))
        tab <- suppressMessages(AnnotationDbi::select(orgdb, keys = valid, keytype = keytype,
                                                      columns = "ENSEMBL"))
        tab <- distinct(tab[!is.na(tab$ENSEMBL), c(keytype, "ENSEMBL")])
        n <- table(tab[[keytype]])
        keep <- tab[[keytype]] %in% names(n)[n == 1]
        list(one = setNames(tab$ENSEMBL[keep], tab[[keytype]][keep]), found = names(n))
      }
      by_symbol <- unique_ensembl(ids, "SYMBOL")
      key <- unname(by_symbol$one[ids])
      status <- ifelse(!is.na(key), "symbol",
                       ifelse(ids %in% by_symbol$found, "ambiguous_symbol", "no_ensembl"))
      # Symbols that are not current (e.g. Ctgf, now Ccn2): try them as aliases, only when unique.
      # The alias must not be a current symbol of another gene.
      old <- status == "no_ensembl" & !(ids %in% AnnotationDbi::keys(orgdb, keytype = "SYMBOL"))
      by_alias <- unique_ensembl(ids[old], "ALIAS")
      alias_key <- unname(by_alias$one[ids[old]])
      key[old] <- alias_key
      status[old] <- ifelse(!is.na(alias_key), "alias",
                            ifelse(ids[old] %in% by_alias$found, "ambiguous_symbol", "no_ensembl"))
    }
  } else {
    hom <- read_homologies(species)
    one2one <- hom %>% filter(homology_type == "ortholog_one2one")
    stopifnot(!anyDuplicated(one2one$gene_id))
    key <- one2one$key[match(ids, one2one$gene_id)]
    status <- ifelse(!is.na(key), "ortholog_1to1",
                     ifelse(ids %in% hom$gene_id, "not_1to1", "no_ortholog"))
  }

  # An alias whose key is already the current symbol of another gene of the dataset is dropped; other
  # shared keys drop every gene that has them.
  alias_dup <- status == "alias" & key %in% key[status != "alias"]
  key[alias_dup] <- NA
  status[alias_dup] <- "duplicate_key"
  dup <- !is.na(key) & (duplicated(key) | duplicated(key, fromLast = TRUE))
  key[dup] <- NA
  status[dup] <- "duplicate_key"

  message(gse, " (", species, "): ", length(ids), " genes; ",
          paste(names(table(status)), table(status), sep = " ", collapse = ", "))
  tibble(dataset = gse, species = d$species, gene_id = tt$gene_id, SYMBOL = tt$SYMBOL,
         key = key, map_status = status)
}

gene_map <- bind_rows(Map(map_one, names(cfg$datasets), cfg$datasets))

# Symbol of the key gene, for display: from the key-species datasets, current symbols before aliases.
key_sym <- gene_map %>%
  filter(compara_name(species) == key_species, !is.na(key), !is.na(SYMBOL), SYMBOL != "") %>%
  arrange(map_status == "alias") %>%
  distinct(key, SYMBOL) %>% group_by(key) %>% summarise(key_symbol = SYMBOL[1], .groups = "drop")
gene_map <- gene_map %>% left_join(key_sym, by = "key") %>%
  mutate(key_symbol = ifelse(is.na(key_symbol) & !is.na(key), SYMBOL, key_symbol))

write_csv(gene_map, snakemake@output$table, na = "")

keys <- split(gene_map$key[!is.na(gene_map$key)], gene_map$dataset[!is.na(gene_map$key)])
for (p in combn(names(keys), 2, simplify = FALSE)) {
  message(p[1], " & ", p[2], ": ", length(intersect(keys[[p[1]]], keys[[p[2]]])), " shared keys")
}
message("Keys in every dataset: ", length(Reduce(intersect, keys)))
