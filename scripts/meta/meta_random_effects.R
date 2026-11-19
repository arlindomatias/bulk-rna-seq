# ============================================================
# Random-effects meta-analysis per gene (meta_analysis/results/meta_random_effects.csv).
# For every gene (key from results/gene_map.csv) tested in at least `meta.min_datasets` datasets, the
# log2FoldChange of each dataset is combined with its standard error (lfcSE) by metafor::rma:
#   - tau² by REML (DL when REML does not converge), z test, BH over all genes -> padj;
#   - Knapp-Hartung test as a sensitivity column (pvalue_knha; very wide with 2-3 datasets);
#   - heterogeneity: Cochran's Q (and p), I², tau²;
#   - leave-one-out: largest p-value and whether the pooled sign flips when each dataset is dropped (with
#     2 datasets this is the p-value of each dataset alone, so loo_max_pvalue is that of the weaker one).
# A meta DEG has padj < meta.padj and |estimate| > meta.lfc. DEG status in each project is added for
# comparison. DESeq2 genes with pvalue NA (Cook's outliers) are left out of that dataset.
# ============================================================
if (!exists("snakemake")) {
  stop("Run this script through Snakemake (`snakemake --use-conda --cores 1`).")
}

log_con <- file(snakemake@log[[1]], open = "wt")
sink(log_con)
sink(log_con, type = "message")
options(warn = 1)  # print warnings to the log as they happen

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
  library(metafor)
})

cfg <- snakemake@config
mcfg <- cfg$meta
datasets <- names(cfg$datasets)

gene_map <- read_csv(snakemake@input$gene_map, show_col_types = FALSE,
                     col_types = cols(.default = col_character()))
degs <- read_csv(snakemake@input$degs, show_col_types = FALSE,
                 col_types = cols(.default = col_character()))

# One row per dataset and key: log2FC, SE and p of the project's own table.
effects <- bind_rows(lapply(datasets, function(gse) {
  d <- cfg$datasets[[gse]]
  tt <- read_csv(file.path(d$project, "results", "top_table.csv"), show_col_types = FALSE,
                 col_select = c("gene_id", "log2FoldChange", "lfcSE", "pvalue"))
  gm <- gene_map[gene_map$dataset == gse & !is.na(gene_map$key), c("gene_id", "key")]
  out <- inner_join(tt, gm, by = "gene_id") %>%
    filter(!is.na(log2FoldChange), !is.na(lfcSE), lfcSE > 0, !is.na(pvalue)) %>%
    transmute(dataset = gse, key, yi = log2FoldChange, sei = lfcSE)
  stopifnot(!anyDuplicated(out$key))
  message(gse, ": ", nrow(out), " genes with a key and a usable log2FC/SE")
  out
}))

k_per_key <- count(effects, key, name = "k")
keep <- k_per_key$key[k_per_key$k >= mcfg$min_datasets]
effects <- effects[effects$key %in% keep, ]
message("Genes in >= ", mcfg$min_datasets, " datasets: ", length(keep),
        " (k = ", paste(names(table(k_per_key$k[k_per_key$k >= mcfg$min_datasets])),
                        table(k_per_key$k[k_per_key$k >= mcfg$min_datasets]), sep = ": ", collapse = ", "),
        ")")

fit_rma <- function(yi, sei, test = "z") {
  fit <- tryCatch(rma(yi = yi, sei = sei, method = mcfg$method, test = test,
                      control = list(stepadj = 0.5, maxiter = 1000)),
                  error = function(e) NULL)
  if (is.null(fit)) fit <- rma(yi = yi, sei = sei, method = "DL", test = test)
  fit
}

meta_one <- function(e) {
  fit <- fit_rma(e$yi, e$sei)
  knha <- tryCatch(fit_rma(e$yi, e$sei, test = "knha")$pval, error = function(err) NA_real_)
  loo_p <- NA_real_
  loo_flip <- NA
  if (nrow(e) >= 2) {
    loo <- vapply(seq_len(nrow(e)), function(i) {
      f <- fit_rma(e$yi[-i], e$sei[-i])
      c(as.numeric(f$beta), f$pval)
    }, numeric(2))
    loo_p <- max(loo[2, ])
    loo_flip <- any(sign(loo[1, ]) != sign(as.numeric(fit$beta)))
  }
  tibble(
    k             = nrow(e),
    datasets      = paste(e$dataset, collapse = ";"),
    estimate      = as.numeric(fit$beta),
    se            = fit$se,
    ci_lb         = fit$ci.lb,
    ci_ub         = fit$ci.ub,
    z             = fit$zval,
    pvalue        = fit$pval,
    pvalue_knha   = knha,
    tau2          = fit$tau2,
    tau2_method   = fit$method,
    I2            = fit$I2,
    Q             = fit$QE,
    Q_pvalue      = fit$QEp,
    n_same_sign   = sum(sign(e$yi) == sign(as.numeric(fit$beta))),
    loo_max_pvalue = loo_p,
    loo_sign_flip = loo_flip
  )
}

t0 <- Sys.time()
groups <- split(effects, effects$key)
res <- bind_rows(lapply(groups, meta_one), .id = "key")
message("Fitted ", nrow(res), " genes in ", round(difftime(Sys.time(), t0, units = "mins"), 1), " min; ",
        sum(res$tau2_method != mcfg$method), " fell back to DL")

# Per-dataset log2FC / SE (wide) and project DEG status.
wide <- effects %>%
  pivot_wider(id_cols = key, names_from = dataset, values_from = c(yi, sei), names_glue = "{.value}_{dataset}") %>%
  rename_with(~ sub("^yi_", "log2FC_", sub("^sei_", "lfcSE_", .x)))
deg_in <- degs %>% filter(!is.na(key)) %>% group_by(key) %>%
  summarise(deg_in = paste(sort(unique(dataset)), collapse = ";"), .groups = "drop")
sym <- gene_map %>% filter(!is.na(key), !is.na(key_symbol)) %>% distinct(key, key_symbol) %>%
  group_by(key) %>% summarise(key_symbol = key_symbol[1], .groups = "drop")

res <- res %>%
  mutate(padj = p.adjust(pvalue, method = "BH"),
         meta_deg = !is.na(padj) & padj < mcfg$padj & abs(estimate) > mcfg$lfc,
         direction = ifelse(estimate > 0, "up", "down")) %>%
  left_join(sym, by = "key") %>%
  left_join(deg_in, by = "key") %>%
  left_join(wide, by = "key") %>%
  relocate(key, key_symbol, k, datasets, estimate, se, ci_lb, ci_ub, z, pvalue, padj, meta_deg, direction) %>%
  arrange(pvalue)

write_csv(res, snakemake@output$table, na = "")
writeLines(capture.output(sessionInfo()), snakemake@output$session)

md <- res[res$meta_deg, ]
message("Meta DEGs (padj < ", mcfg$padj, ", |estimate| > ", mcfg$lfc, "): ", nrow(md), " (",
        sum(md$direction == "up"), " up / ", sum(md$direction == "down"), " down); padj < ", mcfg$padj,
        " any effect: ", sum(res$padj < mcfg$padj, na.rm = TRUE), "; meta DEGs not DEG in any project: ",
        sum(is.na(md$deg_in)))
