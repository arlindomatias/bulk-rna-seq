"""HTML report for a DE (DESeq2 or limma) + enrichment project: input → QC → DE → GSEA → STRING.

Run through Snakemake (`script:` directive, rule `report`). Everything specific to a project comes
from `snakemake.params` (title, design text, cutoffs, notes) and the config (`de.method`, organism),
so the same script is used by bulk_rna_seq_tximport, GSE178812 (DESeq2) and GSE236404 (limma).
It only reads plots/ and results/.
"""
import glob
import os
import sys

import pandas as pd

# Shared report module: <repository>/scripts/html_report.py (this script lives in scripts/de/).
sys.path.insert(0, os.path.join(snakemake.scriptdir, ".."))  # noqa: F821
from html_report import Report  # noqa: E402

inp, out, prm = snakemake.input, snakemake.output, snakemake.params  # noqa: F821
P = "plots"
padj_cut, lfc_cut = float(prm.padj), float(prm.lfc)
trt, ref = prm.treatment, prm.reference
cfg = snakemake.config  # noqa: F821
ecfg = cfg.get("enrichment", {})
n_gsea_plot = ecfg.get("gsea_plot_top_n", 3)
limma = cfg.get("de", {}).get("method", "deseq2") == "limma"
kegg_org = ecfg.get("organism_kegg", "mmu")
species = ecfg.get("species_name", "Mus musculus")


def read_csv(path):
    try:
        return pd.read_csv(path)
    except (pd.errors.EmptyDataError, FileNotFoundError):
        return pd.DataFrame()


coldata = read_csv(inp.coldata)
tt = read_csv(inp.top_table)
go = read_csv(inp.go)
kegg = read_csv(inp.kegg)
string = read_csv(inp.string)
nodes = read_csv(inp.nodes)
edges = read_csv(inp.edges)

tested = len(tt)  # every gene in top_table was tested (padj is NA after independent filtering)
# Same rule as deg_criterion() in 00_common.R: padj, or the unadjusted p-value when no gene has padj < cutoff.
fallback = not bool((tt["padj"] < padj_cut).any())
sig_col = "pvalue" if fallback else "padj"
sig_cut = float(cfg.get("de", {}).get("pvalue_fallback", 0.01)) if fallback else padj_cut
deg = tt[(tt[sig_col] < sig_cut) & (tt["log2FoldChange"].abs() > lfc_cut)]
up = deg[deg["log2FoldChange"] > 0]
down = deg[deg["log2FoldChange"] < 0]
padj_only = int((tt["padj"] < padj_cut).sum())


def sig(df, col="p.adjust"):
    return df[df[col] < 0.05] if len(df) and col in df else df.iloc[0:0]


go_sig, kegg_sig = sig(go), sig(kegg)
string_sig = sig(string, "fdr")
groups = coldata["group"].value_counts() if "group" in coldata else pd.Series(dtype=int)

r = Report(prm.title, subtitle=prm.subtitle, code_files=list(inp.code))
deg_kind = f" (unadjusted p < {sig_cut:g})" if fallback else ""
r.kpis([
    (len(coldata), "samples"),
    (f"{tested:,}", "genes tested"),
    (f"{len(up):,}", f"up in {trt}{deg_kind}"),
    (f"{len(down):,}", f"down in {trt}{deg_kind}"),
    (f"{len(go_sig):,} / {len(kegg_sig):,}", "GSEA terms GO / KEGG (padj < 0.05)"),
    (f"{len(nodes):,}", "DEGs mapped to STRING"),
])

# ---------------------------------------------------------------- design
r.h2("Study design")
r.md(prm.design)
r.md("\n".join(f"* **{g}**: {groups.get(g, 0)} samples." for g in (ref, trt)))
if prm.notes:
    r.md(prm.notes)

if limma:
    model_name, rank_stat = "limma-trend", "moderated t statistic"
    weights = " with sample quality weights (`arrayWeights`)" if cfg.get("de", {}).get("array_weights") else ""
    model_text = (f"`lmFit` on log2(normalized counts + 1){weights} and `eBayes(trend = TRUE, robust = TRUE)` "
                  "(the prior variance follows the mean-variance trend)")
else:
    model_name, rank_stat = "DESeq2", "Wald statistic"
    model_text = '`fitType = "local"`'

r.h2("Pipeline")
r.md(prm.pipeline)
r.md(f"""
Shared settings:

* **Filter:** {prm.filter_text}
* **{model_name}:** {model_text}, {ref} as reference. Positive log2FC = higher in **{trt}**.
  DEGs: `padj < {padj_cut}` and `|log2FC| > {lfc_cut}`; if no gene has `padj < {padj_cut}`, the unadjusted
  `pvalue < {cfg.get("de", {}).get("pvalue_fallback", 0.01):g}` is used instead (exploratory).
* **GSEA:** genes ranked by the {rank_stat}; GO (BP, CC, MF) and KEGG (`{kegg_org}`) with clusterProfiler,
  Benjamini-Hochberg adjustment, fixed seed. Positive NES = enriched in {trt}.
* **STRING:** v12.0, *{species}*, combined score ≥ 700, all tested genes as background.
""")

# ---------------------------------------------------------------- input
r.h2("Input data")
r.md(prm.input_text)
if prm.upstream_qc and os.path.exists(prm.upstream_qc):
    up_qc = pd.read_csv(prm.upstream_qc, sep="\t", dtype={"mouse_id": str})
    r.table(up_qc[["run", "condition", "mouse_id", "num_processed", "percent_mapped", "library_type"]],
            "Read pairs and Salmon mapping rate per sample (from the upstream project).", digits=4)
show_cols = [c for c in coldata.columns if c not in ("fastq_1", "fastq_2", "md5_1", "md5_2")]
r.table(coldata[show_cols], "Samples used in the model (`data/colData.csv`).", collapsed=True)

# ---------------------------------------------------------------- QC
r.h2("Quality control")
r.md("""
These plots check the normalized data before testing. Samples should have similar count
distributions, and replicates of the same group should be close in the PCA and distance heatmap.
""")
r.h3("Count distributions")
r.figure(f"{P}/density_counts.png", "Density of normalized log counts per sample.", optional=True)
r.figure(f"{P}/boxplot.png", "Normalized log counts per sample, colored by group.", optional=True)
r.h3("Sample similarity")
if limma:
    r.figure(f"{P}/pca.png", "PCA of log2(normalized counts + 1), 500 most variable genes.",
             width="80%", optional=True)
    r.figure(f"{P}/sample_distances.png", "Euclidean distances between samples (log2 normalized counts).",
             width="75%", optional=True)
else:
    r.figure(f"{P}/pca_vst.png", "PCA of variance-stabilized counts (VST).", width="80%", optional=True)
    r.figure(f"{P}/sample_distances.png", "Euclidean distances between samples (VST).", width="75%", optional=True)
if snakemake.config.get("plots", {}).get("marker_genes"):  # noqa: F821
    r.h3("Marker genes")
    r.figure(f"{P}/marker_genes.png",
             "Normalized counts of sex and cell-type marker genes per sample.", optional=True)
r.h3("Model fit")
if limma:
    r.figure(f"{P}/sample_weights.png", "Sample weights used in the limma fit (arrayWeights; "
             "1 = average quality, lower = less consistent with its group).", width="75%", optional=True)
    r.figure(f"{P}/mean_variance.png", "limma-trend: residual standard deviation (sqrt sigma) vs average "
             "log2 expression, with the fitted prior trend.", width="70%", optional=True)
else:
    r.figure(f"{P}/dispersion_plot.png", "DESeq2 dispersion estimates and fitted trend.", width="70%", optional=True)
    r.figure(f"{P}/cooks_distance.png", "Cook's distance per sample (outlier counts).", optional=True)

# ---------------------------------------------------------------- DE
r.h2("Differential expression")
if fallback:
    n_p = int((tt["pvalue"] < sig_cut).sum())
    r.md(f"""
> **No gene reaches `padj < {padj_cut}`** (smallest padj: {tt["padj"].min():.3g}). The DEGs below, the
> volcano/MA/heatmap colors and the STRING analysis therefore use the **unadjusted p-value**
> (`pvalue < {sig_cut:g}` and `|log2FC| > {lfc_cut}`): **{len(deg):,} genes**, {len(up):,} up and
> {len(down):,} down in {trt}. With {tested:,} tests, about {round(sig_cut * tested):,} genes are expected
> below p = {sig_cut:g} by chance alone ({n_p:,} observed, without the fold-change cutoff), so this list is
> exploratory and contains many false positives. GSEA does not depend on this cutoff.
""")
else:
    r.md(f"""
{tested:,} genes were tested. **{len(deg):,} genes** pass `padj < {padj_cut}` and `|log2FC| > {lfc_cut}`:
**{len(up):,} up** and **{len(down):,} down** in {trt}. Without the fold-change cutoff,
{padj_only:,} genes have `padj < {padj_cut}`.
""")
r.figure(f"{P}/volcano.png", "Volcano plot: log2 fold change vs "
         + ("unadjusted p-value." if fallback else "adjusted p-value."), optional=True)
if limma:
    r.figure(f"{P}/MAplot.png", "MA plot: average log2 expression vs log2 fold change.", optional=True)
    r.figure(f"{P}/heatmap.png", "Top genes by FDR, row-scaled log2 normalized counts.", optional=True)
else:
    r.figure(f"{P}/MAplot.png", "MA plot: mean expression vs shrunken log2 fold change.", optional=True)
    r.figure(f"{P}/heatmap.png", "Top genes by FDR, scaled VST counts.", optional=True)
r.h3("P-value diagnostics")
r.md("A well-behaved test has a flat p-value histogram with a peak near zero, and a QQ plot that "
     "follows the diagonal for most genes.")
r.figure(f"{P}/fdr_histogram.png", "Histogram of p-values.", width="70%", optional=True)
r.figure(f"{P}/QQPlot.png", "QQ plot of p-values.", width="60%", optional=True)

r.h3("Top genes")
gene_col = "SYMBOL"
cols = [c for c in ["gene_id", gene_col, "baseMean", "log2FoldChange", "log2FoldChange_shrunk",
                    "pvalue" if fallback else None, "padj"] if c in tt]
r.table(up.sort_values(sig_col)[cols], f"Top up-regulated genes in {trt} (by {sig_col}).", max_rows=15)
r.table(down.sort_values(sig_col)[cols], f"Top down-regulated genes in {trt} (by {sig_col}).", max_rows=15)
r.md("The full table is `results/top_table.csv`.")

# ---------------------------------------------------------------- GSEA
r.h2("Gene set enrichment (GSEA)")
if len(go):
    n_go_up = int((go_sig["NES"] > 0).sum())
    r.md(f"""
GSEA uses all tested genes, so it also catches pathways where many genes change a little.
**{len(go_sig):,} GO terms** have `padj < 0.05` ({n_go_up} enriched in {trt}, {len(go_sig) - n_go_up} in {ref}).
""")
r.h3("GO")
r.figure(f"{P}/dotplot_GO.png", "Top 10 GO terms per ontology and direction (activated: NES > 0, higher in the treated group; suppressed: NES < 0).", optional=True)
r.figure(f"{P}/ridge.png", "Distribution of fold changes of the core genes of each term.", optional=True)
r.figure(f"{P}/emaplot.png", "Enrichment map: terms that share genes are linked.", optional=True)
r.figure(f"{P}/cnet.png", "Gene-concept network of the top terms.", optional=True)
r.figure(f"{P}/gsea.png", f"Running enrichment score of the {n_gsea_plot} most significant GO terms.", optional=True)
gsea_cols = [c for c in ["ONTOLOGY", "ID", "Description", "setSize", "NES", "p.adjust"] if c in go]


def top_by_direction(df, cols, n=10):
    if not len(df):
        return df
    pos = df[df["NES"] > 0].sort_values("p.adjust").head(n)
    neg = df[df["NES"] < 0].sort_values("p.adjust").head(n)
    return pd.concat([pos, neg])[cols]


if len(go_sig):
    r.table(top_by_direction(go_sig, gsea_cols),
            f"Top 10 GO terms per direction (NES > 0: enriched in {trt}).")

r.h3("KEGG")
if len(kegg):
    n_k_up = int((kegg_sig["NES"] > 0).sum())
    r.md(f"**{len(kegg_sig):,} KEGG pathways** have `padj < 0.05` ({n_k_up} enriched in {trt}, "
         f"{len(kegg_sig) - n_k_up} in {ref}).")
r.figure(f"{P}/dotplot_kegg.png", "Top KEGG pathways, split by direction.", optional=True)
r.figure(f"{P}/ridge_kegg.png", "Fold-change distribution of the core genes per pathway.", optional=True)
r.figure(f"{P}/ema_kegg.png", "Enrichment map of KEGG pathways.", optional=True)
r.figure(f"{P}/cnet_kegg.png", "Gene-concept network of the top KEGG pathways.", optional=True)
r.figure(f"{P}/gsea_kegg.png", f"Running enrichment score of the {n_gsea_plot} most significant KEGG pathways.", optional=True)
if len(kegg_sig):
    r.table(top_by_direction(kegg_sig, [c for c in ["ID", "Description", "setSize", "NES", "p.adjust"] if c in kegg]),
            f"Top 10 KEGG pathways per direction (NES > 0: enriched in {trt}).")

pv = sorted(glob.glob(f"{P}/pathview/*.pathview.png") + glob.glob(f"{P}/*.pathview.png"))
if pv:
    r.h3("Pathway maps (pathview)")
    r.md("Genes are colored by log2 fold change (red = higher in the treated group, green = lower).")
    names = dict(zip(kegg["ID"], kegg["Description"])) if len(kegg) else {}
    for f in pv:
        pid = os.path.basename(f).split(".")[0]
        r.figure(f, f"KEGG {pid}: {names.get(pid, '')}.")

# ---------------------------------------------------------------- STRING
r.h2("Protein-protein interactions (STRING)")
linked = len(set(edges["from"]).union(edges["to"])) if len(edges) else 0
r.md(f"""
The {len(deg):,} DEGs{deg_kind} were mapped to STRING: **{len(nodes):,}** proteins were found, **{linked:,}** of them with at
least one interaction (**{len(edges):,}** edges, score ≥ 700). Edge and node tables for Cytoscape are
`results/edges_PPI.csv` and `results/nodes_PPI.csv`.
""")
r.figure(f"{P}/rede_PPI.png", "STRING network of the DEGs (image from the STRING server).", width="85%", optional=True)
r.h3("STRING functional enrichment")
r.md("Over-representation of the DEGs against all tested genes (FDR < 0.05).")
r.figure(f"{P}/enrichment_combined.png", "Top terms per category (−log10 FDR).", optional=True)
r.figure(f"{P}/reactome.png", "Reactome pathways.", optional=True, code=["plot_enrichment", "reactome"])
r.figure(f"{P}/wikipathways.png", "WikiPathways.", optional=True, code=["plot_enrichment", "reactome"])
if len(string_sig):
    cnt = string_sig["category"].value_counts().rename_axis("category").reset_index(name="terms (FDR < 0.05)")
    r.table(cnt, "Significant STRING terms per category.")
    keep = ["Process", "KEGG", "RCTM", "WikiPathways", "MPO"]
    top_s = (string_sig[string_sig["category"].isin(keep)].sort_values("fdr")
             .groupby("category", sort=False).head(5))
    r.table(top_s[["category", "term", "description", "number_of_genes", "fdr"]],
            "Top 5 terms per category.")

# ---------------------------------------------------------------- appendix
r.appendix(
    outputs=[
        ("results/top_table.csv", f"{model_name} results for every tested gene"),
        ("results/GSEA_GO_all_full.csv", "Full GSEA table, GO (BP, CC, MF)"),
        ("results/GSEA_KEGG_full.csv", "Full GSEA table, KEGG"),
        ("results/enrichment_STRING_full.csv", "STRING enrichment of the DEGs"),
        ("results/edges_PPI.csv", "PPI edges for Cytoscape"),
        ("results/nodes_PPI.csv", "PPI node attributes (log2FC, padj)"),
        ("data/fit.rds", "Fitted limma model (MArrayLM)") if limma else ("data/dds.rds", "Fitted DESeq2 object"),
    ],
    env_files=sorted(inp.envs),
    config_files=list(inp.config),
    session_file=inp.session,
)

r.write(out.html)
