"""HTML report for a DESeq2 + enrichment project: input → QC → DE → GSEA → STRING.

Run through Snakemake (`script:` directive, rule `report`). Everything specific to a project comes
from `snakemake.params` (title, design text, cutoffs, notes), so the same script is used by
bulk_rna_seq_tximport and GSE178812. It only reads plots/ and results/.
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
n_gsea_plot = snakemake.config.get("enrichment", {}).get("gsea_plot_top_n", 3)  # noqa: F821


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
deg = tt[(tt["padj"] < padj_cut) & (tt["log2FoldChange"].abs() > lfc_cut)]
up = deg[deg["log2FoldChange"] > 0]
down = deg[deg["log2FoldChange"] < 0]
padj_only = int((tt["padj"] < padj_cut).sum())


def sig(df, col="p.adjust"):
    return df[df[col] < 0.05] if len(df) and col in df else df.iloc[0:0]


go_sig, kegg_sig = sig(go), sig(kegg)
string_sig = sig(string, "fdr")
groups = coldata["group"].value_counts() if "group" in coldata else pd.Series(dtype=int)

r = Report(prm.title, subtitle=prm.subtitle, code_files=list(inp.code))
r.kpis([
    (len(coldata), "samples"),
    (f"{tested:,}", "genes tested"),
    (f"{len(up):,}", f"up in {trt}"),
    (f"{len(down):,}", f"down in {trt}"),
    (f"{len(go_sig):,} / {len(kegg_sig):,}", "GSEA terms GO / KEGG (padj < 0.05)"),
    (f"{len(nodes):,}", "DEGs mapped to STRING"),
])

# ---------------------------------------------------------------- design
r.h2("Study design")
r.md(prm.design)
r.md("\n".join(f"* **{g}**: {groups.get(g, 0)} samples." for g in (ref, trt)))
if prm.notes:
    r.md(prm.notes)

r.h2("Pipeline")
r.md(prm.pipeline)
r.md(f"""
Shared settings:

* **Filter:** {prm.filter_text}
* **DESeq2:** `fitType = "local"`, {ref} as reference. Positive log2FC = higher in **{trt}**.
  DEGs: `padj < {padj_cut}` and `|log2FC| > {lfc_cut}`.
* **GSEA:** genes ranked by the Wald statistic; GO (BP, CC, MF) and KEGG (`mmu`) with clusterProfiler,
  Benjamini-Hochberg adjustment, fixed seed. Positive NES = enriched in {trt}.
* **STRING:** v12.0, *Mus musculus*, combined score ≥ 700, all tested genes as background.
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
r.figure(f"{P}/pca_vst.png", "PCA of variance-stabilized counts (VST).", width="80%", optional=True)
r.figure(f"{P}/sample_distances.png", "Euclidean distances between samples (VST).", width="75%", optional=True)
if snakemake.config.get("plots", {}).get("marker_genes"):  # noqa: F821
    r.h3("Marker genes")
    r.figure(f"{P}/marker_genes.png",
             "Normalized counts of sex and cell-type marker genes per sample.", optional=True)
r.h3("Model fit")
r.figure(f"{P}/dispersion_plot.png", "DESeq2 dispersion estimates and fitted trend.", width="70%", optional=True)
r.figure(f"{P}/cooks_distance.png", "Cook's distance per sample (outlier counts).", optional=True)

# ---------------------------------------------------------------- DE
r.h2("Differential expression")
r.md(f"""
{tested:,} genes were tested. **{len(deg):,} genes** pass `padj < {padj_cut}` and `|log2FC| > {lfc_cut}`:
**{len(up):,} up** and **{len(down):,} down** in {trt}. Without the fold-change cutoff,
{padj_only:,} genes have `padj < {padj_cut}`.
""")
r.figure(f"{P}/volcano.png", "Volcano plot: log2 fold change vs adjusted p-value.", optional=True)
r.figure(f"{P}/MAplot.png", "MA plot: mean expression vs shrunken log2 fold change.", optional=True)
r.figure(f"{P}/heatmap.png", "Top genes by FDR, scaled VST counts.", optional=True)
r.h3("P-value diagnostics")
r.md("A well-behaved test has a flat p-value histogram with a peak near zero, and a QQ plot that "
     "follows the diagonal for most genes.")
r.figure(f"{P}/fdr_histogram.png", "Histogram of p-values.", width="70%", optional=True)
r.figure(f"{P}/QQPlot.png", "QQ plot of p-values.", width="60%", optional=True)

r.h3("Top genes")
gene_col = "SYMBOL"
cols = [c for c in ["gene_id", gene_col, "baseMean", "log2FoldChange", "log2FoldChange_shrunk", "padj"] if c in tt]
r.table(up.sort_values("padj")[cols], f"Top up-regulated genes in {trt} (by padj).", max_rows=15)
r.table(down.sort_values("padj")[cols], f"Top down-regulated genes in {trt} (by padj).", max_rows=15)
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
The {len(deg):,} DEGs were mapped to STRING: **{len(nodes):,}** proteins were found, **{linked:,}** of them with at
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
        ("results/top_table.csv", "DESeq2 results for every tested gene"),
        ("results/GSEA_GO_all_full.csv", "Full GSEA table, GO (BP, CC, MF)"),
        ("results/GSEA_KEGG_full.csv", "Full GSEA table, KEGG"),
        ("results/enrichment_STRING_full.csv", "STRING enrichment of the DEGs"),
        ("results/edges_PPI.csv", "PPI edges for Cytoscape"),
        ("results/nodes_PPI.csv", "PPI node attributes (log2FC, padj)"),
        ("data/dds.rds", "Fitted DESeq2 object"),
    ],
    env_files=sorted(inp.envs),
    config_files=list(inp.config),
    session_file=inp.session,
)

r.write(out.html)
