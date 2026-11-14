"""HTML report for raw_reads_processing: from SRA runs to the gene count table.

Run through Snakemake (`script:` directive, rule `report`). Uses the figures in plots/ and the
tables in results/; it does not recompute anything heavy.
"""
import os
import sys

import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from html_report import Report  # noqa: E402

inp, out, cfg = snakemake.input, snakemake.output, snakemake.config  # noqa: F821
P = os.path.dirname(inp.reads_plot)

samples = pd.read_csv(inp.samples, sep="\t", dtype={"mouse_id": str})
qc = pd.read_csv(inp.qc, sep="\t", dtype={"mouse_id": str})
counts = pd.read_csv(inp.counts, sep="\t", index_col="gene_id")
comp = pd.read_csv(inp.composition, sep="\t", index_col=0)
fq_status = pd.read_csv(inp.fastqc_status, sep="\t")
top = pd.read_csv(inp.top_genes, sep="\t")
corr = pd.read_csv(inp.correlation, sep="\t", index_col=0)

n = len(samples)
groups = samples["condition"].value_counts()
mat = counts.drop(columns="gene_name")
detected_all = int((mat >= 10).all(axis=1).sum())
off_diag = corr.where(~pd.DataFrame(
    [[i == j for j in range(len(corr))] for i in range(len(corr))], index=corr.index, columns=corr.columns))
rate_min, rate_max = qc["percent_mapped"].min(), qc["percent_mapped"].max()
salmon_versions = ", ".join(sorted(qc["salmon_version"].astype(str).unique()))
lib_types = ", ".join(sorted(qc["library_type"].unique()))

r = Report(f"{cfg['data_subdir']}: from raw reads to a gene count table",
           subtitle="Bulk RNA-seq · SRA → FastQC → Salmon → tximport", code_files=list(inp.code))

r.kpis([
    (n, "samples"),
    (f"{qc['num_processed'].sum() / 1e6:,.0f} M", "read pairs processed"),
    (f"{rate_min:.1f}–{rate_max:.1f}%", "Salmon mapping rate"),
    (f"{len(mat):,}", "genes in the count table"),
    (f"{detected_all:,}", "genes with ≥ 10 counts in every sample"),
])

# ---------------------------------------------------------------- design
r.h2("Study design")
r.md(f"""
GEO **{cfg['data_subdir']}** (BioProject PRJNA1254109): hippocampus of 12-week-old female C57BL/6J mice,
one month after AC-T chemotherapy (doxorubicin + cyclophosphamide, then paclitaxel) or saline.
Paired-end RNA-seq on an Illumina NextSeq 2000.

* **Saline**: {groups.get('Saline', 0)} mice.
* **AC-T**: {groups.get('AC-T', 0)} mice.

This project stops at the count table. Differential expression is done in
[`bulk_rna_seq_tximport`](../../bulk_rna_seq_tximport/reports/report.html).
""")
r.table(samples[["run", "sample", "condition", "mouse_id", "read_pairs"]]
        .rename(columns={"read_pairs": "read pairs (ENA)"}), "Samples (`data/samples.tsv`).")

r.h2("Pipeline")
r.md(f"""
1. **Download.** `.sra` files come from the NCBI mirror on AWS and are checked with `vdb-validate`.
   `fasterq-dump` writes temporary FASTQ files to a local SSD. They are deleted after FastQC and Salmon.
2. **Read QC.** FastQC runs on each read, and MultiQC collects FastQC and Salmon
   ([MultiQC report](multiqc_report.html)).
3. **Reference.** GENCODE mouse **{cfg['gencode_release']}** transcriptome (GRCm39). The Salmon index uses
   `-k 31 --gencode`, transcriptome only (no decoys).
4. **Quantification.** `salmon quant -l A --gcBias --seqBias`. A sample fails if its mapping rate is below
   {cfg['min_mapping_rate']}%. Every sample uses the same Salmon version and index.
5. **Gene level.** tximport summarizes transcripts to genes and writes `counts_gene.tsv`, `tpm_gene.tsv`,
   `counts_gene_lengthScaledTPM.tsv`, `length_gene.tsv` and `data/txi_gene.rds` (input of the DE project).

Large files (SRA, reference, index, Salmon outputs) are stored under `$DATA/{cfg['data_subdir']}/`.
""")

# ---------------------------------------------------------------- FastQC
r.h2("Raw read quality (FastQC)")
fail = (fq_status.query("status == 'FAIL'").groupby("module").size()
        .reindex(fq_status["module"].unique()).fillna(0).astype(int))
r.md("""
FastQC flags compare each file with a generic genomic library. Several flags are expected for RNA-seq
and do not mean the data are poor.
""")
r.figure(os.path.join(P, "fastqc_module_status.png"),
         "FastQC module status per sample and read (P = pass, W = warning, F = fail).")
r.md(f"""
* **Sequence duplication** fails in {fail.get('Sequence Duplication Levels', 0)} of {len(fq_status['run'].unique()) * 2}
  files. This is expected in RNA-seq: highly expressed transcripts produce many identical reads.
* **Per base sequence content** warns or fails at the read start. Random-hexamer priming causes this
  bias at the first bases of every RNA-seq library.
* **Per tile sequence quality** fails in every file. NextSeq flow cells often trigger this module.
* **Per base sequence quality** fails in {fail.get('Per base sequence quality', 0)} files. In the plot below,
  mean quality stays above Q30 along the read and drops only at the last cycle of read 2. Salmon's selective
  alignment is not affected by one low-quality base, so no trimming was applied.
* **GC content**, **N content** and **overrepresented sequences** fail mostly in read 2.
""")
r.h3("Base quality")
r.figure(os.path.join(P, "fastqc_per_base_quality.png"),
         "Mean Phred quality per position. The shaded band marks good quality (≥ 28). One line per sample.")
r.h3("GC content")
r.figure(os.path.join(P, "fastqc_gc_content.png"),
         "Distribution of mean GC content per read.")
r.h3("Adapter content")
r.figure(os.path.join(P, "fastqc_adapter_content.png"),
         "Share of reads with an adapter sequence at each position (highest adapter type).")

# ---------------------------------------------------------------- Salmon
r.h2("Quantification (Salmon)")
same_pairs = int((qc["num_processed"] == qc["read_pairs_ena"]).sum())
r.md(f"""
All samples were quantified with Salmon **{salmon_versions}**. The library type detected automatically was
**{lib_types}** in every sample (stranded, read 1 from the reverse strand). Mapping rates ranged from
**{rate_min:.1f}%** to **{rate_max:.1f}%**. The number of pairs Salmon processed matches the count
reported by ENA in {"every sample" if same_pairs == n else f"{same_pairs} of {n} samples"}.
""")
r.figure(os.path.join(P, "qc_reads_mapping.png"),
         "Read pairs per sample, split into pairs mapped by Salmon and pairs not mapped. Labels give the mapping rate.")
r.table(qc[["run", "condition", "mouse_id", "read_pairs_ena", "num_processed", "num_mapped",
            "percent_mapped", "library_type", "compatible_fragment_ratio"]],
        "Salmon QC (`results/qc_summary.tsv`).", digits=4)

# ---------------------------------------------------------------- count table
r.h2("Gene count table")
pc = comp["Protein coding"]
mt = comp["Mitochondrial"]
r.md(f"""
tximport gives **{len(mat):,} genes × {n} samples**. Of these, **{detected_all:,}** genes have at least 10
counts in every sample.
""")
r.h3("Read composition")
r.figure(os.path.join(P, "count_composition.png"),
         "Share of assigned reads by gene category (from GENCODE transcript biotypes; `mt-` genes are grouped as mitochondrial).")
r.md(f"""
Protein-coding genes take {pc.min():.0f}–{pc.max():.0f}% of the reads and mitochondrial genes
{mt.min():.0f}–{mt.max():.0f}%. rRNA genes take at most {comp["rRNA"].max():.2f}% of the reads.
""")
r.table(comp.round(2).reset_index().rename(columns={"index": "sample"}),
        "Share of assigned reads (%) by gene category (`results/count_composition.tsv`).", collapsed=True)
r.h3("Detected and most abundant genes")
r.figure(os.path.join(P, "detected_genes.png"), "Genes with at least 10 counts per sample.", width="55%")
r.figure(os.path.join(P, "top_expressed_genes.png"),
         "The 20 genes with the highest mean share of reads.", width="75%")
n_mt_top = int(top["gene_name"].str.startswith("mt-").sum())
other_top = [g for g in top["gene_name"] if not g.startswith("mt-")][:6]
r.md(f"""
{n_mt_top} of the 20 most abundant genes are mitochondrial. The others include
{', '.join('*' + g + '*' for g in other_top)}.
""")
r.h3("Sample similarity")
r.figure(os.path.join(P, "sample_correlation.png"),
         "Pearson correlation of log2 CPM between samples (genes with ≥ 10 counts in every sample).")
r.figure(os.path.join(P, "pca_logcpm.png"),
         "PCA of log2 CPM, 1,000 most variable genes. Labels are mouse IDs.", width="75%")
r.md(f"""
Correlation between samples ranges from **{off_diag.min().min():.3f}** to **{off_diag.max().max():.3f}**.
No sample is a clear outlier. The PCA here is only an overview: the DE project repeats it on
variance-stabilized counts after filtering.
""")

# ---------------------------------------------------------------- appendix
r.appendix(
    outputs=[
        ("results/counts_gene.tsv", "Estimated counts per gene (for DESeq2 use `data/txi_gene.rds`)"),
        ("results/counts_gene_lengthScaledTPM.tsv", "Length-scaled TPM counts (edgeR/limma)"),
        ("results/tpm_gene.tsv", "TPM per gene (not for DE)"),
        ("results/length_gene.tsv", "Effective gene lengths"),
        ("results/qc_summary.tsv", "Reads, mapping rate and library type per sample"),
        ("data/txi_gene.rds", "tximport object for `DESeqDataSetFromTximport`"),
        ("data/tx2gene.tsv", "Transcript → gene table from GENCODE"),
        ("reports/multiqc_report.html", "MultiQC (FastQC + Salmon)"),
    ],
    env_files=sorted(inp.envs),
    config_files=[inp.config],
    session_file=inp.session,
)

r.write(out.html)
