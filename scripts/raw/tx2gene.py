"""Builds tx2gene from the headers of the GENCODE transcript FASTA.

Header: >transcript_id|gene_id|havana_gene|havana_tx|tx_name|gene_name|length|biotype|
With `salmon index --gencode`, the transcript name in quant.sf is the first field.
Usage: python3 tx2gene.py gencode.transcripts.fa.gz tx2gene.tsv
"""
import gzip
import sys

fasta, out = sys.argv[1], sys.argv[2]
n = 0
with gzip.open(fasta, "rt") as fin, open(out, "w") as fout:
    fout.write("transcript_id\tgene_id\tgene_name\ttranscript_biotype\n")
    for line in fin:
        if line.startswith(">"):
            f = line[1:].rstrip("\n").split("|")
            fout.write(f"{f[0]}\t{f[1]}\t{f[5]}\t{f[7]}\n")
            n += 1
print(f"{n} transcripts written to {out}")
