#!/usr/bin/env python3
"""
step0_gene_lengths.py  —  One-time setup: generate dm6_gene_lengths.tsv from GTF.

Run once per genome build.  Output is used by step6_tpm.py.

Usage:
    python step0_gene_lengths.py --gtf /path/to/dm6.ensGene.gtf \
                                 --out  /path/to/dm6_gene_lengths.tsv

NOTE: exon lengths are summed naively (overlapping exons from multiple isoforms
are double-counted).  This is acceptable for exploratory TPM; for publication
figures consider using a tool that computes non-overlapping effective length
(e.g. the tximeta/tximport workflow or the featureCounts effective-length output).
"""

import argparse
import pandas as pd
from collections import defaultdict

def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("--gtf", required=True, help="Path to GTF annotation file")
    p.add_argument("--out", required=True, help="Output TSV path")
    return p.parse_args()

def main():
    args = parse_args()

    gene_lengths = defaultdict(int)
    n_exons = 0

    with open(args.gtf) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            fields = line.strip().split("\t")
            if len(fields) < 9 or fields[2] != "exon":
                continue
            start  = int(fields[3])
            end    = int(fields[4])
            length = end - start + 1
            attr   = fields[8]
            try:
                gene_id = [x for x in attr.split(";") if "gene_id" in x][0].split('"')[1]
            except IndexError:
                continue
            gene_lengths[gene_id] += length
            n_exons += 1

    df = pd.DataFrame.from_dict(gene_lengths, orient="index", columns=["Length"])
    df.index.name = "GeneID"
    df.to_csv(args.out, sep="\t")

    print(f"Parsed {n_exons:,} exons across {len(df):,} genes → {args.out}")

if __name__ == "__main__":
    main()
