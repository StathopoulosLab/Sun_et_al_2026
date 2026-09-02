#!/usr/bin/env python3
"""
step4_umi_collapse.py  —  Collapse UMIs per gene to produce per-embryo count tables.

Reads the gene-tagged BAMs from step3, extracts the XT tag (gene) and the
UMI encoded in the read name, and writes one count file per embryo:
    results/<sample>/counts/<sample>_<embryo>_counts.txt

Usage:
    python step4_umi_collapse.py --sample sample_4embryos

Note on UMI collapse strategy:
    Unique UMI sequences per gene are counted (i.e. set cardinality).
    This corrects for PCR amplification bias but does NOT use read-position
    information.  It is equivalent to the simple UMI counting approach and
    is appropriate for most CEL-Seq2 datasets.  For very highly expressed
    genes at high sequencing depth, consider umi-tools dedup for position-
    aware correction.
"""

import argparse
import glob
import os
import sys
from collections import defaultdict

try:
    import pysam
except ImportError:
    sys.exit("[step4] pysam is required: pip install pysam")

sys.path.insert(0, os.path.dirname(__file__))
import config

def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("--sample", required=True)
    return p.parse_args()

def collapse_one(bam_path: str, counts_path: str):
    bam = pysam.AlignmentFile(bam_path, "rb")
    gene_umis: dict[str, set] = defaultdict(set)
    n_reads = 0

    for read in bam:
        if read.is_unmapped or read.is_secondary or read.is_supplementary:
            continue
        if not read.has_tag("XT"):
            continue                        # not assigned to a gene
        gene = read.get_tag("XT")
        # UMI is appended to read name as "_UMI:XXXXXX" during demux
        try:
            umi = read.query_name.split("_UMI:")[-1]
            if not umi:
                continue
        except Exception:
            continue
        gene_umis[gene].add(umi)
        n_reads += 1

    bam.close()

    with open(counts_path, "w") as out:
        for gene, umis in sorted(gene_umis.items()):
            out.write(f"{gene}\t{len(umis)}\n")

    return n_reads, len(gene_umis)

def main():
    args       = parse_args()
    mapped_dir = os.path.join(config.OUTPUT_ROOT, args.sample, "mapped")
    counts_dir = os.path.join(config.OUTPUT_ROOT, args.sample, "counts")
    os.makedirs(counts_dir, exist_ok=True)

    fc_bams = sorted(glob.glob(
        os.path.join(mapped_dir, f"{args.sample}_*.sorted.featureCounts.bam")))
    if not fc_bams:
        sys.exit(f"[step4] No featureCounts BAMs found in {mapped_dir}\n"
                  "        Did step3_featurecounts.py run successfully?")

    print(f"[step4] Sample  : {args.sample}")
    print(f"[step4] BAMs    : {len(fc_bams)}")

    for bam in fc_bams:
        base   = os.path.basename(bam)   # sample_4embryos_embryo1.sorted.featureCounts.bam
        embryo = (base
                  .replace(f"{args.sample}_", "")
                  .replace(".sorted.featureCounts.bam", ""))
        out_path = os.path.join(counts_dir, f"{args.sample}_{embryo}_counts.txt")

        if os.path.exists(out_path):
            print(f"[step4]   Skipping {embryo} — counts file already exists")
            continue

        print(f"[step4]   Collapsing {embryo} …")
        n_reads, n_genes = collapse_one(bam, out_path)
        print(f"[step4]   {n_reads:,} assigned reads → {n_genes:,} genes → {os.path.basename(out_path)}")

    print(f"[step4] UMI collapse complete for {args.sample}")

if __name__ == "__main__":
    main()
