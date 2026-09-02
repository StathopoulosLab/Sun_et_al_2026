#!/usr/bin/env python3
"""
step3_featurecounts.py  —  Assign reads to genes with featureCounts.

Processes every sorted BAM in results/<sample>/mapped/ and writes
gene-tagged BAMs needed by step4 for UMI collapse.

Usage:
    python step3_featurecounts.py --sample sample_4embryos [--threads 16] [--hpc]

The -R BAM flag is mandatory here: it produces the per-read gene tag (XT)
that step4_umi_collapse.py depends on.  Strandedness is fixed at reverse
(-s 2) which is correct for CEL-Seq2 / CEL-Seq libraries.
"""

import argparse
import glob
import os
import subprocess
import sys

sys.path.insert(0, os.path.dirname(__file__))
import config

def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("--sample",  required=True)
    p.add_argument("--threads", type=int, default=None)
    p.add_argument("--hpc",     action="store_true")
    p.add_argument("--strand",  type=int, default=1,
                   help="featureCounts strandedness: 0=unstranded, 1=forward, "
                        "2=reverse (default 2 for CEL-Seq2)")
    return p.parse_args()

def main():
    args    = parse_args()
    threads = args.threads or (config.THREADS_HPC if args.hpc else config.THREADS_DESKTOP)

    mapped_dir = os.path.join(config.OUTPUT_ROOT, args.sample, "mapped")

    bams = sorted(glob.glob(os.path.join(mapped_dir, f"{args.sample}_*.sorted.bam")))
    if not bams:
        sys.exit(f"[step3] No sorted BAMs found in {mapped_dir}")

    print(f"[step3] Sample   : {args.sample}")
    print(f"[step3] BAMs     : {len(bams)}")
    print(f"[step3] Strand   : {args.strand}")

    for bam in bams:
        fc_bam_out = bam.replace(".sorted.bam", ".sorted.featureCounts.bam")
        if os.path.exists(fc_bam_out):
            print(f"[step3]   Skipping {os.path.basename(bam)} — already processed")
            continue

        fc_counts = bam.replace(".sorted.bam", ".fc.txt")
        print(f"[step3]   featureCounts on {os.path.basename(bam)} …")

        cmd = (
            f"featureCounts "
            f"-a {config.GTF_FILE} "
            f"-o {fc_counts} "
            f"-R BAM "                # REQUIRED: tags each read with its gene (XT)
            f"-T {threads} "
            f"-s {args.strand} "
            f"{bam}"
        )
        subprocess.check_call(cmd, shell=True)

        # featureCounts appends .featureCounts.bam to the input filename
        raw_fc_bam = bam + ".featureCounts.bam"
        os.rename(raw_fc_bam, fc_bam_out)
        print(f"[step3]   → {os.path.basename(fc_bam_out)}")

    print(f"[step3] featureCounts complete for {args.sample}")

if __name__ == "__main__":
    main()
