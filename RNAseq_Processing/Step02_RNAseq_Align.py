#!/usr/bin/env python3
"""
step2_align.py  —  Align all per-embryo R2 FASTQs for one sample using STAR.

Reads the demux output from step1 and writes sorted BAMs under
results/<sample>/mapped/.

Usage:
    python step2_align.py --sample sample_4embryos [--threads 16] [--hpc]

  --hpc flag raises thread count to THREADS_HPC from config.py;
  omit it for desktop runs (uses THREADS_DESKTOP).
"""

import argparse
import glob
import os
import shutil
import subprocess
import sys

sys.path.insert(0, os.path.dirname(__file__))
import config

def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("--sample",  required=True, help="Sample name (e.g. sample_4embryos)")
    p.add_argument("--threads", type=int, default=None,
                   help="Override thread count from config")
    p.add_argument("--hpc",     action="store_true",
                   help="Use THREADS_HPC instead of THREADS_DESKTOP")
    return p.parse_args()

def main():
    args    = parse_args()
    threads = args.threads or (config.THREADS_HPC if args.hpc else config.THREADS_DESKTOP)

    demux_dir  = os.path.join(config.OUTPUT_ROOT, args.sample, "demux")
    mapped_dir = os.path.join(config.OUTPUT_ROOT, args.sample, "mapped")
    os.makedirs(mapped_dir, exist_ok=True)

    # Find all R2 FASTQs for this sample (one per embryo)
    r2_files = sorted(glob.glob(os.path.join(demux_dir, f"{args.sample}_*_R2.fastq")))
    if not r2_files:
        sys.exit(f"[step2] No demuxed R2 FASTQs found in {demux_dir}\n"
                  "        Did step1_demux.py run successfully?")

    print(f"[step2] Sample     : {args.sample}")
    print(f"[step2] Threads    : {threads}")
    print(f"[step2] Embryos    : {len(r2_files)}")

    for r2 in r2_files:
        # embryo name sits between <sample>_ and _R2.fastq
        base   = os.path.basename(r2)                 # sample_4embryos_embryo1_R2.fastq
        embryo = base.replace(f"{args.sample}_", "").replace("_R2.fastq", "")
        prefix = os.path.join(mapped_dir, f"{args.sample}_{embryo}_")

        bam_src = f"{prefix}Aligned.sortedByCoord.out.bam"
        bam_dst = os.path.join(mapped_dir, f"{args.sample}_{embryo}.sorted.bam")

        if os.path.exists(bam_dst):
            print(f"[step2]   Skipping {embryo} — BAM already exists")
            continue

        # If a previous run was interrupted mid-alignment for this embryo
        # (e.g. disk-full crash), leftover partial files under this prefix
        # can make STAR fail on retry rather than cleanly overwrite. Clear
        # them out before trying again.
        stale = glob.glob(f"{prefix}*")
        if stale:
            print(f"[step2]   Found {len(stale)} leftover file(s) from an "
                  f"interrupted run for {embryo} — cleaning up before retry")
            for path in stale:
                if os.path.isdir(path):
                    shutil.rmtree(path)
                else:
                    os.remove(path)

        print(f"[step2]   Aligning {embryo} …")
        cmd = (
            f"STAR --runThreadN {threads} "
            f"--genomeDir {config.STAR_INDEX} "
            f"--readFilesIn {r2} "
            f"--outSAMtype BAM SortedByCoordinate "
            f"--outFileNamePrefix {prefix} "
            f"--outFilterMultimapNmax 1 "
            f"--outSAMattributes NH HI AS nM MD "
        )
        subprocess.check_call(cmd, shell=True)
        os.rename(bam_src, bam_dst)
        print(f"[step2]   Done → {bam_dst}")

    print(f"[step2] Alignment complete for {args.sample}")

if __name__ == "__main__":
    main()
