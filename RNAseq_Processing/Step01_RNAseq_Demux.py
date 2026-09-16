#!/usr/bin/env python3
"""
step1_demux.py  —  Demultiplex a pooled CEL-Seq2 run into per-embryo FASTQs.

Accepts both plain (.fastq) and gzip-compressed (.fastq.gz) input.
Output files are named  <sample>/<embryo_name>_R1.fastq  and  _R2.fastq
so they can never collide across samples.

Usage:
    python step1_demux.py --r1 raw/sample_4embryos_R1.fastq.gz \
                          --r2 raw/sample_4embryos_R2.fastq.gz \
                          --barcodes barcodes/4embryos.tsv

The sample name is derived automatically from the R1 filename:
    sample_4embryos_R1.fastq.gz  →  sample name = "sample_4embryos"
    Output directory             →  results/sample_4embryos/demux/

Barcode TSV format (no header, tab-separated):
    embryo1    AGTGTC
    embryo2    ACCATG
    ...
"""

import argparse
import gzip
import os
import re
import sys

# ── helpers ──────────────────────────────────────────────────────────────────

def sample_name_from_path(r1_path: str) -> str:
    base = os.path.basename(r1_path)
    # strip _R1.fastq.gz / _R1.fastq / _1.fastq.gz / _1.fastq
    name = re.sub(r'[_\.]R?1(\.fastq)?(\.gz)?$', '', base, flags=re.IGNORECASE)
    return name

def smart_open(path: str):
    """Return a text-mode file handle regardless of gzip compression."""
    if path.endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path, "r")

def load_barcodes(tsv_path: str) -> dict:
    barcodes = {}
    with open(tsv_path) as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split("\t")
            if len(parts) != 2:
                sys.exit(f"[step1] Bad barcode line (expected 2 columns): {line!r}")
            embryo, bc = parts[0].strip(), parts[1].strip().upper()
            barcodes[bc] = embryo
    return barcodes

def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("--r1",       required=True, help="R1 FASTQ (plain or .gz)")
    p.add_argument("--r2",       required=True, help="R2 FASTQ (plain or .gz)")
    p.add_argument("--barcodes", required=True, help="Barcode TSV (embryo<TAB>barcode)")
    p.add_argument("--outdir",   default=None,
                   help="Output directory (default: results/<sample>/demux)")
    p.add_argument("--mismatches", type=int, default=0,
                   help="Allowed Hamming-distance mismatches in barcode (default: 0)")
    return p.parse_args()

# ── barcode matching (optional 1-mismatch tolerance) ─────────────────────────

def hamming(a: str, b: str) -> int:
    return sum(x != y for x, y in zip(a, b))

def match_barcode(bc: str, barcodes: dict, max_mm: int) -> str | None:
    """Return embryo name for bc, or None if no match within max_mm mismatches."""
    if bc in barcodes:
        return barcodes[bc]
    if max_mm == 0:
        return None
    hits = [(hamming(bc, ref), name) for ref, name in barcodes.items()
            if len(bc) == len(ref)]
    hits = [(d, n) for d, n in hits if d <= max_mm]
    if len(hits) == 1:
        return hits[0][1]
    return None   # 0 hits or ambiguous

# ── main ─────────────────────────────────────────────────────────────────────

def main():
    args = parse_args()
    barcodes   = load_barcodes(args.barcodes)
    bc_len     = len(next(iter(barcodes)))   # assumes uniform barcode length
    umi_len    = 6                           # CEL-Seq2 fixed UMI length

    sample  = sample_name_from_path(args.r1)
    outdir  = args.outdir or os.path.join("results", sample, "demux")
    os.makedirs(outdir, exist_ok=True)

    print(f"[step1] Sample      : {sample}")
    print(f"[step1] Barcodes    : {len(barcodes)} embryos")
    print(f"[step1] Mismatches  : {args.mismatches}")
    print(f"[step1] Output dir  : {outdir}")

    # Open one output handle per embryo
    out_r1 = {name: open(os.path.join(outdir, f"{sample}_{name}_R1.fastq"), "w")
               for name in barcodes.values()}
    out_r2 = {name: open(os.path.join(outdir, f"{sample}_{name}_R2.fastq"), "w")
               for name in barcodes.values()}

    total = assigned = 0

    with smart_open(args.r1) as r1, smart_open(args.r2) as r2:
        while True:
            # ── read one FASTQ record from each file ──────────────────────
            h1 = r1.readline();  h2 = r2.readline()
            if not h1:
                break
            s1 = r1.readline().strip();  s2 = r2.readline()
            p1 = r1.readline();          p2 = r2.readline()
            q1 = r1.readline().strip();  q2 = r2.readline()

            total += 1

            bc   = s1[:bc_len].upper()
            umi  = s1[bc_len:bc_len + umi_len]
            seq  = s1[bc_len + umi_len:]
            qual = q1[bc_len + umi_len:]

            embryo = match_barcode(bc, barcodes, args.mismatches)
            if embryo is None:
                continue

            assigned += 1
            # Add UMI tag to QNAME portion of R2 header (before first space)
            # so STAR preserves it in the BAM read name
            parts = h2.strip().split(" ", 1)
            parts[0] += f"_UMI:{umi}"
            h2_mod = " ".join(parts) + "\n"

            # Write R1 unchanged (not aligned, UMI not needed here)
            out_r1[embryo].write(h1)
            out_r1[embryo].write(seq + "\n")
            out_r1[embryo].write(p1)
            out_r1[embryo].write(qual + "\n")

            # Write modified R2 with UMI tag
            out_r2[embryo].write(h2_mod)
            out_r2[embryo].write(s2)
            out_r2[embryo].write(p2)
            out_r2[embryo].write(q2)

    for fh in list(out_r1.values()) + list(out_r2.values()):
        fh.close()

    pct = 100 * assigned / total if total else 0
    print(f"[step1] Done. Assigned {assigned:,} / {total:,} reads ({pct:.1f}%)")
    if pct < 50:
        print("[step1] WARNING: <50% assignment rate — check barcodes and R1/R2 orientation")

if __name__ == "__main__":
    main()
