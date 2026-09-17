################################################################################
# generate_toy_dar_cache.R
#
# Generates a small synthetic peak universe + pre-populated counts_cache/,
# structured EXACTLY as 01_dar_calling_seven_genotype_nc14b.r expects, so a
# reviewer can run the real script's full statistical pipeline (DESeq2 +
# limma + edgeR + RUVg, all 13 contrasts) WITHOUT needing real BAM files.
#
# HOW THIS WORKS:
#   01_dar_calling_seven_genotype_nc14b.r always checks its cache directory
#   before ever touching a BAM file. Every count matrix it needs (narrow
#   peaks, per-group peaks, boundary-extended peaks, per-group extended
#   peaks) has a cache-hit path that loads pre-saved data and a cache-miss
#   fallback that calls featureCounts() on the real BAMs. This script
#   pre-populates every one of those cache files, so the real script takes
#   the cache-hit path everywhere and the BAM list at the top of the script
#   (which points at files that don't exist on your machine) is never
#   actually opened.
#
# WHAT YOU GET: real DESeq2/limma/edgeR/RUVg statistics run on synthetic
# data with deliberately engineered signal (see "SIMULATED BIOLOGY" below),
# not just placeholder zeros — so the contrasts should produce a realistic
# mix of DARs and non-DARs, exercising the tiering/confidence-score logic.
#
# USAGE:
#   Place this script inside ATACseq_Processing/02_DAR_Calling/ (alongside
#   01_dar_calling_seven_genotype_nc14b.r itself) and run:
#     Rscript generate_toy_dar_cache.R
#   Then run the real script from the same directory, unmodified:
#     Rscript 01_dar_calling_seven_genotype_nc14b.r
################################################################################

suppressPackageStartupMessages({
  library(GenomicRanges)
})

set.seed(42)

cat("Generating toy peak universe + cache for DAR-calling test...\n\n")

# ---------------------------------------------------------------------------
# Must exactly match the real script's bamfiles/sample_names/coldata (lines
# 121-151 of 01_dar_calling_seven_genotype_nc14b.r).
#
# IMPORTANT: narrow_counts.txt is plain text, loaded via read.table(header=TRUE)
# with default check.names=TRUE -- R silently converts "-" to "." in column
# headers on READ, regardless of how the file was written. sample_names below
# contains real hyphens (e.g. "FoxL1-High_..."), so writing THOSE directly as
# narrow_counts.txt's header would round-trip back mismatched. This is exactly
# why CACHE_KEY_MAP exists in the real script: short, hyphen-free keys (e.g.
# "BOTv_rep1") survive read.table's mangling untouched, then get translated
# back to full BAM names via CACHE_KEY_MAP after loading. So narrow_counts.txt
# below is written with THESE short keys as column headers, not sample_names
# directly -- the RDS-based caches (saveRDS/readRDS) don't have this problem,
# since they preserve strings exactly, so they still use sample_names as-is.
# ---------------------------------------------------------------------------
sample_names <- c(
  "FoxL1-High_Dm_ATAC_Nc14b_rep3_noq_rmdup.noChrM.bam",
  "FoxL1-High_Dm_ATAC_Nc14b_07_noq_rmdup.noChrM.bam",
  "HLH54F-High_Dm_ATAC_Nc14b_rep4_noq_rmdup.noChrM.bam",
  "HLH54F-High_Dm_ATAC_Nc14b_rep5_noq_rmdup.noChrM.bam",
  "BOT_D7_1_nc14b_IR_noq_rmdup.noChrM.bam",
  "BOT_D7_2_nc14b_IR_noq_rmdup.noChrM.bam",
  "Run-D7_Dm_ATAC_Nc14b_01_noq_rmdup.noChrM.bam",
  "Run-D7_Dm_ATAC_Nc14b_02_noq_rmdup.noChrM.bam",
  "Run-D7_Dm_ATAC_Nc14b_03_noq_rmdup.noChrM.bam",
  "Run-D7_Dm_ATAC_Nc14b_07_noq_rmdup.noChrM.bam",
  "HLH_B6_Nc14b_01_noq_rmdup.noChrM.bam",
  "HLH_B6_Nc14b_02_noq_rmdup.noChrM.bam",
  "Mat_Run_B6_3_nc14b_IR_noq_rmdup.noChrM.bam",
  "Mat_Run_B6_5_nc14b_IR_noq_rmdup.noChrM.bam"
)
genotype <- c("BOTv","BOTv","BOTCv","BOTCv","BOT","BOT",
              "BOT_hR","BOT_hR","BOTR","BOTR","BOTC","BOTC","BOTC_oR","BOTC_oR")
names(genotype) <- sample_names
groups <- c("BOTv","BOTCv","BOT","BOT_hR","BOTR","BOTC","BOTC_oR")

# Short, hyphen-free cache keys -- must match CACHE_KEY_MAP's names in the
# real script exactly (same order as sample_names above).
cache_short_keys <- c(
  "BOTv_rep1", "BOTv_rep2", "BOTCv_rep1", "BOTCv_rep2",
  "BOT_nc14b_rep1", "BOT_nc14b_rep2", "BOT_hR_rep1", "BOT_hR_rep2",
  "BOTR_rep1", "BOTR_rep2", "BOTC_rep1", "BOTC_rep2",
  "BOTC_oR_rep1", "BOTC_oR_rep2"
)

# ---------------------------------------------------------------------------
# PEAK UNIVERSE: 300 peaks, 200bp wide, spaced 2000bp apart on chr2L.
# Spacing chosen so boundary-extension (±200bp, ±300bp) never causes
# adjacent peaks to merge -- keeps the extended-peak logic simple and 1:1
# with the narrow peaks for this toy dataset.
# ---------------------------------------------------------------------------
N_PEAKS <- 300
CHROM   <- "chr2L"
WIDTH   <- 200
SPACING <- 2000
BASE    <- 1000000

peak_names <- sprintf("peak_%04d", seq_len(N_PEAKS))
starts     <- BASE + (seq_len(N_PEAKS) - 1) * SPACING
ends       <- starts + WIDTH

peaks_gr <- GRanges(seqnames = CHROM, ranges = IRanges(start = starts, end = ends))
names(peaks_gr) <- peak_names

# ---------------------------------------------------------------------------
# SIMULATED BIOLOGY: most peaks are background noise; a few blocks carry
# deliberate group-specific signal so the contrasts have real DARs to find
# (not just an all-zero / all-identical matrix).
#   peaks 1-15   "ventral_specific"    high in BOTv + BOTCv only
#   peaks 16-30  "runt_dosage_graded"  BOTR(low) < BOT_hR < BOT < BOTv(high);
#                                      BOTC/BOTC_oR held at an intermediate baseline
#   peaks 31-45  "cic_effect"          high when Cic intact (BOTv/BOT/BOT_hR/BOTR),
#                                      low when Cic deleted (BOTCv/BOTC/BOTC_oR)
#   peaks 46-300 background            same mean in every genotype (noise only)
# ---------------------------------------------------------------------------
mean_for <- function(peak_idx, geno) {
  if (peak_idx <= 15) {                      # ventral_specific
    if (geno %in% c("BOTv","BOTCv")) return(150) else return(40)
  } else if (peak_idx <= 30) {                # runt_dosage_graded
    switch(geno,
      BOTR    = 20,  BOT_hR = 50, BOT  = 80, BOTv = 120,
      BOTCv   = 55,  BOTC   = 55, BOTC_oR = 55)
  } else if (peak_idx <= 45) {                # cic_effect (Cic intact vs deleted)
    if (geno %in% c("BOTv","BOT","BOT_hR","BOTR")) return(130) else return(30)
  } else {
    return(40)                                 # background
  }
}

cat("Simulating narrow peak counts (negative binomial, size=8)...\n")
counts_narrow <- matrix(0L, nrow = N_PEAKS, ncol = length(sample_names),
                         dimnames = list(peak_names, sample_names))
for (i in seq_len(N_PEAKS)) {
  for (s in sample_names) {
    mu <- mean_for(i, genotype[[s]])
    counts_narrow[i, s] <- rnbinom(1, mu = mu, size = 8)
  }
}

dir.create("./ATAC_NarrowPeaks", showWarnings = FALSE)
dir.create("./counts_cache/SevenGeno_nc14b", recursive = TRUE, showWarnings = FALSE)

# --- full union peak BED (0-based start, as import_bed() expects) ----------
write.table(
  data.frame(chrom = as.character(seqnames(peaks_gr)),
             start = start(peaks_gr) - 1L, end = end(peaks_gr),
             name  = names(peaks_gr)),
  "./ATAC_NarrowPeaks/FullUniverse_nc14b_AllGeno_union_peaks.bed",
  sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)

# --- per-group peak BEDs ----------------------------------------------------
# Simplified: every group's "own called peaks" = the full universe. Tier 1
# condition-specific detection is therefore driven entirely by the simulated
# read-count differences above (via CPM thresholds), not by BED membership.
# This is a simplification vs. real data (where groups call somewhat
# different peak sets) but keeps the toy dataset's structure easy to reason
# about while still exercising the same statistical code path.
for (g in groups) {
  write.table(
    data.frame(chrom = as.character(seqnames(peaks_gr)),
               start = start(peaks_gr) - 1L, end = end(peaks_gr),
               name  = names(peaks_gr)),
    sprintf("./ATAC_NarrowPeaks/GroupPeaks_%s.bed", g),
    sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
}

# --- narrow_counts.txt cache -------------------------------------------------
# Written with cache_short_keys as headers (see note above) -- NOT
# sample_names -- so it survives read.table()'s hyphen-to-dot mangling on
# the real script's side. realign_cache_cols() then translates these short
# keys back to full BAM names via CACHE_KEY_MAP automatically.
counts_narrow_for_cache <- counts_narrow
colnames(counts_narrow_for_cache) <- cache_short_keys
write.table(counts_narrow_for_cache, "./counts_cache/SevenGeno_nc14b/narrow_counts.txt",
            sep = "\t", quote = FALSE, col.names = NA)

# --- grp_counts_cache.rds ----------------------------------------------------
# Every group's own-peak counts = the same matrix (consistent with group
# peaks = full universe above; a given genomic region has one true read
# count regardless of which "group peak set" it's being counted through).
grp_counts <- setNames(replicate(length(groups), counts_narrow, simplify = FALSE), groups)
grp_peaks  <- setNames(replicate(length(groups), peaks_gr,      simplify = FALSE), groups)
saveRDS(list(counts = grp_counts, peaks = grp_peaks),
        "./counts_cache/SevenGeno_nc14b/grp_counts_cache.rds")

# --- ext{200,300}_counts.rds -------------------------------------------------
# Boundary-extended peaks capture slightly more reads than the narrow call;
# simulate that as narrow counts + a modest Poisson increment. Names follow
# the real script's convention: "bext<ebp>_<original_peak_name>".
build_ext_counts <- function(ebp) {
  ext_gr <- GRanges(seqnames = CHROM,
                     ranges = IRanges(start = pmax(1, start(peaks_gr) - ebp),
                                       end  = end(peaks_gr) + ebp))
  names(ext_gr) <- paste0("bext", ebp, "_", names(peaks_gr))
  ext_counts <- counts_narrow + matrix(
    rpois(length(counts_narrow), lambda = 8),
    nrow = nrow(counts_narrow), dimnames = dimnames(counts_narrow))
  rownames(ext_counts) <- names(ext_gr)
  list(counts = ext_counts, peaks = ext_gr)
}

for (ebp in c(200L, 300L)) {
  saveRDS(build_ext_counts(ebp),
          sprintf("./counts_cache/SevenGeno_nc14b/ext%d_counts.rds", ebp))
}

# --- grp_ext{200,300}_counts.rds --------------------------------------------
build_grp_ext_counts <- function(ebp) {
  ext_obj <- build_ext_counts(ebp)
  counts_list <- setNames(replicate(length(groups), ext_obj$counts, simplify = FALSE), groups)
  peaks_list  <- setNames(replicate(length(groups), ext_obj$peaks,  simplify = FALSE), groups)
  list(counts = counts_list, peaks = peaks_list)
}

for (ebp in c(200L, 300L)) {
  saveRDS(build_grp_ext_counts(ebp),
          sprintf("./counts_cache/SevenGeno_nc14b/grp_ext%d_counts.rds", ebp))
}

cat("\nDone. Created:\n")
cat("  ./ATAC_NarrowPeaks/  (1 union BED + 7 group BEDs, ", N_PEAKS, "peaks each)\n")
cat("  ./counts_cache/SevenGeno_nc14b/  (narrow + grp + ext200/300 + grp_ext200/300)\n\n")
cat("Now run the real script from this same directory:\n")
cat("  Rscript 01_dar_calling_seven_genotype_nc14b.r\n\n")
cat("It should load every count matrix from cache (you'll see \"Loading ... from cache\"\n")
cat("messages, never \"running featureCounts\") and complete all 13 contrasts without\n")
cat("touching any BAM file.\n")
