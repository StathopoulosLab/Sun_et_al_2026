################################################################################
# 04_peak_fate_alluvial_plots_strict.r
#
# STRICT-MODE COPY of 03_peak_fate_alluvial_plots.r. The ONLY change from the
# original is the `out_contrast` path below, which is redirected into the
# STRICT/ subfolder that 02_peak_fate_classification_strict.r writes to. Because every
# read and write in this script goes through `out_contrast` (the fate-table
# read, all PDF ggsave() calls, and the emerged-gene CSV writes), that one
# redirect is enough to make the whole script self-contained to strict-mode
# output -- it reads STRICT/peak_fate_data.csv (produced by
# 02_peak_fate_classification_strict.r, not the default pipeline's peak_fate_data.csv)
# and writes all of its own PDFs/CSVs into STRICT/ as well, so nothing here
# touches or overwrites the default pipeline's alluvial output.
#
# Run 02_peak_fate_classification_strict.r first so STRICT/peak_fate_data.csv exists.
#
# Everything below this CONFIG block is otherwise IDENTICAL to
# 03_peak_fate_alluvial_plots.r -- see that file for plot-by-plot documentation.
# "Deepened" fate labels/colors are retained here only because they're part
# of the original script's palette; the current fate_split scheme (from
# 02_peak_fate_classification_strict.r, matching v16) never produces that label, so
# those bins will simply be empty.
#
# Improved alluvial plots for BOTv vs BOTCv nc14b → nc14late peak fate.
# Produces three versions:
#
#   peak_fate_alluvial_v2_full.pdf
#     — Full alluvial with:
#       • Emerged split by direction (BOTv-open vs BOTCv-open emerged)
#       • Converged greyed out and flows de-emphasized
#       • % of starting direction annotated on left nodes
#       • n + % on all right nodes
#
#   peak_fate_alluvial_v2_active_fates.pdf
#     — Interesting-fates only (excludes Converged peaks entirely)
#       Shows the 292 peaks that did something other than disappear.
#       Greatly improves readability of Deepened/Maintained/Reversed/Emerged.
#
#   peak_fate_alluvial_v2_panel.pdf
#     — Side-by-side: full (left) + active-fates zoom (right) with shared legend
#
# Run after plot_peak_fate_v10_comprehensive.r, or source standalone —
# it will re-load the data from disk using the same paths.
################################################################################

suppressPackageStartupMessages({
  library(GenomicRanges)   # still needed for gr_em (emerged-peak annotation section below)
  library(ggplot2)
  library(dplyr)
  library(patchwork)
  library(scales)
  library(ggalluvial)
})

################################################################################
# CONFIG — must match plot_peak_fate_v11.r (grpA/grpB/labels only; fate
# classification itself is no longer computed here -- see LOAD FATE TABLE
# below, which reads v11's peak_fate_data.csv directly)
################################################################################

out_dir  <- "./Overview_Plots"
out_contrast <- file.path(out_dir, "peak_fate", "BOTv_vs_BOTCv", "STRICT")   # [STRICT] redirected

# Cosmetic only now (used in a subtitle string) -- must match
# deepened_pct_thresh in plot_peak_fate_v11.r to avoid a mismatched label,
# but no longer drives any classification.
deepened_pct_thresh <- 0.2   # keep in sync with plot_peak_fate_v11.r

grpA <- "BOTv";  col_grpA <- "#2ca02c"
grpB <- "BOTCv"; col_grpB <- "#e377c2"
lbl_grpA <- "BOT.v-biased"     # nc14b direction label (PI-requested)
lbl_grpB <- "BOTC.v-biased"    # nc14b direction label (PI-requested)

# PI-requested label scheme -- must exactly match plot_peak_fate_v11.r's
# lbl_* definitions (see that script's comment for why these are literal
# strings rather than paste0(grpA,...)). Reversed labels are keyed by ENDING
# genotype: lbl_rev_A (started grpA-open, ends grpB-open) -> "Reversed to
# BOTC.v"; lbl_rev_B (started grpB-open, ends grpA-open) -> "Reversed to
# BOT.v". This matches the color comments just below (lbl_rev_A is
# BOTCv-pink/dark-magenta because it ENDS BOTCv-open; lbl_rev_B is
# BOTv-green because it ENDS BOTv-open) -- do not swap these.
lbl_emrg_A <- "Emerged BOT.v bias"
lbl_emrg_B <- "Emerged BOTC.v bias"
lbl_dpn_A <- "Deepened BOT.v bias"        # ends grpA-open
lbl_dpn_B <- "Deepened BOTC.v bias"       # ends grpB-open
lbl_mnt_A <- "Maintained BOT.v bias"      # ends grpA-open
lbl_mnt_B <- "Maintained BOTC.v bias"     # ends grpB-open
lbl_rev_A <- "Reversed to BOTC.v"         # started grpA-open, ENDS grpB-open
lbl_rev_B <- "Reversed to BOT.v"          # started grpB-open, ENDS grpA-open

fate_colors <- c(
  "Deepened"   = "#c9579a",   # BOTCv-dominant: medium-deep pink
  "Maintained" = "#e8a8cf",   # BOTCv-dominant: lighter pink
  "Reversed"   = "#8B0057",   # BOTv-dominant: dark magenta/deep pink
  "Converged"   = "#CCCCCC",   # light grey — de-emphasised
  "Emerged"    = "#F39C12"
)
# Direction-split fate colors. Deepened/Maintained-BOTCv keep their original
# pink shades; the BOTv-associated counterparts are new distinct greens
# instead of inheriting the same pink.
fate_colors_split <- c(
  "Converged" = "#CCCCCC"
)
fate_colors_split[lbl_dpn_A]  <- "#2ca02c"   # BOTv Deepened    -- bright green
fate_colors_split[lbl_mnt_A]  <- "#a8d8a8"   # BOTv Maintained  -- light green
fate_colors_split[lbl_dpn_B]  <- "#c9579a"   # BOTCv Deepened   -- medium-deep pink (unchanged)
fate_colors_split[lbl_mnt_B]  <- "#e8a8cf"   # BOTCv Maintained -- light pink (unchanged)
fate_colors_split[lbl_emrg_A] <- "#4aab4a"   # medium green (darker than BOTv, lighter than Deepened)
fate_colors_split[lbl_emrg_B] <- "#d44fa8"   # medium pink (darker than BOTCv, not as deep as Reversed)
fate_colors_split[lbl_rev_A]  <- "#8B0057"   # dark magenta (started BOTv-open, reverses INTO BOTCv-open)
fate_colors_split[lbl_rev_B]  <- "#1a7a3d"   # dark green   (started BOTCv-open, reverses INTO BOTv-open)

fate_order_split <- c(lbl_dpn_A,lbl_dpn_B,lbl_mnt_A,lbl_mnt_B,lbl_rev_A,lbl_rev_B,"Converged",
                      lbl_emrg_A, lbl_emrg_B)
fate_order_nores <- c(lbl_dpn_A,lbl_dpn_B,lbl_mnt_A,lbl_mnt_B,lbl_rev_A,lbl_rev_B,
                      lbl_emrg_A, lbl_emrg_B)

dir_colors <- c(
  setNames(col_grpA, lbl_grpA),
  setNames(col_grpB, lbl_grpB),
  "Not differential at nc14b" = "#F39C12"   # amber for the de novo left node
)

################################################################################
# LOAD FATE TABLE — single source of truth
################################################################################
# This used to recompute fate independently from the raw BEDs, with its own
# deepened_pct_thresh and a single-tier overlap match. That's exactly what let
# it drift out of sync with plot_peak_fate_v11.r (different threshold, and
# missing v11's tier-2/tier-3 rescue matching) -- a peak could legitimately
# come out "Deepened" in v11's peak_fate_data.csv and "Maintained" here.
#
# Instead, read v11's peak_fate_data.csv directly. v11 already applies the
# full three-tier matching (significant DAR / unfiltered trend / trajectory
# imputation) and the current deepened_pct_thresh -- every plot below now
# reflects exactly the same fate calls as the annotation/gene-list exports.
# Run plot_peak_fate_v11.r first so peak_fate_data.csv exists and is current.

cat("Loading fate table from 02_peak_fate_classification_strict.r output...\n")
fate_csv_path <- file.path(out_contrast, "peak_fate_data.csv")
if (!file.exists(fate_csv_path))
  stop("peak_fate_data.csv not found at ", fate_csv_path,
       " -- run 02_peak_fate_classification_strict.r first. It is now the single source of ",
       "truth for strict-mode fate classification; this script only reshapes its output.")

fate_raw <- read.csv(fate_csv_path, stringsAsFactors=FALSE)

is_emerged_row <- fate_raw$direction_nc14b == "Not differential at nc14b"

fate_df <- data.frame(
  direction    = fate_raw$direction_nc14b,
  lfc_nc14b    = fate_raw$lfc_nc14b,
  lfc_nc14late = fate_raw$lfc_nc14late,
  matched      = fate_raw$matched_nc14late,
  # fate_split already carries v11's BOTv/BOTCv-specific Deepened/Maintained/
  # Reversed labels and "Converged" unchanged -- no re-derivation here, so it
  # can never drift from v11's thresholds/matching again. Emerged rows get
  # their fate_split overwritten just below with the direction split.
  fate_base    = fate_raw$fate_split,
  # resolved_breakdown: v11's 6-way Converged mechanism split (BOTv/BOTCv
  # opens/closes, both shift, no temporal signal). NA for non-Converged rows.
  # Used only by the new "detailed converged subcategories" plot below.
  resolved_breakdown = if ("resolved_breakdown" %in% colnames(fate_raw))
                          fate_raw$resolved_breakdown else NA_character_,
  stringsAsFactors = FALSE)

# Emerged peaks — split by direction at nc14late, using v11's own
# direction_nc14late call (already negated-convention correct: lfc>0=grpA).
emrg_dir <- ifelse(fate_raw$direction_nc14late[is_emerged_row] == lbl_grpA,
                   lbl_emrg_A, lbl_emrg_B)
fate_df$fate_base[is_emerged_row] <- emrg_dir
n_emerged <- sum(is_emerged_row)

n_botv   <- sum(fate_df$direction==lbl_grpA)
n_botcv  <- sum(fate_df$direction==lbl_grpB)

# GRanges of emerged peaks, for the gene-annotation section further down.
# Built from the same is_emerged_row subset used for emrg_dir above, so row
# order (and therefore the emrg_dir alignment relied on below) is preserved.
emerged_rows <- fate_raw[is_emerged_row, ]
gr_em <- GRanges(seqnames = emerged_rows$chr,
                 ranges    = IRanges(emerged_rows$start + 1L, emerged_rows$end),
                 name      = emerged_rows$name,
                 log2FC    = emerged_rows$lfc_nc14late)
n_emrg_A <- sum(emrg_dir==lbl_emrg_A)
n_emrg_B <- sum(emrg_dir==lbl_emrg_B)

cat(sprintf("  nc14b: %d BOTv-open, %d BOTCv-open\n", n_botv, n_botcv))
cat(sprintf("  Emerged: %d (%d opens in %s, %d opens in %s)\n",
            n_emerged,
            sum(emrg_dir==lbl_emrg_A), grpA,
            sum(emrg_dir==lbl_emrg_B), grpB))
cat("\n  Fate breakdown (nc14b peaks only) -- matches plot_peak_fate_v11.r's printout exactly:\n")
print(table(fate_df$direction[fate_df$direction %in% c(lbl_grpA,lbl_grpB)],
            fate_df$fate_base[fate_df$direction %in% c(lbl_grpA,lbl_grpB)]))

################################################################################
# COLOR LOOKUP — must be defined before wide_to_lode() calls it
################################################################################

# Master stratum color lookup for all alluvial plots.
# Keys must NOT contain \n — after_stat(stratum) strips \n before lookup,
# causing scale_fill_manual to return NA and assign wrong colors.
stratum_colors <- c(
  setNames(col_grpA, lbl_grpA),
  setNames(col_grpB, lbl_grpB),
  "Not differential at nc14b" = "#F39C12",
  fate_colors_split
)

################################################################################
# BUILD ALLUVIAL DATA FRAMES
################################################################################

# ── Helper: wide → lode form ──────────────────────────────────────────────────
# Wide format (axis1=/axis2=) causes ggalluvial to mis-route ribbons when
# multiple direction groups share a fate stratum, because lode stacking within
# each stratum is determined by row order and factor levels independently on
# each axis — they don't stay in sync. Converting to explicit lode form with
# integer alluvium IDs guarantees each ribbon connects exactly the right pair.
#
# wide_df must have columns: direction (left), fate_base (right), n (weight)
# Returns a long data frame ready for aes(x, stratum, alluvium, y, fill)
wide_to_lode <- function(wide_df, dir_levels, fate_levels, grey_resolved=FALSE) {
  wide_df <- wide_df %>%
    dplyr::mutate(
      alluvium  = dplyr::row_number(),
      direction = factor(as.character(direction), levels=dir_levels),
      fate_base = factor(as.character(fate_base), levels=fate_levels),
      flow_alpha = dplyr::case_when(
        grey_resolved & as.character(fate_base)=="Converged" ~ 0.18,
        as.character(direction)=="Not differential at nc14b"       ~ 0.55,
        TRUE                                                 ~ 0.65
      )
    )
  bind_rows(
    wide_df %>% dplyr::transmute(
      x=1L, stratum=as.character(direction),
      alluvium=alluvium, y=n,
      fill_stratum=stratum_colors[as.character(direction)],
      flow_alpha=flow_alpha),
    wide_df %>% dplyr::transmute(
      x=2L, stratum=as.character(fate_base),
      alluvium=alluvium, y=n,
      fill_stratum=stratum_colors[as.character(fate_base)],
      flow_alpha=flow_alpha)
  ) %>%
    dplyr::mutate(
      stratum = factor(stratum, levels=c(dir_levels, fate_levels)),
      x       = factor(x, levels=c(1L,2L))
    )
}

# ── Build wide tables ─────────────────────────────────────────────────────────

dir_levels_full <- c(lbl_grpA, lbl_grpB, "Not differential at nc14b")

wide_full <- bind_rows(
  fate_df %>%
    dplyr::filter(direction %in% c(lbl_grpA, lbl_grpB)) %>%
    dplyr::count(direction, fate_base, name="n"),
  data.frame(direction="Not differential at nc14b", fate_base=lbl_emrg_A,
             n=sum(emrg_dir==lbl_emrg_A), stringsAsFactors=FALSE),
  data.frame(direction="Not differential at nc14b", fate_base=lbl_emrg_B,
             n=sum(emrg_dir==lbl_emrg_B), stringsAsFactors=FALSE)
)

wide_nores <- bind_rows(
  fate_df %>%
    dplyr::filter(direction %in% c(lbl_grpA, lbl_grpB),
                  fate_base != "Converged") %>%
    dplyr::count(direction, fate_base, name="n"),
  data.frame(direction="Not differential at nc14b", fate_base=lbl_emrg_A,
             n=sum(emrg_dir==lbl_emrg_A), stringsAsFactors=FALSE),
  data.frame(direction="Not differential at nc14b", fate_base=lbl_emrg_B,
             n=sum(emrg_dir==lbl_emrg_B), stringsAsFactors=FALSE)
)

n_interesting <- sum(fate_df$direction %in% c(lbl_grpA,lbl_grpB) &
                       fate_df$fate_base != "Converged")

# Counts still needed for subtitles / summary
n_botv_nores  <- sum(fate_df$direction==lbl_grpA & fate_df$fate_base!="Converged")
n_botcv_nores <- sum(fate_df$direction==lbl_grpB & fate_df$fate_base!="Converged")

lode_full  <- wide_to_lode(wide_full,  dir_levels_full, fate_order_split,
                            grey_resolved=TRUE)
lode_nores <- wide_to_lode(wide_nores, dir_levels_full, fate_order_nores,
                            grey_resolved=FALSE)

# Convenience: n per direction for subtitle counts (still needed)
alluv_full  <- wide_full  # keep for subtitle count code below
alluv_nores <- wide_nores # keep for subtitle count code below

################################################################################
# ALLUVIAL THEME
################################################################################

alluv_theme <- theme_minimal(base_size=11) +
  theme(
    panel.grid        = element_blank(),
    axis.text.x       = element_text(size=10, face="bold", colour="grey20"),
    axis.text.y       = element_blank(),
    axis.ticks        = element_blank(),
    plot.title        = element_text(size=13, face="bold"),
    plot.subtitle     = element_text(size=8.5, colour="grey40", lineheight=1.35),
    plot.caption      = element_text(size=7.5, colour="grey55"),
    legend.position   = "none"
  )

# Label function: n + % for right nodes, n + % for left nodes
make_stratum_label <- function(stratum, count, total_left) {
  # For left axis nodes: show n + what % of starting pop they represent
  # For right axis nodes: show n + % of starting pop that ended here
  paste0(stratum, "\nn=", prettyNum(count, big.mark=","))
}

################################################################################
# PLOT 1: FULL ALLUVIAL (with emerged split, resolved de-emphasised)
################################################################################

n_resolved <- sum(fate_df$fate_base=="Converged" &
                    fate_df$direction %in% c(lbl_grpA,lbl_grpB))

p_full <- ggplot(lode_full,
                 aes(x=x, stratum=stratum, alluvium=alluvium,
                     y=y, fill=stratum)) +
  geom_flow(aes(alpha=flow_alpha), aes.flow="backward",
            width=0.38, knot.anchor=0.45, curve_type="sigmoid") +
  scale_alpha_identity() +
  geom_stratum(width=0.38, colour="white", linewidth=0.5) +
  scale_fill_manual(values=stratum_colors, guide="none") +
  geom_text(stat="stratum",
            aes(label=after_stat(paste0(stratum, "\n",
                                        prettyNum(count, big.mark=",")))),
            size=2.6, lineheight=1.25, colour="grey10") +
  scale_x_discrete(limits=c("1","2"),
                   labels=c("Direction at nc14b","Fate at nc14late"),
                   expand=c(0.20, 0.20)) +
  labs(
    title = sprintf("Peak fate: %s vs %s   nc14b \u2192 nc14late", grpA, grpB),
    subtitle = sprintf(
      "Starting: %s BOTv-open + %s BOTCv-open peaks  |  %s emerged de novo at nc14late (%s opens in %s, %s in %s)\nConverged flows greyed \u2014 %s peaks (%.0f%% of starting pop) no longer a DAR at nc14late",
      prettyNum(n_botv, big.mark=","),
      prettyNum(n_botcv, big.mark=","),
      prettyNum(n_emerged, big.mark=","),
      prettyNum(sum(emrg_dir==lbl_emrg_A), big.mark=","), grpA,
      prettyNum(sum(emrg_dir==lbl_emrg_B), big.mark=","), grpB,
      prettyNum(n_resolved, big.mark=","),
      100 * n_resolved / (n_botv + n_botcv)
    ),
    y="Number of peaks", x=NULL,
    caption=sprintf(
      "Deepened=fold-change magnitude of the genotype-difference grew >%.0f%% (same direction)  |  Reversed=direction flipped  |  Converged=no longer a DAR at nc14late\nEmerged split by direction at nc14late: opens in %s (light green) vs opens in %s (light pink)",
      100*deepened_pct_thresh, grpA, grpB)
  ) +
  alluv_theme

fn_full <- file.path(out_contrast, "peak_fate_alluvial_v2_full.pdf")
ggsave(fn_full, p_full, width=10, height=8)
cat("Saved:", fn_full, "\n")

################################################################################
# PLOT 2: ACTIVE-FATES ZOOM
################################################################################

pct_botv_nores  <- round(100 * n_botv_nores  / n_botv,  1)
pct_botcv_nores <- round(100 * n_botcv_nores / n_botcv, 1)

p_nores <- ggplot(lode_nores,
                  aes(x=x, stratum=stratum, alluvium=alluvium,
                      y=y, fill=stratum)) +
  geom_flow(aes(alpha=flow_alpha), aes.flow="backward",
            width=0.38, knot.anchor=0.45, curve_type="sigmoid") +
  scale_alpha_identity() +
  geom_stratum(width=0.38, colour="white", linewidth=0.5) +
  scale_fill_manual(values=stratum_colors, guide="none") +
  geom_text(stat="stratum",
            aes(label=after_stat(paste0(stratum, "\nn=",
                                        prettyNum(count, big.mark=",")))),
            size=2.8, lineheight=1.25, colour="grey10") +
  scale_x_discrete(limits=c("1","2"),
                   labels=c("Direction at nc14b","Fate at nc14late"),
                   expand=c(0.22, 0.22)) +
  labs(
    title = sprintf("Interesting fates only: %s vs %s   nc14b \u2192 nc14late", grpA, grpB),
    subtitle = sprintf(
      "Converged peaks excluded  |  Showing %s peaks (%s BOTv-open [%.1f%% of BOTv DARs], %s BOTCv-open [%.1f%%])  +  %s emerged",
      prettyNum(n_interesting, big.mark=","),
      prettyNum(n_botv_nores,  big.mark=","), pct_botv_nores,
      prettyNum(n_botcv_nores, big.mark=","), pct_botcv_nores,
      prettyNum(n_emerged, big.mark=",")),
    y="Number of peaks", x=NULL,
    caption=sprintf(
      "Converged = %.0f%% of nc14b starting population \u2014 excluded here to zoom in on active fate transitions\nEmerged split: opens in %s (light green) vs opens in %s (light pink)",
      100 * n_resolved / (n_botv + n_botcv), grpA, grpB)
  ) +
  alluv_theme

fn_nores <- file.path(out_contrast, "peak_fate_alluvial_v2_active_fates.pdf")
ggsave(fn_nores, p_nores, width=9, height=7)
cat("Saved:", fn_nores, "\n")

################################################################################
# PLOT 1b: DETAILED — CONVERGED BROKEN INTO SUBCATEGORIES
################################################################################
# Companion to PLOT 1 (which collapses Converged into one grey bucket). Same
# active-fate splits (Deepened/Maintained/Reversed/Emerged, PI-requested
# labels), but Converged is broken into v11.r's 6 mechanism subcategories
# (resolved_breakdown) instead of staying one node. Two versions were
# requested -- simplified (PLOT 1, single Converged) and this one (Converged
# subcategories) -- so both exist side by side rather than one replacing the
# other.

lbl_res <- c("Converged: BOTv opens","Converged: BOTv closes","Converged: both shift",
            "Converged: BOTCv opens","Converged: BOTCv closes","Converged: no temporal signal")

fate_colors_resbreak <- c(
  fate_colors_split[setdiff(names(fate_colors_split), "Converged")],
  "Converged: BOTv opens"          = "#7fb3a0",
  "Converged: BOTv closes"         = "#a68a5b",
  "Converged: both shift"          = "#8067b7",
  "Converged: BOTCv opens"         = "#d99cc0",
  "Converged: BOTCv closes"        = "#8a6a7a",
  "Converged: no temporal signal"  = "#d9d9d9"
)
fate_order_resbreak <- c(lbl_dpn_A,lbl_dpn_B,lbl_mnt_A,lbl_mnt_B,lbl_rev_A,lbl_rev_B,
                         lbl_res, lbl_emrg_A, lbl_emrg_B)
stratum_colors_resbreak <- c(
  setNames(col_grpA, lbl_grpA),
  setNames(col_grpB, lbl_grpB),
  "Not differential at nc14b" = "#F39C12",
  fate_colors_resbreak
)

if (!("resolved_breakdown" %in% colnames(fate_df)) || all(is.na(fate_df$resolved_breakdown))) {
  cat("  [NOTE] resolved_breakdown column missing/empty in peak_fate_data.csv --\n",
      "  skipping detailed-converged-subcategories plot. Re-run plot_peak_fate_v11.r\n",
      "  (which computes resolved_breakdown unconditionally) to enable it.\n")
} else {
  n_resolved_detail <- sum(fate_df$fate_base=="Converged" &
                           fate_df$direction %in% c(lbl_grpA,lbl_grpB))

  wide_resbreak <- bind_rows(
    fate_df %>%
      dplyr::filter(direction %in% c(lbl_grpA, lbl_grpB)) %>%
      dplyr::mutate(fate_final = ifelse(fate_base=="Converged", resolved_breakdown, fate_base)) %>%
      dplyr::count(direction, fate_final, name="n") %>%
      dplyr::rename(fate_base = fate_final),
    data.frame(direction="Not differential at nc14b", fate_base=lbl_emrg_A,
               n=sum(emrg_dir==lbl_emrg_A), stringsAsFactors=FALSE),
    data.frame(direction="Not differential at nc14b", fate_base=lbl_emrg_B,
               n=sum(emrg_dir==lbl_emrg_B), stringsAsFactors=FALSE)
  )

  # Local wide_to_lode call needs the resbreak-specific color/order maps, not
  # the master stratum_colors/fate_order_split used elsewhere in this script
  # -- build the lode form inline rather than reusing wide_to_lode() as-is.
  lode_resbreak <- wide_resbreak %>%
    dplyr::mutate(
      alluvium  = dplyr::row_number(),
      direction = factor(as.character(direction), levels=dir_levels_full),
      fate_base = factor(as.character(fate_base), levels=fate_order_resbreak),
      flow_alpha = ifelse(as.character(direction)=="Not differential at nc14b", 0.55, 0.65)
    ) %>%
    { bind_rows(
        dplyr::transmute(., x=1L, stratum=as.character(direction), alluvium=alluvium,
                         y=n, flow_alpha=flow_alpha),
        dplyr::transmute(., x=2L, stratum=as.character(fate_base), alluvium=alluvium,
                         y=n, flow_alpha=flow_alpha)
      ) } %>%
    dplyr::mutate(stratum=factor(stratum, levels=c(dir_levels_full, fate_order_resbreak)),
                  x=factor(x, levels=c(1L,2L)))

  p_resbreak <- ggplot(lode_resbreak,
                       aes(x=x, stratum=stratum, alluvium=alluvium, y=y, fill=stratum)) +
    geom_flow(aes(alpha=flow_alpha), aes.flow="backward",
              width=0.38, knot.anchor=0.45, curve_type="sigmoid") +
    scale_alpha_identity() +
    geom_stratum(width=0.38, colour="white", linewidth=0.5) +
    scale_fill_manual(values=stratum_colors_resbreak, guide="none") +
    geom_text(stat="stratum",
              aes(label=after_stat(paste0(stratum,"\nn=",prettyNum(count,big.mark=",")))),
              size=2.3, lineheight=1.2, colour="grey10") +
    scale_x_discrete(limits=c("1","2"),
                     labels=c("Direction at nc14b","Fate at nc14late (Converged detailed)"),
                     expand=c(0.20,0.20)) +
    labs(
      title = sprintf("Peak fate (Converged detailed \u2014 subcategories shown): %s vs %s   nc14b \u2192 nc14late", grpA, grpB),
      subtitle = sprintf(
        "Converged (%s peaks, %.0f%% of starting pop) split by which genotype's own accessibility moved during resolution",
        prettyNum(n_resolved_detail,big.mark=","), 100*n_resolved_detail/(n_botv+n_botcv)),
      y="Number of peaks", x=NULL,
      caption="\"Converged\" alone only says the BOTv-vs-BOTCv difference disappeared -- these sub-categories say why.\nSee peak_fate_alluvial_v2_full.pdf for the simplified (single-Converged-bucket) version of this same plot.") +
    alluv_theme +
    theme(plot.subtitle=element_text(size=8, colour="grey40", lineheight=1.3))

  fn_resbreak <- file.path(out_contrast, "peak_fate_alluvial_v2_detailed_converged_subcategories.pdf")
  ggsave(fn_resbreak, p_resbreak, width=12, height=8.5)
  cat("Saved:", fn_resbreak, "\n")
}

################################################################################
# PLOT 2b: NC14B-ONLY — emerged as right-side context, no spurious flows
################################################################################

# Uses the same stratum_colors master lookup — no separate fill map needed
fate_order_nores_v2 <- c(lbl_dpn_A,lbl_dpn_B,lbl_mnt_A,lbl_mnt_B,lbl_rev_A,lbl_rev_B, lbl_emrg_A, lbl_emrg_B)
dir_order_nores_v2  <- c(lbl_grpA, lbl_grpB, "Not differential at nc14b")

n_interesting <- sum(fate_df$direction %in% c(lbl_grpA,lbl_grpB) &
                     fate_df$fate_base %in% c(lbl_dpn_A,lbl_dpn_B,lbl_mnt_A,lbl_mnt_B,lbl_rev_A,lbl_rev_B))

# ── Build lode-form data for Plot 2b ──────────────────────────────────────────
# Wide format with split geom_flow calls causes ggalluvial to mis-route ribbons
# because it can't reconcile lode stacking across the two separate data frames.
# Using explicit lode (long) format with alluvium IDs guarantees each ribbon
# connects exactly the right left-node to the right-node.
#
# Each "alluvium" is one (direction, fate) combination — an independent ribbon.
# We assign a unique integer ID per row, then pivot to long form with x=1 (left
# axis) and x=2 (right axis), so ggalluvial draws one ribbon per alluvium ID.

wide_v2 <- bind_rows(
  # nc14b peaks → their fates (Deepened/Maintained/Reversed, each split BOTv/BOTCv)
  fate_df %>%
    dplyr::filter(direction %in% c(lbl_grpA, lbl_grpB),
                  fate_base %in% c(lbl_dpn_A,lbl_dpn_B,lbl_mnt_A,lbl_mnt_B,lbl_rev_A,lbl_rev_B)) %>%
    dplyr::count(direction, fate_base, name="n"),
  # Not differential at nc14b left node → its two right-side fate splits
  data.frame(direction="Not differential at nc14b",
             fate_base=c(lbl_emrg_A, lbl_emrg_B),
             n=c(n_emrg_A, n_emrg_B),
             stringsAsFactors=FALSE)
) %>%
  dplyr::mutate(alluvium = row_number(),
                direction = factor(direction, levels=dir_order_nores_v2),
                fate_base = factor(fate_base, levels=fate_order_nores_v2),
                flow_fill = stratum_colors[as.character(direction)],
                flow_alpha = ifelse(as.character(direction)=="Not differential at nc14b",
                                    0.55, 0.65))

# Convert to lode form: x=1 → left stratum (direction), x=2 → right (fate_base)
lode_v2 <- bind_rows(
  wide_v2 %>% dplyr::transmute(
    x        = 1L,
    stratum  = as.character(direction),
    alluvium = alluvium,
    y        = n,
    fill_col = flow_fill,
    flow_alpha = flow_alpha
  ),
  wide_v2 %>% dplyr::transmute(
    x        = 2L,
    stratum  = as.character(fate_base),
    alluvium = alluvium,
    y        = n,
    fill_col = stratum_colors[as.character(fate_base)],
    flow_alpha = flow_alpha
  )
) %>%
  dplyr::mutate(
    stratum = factor(stratum, levels=c(dir_order_nores_v2, fate_order_nores_v2)),
    x       = factor(x, levels=c(1L, 2L))
  )

p_nores_v2 <- ggplot(lode_v2,
                     aes(x=x, stratum=stratum, alluvium=alluvium,
                         y=y, fill=stratum)) +
  geom_flow(aes(alpha=flow_alpha), aes.flow="backward",
            width=0.38, knot.anchor=0.45, curve_type="sigmoid") +
  scale_alpha_identity() +
  geom_stratum(width=0.38, colour="white", linewidth=0.5) +
  scale_fill_manual(values=stratum_colors, guide="none") +
  geom_text(stat="stratum",
            aes(label=after_stat(paste0(stratum, "\nn=",
                                        prettyNum(count, big.mark=",")))),
            size=2.8, lineheight=1.25, colour="grey10") +
  scale_x_discrete(limits=c("1","2"),
                   labels=c("Direction at nc14b","Fate at nc14late"),
                   expand=c(0.22, 0.22)) +
  labs(
    title    = sprintf("Interesting fates: %s vs %s   nc14b \u2192 nc14late", grpA, grpB),
    subtitle = sprintf(
      "Converged excluded  |  %s nc14b peaks with active fates  +  %s de novo emerged (%s %s-open, %s %s-open)",
      prettyNum(n_interesting, big.mark=","),
      prettyNum(n_emerged,  big.mark=","),
      prettyNum(n_emrg_A, big.mark=","), grpA,
      prettyNum(n_emrg_B, big.mark=","), grpB),
    y="Number of peaks", x=NULL,
    caption=sprintf(
      "nc14b peaks (left) flow to fate at nc14late (right)\nEmerged peaks (de novo, no nc14b state) shown as right-side context: %s (green) vs %s (pink)",
      grpA, grpB)
  ) +
  alluv_theme

fn_nores_v2 <- file.path(out_contrast, "peak_fate_alluvial_v2_active_fates_nc14bonly.pdf")
ggsave(fn_nores_v2, p_nores_v2, width=9, height=7)
cat("Saved:", fn_nores_v2, "\n")

################################################################################
# PLOT 2c: EMERGED-AS-LEFT-NODES
#
# Layout: emerged peaks appear as LEFT-axis source nodes (not aggregated into
# a single "De novo" node). Each emerged direction has its own left stratum
# that flows straight across to its matching right stratum.
#
# Left axis (top→bottom):
#   Emerged (opens in BOTCv)   [amber-pink]
#   Emerged (opens in BOTv)    [amber-green]
#   More open in BOTCv         [pink]
#   More open in BOTv          [green]
#
# Right axis (top→bottom):
#   Reversed                   [dark magenta]
#   Maintained                 [light pink]
#   Deepened                   [medium pink]
#   Emerged (opens in BOTCv)   [amber-pink]
#   Emerged (opens in BOTv)    [amber-green]
#
# nc14b peaks (bottom two left nodes) flow to Deepened/Maintained/Reversed.
# Emerged peaks (top two left nodes) flow straight across to matching right node.
################################################################################

# Left-axis level order (bottom to top in ggalluvial = last to first in levels)
dir_order_v3 <- c(lbl_grpA, lbl_grpB, lbl_emrg_A, lbl_emrg_B)

# Right-axis level order
fate_order_v3 <- c(lbl_dpn_A, lbl_dpn_B, lbl_mnt_A, lbl_mnt_B, lbl_rev_A, lbl_rev_B, lbl_emrg_A, lbl_emrg_B)

# Colors for the new left-axis emerged nodes
stratum_colors_v3 <- c(
  stratum_colors,
  # lbl_emrg_A / lbl_emrg_B are already in stratum_colors — just making explicit:
  setNames("#4aab4a", lbl_emrg_A),
  setNames("#d44fa8", lbl_emrg_B)
)
# stratum_colors already has both — deduplicate by just reusing it
stratum_colors_v3 <- stratum_colors_v3[!duplicated(names(stratum_colors_v3))]

# Wide table: one row per (left-node, right-node) pair
wide_v3 <- bind_rows(
  # nc14b peaks → Deepened / Maintained / Reversed, each split BOTv/BOTCv
  fate_df %>%
    dplyr::filter(direction %in% c(lbl_grpA, lbl_grpB),
                  fate_base %in% c(lbl_dpn_A,lbl_dpn_B,lbl_mnt_A,lbl_mnt_B,lbl_rev_A,lbl_rev_B)) %>%
    dplyr::count(direction, fate_base, name="n") %>%
    dplyr::rename(left_node=direction, right_node=fate_base),
  # Emerged peaks: left node = emerged direction, right node = same emerged direction
  data.frame(
    left_node  = c(lbl_emrg_A, lbl_emrg_B),
    right_node = c(lbl_emrg_A, lbl_emrg_B),
    n          = c(n_emrg_A,   n_emrg_B),
    stringsAsFactors = FALSE
  )
) %>%
  dplyr::mutate(
    alluvium   = dplyr::row_number(),
    left_node  = factor(left_node,  levels=dir_order_v3),
    right_node = factor(right_node, levels=fate_order_v3),
    flow_alpha = dplyr::case_when(
      as.character(left_node) %in% c(lbl_emrg_A, lbl_emrg_B) ~ 0.50,
      TRUE ~ 0.65
    )
  )

# Lode form
lode_v3 <- dplyr::bind_rows(
  wide_v3 %>% dplyr::transmute(
    x        = 1L,
    stratum  = as.character(left_node),
    alluvium = alluvium,
    y        = n,
    flow_alpha = flow_alpha
  ),
  wide_v3 %>% dplyr::transmute(
    x        = 2L,
    stratum  = as.character(right_node),
    alluvium = alluvium,
    y        = n,
    flow_alpha = flow_alpha
  )
) %>%
  dplyr::mutate(
    stratum = factor(stratum, levels=unique(c(dir_order_v3, fate_order_v3))),
    x       = factor(x, levels=c(1L, 2L))
  )

p_nores_v3 <- ggplot(lode_v3,
                     aes(x=x, stratum=stratum, alluvium=alluvium,
                         y=y, fill=stratum)) +
  geom_flow(aes(alpha=flow_alpha), aes.flow="backward",
            width=0.38, knot.anchor=0.45, curve_type="sigmoid") +
  scale_alpha_identity() +
  geom_stratum(width=0.38, colour="white", linewidth=0.5) +
  scale_fill_manual(values=stratum_colors, guide="none") +
  geom_text(stat="stratum",
            aes(label=after_stat(paste0(stratum, "\nn=",
                                        prettyNum(count, big.mark=",")))),
            size=2.8, lineheight=1.25, colour="grey10") +
  scale_x_discrete(limits=c("1","2"),
                   labels=c("Direction at nc14b","Fate at nc14late"),
                   expand=c(0.24, 0.24)) +
  labs(
    title    = sprintf("Interesting fates + emerged: %s vs %s   nc14b \u2192 nc14late",
                       grpA, grpB),
    subtitle = sprintf(
      "Converged excluded  |  %s nc14b peaks with active fates  |  %s emerged de novo (%s %s-open, %s %s-open)\nEmerged peaks shown as left-axis source nodes flowing to their direction at nc14late",
      prettyNum(n_interesting, big.mark=","),
      prettyNum(n_emerged,     big.mark=","),
      prettyNum(n_emrg_A, big.mark=","), grpA,
      prettyNum(n_emrg_B, big.mark=","), grpB),
    y="Number of peaks", x=NULL,
    caption=sprintf(
      "nc14b peaks (bottom two left nodes) flow to Deepened/Maintained/Reversed at nc14late\nEmerged peaks (top two left nodes) flow straight to their nc14late direction: %s (green) vs %s (pink)",
      grpA, grpB)
  ) +
  alluv_theme

fn_nores_v3 <- file.path(out_contrast, "peak_fate_alluvial_v2_active_fates_emerged_left.pdf")
ggsave(fn_nores_v3, p_nores_v3, width=9, height=7)
cat("Saved:", fn_nores_v3, "\n")

################################################################################
# PLOT 3: PANEL A+B+C+D
################################################################################

p_full_panel   <- p_full     + labs(tag="A") + theme(plot.tag=element_text(size=14,face="bold"))
p_nores_panel  <- p_nores    + labs(tag="B") + theme(plot.tag=element_text(size=14,face="bold"))
p_nores2_panel <- p_nores_v2 + labs(tag="C") + theme(plot.tag=element_text(size=14,face="bold"))
p_nores3_panel <- p_nores_v3 + labs(tag="D") + theme(plot.tag=element_text(size=14,face="bold"))

p_panel <- p_full_panel + p_nores_panel + p_nores2_panel + p_nores3_panel +
  plot_annotation(
    title    = sprintf("Peak fate: %s vs %s   nc14b \u2192 nc14late", grpA, grpB),
    subtitle = paste0("A: Full alluvial (Converged flows greyed)   |   ",
                      "B: Zoom \u2014 Converged excluded (nc14b + De novo left nodes)   |   ",
                      "C: Zoom \u2014 nc14b peaks only, emerged as right-side context   |   ",
                      "D: Zoom \u2014 emerged as left-axis source nodes"),
    theme    = theme(plot.title    = element_text(size=14, face="bold"),
                     plot.subtitle = element_text(size=9,  colour="grey40"))
  )

fn_panel <- file.path(out_contrast, "peak_fate_alluvial_v2_panel.pdf")
ggsave(fn_panel, p_panel, width=36, height=8)
cat("Saved:", fn_panel, "\n")

################################################################################
# GENE ANNOTATION FOR EMERGED PEAKS
#
# Annotates emerged peaks (those appearing only at nc14late, split by direction)
# with the nearest gene, using one of three strategies in order of preference:
#
#   1. ChIPseeker::annotatePeak  (full annotation: nearest gene, feature type,
#                                 distance to TSS) — preferred
#   2. Nearest-gene via TxDb + distanceToNearest  (lightweight fallback)
#   3. Peak name / coordinates only (last resort if no annotation DB available)
#
# OUTPUTS (written to out_contrast/):
#   emerged_BOTv_open_genes.csv     — peaks more open in BOTv at nc14late
#   emerged_BOTCv_open_genes.csv    — peaks more open in BOTCv at nc14late
#   emerged_all_annotated.csv       — combined table with group column
#
# Column key in output CSVs:
#   peak_name        — peak ID from BED col 4 (or auto-generated)
#   chr / start / end — genomic coordinates (0-based BED convention)
#   lfc_nc14late     — log2FC at nc14late (negated/BOTv-numerator convention: >0 = BOTv-open)
#   emerged_group    — "BOTv_open" or "BOTCv_open"
#   gene_name        — nearest gene symbol (SYMBOL column from annotation)
#   gene_id          — Ensembl / FlyBase gene ID
#   annotation       — genomic feature (Promoter, Intron, Exon, Distal IGR …)
#   distance_to_tss  — bp distance to nearest TSS (negative = upstream)
#   distanceToTSS    — alias, kept for ChIPseeker compatibility
################################################################################

cat("\n=== EMERGED PEAK GENE ANNOTATION ===\n")

dir.create(out_contrast, recursive = TRUE, showWarnings = FALSE)

# ── Build coordinate data frame for emerged peaks ─────────────────────────────

emerged_coord_df <- data.frame(
  peak_name    = if (!is.null(mcols(gr_em)$name)) as.character(mcols(gr_em)$name)
                 else paste0("emerged_", seq_len(length(gr_em))),
  chr          = as.character(seqnames(gr_em)),
  start        = start(gr_em) - 1L,   # back to 0-based BED
  end          = end(gr_em),
  lfc_nc14late = mcols(gr_em)$log2FC,
  emerged_group = ifelse(emrg_dir == lbl_emrg_A, "BOTv_open", "BOTCv_open"),
  stringsAsFactors = FALSE
)

# ── Annotation strategy 1: ChIPseeker ─────────────────────────────────────────

has_chipseeker <- requireNamespace("ChIPseeker",   quietly = TRUE)
has_txdb_dm6   <- requireNamespace("TxDb.Dmelanogaster.UCSC.dm6.ensGene", quietly = TRUE)
has_txdb_dm3   <- requireNamespace("TxDb.Dmelanogaster.UCSC.dm3.ensGene", quietly = TRUE)
has_orgdb      <- requireNamespace("org.Dm.eg.db",  quietly = TRUE)

txdb <- NULL
if (has_txdb_dm6) {
  txdb <- TxDb.Dmelanogaster.UCSC.dm6.ensGene::TxDb.Dmelanogaster.UCSC.dm6.ensGene
  cat("  Using TxDb: dm6\n")
} else if (has_txdb_dm3) {
  txdb <- TxDb.Dmelanogaster.UCSC.dm3.ensGene::TxDb.Dmelanogaster.UCSC.dm3.ensGene
  cat("  Using TxDb: dm3 (dm6 not found)\n")
} else {
  cat("  No TxDb found — annotation will fall back to coordinates only.\n")
}

annot_df <- NULL   # will be populated by whichever strategy succeeds

if (has_chipseeker && !is.null(txdb)) {

  cat("  Strategy 1: ChIPseeker::annotatePeak\n")

  # ChIPseeker can optionally add gene symbols via org.Db
  anno_args <- list(
    peak          = gr_em,
    TxDb          = txdb,
    tssRegion     = c(-3000, 3000),
    verbose       = FALSE
  )
  if (has_orgdb) {
    anno_args$annoDb <- "org.Dm.eg.db"
    cat("  Adding gene symbols via org.Dm.eg.db\n")
  }

  peak_anno <- tryCatch(
    do.call(ChIPseeker::annotatePeak, anno_args),
    error = function(e) {
      message("  ChIPseeker::annotatePeak failed: ", conditionMessage(e))
      NULL
    }
  )

  if (!is.null(peak_anno)) {
    cs_df <- as.data.frame(peak_anno)

    # Normalise column names across ChIPseeker versions
    if ("SYMBOL"   %in% colnames(cs_df)) cs_df$gene_name <- cs_df$SYMBOL
    if ("geneId"   %in% colnames(cs_df)) cs_df$gene_id   <- cs_df$geneId
    if ("annotation" %in% colnames(cs_df)) cs_df$annotation <- cs_df$annotation
    if ("distanceToTSS" %in% colnames(cs_df)) cs_df$distance_to_tss <- cs_df$distanceToTSS

    # Merge with our coord frame (match by row order — gr_em order preserved)
    annot_df <- cbind(emerged_coord_df, cs_df[, intersect(
      c("gene_name","gene_id","annotation","distance_to_tss","distanceToTSS"),
      colnames(cs_df)), drop = FALSE])

    cat(sprintf("  ChIPseeker annotated %d emerged peaks.\n", nrow(annot_df)))
  }
}

# ── Annotation strategy 2: TxDb nearest-gene (no ChIPseeker) ─────────────────

if (is.null(annot_df) && !is.null(txdb)) {

  cat("  Strategy 2: distanceToNearest via TxDb genes\n")

  suppressPackageStartupMessages(library(GenomicFeatures))
  genes_gr <- tryCatch(genes(txdb), error = function(e) NULL)

  if (!is.null(genes_gr)) {
    d2n   <- distanceToNearest(gr_em, genes_gr, ignore.strand = TRUE)
    hit_s <- subjectHits(d2n)
    dist_v <- mcols(d2n)$distance

    gene_id_v   <- rep(NA_character_, length(gr_em))
    gene_name_v <- rep(NA_character_, length(gr_em))

    gene_id_v[queryHits(d2n)] <- as.character(genes_gr$gene_id[hit_s])

    # Map gene IDs to symbols if org.Db available
    if (has_orgdb) {
      suppressPackageStartupMessages(library(org.Dm.eg.db))
      sym_map <- tryCatch(
        AnnotationDbi::select(org.Dm.eg.db,
                              keys    = na.omit(unique(gene_id_v)),
                              columns = "SYMBOL",
                              keytype = "ENSEMBL"),
        error = function(e) NULL
      )
      if (!is.null(sym_map)) {
        sym_map <- sym_map[!duplicated(sym_map$ENSEMBL), ]
        gene_name_v[queryHits(d2n)] <-
          sym_map$SYMBOL[match(gene_id_v[queryHits(d2n)], sym_map$ENSEMBL)]
      }
    }

    dist_full <- rep(NA_integer_, length(gr_em))
    dist_full[queryHits(d2n)] <- dist_v

    annot_df <- emerged_coord_df
    annot_df$gene_id        <- gene_id_v
    annot_df$gene_name      <- gene_name_v
    annot_df$distance_to_tss <- dist_full
    annot_df$annotation     <- NA_character_

    cat(sprintf("  Nearest-gene annotated %d emerged peaks.\n", nrow(annot_df)))
  }
}

# ── Annotation strategy 3: coordinates only ───────────────────────────────────

if (is.null(annot_df)) {
  cat("  Strategy 3: No annotation DB available — outputting coordinates only.\n")
  annot_df <- emerged_coord_df
  annot_df$gene_name      <- NA_character_
  annot_df$gene_id        <- NA_character_
  annot_df$annotation     <- NA_character_
  annot_df$distance_to_tss <- NA_integer_
}

# ── Write per-group and combined CSVs ─────────────────────────────────────────

# Ensure consistent column order
core_cols   <- c("peak_name","chr","start","end","lfc_nc14late","emerged_group")
annot_cols  <- c("gene_name","gene_id","annotation","distance_to_tss")
# Only keep annot_cols that actually exist
annot_cols  <- intersect(annot_cols, colnames(annot_df))
out_cols    <- c(core_cols, annot_cols)
annot_df    <- annot_df[, out_cols, drop = FALSE]

# Sort: BOTv-open by lfc ascending (most open first), BOTCv-open by lfc descending
annot_botv  <- annot_df[annot_df$emerged_group == "BOTv_open",  ]
annot_botcv <- annot_df[annot_df$emerged_group == "BOTCv_open", ]
annot_botv  <- annot_botv[order(annot_botv$lfc_nc14late,  na.last = TRUE), ]
annot_botcv <- annot_botcv[order(-annot_botcv$lfc_nc14late, na.last = TRUE), ]

fn_botv  <- file.path(out_contrast, "emerged_BOTv_open_genes.csv")
fn_botcv <- file.path(out_contrast, "emerged_BOTCv_open_genes.csv")
fn_all   <- file.path(out_contrast, "emerged_all_annotated.csv")

write.csv(annot_botv,  fn_botv,  row.names = FALSE)
write.csv(annot_botcv, fn_botcv, row.names = FALSE)
write.csv(rbind(annot_botv, annot_botcv), fn_all, row.names = FALSE)

cat(sprintf("\n  Saved emerged gene lists:\n"))
cat(sprintf("    BOTv-open   (%3d peaks) : %s\n",  nrow(annot_botv),  fn_botv))
cat(sprintf("    BOTCv-open  (%3d peaks) : %s\n",  nrow(annot_botcv), fn_botcv))
cat(sprintf("    Combined    (%3d peaks) : %s\n",
            nrow(annot_botv) + nrow(annot_botcv), fn_all))

# ── Console preview ────────────────────────────────────────────────────────────

preview_cols <- intersect(c("peak_name","lfc_nc14late","gene_name","annotation",
                             "distance_to_tss"), colnames(annot_df))

cat("\n  BOTv-open emerged (top 10 by |LFC|):\n")
top_botv <- head(annot_botv[order(abs(annot_botv$lfc_nc14late), decreasing = TRUE), preview_cols], 10)
print(as.data.frame(top_botv), row.names = FALSE)

cat("\n  BOTCv-open emerged (top 10 by |LFC|):\n")
top_botcv <- head(annot_botcv[order(abs(annot_botcv$lfc_nc14late), decreasing = TRUE), preview_cols], 10)
print(as.data.frame(top_botcv), row.names = FALSE)

################################################################################
# SUMMARY
################################################################################

cat("\n")
cat(strrep("=",65), "\n")
cat("ALLUVIAL SUMMARY\n")
cat(strrep("=",65), "\n")
cat(sprintf("  BOTv-open peaks       : %d\n", n_botv))
cat(sprintf("  BOTCv-open peaks      : %d\n", n_botcv))
cat(sprintf("  Total starting pop    : %d\n", n_botv+n_botcv))
cat(sprintf("  Converged              : %d (%.1f%%)\n",
            sum(fate_df$fate_base=="Converged" & fate_df$direction %in% c(lbl_grpA,lbl_grpB)),
            100*sum(fate_df$fate_base=="Converged" & fate_df$direction %in% c(lbl_grpA,lbl_grpB)) /
              (n_botv+n_botcv)))
cat(sprintf("  Interesting fates     : %d (%.1f%%)\n",
            n_interesting,
            100*n_interesting/(n_botv+n_botcv)))
cat(sprintf("  Emerged (BOTv-open)   : %d\n", sum(emrg_dir==lbl_emrg_A)))
cat(sprintf("  Emerged (BOTCv-open)  : %d\n", sum(emrg_dir==lbl_emrg_B)))

cat("\n  Interesting fate breakdown:\n")
ft_summ <- fate_df %>%
  filter(direction %in% c(lbl_grpA,lbl_grpB), fate_base!="Converged") %>%
  dplyr::count(direction, fate_base) %>%
  mutate(pct_of_direction=round(100*n / ifelse(direction==lbl_grpA, n_botv, n_botcv), 2))
print(as.data.frame(ft_summ), row.names=FALSE)
cat(strrep("=",65), "\n\n")
