"""
Config_SampleMetadata.py  —  Single source of truth for all validated sample sets.

Canonical genotype vocabulary (use these labels everywhere downstream)
----------------------------------------------------------------------
  WT          — wild type (no BOT, no toll)
  runt        — runt truncation mutant only (no BOT background, no toll)
  BOTR   — BOT background (bcd- osk- tsl-) + runt mutation
                (was "BOTR" / "runtD7" in old naming)
  FoxL1_BOTv  — BOT background + toll ventralization
                (was "foxl1" / "foxl1-high" in old naming)
  B6_BOTC    — BOTC (BOT+cic) background, non-ventralized
                (was "B6" in old naming)
  HLH54F_BOTCv— BOTC background + toll ventralization
                (was "HLH54F" in old naming)

The _BOTv / _BOTCv suffixes encode the maternal background:
  BOT   = bcd- osk- tsl-  (Bicoid, Oskar, Torso-like maternal triple mutant)
  BOTC  = bcd- osk- tsl- cic-  (BOT + capicua)

D/V classification
------------------
  Only FoxL1_BOTv and HLH54F_BOTCv carry the toll allele and are classified
  for ventralization. WT, runt, BOTR, and B6_BOTC lack toll — their
  DV score is suppressed.

Runt expression logic (genotype-aware)
--------------------------------------
  WT / B6_BOTC    : runt at baseline (~WT level)
  runt / BOTR : runt LOW  → true runt mutant allele
  FoxL1_BOTv       : runt HIGH → well-ventralized (ubiquitous)
                     runt LOW  → incomplete ventralization
  HLH54F_BOTCv     : runt LOW  → well-ventralized
                     runt HIGH → NOT fully ventralized

NOTE: BOTR classification uses run TPM only (not RAI) because the BOT
maternal background independently abolishes AP patterning genes (ftz, eve,
slp1 etc.) making RAI uninformative for this genotype.

SET 1  (WT_seRNAseq_4embryos_set1 / Mutant_seRNAseq_19embryos_set2)
---------------------------------------------------------------------
  WT          : WT_embryo_bc2, bc3, bc4
  FoxL1_BOTv  : FoxL1_BOTv_embryo_bc1/bc2/bc8  (nc14b; formerly Mutant_embryo_bc*)
  runt        : runt_embryo_bc13/bc14/bc18      (nc14b; formerly Mutant_embryo_bc*)

SET 2  (28232_28233_28271_28272_28273_S1)
-----------------------------------------
  Plates    Barcodes   Genotype       Stage    N
  P1-9,P16  bc01-09,16 runt           nc14b   10
  P10-15    bc10-15    BOTR      nc14b    6
  P19-21    bc19-21    B6_BOTC       nc14d    3
  P22-24    bc22-24    B6_BOTC       gastr    3
  P25-28    bc25-28    FoxL1_BOTv     nc14d    4
  P29-30    bc29-30    FoxL1_BOTv     gastr    2
  P32-35    bc32-35    HLH54F_BOTCv   nc14b    4
  P40-41    bc40-41    HLH54F_BOTCv   nc14b    2
  P36-37    bc36-37    HLH54F_BOTCv   nc14d    2
  P38-39    bc38-39    HLH54F_BOTCv   gastr    2

SET 3  (28348_28349_28350_28351_S1)
------------------------------------
  IMPORTANT: this pooled run REUSES barcode numbers already assigned to
  specific genotypes in SET2_BARCODE_MAP (bc1-6, bc13-16). Genotype for
  Set3 columns must come ONLY from the "<genotype>_<stage>_bc<N>" name
  written by step1_demux.py / barcodes_set3.tsv — NEVER from a bare
  barcode-number lookup against SET2_BARCODE_MAP. See build_group_map().

  Plates  Barcodes  Genotype  Stage   N
  P1-3    bc1-3     yw        nc14b   3   -> WT
  P4-6    bc4-6     BOTv      nc14b   3   -> FoxL1_BOTv
  P13-16  bc13-16   run       nc14b   4   -> runt
  P42-45  bc42-45   BOT       nc14b   4   -> WT_D7 (new group, distinct
  P46-47  bc46-47   BOT       nc14d   2      from all prior groups)
  P54-56  bc54-56   BOTC      nc14d   3   -> B6_BOTC
"""

# ── Sequencing run / sample names ─────────────────────────────────────────────
WT_SAMPLE     = "WT_seRNAseq_4embryos_set1"
MUTANT_SAMPLE = "Mutant_seRNAseq_19embryos_set2"
SET2_SAMPLE   = "28232_28233_28271_28272_28273_S1"
SET3_SAMPLE   = "28348_28349_28350_28351_S1"

ALL_SAMPLES = [WT_SAMPLE, MUTANT_SAMPLE, SET2_SAMPLE, SET3_SAMPLE]

# ── Canonical group labels ─────────────────────────────────────────────────────
GROUPS = ["WT", "runt", "BOTR", "FoxL1_BOTv", "B6_BOTC", "HLH54F_BOTCv", "WT_D7"]

# ── Genotype prefix → canonical group label ────────────────────────────────────
# Matrix column names still use old demux labels (foxl1_Dm_, runtD7_Dm_ etc.)
# This map translates them to the new canonical names.
GENOTYPE_PREFIX_TO_GROUP: dict[str, str] = {
    # ── Set2 canonical prefixes (post-rename) ─────────────────────────────────
    "BOTR_Dm"          : "BOTR",
    "runt_Dm"          : "runt",
    "B6_BOTC_Dm"       : "B6_BOTC",
    "FoxL1_BOTv_Dm"    : "FoxL1_BOTv",
    "HLH54F_BOTCv_Dm"  : "HLH54F_BOTCv",
    # ── Set2 legacy demux labels (pre-rename — present in current matrix) ─────
    # These are the labels the demultiplexer used before rename_matrix_columns.py
    # was run.  Keep until all matrices are regenerated with canonical names.
    "foxl1_Dm"         : "FoxL1_BOTv",   # foxl1_Dm_nc14_bcN  / foxl1_Dm_nc14d_bcN / foxl1_Dm_gastr_bcN
    "HLH54F_Dm"        : "HLH54F_BOTCv", # HLH54F_Dm_nc14b_bcN / HLH54F_Dm_nc14d_bcN / HLH54F_Dm_gastr_bcN
    "runtD7_Dm"        : "BOTR",         # runtD7_Dm_nc14b_bcN  (runtD7 = BOTR in canonical vocabulary)
    "B6_Dm"            : "B6_BOTC",      # B6_Dm_nc14d_bcN / B6_Dm_gastr_bcN
    # ── Set1 prefixes ─────────────────────────────────────────────────────────
    "WT_WT_embryo"     : "WT",
    "WT_embryo"        : "WT",
    "FoxL1_BOTv_embryo": "FoxL1_BOTv",
    "runt_embryo"      : "runt",
    # Legacy set1 bare Mutant_embryo_ fallback (should not appear after rename)
    "Mutant_embryo"    : "unknown",
    # ── Set3 prefixes (28348_28349_28350_28351_S1) ────────────────────────────
    # NOTE: Set3 REUSES barcode numbers already assigned in SET2_BARCODE_MAP
    # (bc1-6, bc13-16). Resolution for these MUST come from this prefix map,
    # never from the barcode-number lookup — see build_group_map()'s scoping
    # of _barcode_to_genotype() to SET2_SAMPLE columns only.
    "yw_"              : "WT",
    "BOTv_"            : "FoxL1_BOTv",
    "run_"             : "runt",
    "BOT_"             : "WT_D7",   # new group — distinct from all prior BOT* groups
    "BOTC_"            : "B6_BOTC",
}

# ── Set1 embryo → stage ────────────────────────────────────────────────────────
SET1_EMBRYO_STAGE: dict[str, str] = {
    "WT_embryo_bc2":       "nc14b",
    "WT_embryo_bc3":       "nc14b",
    "WT_embryo_bc4":       "nc14b",
    "FoxL1_BOTv_embryo_bc1":   "nc14b",   # FoxL1_BOTv
    "FoxL1_BOTv_embryo_bc2":   "nc14b",   # FoxL1_BOTv
    "FoxL1_BOTv_embryo_bc8":   "nc14b",   # FoxL1_BOTv
    "runt_embryo_bc13":  "nc14b",   # runt
    "runt_embryo_bc14":  "nc14b",   # runt
    "runt_embryo_bc18":  "nc14b",   # runt
}

# ── Set2 barcode map ───────────────────────────────────────────────────────────
SET2_BARCODE_MAP: dict[str, tuple[str, str, str]] = {
    # runt nc14b
    "AGTGTC": ("bc01", "runt",          "nc14b"), "ACCATG": ("bc02", "runt",          "nc14b"),
    "GAGTGA": ("bc03", "runt",          "nc14b"), "CACTCA": ("bc04", "runt",          "nc14b"),
    "CATGTC": ("bc05", "runt",          "nc14b"), "ACAGGA": ("bc06", "runt",          "nc14b"),
    "GTACCA": ("bc07", "runt",          "nc14b"), "ACAGAC": ("bc08", "runt",          "nc14b"),
    "ACGTTG": ("bc09", "runt",          "nc14b"), "CTAGGA": ("bc16", "runt",          "nc14b"),
    # BOTR nc14b
    "ACCAAC": ("bc10", "BOTR",     "nc14b"), "GTGAAG": ("bc11", "BOTR",     "nc14b"),
    "CACTTC": ("bc12", "BOTR",     "nc14b"), "GAGTTG": ("bc13", "BOTR",     "nc14b"),
    "GAAGAC": ("bc14", "BOTR",     "nc14b"), "TGCAGA": ("bc15", "BOTR",     "nc14b"),
    # B6_BOTC nc14d
    "CTAGAC": ("bc19", "B6_BOTC",      "nc14d"), "AGCTCA": ("bc20", "B6_BOTC",      "nc14d"),
    "ACTCGA": ("bc21", "B6_BOTC",      "nc14d"),
    # B6_BOTC gastr
    "CTGTTG": ("bc22", "B6_BOTC",      "gastr"), "CATGCA": ("bc23", "B6_BOTC",      "gastr"),
    "CAGAAG": ("bc24", "B6_BOTC",      "gastr"),
    # FoxL1_BOTv nc14d
    "GTCTCA": ("bc25", "FoxL1_BOTv",    "nc14d"), "GTGATC": ("bc26", "FoxL1_BOTv",    "nc14d"),
    "TGTCTG": ("bc27", "FoxL1_BOTv",    "nc14d"), "GACAGA": ("bc28", "FoxL1_BOTv",    "nc14d"),
    # FoxL1_BOTv gastr
    "ACTCTG": ("bc29", "FoxL1_BOTv",    "gastr"), "TGCAAC": ("bc30", "FoxL1_BOTv",    "gastr"),
    # HLH54F_BOTCv nc14b
    "GTTGAG": ("bc32", "HLH54F_BOTCv",  "nc14b"), "AGACCA": ("bc33", "HLH54F_BOTCv",  "nc14b"),
    "TGGTTG": ("bc34", "HLH54F_BOTCv",  "nc14b"), "GATCTG": ("bc35", "HLH54F_BOTCv",  "nc14b"),
    "GATCGA": ("bc40", "HLH54F_BOTCv",  "nc14b"), "GTACTC": ("bc41", "HLH54F_BOTCv",  "nc14b"),
    # HLH54F_BOTCv nc14d
    "CTAGTG": ("bc36", "HLH54F_BOTCv",  "nc14d"), "CTCAGA": ("bc37", "HLH54F_BOTCv",  "nc14d"),
    # HLH54F_BOTCv gastr
    "CTTCGA": ("bc38", "HLH54F_BOTCv",  "gastr"), "AGCTAG": ("bc39", "HLH54F_BOTCv",  "gastr"),
}

SUBLIBRARY_GENOTYPE: dict[str, str] = {
    "28271": "runt", "28272": "runt", "28232": "BOTR",
    "28233": "HLH54F_BOTCv", "28273": "FoxL1_BOTv_B6_BOTC",
}

# ── Validated embryo lists ─────────────────────────────────────────────────────
VALIDATED_WT: list[str] = [
    "WT_embryo_bc2",
    "WT_embryo_bc3",
    "WT_embryo_bc4",
]

VALIDATED_RUNT: list[str] = [
    # Set1
    "runt_embryo_bc13",
    "runt_embryo_bc14",
    "runt_embryo_bc18",
    # Set2 canonical
    "runt_Dm_nc14b_bc07",
    "runt_Dm_nc14b_bc09",
    # Set2 legacy (bare nc14 stage tag — same barcodes, different stage string)
    "runt_Dm_nc14_bc07",
    "runt_Dm_nc14_bc09",
]

VALIDATED_BOTR: list[str] = [
    # Canonical names (post rename_matrix_columns.py)
    "BOTR_Dm_nc14b_bc11",
    "BOTR_Dm_nc14b_bc12",
    "BOTR_Dm_nc14b_bc13",
    "BOTR_Dm_nc14b_bc15",
    # Legacy demux names (pre-rename — still present in current matrix)
    "runtD7_Dm_nc14b_bc11",
    "runtD7_Dm_nc14b_bc12",
    "runtD7_Dm_nc14b_bc13",
    "runtD7_Dm_nc14b_bc15",
]

VALIDATED_FOXL1: list[str] = [
    # Set1 nc14b
    "FoxL1_BOTv_embryo_bc1",
    "FoxL1_BOTv_embryo_bc2",
    "FoxL1_BOTv_embryo_bc8",
    # Set2 nc14b — canonical names
    "FoxL1_BOTv_Dm_nc14b_bc20",
    "FoxL1_BOTv_Dm_nc14b_bc21",
    "FoxL1_BOTv_Dm_nc14b_bc22",
    "FoxL1_BOTv_Dm_nc14b_bc23",
    # NOTE: foxl1_Dm_nc14_bc20-23 are B6_BOTC (not FoxL1_BOTv) per SET2_BARCODE_MAP
    # They are listed in VALIDATED_B6 and must NOT appear here.
    # Set2 nc14d — canonical
    "FoxL1_BOTv_Dm_nc14d_bc25",
    "FoxL1_BOTv_Dm_nc14d_bc27",
    "FoxL1_BOTv_Dm_nc14d_bc28",
    # Set2 nc14d — legacy
    "foxl1_Dm_nc14d_bc25",
    "foxl1_Dm_nc14d_bc27",
    "foxl1_Dm_nc14d_bc28",
    # Set2 gastr — canonical
    "FoxL1_BOTv_Dm_gastr_bc30",
    # Set2 gastr — legacy
    "foxl1_Dm_gastr_bc30",
]

VALIDATED_B6: list[str] = [
    # nc14d — confirmed real libraries (were mislabelled foxl1_Dm_ by demux)
    "B6_BOTC_Dm_nc14d_bc19",
    "B6_BOTC_Dm_nc14d_bc20",
    "B6_BOTC_Dm_nc14d_bc21",
    # Legacy demux names (foxl1_Dm_ prefix used for whole sublibrary 28273)
    "foxl1_Dm_nc14_bc19",
    "foxl1_Dm_nc14_bc20",
    "foxl1_Dm_nc14_bc21",
    # gastr — confirmed
    "B6_BOTC_Dm_gastr_bc22",
    "B6_BOTC_Dm_gastr_bc23",
    "B6_BOTC_Dm_gastr_bc24",
    "foxl1_Dm_nc14_bc22",
    "foxl1_Dm_nc14_bc23",
    "foxl1_Dm_nc14_bc24",
]

VALIDATED_HLH54F: list[str] = [
    # Canonical names (post rename_matrix_columns.py)
    "HLH54F_BOTCv_Dm_nc14b_bc32",
    "HLH54F_BOTCv_Dm_nc14b_bc33",
    "HLH54F_BOTCv_Dm_nc14b_bc34",
    "HLH54F_BOTCv_Dm_nc14b_bc35",
    "HLH54F_BOTCv_Dm_nc14d_bc36",
    "HLH54F_BOTCv_Dm_gastr_bc39",
    # Legacy demux names (HLH54F_Dm_ prefix — present in current matrix)
    "HLH54F_Dm_nc14b_bc32",
    "HLH54F_Dm_nc14b_bc33",
    "HLH54F_Dm_nc14b_bc34",
    "HLH54F_Dm_nc14b_bc35",
    "HLH54F_Dm_nc14d_bc36",
    "HLH54F_Dm_gastr_bc39",
]

VALIDATED_BY_GROUP: dict[str, list[str]] = {
    "WT":             VALIDATED_WT,
    "runt":           VALIDATED_RUNT,
    "BOTR":      VALIDATED_BOTR,
    "FoxL1_BOTv":     VALIDATED_FOXL1,
    "B6_BOTC":       VALIDATED_B6,
    "HLH54F_BOTCv":   VALIDATED_HLH54F,
}

# ── Manual exclusions ─────────────────────────────────────────────────────────
# Embryos excluded from all downstream analyses (failed demux, sub-threshold
# UMI capture, or duplicate/ambiguous barcode assignment). See step5c QC
# filter output (*_excluded_embryos.tsv) for per-embryo QC metrics.
MANUAL_EXCLUDE: set[str] = {
    # ── Set1 ───────────────────────────────────────────────────────────────
    "WT_WT_embryo_bc1",
    "WT_embryo_bc1",
    "FoxL1_BOTv_embryo_bc4",
    # Full-prefix form needed to catch this one after sample-prefixing;
    # the bare name is a valid embryo elsewhere and must stay excluded here only.
    "WT_seRNAseq_4embryos_set1__FoxL1_BOTv_embryo_bc1",
    # ── Set2 ───────────────────────────────────────────────────────────────
    "runt_Dm_nc14b_bc16",
    "runtD7_Dm_nc14b_bc10",
    "runtD7_Dm_nc14b_bc14",
    "foxl1_Dm_nc14_bc19",
    "foxl1_Dm_nc14_bc24",
    "foxl1_Dm_nc14d_bc26",
    "foxl1_Dm_gastr_bc29",
    "HLH54F_Dm_nc14d_bc37",
    "HLH54F_Dm_gastr_bc38",
    "runt_Dm_nc14_bc41",
    "runt_Dm_nc14_bc40",
    "runt_Dm_nc14b_bc01",
    "runt_Dm_nc14b_bc02",
    "runt_Dm_nc14b_bc03",
    "runt_Dm_nc14b_bc04",
    "runt_Dm_nc14b_bc05",
    "runt_Dm_nc14b_bc06",
}

# ── DV classification (toll allele present) ───────────────────────────────────
DV_CLASSIFIED_GROUPS: set[str] = {"FoxL1_BOTv", "HLH54F_BOTCv"}

# ── Group display colours ─────────────────────────────────────────────────────
GROUP_COLOURS: dict[str, str] = {
    "WT":             "#4878CF",
    "runt":           "#E53935",
    "BOTR":      "#7B1FA2",
    "FoxL1_BOTv":     "#2E7D32",
    "B6_BOTC":       "#00838F",
    "HLH54F_BOTCv":   "#E65100",
    "WT_D7":          "#F9A825",
    "unknown":        "#AAAAAA",
}

STAGE_ALPHA: dict[str, float] = {
    "nc14b": 1.0,
    "nc14d": 0.7,
    "gastr": 0.45,
}

GROUP_ORDER: list[str] = [
    "WT", "runt", "BOTR", "FoxL1_BOTv", "B6_BOTC", "HLH54F_BOTCv", "WT_D7"
]

# ── Library QC thresholds ─────────────────────────────────────────────────────
QC_MIN_TOTAL_UMIS     = 8_000
QC_MIN_DETECTED_GENES = 5_000
QC_MIN_MEDIAN_LOG2TPM = 4.0
QC_EXPANDED_MIN_UMIS  = 3_000
QC_EXPANDED_MIN_GENES = 2_000


# ── Helpers ───────────────────────────────────────────────────────────────────

def _barcode_to_genotype(bare: str) -> "str | None":
    """
    For Set2 columns, look up the exact genotype from SET2_BARCODE_MAP
    using the bcNN suffix.  This is the authoritative source for columns
    from sublibrary 28273 which contains both FoxL1_BOTv (bc25-30) AND
    B6_BOTC (bc19-24) — the prefix-based assignment cannot distinguish them
    because the demux labelled all of them with the foxl1_Dm_ prefix.

    Returns the canonical group string, or None if the barcode is not in the map.
    """
    import re as _re
    m = _re.search(r"bc(\d+)$", bare)
    if m is None:
        return None
    bc_num = int(m.group(1))
    # SET2_BARCODE_MAP keys are hex sequences; values are (bcNN, genotype, stage)
    # Build a bcNN → genotype lookup (lazy, cached on first call)
    if not hasattr(_barcode_to_genotype, "_cache"):
        _barcode_to_genotype._cache = {}
        for _hex, (_bc, _gt, _st) in SET2_BARCODE_MAP.items():
            _n = int(_bc.replace("bc", ""))
            _barcode_to_genotype._cache[_n] = _gt
    return _barcode_to_genotype._cache.get(bc_num)


def build_group_map(columns: list[str]) -> dict[str, str]:
    """Return {col: canonical_group} for every column.

    Resolution order (first match wins):
      1. Explicit validated list — handles both canonical and legacy bare IDs.
      2. Barcode-number lookup via SET2_BARCODE_MAP — authoritative ONLY for
         columns from SET2_SAMPLE, where the demux prefix is ambiguous (e.g.
         foxl1_Dm_ prefix used for both FoxL1_BOTv bc25-30 AND B6_BOTC bc19-24
         in sublibrary 28273). Never applied to other samples, since other
         pooled runs (e.g. Set3) independently reuse the same barcode numbers.
      3. Prefix map sorted longest-first.
      4. 'unknown' fallback.
    """
    # Build bare-id → group lookup from validated lists
    col_to_group: dict[str, str] = {}
    for group, validated in VALIDATED_BY_GROUP.items():
        for bare_id in validated:
            col_to_group[bare_id] = group

    # Sort prefixes by length descending so more-specific entries match first
    sorted_prefixes = sorted(GENOTYPE_PREFIX_TO_GROUP.items(),
                             key=lambda kv: -len(kv[0]))

    result: dict[str, str] = {}
    for col in columns:
        bare = col.split("__")[-1] if "__" in col else col
        sample_prefix = col.split("__")[0] if "__" in col else None

        # 1. Explicit validated list
        if bare in col_to_group:
            result[col] = col_to_group[bare]
            continue

        # 2. Barcode-number lookup — ONLY authoritative for genuine Set2
        #    columns. Other pooled runs (e.g. Set3) independently reuse the
        #    same barcode numbers for entirely different genotypes, so this
        #    lookup must not apply to them or it silently mis-assigns groups.
        if sample_prefix == SET2_SAMPLE:
            bc_genotype = _barcode_to_genotype(bare)
            if bc_genotype is not None:
                result[col] = bc_genotype
                continue

        # 3. Prefix map
        matched = False
        for prefix, group in sorted_prefixes:
            if bare.startswith(prefix):
                result[col] = group
                matched = True
                break

        if not matched:
            result[col] = "unknown"

    return result


def get_stage(col: str) -> str:
    import re as _re
    bare = col.split("__")[-1] if "__" in col else col
    if "gastr" in bare:  return "gastr"
    if "nc14d" in bare:  return "nc14d"
    if "nc14b" in bare:  return "nc14b"
    # Bare "nc14" without b/d suffix: check SET2_BARCODE_MAP for authoritative stage
    if "nc14" in bare:
        m = _re.search(r"bc(\d+)$", bare)
        if m:
            bc_num = int(m.group(1))
            for _hex, (_bc, _gt, _st) in SET2_BARCODE_MAP.items():
                if int(_bc.replace("bc","")) == bc_num:
                    return _st
        return "nc14b"  # fallback
    return SET1_EMBRYO_STAGE.get(bare, "nc14b")


def apply_manual_exclude(columns: list[str], verbose: bool = True) -> list[str]:
    """Remove manually excluded embryos from a column list.

    Matching strategy:
      1. Exact bare-id match (col without sample prefix).
      2. Exact full-column match.
      3. Barcode-number match: if both the column bare-id and an exclusion
         entry share the same bcNN suffix AND the same barcode number, they
         are treated as the same embryo regardless of prefix differences
         (e.g. 'runtD7_Dm_nc14b_bc11' matches exclusion 'BOTR_Dm_nc14b_bc11').
         This is intentionally NOT applied when the barcode number is ambiguous
         (i.e. the same bcNN exists in multiple genotype groups).
    """
    import re as _re

    def _bc_num(s: str) -> str | None:
        m = _re.search(r"bc(\d+)$", s)
        return m.group(1) if m else None

    # Pre-build set of excluded barcode numbers for fast lookup
    excl_bc_nums: set[str] = set()
    for excl in MANUAL_EXCLUDE:
        bc = _bc_num(excl)
        if bc:
            excl_bc_nums.add(bc)

    def _is_excluded(col: str) -> bool:
        bare = col.split("__")[-1] if "__" in col else col
        # 1. Exact match
        if bare in MANUAL_EXCLUDE or col in MANUAL_EXCLUDE:
            return True
        # 2. Legacy _Dm_ suffix match
        for excl in MANUAL_EXCLUDE:
            if "_Dm_" in excl and col.endswith(f"__{excl}"):
                return True
        # NOTE: barcode-number matching disabled — caused false positives for
        # low bc numbers (bc1, bc2) shared across groups. Exact name matching
        # above is sufficient for all current exclusion cases.
        return False

    kept    = [c for c in columns if not _is_excluded(c)]
    removed = [c for c in columns if     _is_excluded(c)]
    if verbose and removed:
        print(f"[sample_config] Excluded: {[c.split('__')[-1] for c in removed]}")
    return kept


def get_all_validated() -> list[str]:
    out = []
    for lst in VALIDATED_BY_GROUP.values():
        out.extend(lst)
    return out


def validate_set2_barcodes(tsv_path: str) -> bool:
    import csv as _csv
    mismatches = 0
    with open(tsv_path) as fh:
        for lineno, row in enumerate(_csv.reader(fh, delimiter="\t"), 1):
            if not row or row[0].startswith("#"):
                continue
            bc  = row[0].upper()
            gt  = row[2] if len(row) > 2 else ""
            st  = row[3] if len(row) > 3 else ""
            expected = SET2_BARCODE_MAP.get(bc)
            if expected is None:
                print(f"[validate] Line {lineno}: {bc!r} not in SET2_BARCODE_MAP")
                mismatches += 1
            else:
                exp_gt, exp_st = expected[1], expected[2]
                if gt and gt != exp_gt:
                    print(f"[validate] Line {lineno}: {bc} genotype mismatch "
                          f"(TSV={gt!r}, expected={exp_gt!r})")
                    mismatches += 1
                if st and st != exp_st:
                    print(f"[validate] Line {lineno}: {bc} stage mismatch "
                          f"(TSV={st!r}, expected={exp_st!r})")
                    mismatches += 1
    if mismatches == 0:
        print(f"[validate] {tsv_path}: all barcodes OK")
    return mismatches == 0
