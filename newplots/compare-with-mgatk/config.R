#!/usr/bin/env Rscript
# Metainfo ----------------------------------------------------------------
# @AUTHOR: Chun-Jie Liu
# @CONTACT: chunjie.sam.liu.at.gmail.com
# @DATE: 2026-09-09
# @DESCRIPTION: Stage constants, paths and shared helpers for the comparison
#               between scMOCHA-updated mgatk and original mgatk variant
#               calling. Sourced by every step in this stage.
# @VERSION: v0.1.0

# Stage identity -----------------------------------------------------------
STAGE <- "compare-with-mgatk"

# Set once the GSE/GSM/SRR identifier of the compared sample is known; every
# figure subtitle reads it from here.
SAMPLE_LABEL <- "<pending>"

# Paths --------------------------------------------------------------------
# .env stores REPODIR and HIGHRESDIR with a leading "~", so every env-derived
# path is expanded before use.
stage_paths <- function() {
  repodir <- fs::path_expand(Sys.getenv("REPODIR"))
  indir <- fs::path_expand(fs::path(Sys.getenv("ISILON_BASE"), STAGE))
  stagedir <- fs::path(repodir, "newplots", STAGE)

  list(
    repodir = repodir,
    indir = indir,
    cachedir = fs::path(indir, "derived"),
    stagedir = stagedir,
    figdir = fs::path(stagedir, "figures"),
    tabdir = fs::path(stagedir, "tables"),
    colorfile = fs::path(
      fs::path_expand(Sys.getenv("HIGHRESDIR")),
      "00-colors.R"
    )
  )
}

# Input files --------------------------------------------------------------
stage_inputs <- function(paths = stage_paths(), prefix = "cell") {
  list(
    stats_scmocha = fs::path(
      paths$indir,
      "{prefix}.variant_stats.tsv.gz" |>
        glue::glue()
    ),
    stats_mgatk = fs::path(
      paths$indir,
      "{prefix}.variant_stats_mgatk_original.tsv.gz" |> glue::glue()
    ),
    af_scmocha = fs::path(
      paths$indir,
      "{prefix}.cell_heteroplasmic_df.tsv.gz" |> glue::glue()
    ),
    af_mgatk = fs::path(
      paths$indir,
      "{prefix}.cell_heteroplasmic_df_mgatk_original.tsv.gz" |> glue::glue()
    ),
    af_raw = fs::path(
      paths$indir,
      "{prefix}.cell_heteroplasmic_df_raw.tsv.gz" |> glue::glue()
    ),
    depth = fs::path(paths$indir, "{prefix}.depthTable.txt" |> glue::glue()),
    coverage = fs::path(paths$indir, "{prefix}.coverage.txt.gz" |> glue::glue())
  )
}

# Criterion constants ------------------------------------------------------
# Original mgatk, from newplots/thecode/original-mgatk-variant-calling.py and
# the canonical mgatk/Signac downstream gate.
CUTOFF_CELL_MEANCOV <- 10 # cells kept when mean MT coverage > this
CUTOFF_ALT_STRAND <- 2 # alt reads required on each strand
CUTOFF_NCELLS_CONF <- 3 # cells required per retained variant
CUTOFF_VMR_MGATK <- 0.01
CUTOFF_STRAND_MGATK <- 0.65

# scMOCHA, from newplots/thecode/scmocha-mgatk-variant-calling.py.
CUTOFF_ALT_READS <- 10 # fwd + rev alt reads required in a confident cell

# scMOCHA reliability gate, copied verbatim from
# src/06.1-collect-variants-new.R. Keep in step with that file; it is the
# production classification script and this stage must not drift from it.
CUTOFF_HETEROPLASMIC <- 0.05
CUTOFF_HOMOPLASMIC <- 0.95
CUTOFF_MIN_READS <- 10
CUTOFF_NOTRELIABLE <- 10

# Read-support floor for the loose carrier definition, kept only so the
# sensitivity table can show what that definition does to the answer.
CUTOFF_MIN_ALT_OBS <- 2

POS_RNA_EDITING <- c(
  585,
  1610,
  3238,
  4271,
  5520,
  7526,
  8303,
  9999,
  10413,
  12146,
  12274,
  14734,
  15896,
  295,
  2617,
  13710
)

POS_MISSALIGNMENT_ERROR <- c(
  66:71,
  300:316,
  513:525,
  3106:3107,
  12418:12425,
  16182:16194
)

# Display ------------------------------------------------------------------
AF_BIN_BREAKS <- c(0, 0.01, 0.05, 0.10, 0.20, 0.50, 1)
AF_BIN_LABELS <- c(
  "<0.01",
  "0.01-0.05",
  "0.05-0.10",
  "0.10-0.20",
  "0.20-0.50",
  "\u22650.50"
)
AF_BIN_NONE <- "No detected cell"

CALLER_LEVELS <- c("Original mgatk", "scMOCHA")
VARIANT_SET_LEVELS <- c("Both", "Original mgatk only", "scMOCHA only")

# Three calling arms. scMOCHA's variant call applies no AF filter; the 0.05 AF
# rule belongs to the downstream analysis in src/06.1-collect-variants-new.R,
# so it is a separate arm rather than part of the call.
ARM_MGATK <- "Original mgatk"
ARM_SCMOCHA <- "scMOCHA variant call"
ARM_SCMOCHA_AF5 <- "scMOCHA AF>5%"
ARM_LEVELS <- c(ARM_MGATK, ARM_SCMOCHA, ARM_SCMOCHA_AF5)

EXCLUSION_LEVELS <- c(
  "Retained",
  "Not proposed as candidate",
  "< 3 confident cells",
  "VMR and strand r",
  "VMR only",
  "strand r only",
  "Blacklisted position",
  "< 10 cells at AF >= 0.05"
)

# Helpers ------------------------------------------------------------------
fn_theme <- function() {
  ggplot2::theme_bw(base_size = 12) +
    ggplot2::theme(
      panel.grid = ggplot2::element_blank(),
      strip.background = ggplot2::element_rect(fill = "grey92"),
      plot.title = ggplot2::element_text(face = "bold", hjust = 0.5),
      plot.subtitle = ggplot2::element_text(
        hjust = 0.5,
        color = "grey40",
        size = 10
      )
    )
}

fn_af_bin <- function(af) {
  out <- as.character(cut(
    af,
    breaks = AF_BIN_BREAKS,
    labels = AF_BIN_LABELS,
    include.lowest = TRUE,
    right = FALSE
  ))
  out[is.na(af)] <- AF_BIN_NONE
  factor(out, levels = c(AF_BIN_LABELS, AF_BIN_NONE))
}

fn_sample_note <- function() {
  glue::glue("sample {SAMPLE_LABEL}")
}

# Subtitles here carry N and every threshold, so they routinely run past the
# panel width. ggplot2 does not wrap them; applied at save time.
fn_wrap_labs <- function(p, width = 95) {
  if (!is.null(p$labels$subtitle)) {
    p$labels$subtitle <- paste(
      strwrap(p$labels$subtitle, width = width),
      collapse = "\n"
    )
  }
  p
}

# The first failed criterion in each arm's own order, so every excluded variant
# is attributed to exactly one cutoff. Order matters: a variant that is never
# proposed cannot also fail a later gate, and a variant failing both VMR and
# strand correlation is reported as failing both rather than arbitrarily one.
fn_exclusion_mgatk <- function(d) {
  factor(
    data.table::fcase(
      d$final_mgatk                                                           , "Retained"                  ,
      !d$candidate_mgatk                                                      , "Not proposed as candidate" ,
      !d$s1_mgatk                                                             , "< 3 confident cells"       ,
      d$vmr_mgatk <= CUTOFF_VMR_MGATK & d$strand_mgatk <= CUTOFF_STRAND_MGATK , "VMR and strand r"          ,
      d$vmr_mgatk <= CUTOFF_VMR_MGATK                                         , "VMR only"                  ,
      default = "strand r only"
    ),
    levels = EXCLUSION_LEVELS
  )
}

fn_exclusion_scmocha_af5 <- function(d) {
  factor(
    data.table::fcase(
      d$final_scmocha      , "Retained"                  ,
      !d$candidate_scmocha , "Not proposed as candidate" ,
      !d$s1_scmocha        , "< 3 confident cells"       ,
      d$blacklisted        , "Blacklisted position"      ,
      default = "< 10 cells at AF >= 0.05"
    ),
    levels = EXCLUSION_LEVELS
  )
}

fn_exclusion_scmocha_call <- function(d) {
  factor(
    data.table::fcase(
      d$s1_scmocha         , "Retained"                  ,
      !d$candidate_scmocha , "Not proposed as candidate" ,
      default = "< 3 confident cells"
    ),
    levels = EXCLUSION_LEVELS
  )
}

# TRUE for a variant position excluded by the scMOCHA position blacklist.
fn_is_blacklisted <- function(position) {
  position %in% POS_RNA_EDITING | position %in% POS_MISSALIGNMENT_ERROR
}

# Long cell-by-variant table carrying AF, the cell's depth at that variant's
# position, and the scMOCHA per-cell detection flag. A (position, barcode)
# pair absent from the coverage file means no reads, so depth is 0.
fn_detection_long <- function(cell_af, cell_pos_coverage, variant_position) {
  long <- data.table::melt(
    data.table::as.data.table(cell_af),
    id.vars = "barcode",
    variable.name = "variant",
    value.name = "af",
    variable.factor = FALSE
  )
  long[
    data.table::as.data.table(variant_position),
    on = "variant",
    position := i.position
  ]
  long[
    data.table::as.data.table(cell_pos_coverage),
    on = c("position", "barcode"),
    depth := i.depth
  ]
  long[is.na(depth), depth := 0]
  long[,
    detected := !is.na(af) &
      af >= CUTOFF_HETEROPLASMIC &
      depth >= CUTOFF_MIN_READS
  ]
  long[]
}

# scMOCHA reliability gate: not blacklisted, and detected in enough cells.
fn_scmocha_gate <- function(detection_long, blacklisted_variants) {
  n_cells <- detection_long[detected == TRUE, .(n_cells = .N), by = variant]
  n_cells[
    n_cells >= CUTOFF_NOTRELIABLE & !variant %in% blacklisted_variants,
    variant
  ]
}

# Background alt-read rate: alt reads over total reads across every covered
# cell that is nowhere near carrying the variant, including the zero-alt cells.
# Excluding them would bias the estimate upward.
fn_background_rate <- function(detection_long) {
  detection_long[
    depth >= CUTOFF_MIN_READS & (is.na(af) | af < 0.01),
    sum(round(af * depth), na.rm = TRUE) / sum(depth, na.rm = TRUE)
  ]
}

# Heteroplasmy among the cells that carry the variant, on three definitions.
# Which one is used changes the answer, so the choice is stated everywhere.
#
#   af_gate_*    cells passing the scMOCHA gate (AF >= 0.05, depth >= 10).
#                Censored at 0.05 by construction: it cannot see the
#                low-heteroplasmy region this stage exists to examine.
#   af_loose_*   cells with >= 2 alt reads. Uncensored but noise-dominated:
#                2 alt reads at depth 1000 is expected by chance at a
#                background rate of 0.03%, so the median over these cells
#                describes the background, not the variant.
#   af_carrier_* cells whose alt count is inconsistent with the background
#                error rate under a binomial test, Bonferroni-corrected over
#                every cell-by-variant observation. Uncensored AND noise-
#                resistant. This is the primary measure.
#
# The mgatk `mean` column is not a fourth option: it is total alt reads over
# total coverage across all cells, so it tracks prevalence, not heteroplasmy.
fn_carrier_stats <- function(detection_long, background_rate) {
  d <- detection_long
  alpha <- 0.05 / nrow(d)

  d[, alt_reads := round(af * depth)]
  d[,
    is_carrier := depth >= CUTOFF_MIN_READS &
      !is.na(alt_reads) &
      alt_reads >= 1 &
      stats::pbinom(
        alt_reads - 1,
        depth,
        background_rate,
        lower.tail = FALSE
      ) <
        alpha
  ]

  gated <- d[
    detected == TRUE,
    .(
      n_cells_gate = .N,
      af_gate_median = stats::median(af),
      af_gate_max = max(af)
    ),
    by = variant
  ]

  loose <- d[
    depth >= CUTOFF_MIN_READS & !is.na(af) & alt_reads >= CUTOFF_MIN_ALT_OBS,
    .(
      n_cells_loose = .N,
      af_loose_median = stats::median(af)
    ),
    by = variant
  ]

  carrier <- d[
    is_carrier == TRUE,
    .(
      n_cells_carrier = .N,
      af_carrier_median = stats::median(af),
      af_carrier_min = min(af),
      af_carrier_max = max(af)
    ),
    by = variant
  ]

  Reduce(
    function(x, y) merge(x, y, by = "variant", all = TRUE),
    list(gated, loose, carrier)
  )
}
