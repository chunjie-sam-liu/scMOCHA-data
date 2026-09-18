#!/usr/bin/env Rscript
# Metainfo ----------------------------------------------------------------
# @AUTHOR: Chun-Jie Liu
# @CONTACT: chunjie.sam.liu.at.gmail.com
# @DATE: 2026-09-09
# @DESCRIPTION: Stage constants, paths and shared helpers for the comparison
#               between scMOCHA-updated mgatk and original mgatk variant
#               calling. Sourced by every step in this stage.
# @VERSION: v0.2.0

# Stage identity -----------------------------------------------------------
STAGE <- "compare-with-mgatk"

# Every sample compared by this stage. The single source of truth for which
# samples exist: the extraction step, the per-sample steps and the
# cross-sample step all read this table. sample_id normalises the "-" that two
# archive names carry, so one token is safe as a path, a factor level and a
# file name; archive keeps the real filename.
SAMPLES <- data.table::data.table(
  sample_id = c(
    "GSE149689_GSM4509019_3PV3",
    "GSE163314_GSM4976997_3PV2",
    "GSE163668_GSM4995445_5PR2",
    "GSE181279_GSM5494116_5PPE",
    "GSE271107_GSM8369876_3PV3"
  ),
  archive = c(
    "GSE149689_GSM4509019_3PV3.zip",
    "GSE163314_GSM4976997_3PV2.zip",
    "GSE163668-GSM4995445_5PR2.zip",
    "GSE181279-GSM5494116_5PPE.zip",
    "GSE271107_GSM8369876_3PV3.zip"
  ),
  gse = c("GSE149689", "GSE163314", "GSE163668", "GSE181279", "GSE271107"),
  gsm = c(
    "GSM4509019",
    "GSM4976997",
    "GSM4995445",
    "GSM5494116",
    "GSM8369876"
  ),
  chemistry = c("SC3Pv3", "SC3Pv2", "SC5P-R2", "SC5P-PE", "SC3Pv3")
)

SAMPLE_IDS <- SAMPLES$sample_id

# Output leaf for everything that spans samples rather than describing one.
CROSS_SAMPLE <- "cross-sample"

# Set per run by fn_sample_label(); every figure subtitle reads it through
# fn_sample_note().
SAMPLE_LABEL <- "<all samples>"

# The seven files this stage reads out of each archive. The allele-count
# matrices, the two RDS objects and the PNGs are never opened here.
ARCHIVE_MEMBERS <- c(
  "cell.variant_stats.tsv.gz",
  "cell.variant_stats_mgatk_original.tsv.gz",
  "cell.cell_heteroplasmic_df.tsv.gz",
  "cell.cell_heteroplasmic_df_mgatk_original.tsv.gz",
  "cell.cell_heteroplasmic_df_raw.tsv.gz",
  "cell.depthTable.txt",
  "cell.coverage.txt.gz"
)

# Paths --------------------------------------------------------------------
# .env stores REPODIR and HIGHRESDIR with a leading "~", so every env-derived
# path is expanded before use. sample_id = NULL selects the cross-sample leaf.
stage_paths <- function(sample_id = NULL) {
  repodir <- fs::path_expand(Sys.getenv("REPODIR"))
  root <- fs::path_expand(fs::path(Sys.getenv("ISILON_BASE"), STAGE))
  stagedir <- fs::path(repodir, "newplots", STAGE)
  leaf <- if (is.null(sample_id)) CROSS_SAMPLE else sample_id

  list(
    repodir = repodir,
    root = root,
    indir = fs::path(root, "samples", leaf),
    cachedir = fs::path(root, "derived", leaf),
    stagedir = stagedir,
    figroot = fs::path(stagedir, "figures"),
    tabroot = fs::path(stagedir, "tables"),
    figdir = fs::path(stagedir, "figures", leaf),
    tabdir = fs::path(stagedir, "tables", leaf),
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

# Stops with the valid list rather than letting a typo build paths under a
# directory that will never exist.
fn_check_sample <- function(sample_id) {
  if (length(sample_id) != 1L || !nzchar(sample_id)) {
    stop(
      "--sample is required. One of: ",
      paste(SAMPLE_IDS, collapse = ", ")
    )
  }
  if (!sample_id %in% SAMPLE_IDS) {
    stop(
      "unknown sample '",
      sample_id,
      "'. One of: ",
      paste(SAMPLE_IDS, collapse = ", ")
    )
  }
  sample_id
}

fn_sample_label <- function(sample_id) {
  i <- match(sample_id, SAMPLES$sample_id)
  as.character(glue::glue(
    "{SAMPLES$gse[i]} {SAMPLES$gsm[i]} ({SAMPLES$chemistry[i]})"
  ))
}

# A panel for an arm or a group that turned out to be empty. Four of the five
# samples have zero variants in the original-mgatk arm, and a zero there is a
# result: the figure has to be written and say so, not be skipped.
fn_empty_panel <- function(title, subtitle, note) {
  ggplot2::ggplot() +
    ggplot2::annotate(
      "text",
      x = 0,
      y = 0,
      label = paste(strwrap(note, width = 46), collapse = "\n"),
      size = 4,
      color = "grey35"
    ) +
    ggplot2::xlim(-1, 1) +
    ggplot2::ylim(-1, 1) +
    fn_theme() +
    ggplot2::theme(
      axis.text = ggplot2::element_blank(),
      axis.ticks = ggplot2::element_blank(),
      panel.border = ggplot2::element_blank()
    ) +
    ggplot2::labs(title = title, subtitle = subtitle, x = NULL, y = NULL)
}

# Returns `plot` when `ok`, otherwise the placeholder. `plot` is a promise, so
# a panel that would error on an empty group is never evaluated.
fn_or_empty <- function(ok, plot, title, note) {
  if (isTRUE(ok)) plot else fn_empty_panel(title, fn_sample_note(), note)
}

# TRUE when a two-group comparison has enough variants on both sides to be
# worth testing. Below this the medians are still reported, with P = NA.
CUTOFF_MIN_GROUP <- 3

fn_testable <- function(g) {
  n <- table(droplevels(as.factor(g)))
  length(n) == 2L && all(n >= CUTOFF_MIN_GROUP)
}

# Every per-sample table carries its sample id in the first column, so the
# cross-sample step can bind them without re-deriving provenance.
fn_export_tab <- function(d, paths, name, sample_id) {
  out <- data.table::as.data.table(d)
  out <- data.table::copy(out)
  data.table::set(out, j = "sample", value = sample_id)
  data.table::setcolorder(out, "sample")
  export(out, as.character(fs::path(paths$tabdir, name)))
  invisible(out)
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
  d <- detection_long[depth >= CUTOFF_MIN_READS & (is.na(af) | af < 0.01)]
  total_depth <- sum(d$depth, na.rm = TRUE)
  rate <- if (total_depth > 0) {
    sum(round(d$af * d$depth), na.rm = TRUE) / total_depth
  } else {
    NA_real_
  }

  # A shallow sample can carry no background alt read, or no background cell at
  # all. Either way the estimate is 0 or undefined, which makes the binomial
  # carrier test vacuous: with rate 0 every cell holding one alt read passes,
  # and with the degenerate rate 1 no cell ever does. Fall back to the smallest
  # rate the reads in hand could have resolved.
  if (!is.finite(rate) || rate <= 0) {
    resolvable <- max(total_depth, sum(detection_long$depth, na.rm = TRUE), 1)
    rate <- 1 / resolvable
    log_warn(
      "background alt rate not estimable from {total_depth} background ",
      "reads; floored at 1/{resolvable} = {signif(rate, 3)}"
    )
  }
  rate
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
