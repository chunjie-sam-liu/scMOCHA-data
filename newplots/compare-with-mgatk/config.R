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
    "GSE155673_GSM4712895_3PV3",
    "GSE163314_GSM4976997_3PV2",
    "GSE163668_GSM4995445_5PR2",
    "GSE175499_GSM5335510_3PV3",
    "GSE181279_GSM5494116_5PPE",
    "GSE188632_GSM5687372_3PV3",
    "GSE220189_GSM6793474_3PV3",
    "GSE271107_GSM8369876_3PV3",
    "GSE279945_GSM8583916_3PV3"
  ),
  archive = c(
    "GSE149689_GSM4509019_3PV3.zip",
    "GSE155673_GSM4712895_3PV3.zip",
    "GSE163314_GSM4976997_3PV2.zip",
    "GSE163668-GSM4995445_5PR2.zip",
    "GSE175499_GSM5335510_3PV3.zip",
    "GSE181279-GSM5494116_5PPE.zip",
    "GSE188632_GSM5687372_3PV3.zip",
    "GSE220189_GSM6793474_3PV3.zip",
    "GSE271107_GSM8369876_3PV3.zip",
    "GSE279945_GSM8583916_3PV3.zip"
  ),
  gse = c(
    "GSE149689",
    "GSE155673",
    "GSE163314",
    "GSE163668",
    "GSE175499",
    "GSE181279",
    "GSE188632",
    "GSE220189",
    "GSE271107",
    "GSE279945"
  ),
  gsm = c(
    "GSM4509019",
    "GSM4712895",
    "GSM4976997",
    "GSM4995445",
    "GSM5335510",
    "GSM5494116",
    "GSM5687372",
    "GSM6793474",
    "GSM8369876",
    "GSM8583916"
  ),
  chemistry = c(
    "SC3Pv3",
    "SC3Pv3",
    "SC3Pv2",
    "SC5P-R2",
    "SC3Pv3",
    "SC5P-PE",
    "SC3Pv3",
    "SC3Pv3",
    "SC3Pv3",
    "SC3Pv3"
  )
)

SAMPLE_IDS <- SAMPLES$sample_id

# Samples kept out of the cross-sample FIGURES only. The 5' libraries are the
# two odd ones out: restricting the panels to the 3' samples makes them a
# comparison within one library family rather than across three. Every
# cross-sample TABLE still carries all of SAMPLE_IDS, so nothing is hidden from
# the record -- and GSE181279, the only sample where mgatk retains anything at
# all, is one of the two excluded, so the panels alone understate what mgatk
# does on deep 5' data. Read the tables alongside them.
CROSS_SAMPLE_FIG_EXCLUDE <- c(
  "GSE163668_GSM4995445_5PR2",
  "GSE181279_GSM5494116_5PPE"
)

CROSS_SAMPLE_FIG_IDS <- setdiff(SAMPLE_IDS, CROSS_SAMPLE_FIG_EXCLUDE)

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

# Bands for the call-level spectrum and the cell-level depth panel, steps 07
# and 08. Separate from AF_BIN_* on purpose: those describe carrier AF of a
# final variant set, these run two decades lower because a cell-level AF and a
# prevalence measure both reach far below 1%.
CALL_AF_BIN_BREAKS <- c(0, 0.001, 0.01, 0.05, 0.20, 1)
CALL_AF_BIN_LABELS <- c("<0.1%", "0.1-1%", "1-5%", "5-20%", "20-100%")

# Fixed log10 grid for stacking the call-level dots, so bin width is identical
# across samples and across the two AF measures and the panels stay comparable.
# Values below 1e-5 are clamped into the first bin and counted in the subtitle.
CALL_AF_LOG_MIN <- -5
CALL_AF_STACK_STEP <- 0.125

# A cell counts as carrying the variant when it clears scMOCHA's per-cell rule:
# CUTOFF_ALT_STRAND alt reads on each strand and CUTOFF_MIN_READS reads of
# depth at the position.
#
# The cached data has AF and depth but not the forward/reverse split, so the
# strand pair is approximated by its sum, 2 x CUTOFF_ALT_STRAND alt reads. That
# is permissive: 4 alt reads all on one strand would pass here and fail in the
# caller.
#
# Note that newplots/thecode/scmocha-mgatk-variant-calling.py disagrees with
# itself on the third condition. The comment reads "minimum total coverage
# >=10" but the code is
#   (fwd >= 2) & (rev >= 2) & ((fwd + rev) >= low_coverage_threshold)
# on the per-base alt matrices, which is 10 alt reads, not 10 reads of depth.
# The depth reading is used here; see M19.
CUTOFF_CELL_ALT_PAIR <- 2 * CUTOFF_ALT_STRAND

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
      "background alt rate not estimable: {total_depth} background depth ",
      "carried 0 alt reads; floored at 1/{resolvable} = {signif(rate, 3)}"
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

# Call-level AF spectrum and cell-level read support -----------------------
# Shared by step 08 (one sample) and step 07 (the samples pooled), so the two
# never drift into drawing the same quantity two ways.

fn_call_af_bin <- function(af) {
  factor(
    as.character(cut(
      af,
      breaks = CALL_AF_BIN_BREAKS,
      labels = CALL_AF_BIN_LABELS,
      include.lowest = TRUE,
      right = FALSE
    )),
    levels = CALL_AF_BIN_LABELS
  )
}

# The two per-call AF measures are kept side by side rather than one being
# chosen, because they answer different questions and the stage has twice
# reported the wrong one; see D12, D17 and M14.
CALL_AF_MEASURES <- c(
  mean_scmocha = paste(
    "Prevalence: alt reads / coverage,",
    "pooled over every cell"
  ),
  af_carrier_median = paste(
    "Heteroplasmy: median AF across",
    "carrier cells"
  )
)

# Long table of one row per (sample, variant, measure) for the call arm.
fn_call_af_long <- function(variants) {
  d <- data.table::as.data.table(variants)[arm_scmocha_call == TRUE]
  out <- data.table::melt(
    d[, .(sample, variant, mean_scmocha, af_carrier_median)],
    id.vars = c("sample", "variant"),
    variable.name = "measure",
    value.name = "af",
    variable.factor = FALSE
  )
  out <- out[!is.na(af)]
  out[, measure := factor(CALL_AF_MEASURES[measure], levels = CALL_AF_MEASURES)]
  out[, af_bin := fn_call_af_bin(af)]
  out[]
}

# Stack position for each dot on a fixed log10 grid. Returns x at the bin
# centre and y as the rank within the bin, which is what makes the panel a
# dot histogram rather than a scatter.
fn_call_af_stack <- function(d) {
  out <- data.table::copy(data.table::as.data.table(d))
  out[, af_plot := pmax(af, 10^CALL_AF_LOG_MIN)]
  out[, lg := log10(af_plot)]
  out[, bin := floor((lg - CALL_AF_LOG_MIN) / CALL_AF_STACK_STEP)]
  out[,
    x := 10^(CALL_AF_LOG_MIN + (bin + 0.5) * CALL_AF_STACK_STEP) * 100
  ]
  data.table::setorder(out, measure, bin, af)
  out[, y := seq_len(.N), by = .(measure, bin)]
  out[]
}

# Panel C. One dot per sample-variant call, stacked by AF.
fn_plot_call_af_spectrum <- function(call_long, note) {
  d <- fn_call_af_stack(call_long)
  n_clamped <- d[af < 10^CALL_AF_LOG_MIN, .N]
  n_low <- d[
    measure == levels(d$measure)[1] & af < CUTOFF_HETEROPLASMIC,
    .N
  ]
  n_calls <- d[measure == levels(d$measure)[1], .N]

  clamp_note <- if (n_clamped > 0) {
    glue::glue(
      " \u00b7 {n_clamped} call(s) below ",
      "{scales::percent(10^CALL_AF_LOG_MIN, accuracy = 0.001)} drawn in the ",
      "first bin"
    )
  } else {
    ""
  }

  ggplot2::ggplot(
    d,
    ggplot2::aes(x = x, y = y, color = af_bin)
  ) +
    ggplot2::geom_point(size = 1.6) +
    ggplot2::geom_vline(
      xintercept = CUTOFF_HETEROPLASMIC * 100,
      linetype = "dashed",
      color = "grey30"
    ) +
    ggplot2::facet_wrap(
      ~measure,
      ncol = 1,
      scales = "free_y",
      labeller = ggplot2::labeller(measure = scales::label_wrap(60))
    ) +
    ggplot2::scale_x_log10(
      labels = scales::label_number(big.mark = ",", drop0trailing = TRUE)
    ) +
    ggplot2::scale_y_continuous(
      expand = ggplot2::expansion(mult = c(0, 0.12))
    ) +
    ggplot2::scale_color_manual(values = color_call_af_bin, drop = FALSE) +
    ggplot2::guides(
      color = ggplot2::guide_legend(override.aes = list(size = 2.6))
    ) +
    fn_theme() +
    ggplot2::labs(
      title = "Allele frequency of every scMOCHA variant call",
      subtitle = glue::glue(
        "one dot per sample-variant call, not per distinct variant \u00b7 ",
        "{n_calls} calls \u00b7 dashed line at the ",
        "{scales::percent(CUTOFF_HETEROPLASMIC)} downstream cutoff \u00b7 ",
        "{n_low} call(s) fall below it on the prevalence measure \u00b7 the ",
        "two panels are different quantities and are not interchangeable",
        "{clamp_note} \u00b7 {note}"
      ),
      x = "Variant AF (%, log scale)",
      y = "Number of calls",
      color = "AF range"
    )
}

# Panel D. One dot per cell in which scMOCHA would call the variant, so the
# axis pair shows what read support a given AF needs to be callable at all.
# Cells below either scMOCHA floor are excluded upstream.
fn_plot_cell_af_depth <- function(cells, note, n_cells_total = NA_integer_) {
  d <- data.table::as.data.table(cells)[
    !is.na(af) & af > 0 & !is.na(depth) & depth > 0
  ]
  d[, af_pct := af * 100]
  d[, af_bin := fn_call_af_bin(af)]
  n_low <- d[af < CUTOFF_HETEROPLASMIC, .N]

  # A cell carries several called variants, so rows outnumber cells about four
  # to one. Reporting rows as "cells" reads as a cell count and is wrong.
  n_cells_seen <- nrow(unique(d[, .(sample, barcode)]))
  n_calls <- nrow(unique(d[, .(sample, variant)]))
  cells_note <- if (is.na(n_cells_total)) {
    glue::glue("{format(n_cells_seen, big.mark = ',')} cells")
  } else {
    glue::glue(
      "{format(n_cells_seen, big.mark = ',')} of ",
      "{format(n_cells_total, big.mark = ',')} cells"
    )
  }

  # ggplot only draws a key glyph for a level that appears in the data, so an
  # AF band with no cells would render as a bare label. One invisible point
  # per band forces every key, and override.aes makes the keys opaque.
  keys <- data.table::data.table(
    depth = 1,
    af_pct = 0.01,
    af_bin = factor(CALL_AF_BIN_LABELS, levels = CALL_AF_BIN_LABELS)
  )

  # Show at least 0.01%-100% so the requested decades are always on the axis,
  # but never crop: a deep cell can report an AF below 0.01%, and a fixed floor
  # silently dropped 126 of them in GSE181279.
  y_lo <- min(0.01, min(d$af_pct, na.rm = TRUE)) * 0.9
  x_lo <- min(1, min(d$depth, na.rm = TRUE))

  ggplot2::ggplot(
    d,
    ggplot2::aes(x = depth, y = af_pct, color = af_bin)
  ) +
    ggplot2::geom_point(data = keys, alpha = 0, show.legend = TRUE) +
    ggplot2::geom_point(size = 0.8, alpha = 0.65) +
    ggplot2::geom_hline(
      yintercept = CUTOFF_HETEROPLASMIC * 100,
      linetype = "dashed",
      color = color_cutoff_line
    ) +
    ggplot2::scale_x_log10(
      limits = c(x_lo, NA),
      breaks = 10^(0:6),
      labels = scales::label_log()
    ) +
    ggplot2::scale_y_log10(
      limits = c(y_lo, 100),
      breaks = 10^(floor(log10(y_lo)):2),
      labels = scales::label_log()
    ) +
    ggplot2::scale_color_manual(
      values = color_call_af_bin,
      limits = CALL_AF_BIN_LABELS,
      drop = FALSE,
      na.translate = FALSE
    ) +
    ggplot2::guides(
      color = ggplot2::guide_legend(
        override.aes = list(size = 2.6, alpha = 1)
      )
    ) +
    fn_theme() +
    ggplot2::labs(
      title = "Read support behind each cell-level allele frequency",
      subtitle = glue::glue(
        "one dot per cell-variant pair, not per cell \u00b7 ",
        "{format(nrow(d), big.mark = ',')} pairs from {cells_note}, over ",
        "{n_calls} sample-variant calls ",
        "({data.table::uniqueN(d$variant)} distinct variants) \u00b7 a pair ",
        "is kept when the cell clears scMOCHA's per-cell rule: \u2265 ",
        "{CUTOFF_ALT_STRAND} alt reads on each strand, approximated by ",
        "\u2265 {CUTOFF_CELL_ALT_PAIR} in total, at \u2265 ",
        "{CUTOFF_MIN_READS} reads of depth \u00b7 {n_low} pairs sit below the ",
        "{scales::percent(CUTOFF_HETEROPLASMIC)} cutoff \u00b7 a cell can only ",
        "report an AF as low as its depth allows, so the lower-left edge is a ",
        "detection limit, not a biological floor \u00b7 {note}"
      ),
      x = "Cell mtDNA depth at the variant position (reads, log scale)",
      y = "Variant AF in the cell (%, log scale)",
      color = "AF range"
    )
}
