#!/usr/bin/env Rscript
# Metainfo ----------------------------------------------------------------
# @AUTHOR: Chun-Jie Liu
# @CONTACT: chunjie.sam.liu.at.gmail.com
# @DATE: 2026-09-23
# @DESCRIPTION: Allele-frequency spectrum of the scMOCHA variant calls and the
#               read support behind them, for one sample. Panel a, one dot per
#               sample-variant call stacked by AF, drawn on both the
#               prevalence and the carrier-heteroplasmy measure; panel b, one
#               dot per cell holding an alt read, AF against that cell's depth.
#               Usage: pixi run Rscript newplots/compare-with-mgatk/08-call-af-depth.R --sample=<id>
# @VERSION: v0.1.0

# Reproducibility ----------------------------------------------------------
set.seed(9527)

# Library ------------------------------------------------------------------
suppressMessages({
  library(jutils)
  load_pkg(scales)
})

# Logger -------------------------------------------------------------------
log_layout(layout_glue_colors)
log_threshold(INFO)

# Load data ----------------------------------------------------------------
dotenv()
suppressMessages({
  conflicted::conflicts_prefer(dplyr::filter, fs::path)
})

# Source -------------------------------------------------------------------
stagedir <- fs::path(
  fs::path_expand(Sys.getenv("REPODIR")),
  "newplots",
  "compare-with-mgatk"
)
source(fs::path(stagedir, "config.R"))

# args --------------------------------------------------------------------
GetoptLong.options(help_style = "two-column")
sample <- ""

GetoptLong(
  "sample=s",
  "sample id; one of the ids in SAMPLES in config.R"
)

sample_id <- fn_check_sample(sample)
rm(sample)
SAMPLE_LABEL <- fn_sample_label(sample_id)
log_info("sample {sample_id} ({SAMPLE_LABEL})")

paths <- stage_paths(sample_id)
source(paths$colorfile)
fs::dir_create(c(paths$figdir, paths$tabdir))

stats_gated <- import(as.character(fs::path(
  paths$cachedir,
  "03-variant-gated.qs"
))) |>
  data.table::as.data.table()
cell_af <- import(as.character(fs::path(paths$cachedir, "01-cell-af-s1.qs")))
cell_pos_coverage <- import(as.character(fs::path(
  paths$cachedir,
  "01-cell-pos-coverage.qs"
)))

# body ---------------------------------------------------------------------
stats_gated[, sample := sample_id]

call_long <- fn_call_af_long(stats_gated)
n_calls <- data.table::uniqueN(call_long$variant)

log_info(
  "call arm: {n_calls} variants; measures ",
  "{paste(levels(call_long$measure), collapse = ' | ')}"
)

# Cell-level support. scMOCHA's per-cell rule is CUTOFF_ALT_STRAND alt reads on
# each strand at CUTOFF_MIN_READS reads of depth; the cache has no strand split,
# so the pair is approximated by its sum.
detection <- fn_detection_long(
  cell_af,
  cell_pos_coverage,
  stats_gated[, .(variant, position)]
)
called <- stats_gated[arm_scmocha_call == TRUE | arm_mgatk == TRUE, variant]

cells <- detection[
  variant %in% called & !is.na(af) & af > 0 & depth >= CUTOFF_MIN_READS
][
  round(af * depth) >= CUTOFF_CELL_ALT_PAIR
][,
  .(sample = sample_id, variant, barcode, af, depth)
]

n_cells_total <- data.table::uniqueN(detection$barcode)

log_info(
  "cells clearing scMOCHA's per-cell rule (>= {CUTOFF_CELL_ALT_PAIR} alt ",
  "reads, >= {CUTOFF_MIN_READS} depth): {nrow(cells)} cell-variant pairs ",
  "from {nrow(unique(cells[, .(sample, barcode)]))} of {n_cells_total} cells, ",
  "over {nrow(unique(cells[, .(sample, variant)]))} calls; ",
  "{cells[af < CUTOFF_HETEROPLASMIC, .N]} below the ",
  "{CUTOFF_HETEROPLASMIC} cutoff"
)

p_a <- fn_or_empty(
  nrow(call_long) > 0L,
  fn_plot_call_af_spectrum(call_long, fn_sample_note()),
  "Allele frequency of every scMOCHA variant call",
  "The scMOCHA call arm is empty in this sample."
)

p_b <- fn_or_empty(
  nrow(cells) > 0L,
  fn_plot_cell_af_depth(cells, fn_sample_note(), n_cells_total),
  "Read support behind each cell-level allele frequency",
  glue::glue(
    "No cell reaches {CUTOFF_CELL_ALT_PAIR} alt reads at {CUTOFF_MIN_READS} ",
    "reads of depth for a called variant in this sample."
  )
)

# save ---------------------------------------------------------------------
fn_export_tab(
  call_long[, .(sample, variant, measure = as.character(measure), af, af_bin)],
  paths,
  "08-call-af.tsv",
  sample_id
)
fn_export_tab(cells, paths, "08-cell-af-depth.tsv", sample_id)

saveplot(
  as.character(fs::path(paths$figdir, "08a-call-af-spectrum.pdf")),
  fn_wrap_labs(p_a, width = 95),
  width = 8.5,
  height = 7,
  device = cairo_pdf
)
saveplot(
  as.character(fs::path(paths$figdir, "08b-cell-af-depth.pdf")),
  fn_wrap_labs(p_b, width = 95),
  width = 8.5,
  height = 5.5,
  device = cairo_pdf
)

log_info("saved 2 figures to {paths$figdir} and 2 tables to {paths$tabdir}")
