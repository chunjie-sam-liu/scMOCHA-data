#!/usr/bin/env Rscript
# Metainfo ----------------------------------------------------------------
# @AUTHOR: Chun-Jie Liu
# @CONTACT: chunjie.sam.liu.at.gmail.com
# @DATE: 2026-09-09
# @DESCRIPTION: Criterion C1 and C2. Panel a, per-cell MT coverage with the
#               mgatk cell cutoff; panel b, cells retained by each caller;
#               panel c, minimum detectable AF per cell under each confident
#               detection rule; panel d, variant detections lost with the
#               cells mgatk discards.
#               Usage: pixi run Rscript newplots/compare-with-mgatk/02-cell-inclusion.R
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

paths <- stage_paths()
source(paths$colorfile)
fs::dir_create(c(paths$figdir, paths$tabdir))

stats_joined <- import(as.character(fs::path(
  paths$cachedir,
  "01-variant-joined.qs"
))) |>
  data.table::as.data.table()
cell_depth <- import(as.character(fs::path(
  paths$cachedir,
  "01-cell-depth.qs"
))) |>
  data.table::as.data.table()
cell_af <- import(as.character(fs::path(paths$cachedir, "01-cell-af-s1.qs")))
cell_pos_coverage <- import(as.character(fs::path(
  paths$cachedir,
  "01-cell-pos-coverage.qs"
)))

# body ---------------------------------------------------------------------
cell_depth[,
  inclusion := data.table::fifelse(
    kept_by_mgatk,
    "Kept by both",
    "Dropped by original mgatk"
  )
]
cell_depth[,
  inclusion := factor(inclusion, levels = names(color_cell_inclusion))
]

n_total <- nrow(cell_depth)
n_dropped <- cell_depth[kept_by_mgatk == FALSE, .N]

# a: per-cell coverage ------------------------------------------------------
p_a <- cell_depth |>
  ggplot(aes(x = mean_coverage, fill = inclusion)) +
  geom_histogram(bins = 60, color = NA) +
  geom_vline(
    xintercept = CUTOFF_CELL_MEANCOV,
    linetype = "dashed",
    color = "grey30"
  ) +
  scale_x_log10(labels = scales::label_number()) +
  scale_fill_manual(values = color_cell_inclusion) +
  fn_theme() +
  labs(
    title = "Per-cell mitochondrial coverage and the mgatk cell filter",
    subtitle = glue::glue(
      "{n_total} cells \u00b7 {n_dropped} ",
      "({scales::percent(n_dropped / n_total, accuracy = 0.1)}) dropped by ",
      "mean coverage \u2264 {CUTOFF_CELL_MEANCOV} \u00b7 {fn_sample_note()}"
    ),
    x = "Mean MT coverage per cell (log scale)",
    y = "Cells",
    fill = NULL
  )

# b: cells retained ---------------------------------------------------------
d_b <- data.table::data.table(
  caller = factor(CALLER_LEVELS, levels = CALLER_LEVELS),
  n_cells = c(n_total - n_dropped, n_total)
)

p_b <- d_b |>
  ggplot(aes(x = caller, y = n_cells, fill = caller)) +
  geom_col(width = 0.6) +
  geom_text(aes(label = scales::comma(n_cells)), vjust = -0.4, size = 3.5) +
  scale_y_continuous(
    labels = scales::comma,
    expand = expansion(mult = c(0, 0.12))
  ) +
  scale_fill_manual(values = color_caller) +
  fn_theme() +
  theme(legend.position = "none") +
  labs(
    title = "Cells entering variant calling",
    subtitle = glue::glue(
      "original mgatk requires mean coverage > {CUTOFF_CELL_MEANCOV}; ",
      "scMOCHA applies no cell filter \u00b7 {fn_sample_note()}"
    ),
    x = NULL,
    y = "Cells"
  )

# c: minimum detectable AF --------------------------------------------------
# mgatk needs 2 alt reads on each strand; scMOCHA needs that plus 10 alt reads
# in total, so the floors are 4/depth and 10/depth.
cov_grid <- 10^seq(
  log10(max(1, min(cell_depth$mean_coverage))),
  log10(max(cell_depth$mean_coverage)),
  length.out = 400
)

d_c <- data.table::rbindlist(list(
  data.table::data.table(
    caller = "Original mgatk",
    coverage = cov_grid,
    min_af = pmin(1, (2 * CUTOFF_ALT_STRAND) / cov_grid)
  ),
  data.table::data.table(
    caller = "scMOCHA",
    coverage = cov_grid,
    min_af = pmin(1, CUTOFF_ALT_READS / cov_grid)
  )
))
d_c[, caller := factor(caller, levels = CALLER_LEVELS)]

p_c <- ggplot() +
  geom_rug(
    data = cell_depth[sample(.N, min(.N, 2000))],
    aes(x = mean_coverage),
    sides = "b",
    alpha = 0.15,
    color = "grey50"
  ) +
  geom_line(
    data = d_c,
    aes(x = coverage, y = min_af, color = caller),
    linewidth = 0.9
  ) +
  geom_hline(
    yintercept = CUTOFF_HETEROPLASMIC,
    linetype = "dotted",
    color = "grey30"
  ) +
  scale_x_log10(labels = scales::label_number()) +
  scale_y_log10(labels = scales::label_number()) +
  scale_color_manual(values = color_caller) +
  fn_theme() +
  labs(
    title = "Lowest heteroplasmy a single cell can support",
    subtitle = glue::glue(
      "alt reads required per cell: mgatk {2 * CUTOFF_ALT_STRAND}, ",
      "scMOCHA {CUTOFF_ALT_READS} \u00b7 dotted line at the ",
      "{CUTOFF_HETEROPLASMIC} heteroplasmy threshold \u00b7 rug shows observed ",
      "per-cell coverage"
    ),
    x = "Mean MT coverage per cell (log scale)",
    y = "Minimum detectable AF (log scale)",
    color = NULL
  )

# d: detections lost with the dropped cells ---------------------------------
detection <- fn_detection_long(
  cell_af,
  cell_pos_coverage,
  stats_joined[, .(variant, position)]
)
detection[
  cell_depth[, .(barcode, kept_by_mgatk)],
  on = "barcode",
  kept_by_mgatk := i.kept_by_mgatk
]

n_det_total <- detection[detected == TRUE, .N]
n_det_lost <- detection[detected == TRUE & kept_by_mgatk == FALSE, .N]

per_variant <- detection[
  detected == TRUE,
  .(
    n_cells_all = .N,
    n_cells_mgatk = sum(kept_by_mgatk)
  ),
  by = variant
]
per_variant[, `:=`(
  pass_all = n_cells_all >= CUTOFF_NOTRELIABLE,
  pass_mgatk_cells = n_cells_mgatk >= CUTOFF_NOTRELIABLE
)]
n_var_lost <- per_variant[pass_all == TRUE & pass_mgatk_cells == FALSE, .N]

d_d <- data.table::data.table(
  metric = factor(
    c("Cell-level detections", "Variants reaching 10 cells"),
    levels = c("Cell-level detections", "Variants reaching 10 cells")
  ),
  inclusion = factor(
    rep(names(color_cell_inclusion), each = 2),
    levels = names(color_cell_inclusion)
  ),
  n = c(
    n_det_total - n_det_lost,
    per_variant[pass_mgatk_cells == TRUE, .N],
    n_det_lost,
    n_var_lost
  )
)

p_d <- d_d |>
  ggplot(aes(x = metric, y = n, fill = inclusion)) +
  geom_col(width = 0.6) +
  facet_wrap(~metric, scales = "free", nrow = 1) +
  scale_y_continuous(
    labels = scales::comma,
    expand = expansion(mult = c(0, 0.12))
  ) +
  scale_fill_manual(values = color_cell_inclusion) +
  fn_theme() +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank()
  ) +
  labs(
    title = "What the mgatk cell filter removes",
    subtitle = glue::glue(
      "detection = cell AF \u2265 {CUTOFF_HETEROPLASMIC} and cell depth ",
      "\u2265 {CUTOFF_MIN_READS} \u00b7 {n_det_lost} detections and ",
      "{n_var_lost} variants lost with the {n_dropped} discarded cells"
    ),
    x = NULL,
    y = "Count",
    fill = NULL
  )

# save ---------------------------------------------------------------------
d_out <- data.table::data.table(
  metric = c(
    "cells_total",
    "cells_kept_mgatk",
    "cells_dropped_mgatk",
    "detections_total",
    "detections_in_dropped_cells",
    "variants_ge10cells_all_cells",
    "variants_ge10cells_mgatk_cells",
    "variants_lost_by_cell_filter"
  ),
  value = c(
    n_total,
    n_total - n_dropped,
    n_dropped,
    n_det_total,
    n_det_lost,
    per_variant[pass_all == TRUE, .N],
    per_variant[pass_mgatk_cells == TRUE, .N],
    n_var_lost
  )
)
export(d_out, as.character(fs::path(paths$tabdir, "02-cell-inclusion.tsv")))

saveplot(
  as.character(fs::path(paths$figdir, "02a-cell-coverage.pdf")),
  fn_wrap_labs(p_a),
  width = 8,
  height = 5,
  device = cairo_pdf
)
saveplot(
  as.character(fs::path(paths$figdir, "02b-cells-retained.pdf")),
  fn_wrap_labs(p_b, width = 45),
  width = 5.5,
  height = 5,
  device = cairo_pdf
)
saveplot(
  as.character(fs::path(paths$figdir, "02c-min-detectable-af.pdf")),
  fn_wrap_labs(p_c),
  width = 8,
  height = 5,
  device = cairo_pdf
)
saveplot(
  as.character(fs::path(paths$figdir, "02d-detections-lost.pdf")),
  fn_wrap_labs(p_d),
  width = 8,
  height = 5,
  device = cairo_pdf
)

log_info(
  "02: {n_total} cells, {n_dropped} dropped by mgatk; {n_det_lost} of ",
  "{n_det_total} detections and {n_var_lost} variants lost with them"
)
log_info("saved 4 figures to {paths$figdir} and 02-cell-inclusion.tsv")
