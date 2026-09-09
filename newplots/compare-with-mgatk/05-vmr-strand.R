#!/usr/bin/env Rscript
# Metainfo ----------------------------------------------------------------
# @AUTHOR: Chun-Jie Liu
# @CONTACT: chunjie.sam.liu.at.gmail.com
# @DATE: 2026-09-09
# @DESCRIPTION: Criterion C4 in detail. Panel a, the VMR / strand-correlation
#               plane for each rule set with mgatk's cutoffs drawn; panel b,
#               AF-bin composition of the region mgatk's gate rejects versus
#               the region it accepts; panel c, per-cell alt-read support of
#               the variants only mgatk reports.
#               Usage: pixi run Rscript newplots/compare-with-mgatk/05-vmr-strand.R
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

af_bin_colors <- stats::setNames(
  c(color_af_bin, "grey80"),
  c(AF_BIN_LABELS, AF_BIN_NONE)
)

# body ---------------------------------------------------------------------
# Carrier AF from the background error model; see config.R.md for why the
# other two definitions cannot answer a question about low heteroplasmy.
d <- stats_gated[s1_mgatk == TRUE | s1_scmocha == TRUE]
d[, af_bin := fn_af_bin(af_carrier_median)]

# a: the VMR / strand plane -------------------------------------------------
# Both callers compute vmr and strand_correlation; only mgatk acts on them.
d_a <- data.table::rbindlist(list(
  d[
    s1_mgatk == TRUE,
    .(
      caller = "Original mgatk",
      vmr = vmr_mgatk,
      strand = strand_mgatk,
      af_bin
    )
  ],
  d[
    s1_scmocha == TRUE,
    .(
      caller = "scMOCHA",
      vmr = vmr_scmocha,
      strand = strand_scmocha,
      af_bin
    )
  ]
))[!is.na(vmr) & !is.na(strand) & vmr > 0]
d_a[, caller := factor(caller, levels = CALLER_LEVELS)]

panel_note <- data.table::data.table(
  caller = factor(CALLER_LEVELS, levels = CALLER_LEVELS),
  label = c("cutoffs applied", "cutoffs shown, not applied")
)

p_a <- d_a |>
  ggplot(aes(x = strand, y = vmr, color = af_bin)) +
  geom_point(size = 0.9, alpha = 0.7) +
  geom_hline(
    yintercept = CUTOFF_VMR_MGATK,
    linetype = "dashed",
    color = "grey30"
  ) +
  geom_vline(
    xintercept = CUTOFF_STRAND_MGATK,
    linetype = "dashed",
    color = "grey30"
  ) +
  geom_text(
    data = panel_note,
    aes(x = -Inf, y = Inf, label = label),
    inherit.aes = FALSE,
    hjust = -0.05,
    vjust = 1.4,
    size = 3,
    color = "grey40"
  ) +
  facet_wrap(~caller) +
  scale_y_log10(labels = scales::label_number()) +
  scale_color_manual(values = af_bin_colors) +
  guides(color = guide_legend(override.aes = list(size = 2.5))) +
  fn_theme() +
  labs(
    title = "The VMR and strand-correlation plane",
    subtitle = glue::glue(
      "S1 variants \u00b7 dashed lines are mgatk's vmr > {CUTOFF_VMR_MGATK} ",
      "and strand r > {CUTOFF_STRAND_MGATK} \u00b7 scMOCHA applies neither ",
      "\u00b7 colour is median AF across carrier cells \u00b7 ",
      "{fn_sample_note()}"
    ),
    x = "Strand correlation",
    y = "VMR (log scale)",
    color = "Carrier AF"
  )

# b: what the gate rejects, by AF bin ---------------------------------------
d_b <- d[gate_scmocha == TRUE, .N, by = .(gate_mgatk, af_bin)]
d_b[,
  gate_call := factor(
    data.table::fifelse(
      gate_mgatk,
      "Passes mgatk VMR/strand gate",
      "Rejected by mgatk VMR/strand gate"
    ),
    levels = names(color_mgatk_gate)
  )
]
d_b[, frac := N / sum(N), by = gate_call]

tab_b <- data.table::dcast(
  d_b,
  af_bin ~ gate_call,
  value.var = "N",
  fill = 0L
)

p_b <- d_b |>
  ggplot(aes(x = gate_call, y = frac, fill = af_bin)) +
  geom_col(width = 0.6) +
  scale_x_discrete(labels = scales::label_wrap(18)) +
  scale_y_continuous(labels = scales::percent, expand = c(0, 0)) +
  scale_fill_manual(values = af_bin_colors) +
  fn_theme() +
  labs(
    title = "Heteroplasmy of what mgatk's gate keeps and discards",
    subtitle = glue::glue(
      "variants passing the scMOCHA reliability gate \u00b7 ",
      "{d_b[gate_mgatk == TRUE, sum(N)]} pass and ",
      "{d_b[gate_mgatk == FALSE, sum(N)]} are rejected \u00b7 ",
      "{fn_sample_note()}"
    ),
    x = NULL,
    y = "Fraction of variants",
    fill = "Carrier AF"
  )

# c: read support of the mgatk-only variants --------------------------------
# Alt reads per cell are reconstructed as AF x depth; the callers store AF and
# total coverage, not the alt count itself.
detection <- fn_detection_long(
  cell_af,
  cell_pos_coverage,
  stats_gated[, .(variant, position)]
)
detection[, alt_reads := af * depth]

set_class <- d[, .(
  variant,
  variant_set = data.table::fcase(
    final_mgatk & final_scmocha  , "Both"                ,
    final_mgatk & !final_scmocha , "Original mgatk only" ,
    !final_mgatk & final_scmocha , "scMOCHA only"        ,
    default = NA_character_
  )
)][!is.na(variant_set)]

d_c <- detection[
  detected == TRUE
][
  set_class,
  on = "variant",
  variant_set := i.variant_set
][!is.na(variant_set)]
d_c[, variant_set := factor(variant_set, levels = VARIANT_SET_LEVELS)]

med_c <- d_c[,
  .(
    median_alt_reads = stats::median(alt_reads),
    n_detections = .N
  ),
  by = variant_set
]

p_c <- d_c |>
  ggplot(aes(x = variant_set, y = alt_reads, fill = variant_set)) +
  geom_violin(color = NA, alpha = 0.6, scale = "width") +
  geom_boxplot(width = 0.15, outlier.size = 0.3, fill = "white") +
  geom_hline(
    yintercept = CUTOFF_ALT_READS,
    linetype = "dashed",
    color = "grey30"
  ) +
  scale_x_discrete(labels = scales::label_wrap(12)) +
  scale_y_log10(labels = scales::label_number()) +
  scale_fill_manual(values = color_variant_set) +
  fn_theme() +
  theme(legend.position = "none") +
  labs(
    title = "Read support behind each variant class",
    subtitle = glue::glue(
      "alt reads per detected cell, reconstructed as AF \u00d7 depth \u00b7 ",
      "dashed line at the {CUTOFF_ALT_READS}-read scMOCHA requirement \u00b7 ",
      "{fn_sample_note()}"
    ),
    x = NULL,
    y = "Alt reads in the cell (log scale)"
  )

log_info(
  "median alt reads per detected cell: ",
  "{paste(med_c$variant_set, signif(med_c$median_alt_reads, 3), sep = '=', collapse = ', ')}"
)

# d: the same plane, coloured by which arm reports the variant ---------------
# mgatk's own vmr and strand correlation are what its gate acts on, so the
# question "where do the variants only scMOCHA reports sit in mgatk's plane"
# is asked on mgatk's coordinates.
d_d <- d[!is.na(vmr_mgatk) & !is.na(strand_mgatk) & vmr_mgatk > 0]
d_d[,
  arm_label := data.table::fcase(
    arm_mgatk & arm_scmocha_af5  , "Both arms"                 ,
    arm_mgatk & !arm_scmocha_af5 , "Original mgatk only"       ,
    !arm_mgatk & arm_scmocha_af5 , "scMOCHA AF>5% only"        ,
    arm_scmocha_call             , "scMOCHA variant call only" ,
    default = "Neither arm"
  )
]
d_d[,
  arm_label := factor(
    arm_label,
    levels = c(
      "Both arms",
      "Original mgatk only",
      "scMOCHA AF>5% only",
      "scMOCHA variant call only",
      "Neither arm"
    )
  )
]

arm_label_colors <- c(
  "Both arms" = unname(color_variant_set["Both"]),
  "Original mgatk only" = unname(color_arm[ARM_MGATK]),
  "scMOCHA AF>5% only" = unname(color_arm[ARM_SCMOCHA_AF5]),
  "scMOCHA variant call only" = unname(color_arm[ARM_SCMOCHA]),
  "Neither arm" = "grey85"
)

p_d <- d_d[order(-as.integer(arm_label))] |>
  ggplot(aes(x = strand_mgatk, y = vmr_mgatk, color = arm_label)) +
  geom_point(size = 1.1, alpha = 0.8) +
  geom_hline(
    yintercept = CUTOFF_VMR_MGATK,
    linetype = "dashed",
    color = "grey30"
  ) +
  geom_vline(
    xintercept = CUTOFF_STRAND_MGATK,
    linetype = "dashed",
    color = "grey30"
  ) +
  scale_y_log10(labels = scales::label_number()) +
  scale_color_manual(values = arm_label_colors) +
  guides(color = guide_legend(override.aes = list(size = 2.5))) +
  fn_theme() +
  labs(
    title = "Which arm reports each variant, in mgatk's decision plane",
    subtitle = glue::glue(
      "mgatk's own vmr and strand correlation \u00b7 dashed lines are its ",
      "gate \u00b7 everything mgatk reports lies in the upper right by ",
      "construction, so the informative points are the scMOCHA AF>5% ",
      "variants outside it \u00b7 {fn_sample_note()}"
    ),
    x = "Strand correlation (original mgatk)",
    y = "VMR (original mgatk, log scale)",
    color = NULL
  )

# e: the same plane, faceted by arm membership ------------------------------
p_e <- d_d |>
  ggplot(aes(x = strand_mgatk, y = vmr_mgatk)) +
  geom_point(
    data = d_d[, .(strand_mgatk, vmr_mgatk)],
    color = "grey88",
    size = 0.7
  ) +
  geom_point(aes(color = arm_label), size = 1.1, alpha = 0.85) +
  geom_hline(
    yintercept = CUTOFF_VMR_MGATK,
    linetype = "dashed",
    color = "grey30"
  ) +
  geom_vline(
    xintercept = CUTOFF_STRAND_MGATK,
    linetype = "dashed",
    color = "grey30"
  ) +
  facet_wrap(~arm_label, nrow = 1) +
  scale_y_log10(labels = scales::label_number()) +
  scale_color_manual(values = arm_label_colors) +
  fn_theme() +
  theme(legend.position = "none") +
  labs(
    title = "Arm membership across mgatk's decision plane",
    subtitle = glue::glue(
      "grey points are all S1 variants, repeated in every panel \u00b7 dashed ",
      "lines are mgatk's vmr > {CUTOFF_VMR_MGATK} and strand r > ",
      "{CUTOFF_STRAND_MGATK} \u00b7 {fn_sample_note()}"
    ),
    x = "Strand correlation (original mgatk)",
    y = "VMR (original mgatk, log scale)"
  )

log_info(
  "arm membership in mgatk plane: ",
  "{paste(names(table(d_d$arm_label)), table(d_d$arm_label), sep = '=', collapse = ', ')}"
)

# save ---------------------------------------------------------------------
export(
  tab_b,
  as.character(fs::path(
    paths$tabdir,
    "05-vmr-strand-rejected.tsv"
  ))
)
export(med_c, as.character(fs::path(paths$tabdir, "05-read-support.tsv")))
export(
  d_d[, .N, by = arm_label][order(-N)][, .(arm_label, n_variants = N)],
  as.character(fs::path(paths$tabdir, "05-arm-in-mgatk-plane.tsv"))
)

saveplot(
  as.character(fs::path(paths$figdir, "05a-vmr-strand-plane.pdf")),
  fn_wrap_labs(p_a, width = 110),
  width = 10,
  height = 5,
  device = cairo_pdf
)
saveplot(
  as.character(fs::path(paths$figdir, "05b-gate-rejected-af-bins.pdf")),
  fn_wrap_labs(p_b, width = 58),
  width = 6.5,
  height = 5.5,
  device = cairo_pdf
)
saveplot(
  as.character(fs::path(paths$figdir, "05c-read-support.pdf")),
  fn_wrap_labs(p_c, width = 60),
  width = 7,
  height = 5.5,
  device = cairo_pdf
)
saveplot(
  as.character(fs::path(paths$figdir, "05d-arm-in-mgatk-plane.pdf")),
  fn_wrap_labs(p_d, width = 85),
  width = 8.5,
  height = 5.5,
  device = cairo_pdf
)
saveplot(
  as.character(fs::path(paths$figdir, "05e-arm-plane-facets.pdf")),
  fn_wrap_labs(p_e, width = 115),
  width = 13,
  height = 4.2,
  device = cairo_pdf
)

log_info("saved 5 figures to {paths$figdir} and 3 tables to {paths$tabdir}")
