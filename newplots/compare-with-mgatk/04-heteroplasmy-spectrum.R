#!/usr/bin/env Rscript
# Metainfo ----------------------------------------------------------------
# @AUTHOR: Chun-Jie Liu
# @CONTACT: chunjie.sam.liu.at.gmail.com
# @DATE: 2026-09-09
# @DESCRIPTION: The editor's question. Panel a, ECDF of population AF by
#               variant class; panel b, density of maximum heteroplasmy;
#               panel c, AF-bin composition of each final set; panel d,
#               cell-level AF distribution; panel e, AF of the variants
#               mgatk's VMR/strand gate rejects but the scMOCHA gate keeps.
#               Usage: pixi run Rscript newplots/compare-with-mgatk/04-heteroplasmy-spectrum.R
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
# Heteroplasmy is measured on carrier cells defined by the background error
# model, not by the mgatk `mean` column and not by a fixed AF or read floor.
# See config.R.md: `mean` tracks prevalence, an AF floor censors the region
# under study, and a 2-alt-read floor lets background cells dominate the
# median. Only the error-model definition is both uncensored and noise-
# resistant. Panel f shows what the other choices would have done.
d <- stats_gated[s1_mgatk == TRUE | s1_scmocha == TRUE]
d[,
  variant_set := data.table::fcase(
    final_mgatk & final_scmocha  , "Both"                ,
    final_mgatk & !final_scmocha , "Original mgatk only" ,
    !final_mgatk & final_scmocha , "scMOCHA only"        ,
    default = NA_character_
  )
]
d[, variant_set := factor(variant_set, levels = VARIANT_SET_LEVELS)]
d[, af_bin := fn_af_bin(af_carrier_median)]

d_final <- d[!is.na(variant_set)]
d_af <- d_final[n_cells_carrier > 0]

log_info(
  "variant classes: {paste(names(table(d_final$variant_set)), ",
  "table(d_final$variant_set), sep = '=', collapse = ', ')}"
)
log_info(
  "variants with no carrier cell: ",
  "{d_final[n_cells_carrier == 0, .N]} of {nrow(d_final)}"
)

# a: ECDF of carrier heteroplasmy -------------------------------------------
p_a <- d_af |>
  ggplot(aes(x = af_carrier_median, color = variant_set)) +
  stat_ecdf(linewidth = 0.9) +
  geom_vline(
    xintercept = CUTOFF_HETEROPLASMIC,
    linetype = "dotted",
    color = "grey30"
  ) +
  scale_x_log10(labels = scales::label_number()) +
  scale_color_manual(values = color_variant_set) +
  fn_theme() +
  labs(
    title = "Heteroplasmy of each variant class in the cells that carry it",
    subtitle = glue::glue(
      "{nrow(d_af)} of {nrow(d_final)} variants have a carrier cell \u00b7 ",
      "carrier = alt count inconsistent with the background error rate, ",
      "no AF floor \u00b7 dotted line at {CUTOFF_HETEROPLASMIC} \u00b7 ",
      "{fn_sample_note()}"
    ),
    x = "Median allele frequency across carrier cells (log scale)",
    y = "Cumulative fraction of variants",
    color = NULL
  )

# b: maximum heteroplasmy ---------------------------------------------------
p_b <- d_af |>
  ggplot(aes(x = af_carrier_max, fill = variant_set)) +
  geom_density(alpha = 0.4, color = NA) +
  scale_fill_manual(values = color_variant_set) +
  fn_theme() +
  labs(
    title = "Maximum per-cell heteroplasmy of each variant class",
    subtitle = glue::glue(
      "{nrow(d_af)} variants with at least one carrier cell \u00b7 ",
      "{fn_sample_note()}"
    ),
    x = "Maximum allele frequency across carrier cells",
    y = "Density",
    fill = NULL
  )

# c: AF-bin composition of each arm -----------------------------------------
d_c <- data.table::rbindlist(list(
  d[arm_mgatk == TRUE, .(arm = ARM_MGATK, af_bin)],
  d[arm_scmocha_call == TRUE, .(arm = ARM_SCMOCHA, af_bin)],
  d[arm_scmocha_af5 == TRUE, .(arm = ARM_SCMOCHA_AF5, af_bin)]
))[, .N, by = .(arm, af_bin)]
d_c[, arm := factor(arm, levels = ARM_LEVELS)]
d_c[, frac := N / sum(N), by = arm]

p_c <- d_c |>
  ggplot(aes(x = arm, y = frac, fill = af_bin)) +
  geom_col(width = 0.6) +
  scale_x_discrete(labels = scales::label_wrap(14)) +
  scale_y_continuous(labels = scales::percent, expand = c(0, 0)) +
  scale_fill_manual(values = af_bin_colors) +
  fn_theme() +
  labs(
    title = "Heteroplasmy composition of each arm",
    subtitle = glue::glue(
      "{d_c[arm == ARM_MGATK, sum(N)]} mgatk, ",
      "{d_c[arm == ARM_SCMOCHA, sum(N)]} scMOCHA call and ",
      "{d_c[arm == ARM_SCMOCHA_AF5, sum(N)]} scMOCHA AF>5% variants \u00b7 ",
      "median AF over every carrier cell, carrier set by the background error ",
      "model, one measure for all three arms \u00b7 {fn_sample_note()}"
    ),
    x = NULL,
    y = "Fraction of variants",
    fill = "Carrier AF"
  )

# g: each arm measured on its own terms -------------------------------------
# The AF>5% arm is defined by cells at AF >= 0.05, so it is summarised on those
# cells. The other two arms apply no AF rule, so they keep the error-model
# carrier measure. The three bars therefore do NOT share a measure: the AF>5%
# bar cannot fall below 0.05 while the other two can, so part of its apparent
# shift is definitional. 04c is the like-for-like comparison.
d[, af_bin_gate := fn_af_bin(af_gate_median)]

d_g <- data.table::rbindlist(list(
  d[
    arm_mgatk == TRUE,
    .(arm = ARM_MGATK, measured = "all carrier cells", af_bin_own = af_bin)
  ],
  d[
    arm_scmocha_call == TRUE,
    .(arm = ARM_SCMOCHA, measured = "all carrier cells", af_bin_own = af_bin)
  ],
  d[
    arm_scmocha_af5 == TRUE,
    .(
      arm = ARM_SCMOCHA_AF5,
      measured = "cells at AF >= 0.05 only",
      af_bin_own = af_bin_gate
    )
  ]
))
d_g[, arm := factor(arm, levels = ARM_LEVELS)]
d_g <- d_g[, .N, by = .(arm, measured, af_bin_own)]
d_g[, frac := N / sum(N), by = arm]

p_g <- d_g |>
  ggplot(aes(x = arm, y = frac, fill = af_bin_own)) +
  geom_col(width = 0.6) +
  geom_text(
    data = unique(d_g[, .(arm, measured)]),
    aes(x = arm, y = 1.03, label = scales::label_wrap(18)(measured)),
    inherit.aes = FALSE,
    size = 2.6,
    color = "grey35",
    vjust = 0
  ) +
  scale_x_discrete(labels = scales::label_wrap(14)) +
  scale_y_continuous(
    labels = scales::percent,
    expand = expansion(mult = c(0, 0.16))
  ) +
  scale_fill_manual(values = af_bin_colors) +
  fn_theme() +
  labs(
    title = "Each arm on its own AF definition",
    subtitle = glue::glue(
      "the {CUTOFF_HETEROPLASMIC} AF rule belongs to the scMOCHA AF>5% arm ",
      "alone, so only that arm is measured on cells at AF \u2265 ",
      "{CUTOFF_HETEROPLASMIC} \u00b7 the three bars do not share a measure ",
      "and are not a like-for-like comparison; see panel 04c for that \u00b7 ",
      "{fn_sample_note()}"
    ),
    x = NULL,
    y = "Fraction of variants",
    fill = "AF over the arm's\nown cells"
  )

log_info(
  "own-definition view: {paste(unique(d_g[, .(arm, measured)])$arm, ",
  "unique(d_g[, .(arm, measured)])$measured, sep = ' -> ', collapse = '; ')}"
)

# d: cell-level AF ----------------------------------------------------------
detection <- fn_detection_long(
  cell_af,
  cell_pos_coverage,
  stats_gated[, .(variant, position)]
)
background_rate <- fn_background_rate(detection)
detection[, alt_reads := round(af * depth)]

# Same carrier rule as the per-variant summaries: no AF floor, so the sub-0.05
# tail stays visible, but background cells are excluded by the error model.
alpha <- 0.05 / nrow(detection)
obs <- detection[
  depth >= CUTOFF_MIN_READS &
    !is.na(alt_reads) &
    alt_reads >= 1 &
    stats::pbinom(alt_reads - 1, depth, background_rate, lower.tail = FALSE) <
      alpha
]

d_d <- data.table::rbindlist(list(
  obs[
    variant %in% d[final_mgatk == TRUE, variant],
    .(caller = "Original mgatk", af)
  ],
  obs[
    variant %in% d[final_scmocha == TRUE, variant],
    .(caller = "scMOCHA", af)
  ]
))
d_d[, caller := factor(caller, levels = CALLER_LEVELS)]

# Density rather than counts: the two sets carry very different numbers of
# observations, and empty count bins are -Inf on a log axis.
p_d <- d_d |>
  ggplot(aes(x = af, color = caller)) +
  geom_density(linewidth = 0.8) +
  geom_vline(
    xintercept = CUTOFF_HETEROPLASMIC,
    linetype = "dotted",
    color = "grey30"
  ) +
  scale_x_log10(labels = scales::label_number()) +
  scale_color_manual(values = color_caller) +
  fn_theme() +
  labs(
    title = "Cell-level heteroplasmy of the detections each set carries",
    subtitle = glue::glue(
      "carrier cells only, alt count inconsistent with the ",
      "{signif(background_rate, 3)} background rate \u00b7 ",
      "{d_d[caller == 'Original mgatk', .N]} mgatk and ",
      "{d_d[caller == 'scMOCHA', .N]} scMOCHA cell observations \u00b7 ",
      "dotted line at {CUTOFF_HETEROPLASMIC} \u00b7 {fn_sample_note()}"
    ),
    x = "Cell-level allele frequency (log scale)",
    y = "Density",
    color = NULL
  )

# e: what the mgatk gate rejects --------------------------------------------
# Restricted to variants that reach the scMOCHA reliability gate, so the
# comparison is between two decisions about the same reliable variants. That
# gate already requires CUTOFF_NOTRELIABLE carrier cells, so carrier AF is
# defined for every variant here.
d_e <- d[gate_scmocha == TRUE]
d_e[,
  gate_call := data.table::fifelse(
    gate_mgatk,
    "Passes mgatk VMR/strand gate",
    "Rejected by mgatk VMR/strand gate"
  )
]
d_e[, gate_call := factor(gate_call, levels = names(color_mgatk_gate))]

wt <- stats::wilcox.test(
  af_carrier_median ~ gate_call,
  data = d_e,
  conf.int = TRUE
)
med <- d_e[,
  .(
    median_af = stats::median(af_carrier_median),
    n = .N
  ),
  keyby = gate_call
]

af_pass <- med[gate_call == "Passes mgatk VMR/strand gate", median_af]
af_rej <- med[gate_call == "Rejected by mgatk VMR/strand gate", median_af]
n_pass <- med[gate_call == "Passes mgatk VMR/strand gate", n]
n_rej <- med[gate_call == "Rejected by mgatk VMR/strand gate", n]

log_info(
  "gate test on error-model carrier AF: rejected n={n_rej} ",
  "median={signif(af_rej, 3)}; passed n={n_pass} median={signif(af_pass, 3)}; ",
  "P={signif(wt$p.value, 3)}"
)

p_e <- d_e |>
  ggplot(aes(x = gate_call, y = af_carrier_median, fill = gate_call)) +
  geom_violin(color = NA, alpha = 0.6, scale = "width") +
  geom_boxplot(width = 0.15, outlier.size = 0.4, fill = "white") +
  scale_x_discrete(labels = scales::label_wrap(20)) +
  scale_y_log10(labels = scales::label_number()) +
  scale_fill_manual(values = color_mgatk_gate) +
  fn_theme() +
  theme(legend.position = "none") +
  labs(
    title = "Heteroplasmy of the variants mgatk's gate discards",
    subtitle = glue::glue(
      "variants passing the scMOCHA reliability gate, split by mgatk's ",
      "vmr > {CUTOFF_VMR_MGATK} and strand r > {CUTOFF_STRAND_MGATK} \u00b7 ",
      "median carrier AF {signif(af_pass, 3)} (n = {n_pass}) vs ",
      "{signif(af_rej, 3)} (n = {n_rej}) \u00b7 Wilcoxon P = ",
      "{format.pval(wt$p.value, digits = 3)}"
    ),
    x = NULL,
    y = "Median allele frequency across carrier cells (log scale)"
  )

# f: does the conclusion survive a different carrier definition? ------------
# The choice of carrier definition changed the answer once already, so all
# four are reported side by side rather than only the one that was adopted.
measures <- list(
  list(
    label = "AF >= 0.05 cells (censored at 0.05)",
    col = "af_gate_median"
  ),
  list(
    label = ">= 2 alt reads (background-dominated)",
    col = "af_loose_median"
  ),
  list(
    label = "error-model carriers (adopted)",
    col = "af_carrier_median"
  ),
  list(
    label = "max across error-model carriers",
    col = "af_carrier_max"
  )
)

d_f <- data.table::rbindlist(lapply(measures, function(m) {
  x <- d_e[!is.na(get(m$col))]
  w <- stats::wilcox.test(x[[m$col]] ~ x$gate_call)
  data.table::data.table(
    measure = m$label,
    n = nrow(x),
    median_rejected = stats::median(
      x[gate_call == "Rejected by mgatk VMR/strand gate", get(m$col)]
    ),
    median_passed = stats::median(
      x[gate_call == "Passes mgatk VMR/strand gate", get(m$col)]
    ),
    p_value = w$p.value
  )
}))
d_f[,
  measure := factor(
    measure,
    levels = rev(vapply(
      measures,
      \(m) m$label,
      character(1)
    ))
  )
]

p_f <- data.table::melt(
  d_f,
  id.vars = c("measure", "p_value"),
  measure.vars = c("median_rejected", "median_passed"),
  variable.name = "gate_call",
  value.name = "median_af"
) |>
  ggplot(aes(x = median_af, y = measure, color = gate_call)) +
  geom_line(aes(group = measure), color = "grey70", linewidth = 0.6) +
  geom_point(size = 3) +
  geom_text(
    data = d_f,
    aes(
      x = Inf,
      y = measure,
      label = glue::glue("P = {format.pval(p_value, digits = 2)}")
    ),
    inherit.aes = FALSE,
    hjust = 1.05,
    size = 3,
    color = "grey35"
  ) +
  scale_x_log10(
    labels = scales::label_number(),
    expand = expansion(mult = c(0.08, 0.35))
  ) +
  scale_y_discrete(labels = scales::label_wrap(28)) +
  scale_color_manual(
    values = stats::setNames(
      unname(color_mgatk_gate[c(
        "Rejected by mgatk VMR/strand gate",
        "Passes mgatk VMR/strand gate"
      )]),
      c("median_rejected", "median_passed")
    ),
    labels = c("rejected by mgatk gate", "passes mgatk gate")
  ) +
  fn_theme() +
  labs(
    title = "The conclusion under four carrier definitions",
    subtitle = glue::glue(
      "same {nrow(d_e)} variants throughout \u00b7 the direction holds under ",
      "every definition; the censored one cannot resolve it because it has ",
      "no values below {CUTOFF_HETEROPLASMIC} \u00b7 {fn_sample_note()}"
    ),
    x = "Median allele frequency across carrier cells (log scale)",
    y = NULL,
    color = NULL
  )

log_info(
  "measure sensitivity: ",
  "{paste(d_f$measure, signif(d_f$p_value, 2), sep = ' P=', collapse = '; ')}"
)

# save ---------------------------------------------------------------------
export(
  d_c[, .(arm, af_bin, n_variants = N, fraction = frac)],
  as.character(fs::path(paths$tabdir, "04-af-bins.tsv"))
)
export(
  d_g[, .(arm, measured, af_bin = af_bin_own, n_variants = N, fraction = frac)],
  as.character(fs::path(paths$tabdir, "04-af-bins-own-definition.tsv"))
)
export(
  data.table::data.table(
    test = "Wilcoxon rank sum, median error-model carrier AF by mgatk gate call",
    subset = "variants passing the scMOCHA reliability gate",
    background_alt_rate = background_rate,
    n_passed = n_pass,
    n_rejected = n_rej,
    median_af_passed = af_pass,
    median_af_rejected = af_rej,
    location_shift = unname(wt$estimate),
    ci_low = wt$conf.int[1],
    ci_high = wt$conf.int[2],
    p_value = wt$p.value
  ),
  as.character(fs::path(paths$tabdir, "04-tests.tsv"))
)
export(
  d_f[, .(measure, n, median_rejected, median_passed, p_value)],
  as.character(fs::path(paths$tabdir, "04-measure-sensitivity.tsv"))
)

saveplot(
  as.character(fs::path(paths$figdir, "04a-af-ecdf.pdf")),
  fn_wrap_labs(p_a),
  width = 8,
  height = 5,
  device = cairo_pdf
)
saveplot(
  as.character(fs::path(paths$figdir, "04b-max-heteroplasmy.pdf")),
  fn_wrap_labs(p_b),
  width = 8,
  height = 5,
  device = cairo_pdf
)
saveplot(
  as.character(fs::path(paths$figdir, "04c-af-bins.pdf")),
  fn_wrap_labs(p_c, width = 55),
  width = 6.5,
  height = 5.5,
  device = cairo_pdf
)
saveplot(
  as.character(fs::path(paths$figdir, "04d-cell-level-af.pdf")),
  fn_wrap_labs(p_d),
  width = 8,
  height = 5,
  device = cairo_pdf
)
saveplot(
  as.character(fs::path(paths$figdir, "04e-gate-rejected-af.pdf")),
  fn_wrap_labs(p_e, width = 60),
  width = 7,
  height = 5.5,
  device = cairo_pdf
)
saveplot(
  as.character(fs::path(paths$figdir, "04f-measure-sensitivity.pdf")),
  fn_wrap_labs(p_f, width = 85),
  width = 8.5,
  height = 4.5,
  device = cairo_pdf
)
saveplot(
  as.character(fs::path(paths$figdir, "04g-af-bins-own-definition.pdf")),
  fn_wrap_labs(p_g, width = 62),
  width = 7,
  height = 5.5,
  device = cairo_pdf
)

log_info("saved 7 figures to {paths$figdir} and 4 tables to {paths$tabdir}")
