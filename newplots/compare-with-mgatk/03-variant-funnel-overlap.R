#!/usr/bin/env Rscript
# Metainfo ----------------------------------------------------------------
# @AUTHOR: Chun-Jie Liu
# @CONTACT: chunjie.sam.liu.at.gmail.com
# @DATE: 2026-09-09
# @DESCRIPTION: Criteria C3, C4 and C6 across three calling arms: original
#               mgatk, the scMOCHA variant call, and scMOCHA after its
#               downstream AF>5% gate. Panel a, the funnel; panel b, the
#               three-way overlap; panel c, each reliability gate applied to
#               each S1 set; panels d and e, the cutoff that excludes each
#               arm-specific variant from the other arm; panel f, every arm
#               membership combination.
#               Usage: pixi run Rscript newplots/compare-with-mgatk/03-variant-funnel-overlap.R
# @VERSION: v0.2.0

# Reproducibility ----------------------------------------------------------
set.seed(9527)

# Library ------------------------------------------------------------------
suppressMessages({
  library(jutils)
  load_pkg(scales, ggvenn)
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
cell_af <- import(as.character(fs::path(paths$cachedir, "01-cell-af-s1.qs")))
cell_pos_coverage <- import(as.character(fs::path(
  paths$cachedir,
  "01-cell-pos-coverage.qs"
)))

# body ---------------------------------------------------------------------
detection <- fn_detection_long(
  cell_af,
  cell_pos_coverage,
  stats_joined[, .(variant, position)]
)

blacklisted <- stats_joined[blacklisted == TRUE, variant]

background_rate <- fn_background_rate(detection)
log_info(
  "background alt rate {signif(background_rate, 3)} from ",
  "{format(nrow(detection), big.mark = ',')} cell-by-variant observations"
)

# The scMOCHA gate is evaluated on all cells for both rule sets, so this panel
# isolates the gate; the cell filter is quantified separately in step 02.
gate_scmocha_pass <- fn_scmocha_gate(detection, blacklisted)

stats_joined[, `:=`(
  gate_scmocha = variant %in% gate_scmocha_pass,
  gate_mgatk = !is.na(vmr_mgatk) &
    !is.na(strand_mgatk) &
    vmr_mgatk > CUTOFF_VMR_MGATK &
    strand_mgatk > CUTOFF_STRAND_MGATK
)]

stats_joined[, `:=`(
  final_mgatk = s1_mgatk & gate_mgatk,
  final_scmocha = s1_scmocha & gate_scmocha
)]

# Carrier-level heteroplasmy travels with the variant table so steps 04 and 05
# never fall back on the population AF column.
carrier <- fn_carrier_stats(detection, background_rate)
stats_joined[
  carrier,
  on = "variant",
  `:=`(
    n_cells_gate = i.n_cells_gate,
    af_gate_median = i.af_gate_median,
    af_gate_max = i.af_gate_max,
    n_cells_loose = i.n_cells_loose,
    af_loose_median = i.af_loose_median,
    n_cells_carrier = i.n_cells_carrier,
    af_carrier_median = i.af_carrier_median,
    af_carrier_min = i.af_carrier_min,
    af_carrier_max = i.af_carrier_max
  )
]
stats_joined[is.na(n_cells_gate), n_cells_gate := 0L]
stats_joined[is.na(n_cells_loose), n_cells_loose := 0L]
stats_joined[is.na(n_cells_carrier), n_cells_carrier := 0L]

# Arm membership. The scMOCHA variant call applies no AF filter; AF>5% is the
# downstream gate, so it is a third arm rather than part of the call.
stats_joined[, `:=`(
  arm_mgatk = final_mgatk,
  arm_scmocha_call = s1_scmocha,
  arm_scmocha_af5 = final_scmocha
)]

set_mgatk <- stats_joined[arm_mgatk == TRUE, variant]
set_call <- stats_joined[arm_scmocha_call == TRUE, variant]
set_af5 <- stats_joined[arm_scmocha_af5 == TRUE, variant]

log_info(
  "arms: mgatk {length(set_mgatk)}, scMOCHA call {length(set_call)}, ",
  "scMOCHA AF>5% {length(set_af5)}"
)
log_info(
  "mgatk vs AF>5%: shared {length(intersect(set_mgatk, set_af5))}, ",
  "mgatk only {length(setdiff(set_mgatk, set_af5))}, ",
  "AF>5% only {length(setdiff(set_af5, set_mgatk))}"
)

# a: funnel -----------------------------------------------------------------
d_a <- data.table::data.table(
  arm = factor(
    c(rep(ARM_MGATK, 3), rep(ARM_SCMOCHA, 2), rep(ARM_SCMOCHA_AF5, 3)),
    levels = ARM_LEVELS
  ),
  stage = factor(
    c(
      "S0 candidates",
      "S1 confident cells",
      "S2 reliability gate",
      "S0 candidates",
      "S1 confident cells",
      "S0 candidates",
      "S1 confident cells",
      "S2 reliability gate"
    ),
    levels = c("S0 candidates", "S1 confident cells", "S2 reliability gate")
  ),
  n = c(
    stats_joined[candidate_mgatk == TRUE, .N],
    stats_joined[s1_mgatk == TRUE, .N],
    length(set_mgatk),
    stats_joined[candidate_scmocha == TRUE, .N],
    length(set_call),
    stats_joined[candidate_scmocha == TRUE, .N],
    stats_joined[s1_scmocha == TRUE, .N],
    length(set_af5)
  )
)

p_a <- d_a |>
  ggplot(aes(x = stage, y = n, fill = arm)) +
  geom_col(position = position_dodge2(width = 0.8, preserve = "single")) +
  geom_text(
    aes(label = scales::comma(n)),
    position = position_dodge2(width = 0.8, preserve = "single"),
    vjust = -0.4,
    size = 2.8
  ) +
  scale_y_log10(
    labels = scales::comma,
    expand = expansion(mult = c(0, 0.15))
  ) +
  scale_fill_manual(values = color_arm) +
  fn_theme() +
  labs(
    title = "Variants surviving each filtering stage",
    subtitle = glue::glue(
      "S1 is n_cells_conf_detected \u2265 {CUTOFF_NCELLS_CONF} in all arms, ",
      "from different per-cell rules \u00b7 S2 is vmr > {CUTOFF_VMR_MGATK} and ",
      "strand r > {CUTOFF_STRAND_MGATK} for mgatk, blacklist plus \u2265 ",
      "{CUTOFF_NOTRELIABLE} cells at AF \u2265 {CUTOFF_HETEROPLASMIC} for ",
      "scMOCHA AF>5% \u00b7 the scMOCHA variant call has no S2 \u00b7 ",
      "{fn_sample_note()}"
    ),
    x = NULL,
    y = "Variants (log scale)",
    fill = NULL
  )

# b: overlap ----------------------------------------------------------------
venn_input <- stats::setNames(
  list(set_mgatk, set_call, set_af5),
  ARM_LEVELS
)

p_b <- ggvenn::ggvenn(
  venn_input,
  fill_color = unname(color_arm[ARM_LEVELS]),
  fill_alpha = 0.4,
  stroke_size = 0.4,
  set_name_size = 3.6,
  text_size = 3.4
) +
  labs(
    title = "Variant sets of the three arms",
    subtitle = glue::glue(
      "{length(set_mgatk)} mgatk \u00b7 {length(set_call)} scMOCHA call ",
      "\u00b7 {length(set_af5)} scMOCHA AF>5% \u00b7 {fn_sample_note()}"
    )
  ) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5, color = "grey40", size = 10)
  )

# c: gates cross-applied ----------------------------------------------------
d_c <- data.table::CJ(
  s1 = CALLER_LEVELS,
  gate = CALLER_LEVELS
)[,
  n := c(
    stats_joined[s1_mgatk == TRUE & gate_mgatk == TRUE, .N],
    stats_joined[s1_mgatk == TRUE & gate_scmocha == TRUE, .N],
    stats_joined[s1_scmocha == TRUE & gate_mgatk == TRUE, .N],
    stats_joined[s1_scmocha == TRUE & gate_scmocha == TRUE, .N]
  )
]
d_c[, `:=`(
  s1 = factor(s1, levels = CALLER_LEVELS),
  gate = factor(gate, levels = CALLER_LEVELS)
)]

p_c <- d_c |>
  ggplot(aes(x = gate, y = n, fill = gate)) +
  geom_col(width = 0.6) +
  geom_text(aes(label = scales::comma(n)), vjust = -0.4, size = 3.2) +
  facet_wrap(~s1, labeller = labeller(s1 = \(x) glue::glue("S1 set: {x}"))) +
  scale_y_continuous(
    labels = scales::comma,
    expand = expansion(mult = c(0, 0.15))
  ) +
  scale_fill_manual(values = color_caller) +
  fn_theme() +
  theme(legend.position = "none") +
  labs(
    title = "Each reliability gate applied to each S1 set",
    subtitle = glue::glue(
      "left panel holds the confident-detection rule at mgatk's, right at ",
      "scMOCHA's \u00b7 within a panel the difference is the gate alone \u00b7 ",
      "gate evaluated on all {data.table::uniqueN(detection$barcode)} cells"
    ),
    x = "Reliability gate applied",
    y = "Variants"
  )

# d, e: which cutoff excludes each arm-specific variant ---------------------
stats_joined[, `:=`(
  excl_mgatk = fn_exclusion_mgatk(stats_joined),
  excl_af5 = fn_exclusion_scmocha_af5(stats_joined),
  excl_call = fn_exclusion_scmocha_call(stats_joined)
)]

# scMOCHA AF>5% variants that mgatk does not report, and why.
d_d <- stats_joined[
  arm_scmocha_af5 == TRUE & arm_mgatk == FALSE,
  .N,
  by = excl_mgatk
][order(-N)]

p_d <- d_d |>
  ggplot(aes(x = stats::reorder(excl_mgatk, N), y = N, fill = excl_mgatk)) +
  geom_col(width = 0.65) +
  geom_text(aes(label = N), hjust = -0.25, size = 3.2) +
  coord_flip() +
  scale_y_continuous(
    labels = scales::comma,
    expand = expansion(mult = c(0, 0.18))
  ) +
  scale_fill_manual(values = color_exclusion) +
  fn_theme() +
  theme(legend.position = "none") +
  labs(
    title = "Why original mgatk misses the scMOCHA AF>5% variants",
    subtitle = glue::glue(
      "{sum(d_d$N)} variants in scMOCHA AF>5% but not in mgatk \u00b7 each is ",
      "attributed to the first mgatk criterion it fails \u00b7 {fn_sample_note()}"
    ),
    x = NULL,
    y = "Variants"
  )

# mgatk variants that scMOCHA AF>5% does not report, and why.
d_e <- stats_joined[
  arm_mgatk == TRUE & arm_scmocha_af5 == FALSE,
  .N,
  by = excl_af5
][order(-N)]

p_e <- d_e |>
  ggplot(aes(x = stats::reorder(excl_af5, N), y = N, fill = excl_af5)) +
  geom_col(width = 0.65) +
  geom_text(aes(label = N), hjust = -0.25, size = 3.2) +
  coord_flip() +
  scale_y_continuous(
    labels = scales::comma,
    expand = expansion(mult = c(0, 0.18))
  ) +
  scale_fill_manual(values = color_exclusion) +
  fn_theme() +
  theme(legend.position = "none") +
  labs(
    title = "Why scMOCHA AF>5% misses the original mgatk variants",
    subtitle = glue::glue(
      "{sum(d_e$N)} variants in mgatk but not in scMOCHA AF>5% \u00b7 each is ",
      "attributed to the first scMOCHA criterion it fails \u00b7 ",
      "{fn_sample_note()}"
    ),
    x = NULL,
    y = "Variants"
  )

# f: arm membership combinations --------------------------------------------
d_f <- stats_joined[
  arm_mgatk | arm_scmocha_call | arm_scmocha_af5,
  .N,
  by = .(arm_mgatk, arm_scmocha_call, arm_scmocha_af5)
]
d_f[,
  combination := paste(
    data.table::fifelse(arm_mgatk, "mgatk", ""),
    data.table::fifelse(arm_scmocha_call, "call", ""),
    data.table::fifelse(arm_scmocha_af5, "AF>5%", ""),
    sep = " "
  )
]
d_f[, combination := trimws(gsub(" +", " + ", trimws(combination)))]
d_f <- d_f[order(-N)]

p_f <- d_f |>
  ggplot(aes(x = stats::reorder(combination, N), y = N)) +
  geom_col(width = 0.65, fill = "grey45") +
  geom_text(aes(label = scales::comma(N)), hjust = -0.25, size = 3.2) +
  coord_flip() +
  scale_y_continuous(
    labels = scales::comma,
    expand = expansion(mult = c(0, 0.18))
  ) +
  fn_theme() +
  labs(
    title = "Every arm membership combination",
    subtitle = glue::glue(
      "{d_f[, sum(N)]} variants reported by at least one arm \u00b7 ",
      "scMOCHA AF>5% is a subset of the scMOCHA variant call by construction ",
      "\u00b7 {fn_sample_note()}"
    ),
    x = NULL,
    y = "Variants"
  )

# save ---------------------------------------------------------------------
export(d_a, as.character(fs::path(paths$tabdir, "03-funnel-counts.tsv")))
export(
  d_c[, .(s1_set = s1, gate_applied = gate, n_variants = n)],
  as.character(fs::path(paths$tabdir, "03-gate-crossapplied.tsv"))
)
export(
  stats_joined[
    s1_scmocha == TRUE | s1_mgatk == TRUE,
    .(
      variant,
      position,
      nucleotide,
      s1_mgatk,
      s1_scmocha,
      gate_mgatk,
      gate_scmocha,
      arm_mgatk,
      arm_scmocha_call,
      arm_scmocha_af5,
      excl_mgatk,
      excl_call,
      excl_af5,
      blacklisted,
      n_cells_gate,
      af_gate_median,
      af_gate_max,
      n_cells_loose,
      af_loose_median,
      n_cells_carrier,
      af_carrier_median,
      af_carrier_min,
      af_carrier_max,
      vmr_mgatk,
      strand_mgatk
    )
  ],
  as.character(fs::path(paths$tabdir, "03-variant-membership.tsv"))
)
export(
  data.table::rbindlist(list(
    d_d[, .(
      direction = "in scMOCHA AF>5%, not in mgatk",
      excluded_by = as.character(excl_mgatk),
      n_variants = N
    )],
    d_e[, .(
      direction = "in mgatk, not in scMOCHA AF>5%",
      excluded_by = as.character(excl_af5),
      n_variants = N
    )]
  )),
  as.character(fs::path(paths$tabdir, "03-exclusion-reasons.tsv"))
)
export(
  d_f[, .(combination, n_variants = N)],
  as.character(fs::path(paths$tabdir, "03-arm-combinations.tsv"))
)
export(
  stats_joined,
  as.character(fs::path(paths$cachedir, "03-variant-gated.qs"))
)

saveplot(
  as.character(fs::path(paths$figdir, "03a-variant-funnel.pdf")),
  fn_wrap_labs(p_a),
  width = 9,
  height = 5.5,
  device = cairo_pdf
)
saveplot(
  as.character(fs::path(paths$figdir, "03b-variant-overlap.pdf")),
  fn_wrap_labs(p_b, width = 55),
  width = 6.5,
  height = 6,
  device = cairo_pdf
)
saveplot(
  as.character(fs::path(paths$figdir, "03c-gate-crossapplied.pdf")),
  fn_wrap_labs(p_c),
  width = 8,
  height = 5,
  device = cairo_pdf
)
saveplot(
  as.character(fs::path(paths$figdir, "03d-excluded-from-mgatk.pdf")),
  fn_wrap_labs(p_d, width = 80),
  width = 8,
  height = 4.5,
  device = cairo_pdf
)
saveplot(
  as.character(fs::path(paths$figdir, "03e-excluded-from-scmocha-af5.pdf")),
  fn_wrap_labs(p_e, width = 80),
  width = 8,
  height = 4.5,
  device = cairo_pdf
)
saveplot(
  as.character(fs::path(paths$figdir, "03f-arm-combinations.pdf")),
  fn_wrap_labs(p_f, width = 80),
  width = 8,
  height = 4.5,
  device = cairo_pdf
)

log_info(
  "03: excluded from mgatk by ",
  "{paste(d_d$excl_mgatk, d_d$N, sep = '=', collapse = ', ')}"
)
log_info(
  "03: excluded from scMOCHA AF>5% by ",
  "{paste(d_e$excl_af5, d_e$N, sep = '=', collapse = ', ')}"
)
log_info("saved 6 figures to {paths$figdir} and 5 tables to {paths$tabdir}")
