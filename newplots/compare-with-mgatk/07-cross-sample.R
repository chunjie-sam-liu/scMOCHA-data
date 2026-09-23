#!/usr/bin/env Rscript
# Metainfo ----------------------------------------------------------------
# @AUTHOR: Chun-Jie Liu
# @CONTACT: chunjie.sam.liu.at.gmail.com
# @DATE: 2026-09-17
# @DESCRIPTION: What holds across all five samples. Panel a, variants retained
#               by each arm in each sample; panel b, what the mgatk cell filter
#               costs per sample; panel c, strand correlation of mgatk's own S1
#               variants against its floor; panel d, which mgatk cutoff rejects
#               each S1 variant; panel e, heteroplasmy of what the gate
#               rejects, per sample.
#               Usage: pixi run Rscript newplots/compare-with-mgatk/07-cross-sample.R
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

# function ----------------------------------------------------------------
fn_sample_cache <- function(sample_id) {
  f <- fs::path(stage_paths(sample_id)$cachedir, "03-variant-gated.qs")
  if (!fs::file_exists(f)) {
    stop(
      "missing ",
      f,
      "; run 03-variant-funnel-overlap.R --sample=",
      sample_id
    )
  }
  d <- data.table::as.data.table(import(as.character(f)))
  d[, sample := sample_id]
  d[]
}

fn_sample_tab <- function(sample_id, name) {
  f <- fs::path(stage_paths(sample_id)$tabdir, name)
  if (!fs::file_exists(f)) {
    stop("missing ", f, "; run the step that writes it for ", sample_id)
  }
  data.table::as.data.table(import(as.character(f)))
}

# body ---------------------------------------------------------------------
variants <- data.table::rbindlist(
  lapply(SAMPLE_IDS, fn_sample_cache),
  fill = TRUE
)
variants[, sample := factor(sample, levels = SAMPLE_IDS)]

cell_inclusion <- data.table::dcast(
  data.table::rbindlist(
    lapply(SAMPLE_IDS, fn_sample_tab, name = "02-cell-inclusion.tsv")
  ),
  sample ~ metric,
  value.var = "value"
)
cell_inclusion[, sample := factor(sample, levels = SAMPLE_IDS)]

sample_lab <- stats::setNames(
  as.character(glue::glue_data(SAMPLES, "{gse}\n{gsm}")),
  SAMPLES$sample_id
)

# Two samples share a chemistry, so it cannot identify a sample; GSE plus GSM
# can. Chemistry stays in the SAMPLES registry and the workbook's 00_Samples.
SAMPLE_DISPLAY <- stats::setNames(
  as.character(glue::glue_data(SAMPLES, "{gse}_{gsm}")),
  SAMPLES$sample_id
)

# Tables are read next to the figures, so their sample column carries the same
# label the panels show; sample_id keeps the join back to the per-sample dirs.
fn_relabel <- function(d) {
  out <- data.table::copy(data.table::as.data.table(d))
  out[, sample_id := as.character(sample)]
  out[, sample := SAMPLE_DISPLAY[sample_id]]
  data.table::setcolorder(out, c("sample", "sample_id"))
  out[]
}

# Applied to the panels only. Every table below is still built from the
# unfiltered object, so the excluded samples keep their rows in the record.
# Levels are dropped so an excluded sample leaves no empty facet or legend key.
fn_fig_subset <- function(d) {
  out <- data.table::as.data.table(d)[sample %in% CROSS_SAMPLE_FIG_IDS]
  if (is.factor(out$sample)) {
    out[, sample := droplevels(sample)]
  }
  out[]
}

log_info(
  "loaded {nrow(variants)} variant rows from {length(SAMPLE_IDS)} samples"
)
log_info(
  "cross-sample figures show {length(CROSS_SAMPLE_FIG_IDS)} samples; ",
  "excluded from panels only: ",
  "{paste(CROSS_SAMPLE_FIG_EXCLUDE, collapse = ', ')}"
)

# a: variants retained by each arm ------------------------------------------
d_a <- data.table::melt(
  variants[,
    .(
      n_mgatk = sum(arm_mgatk),
      n_call = sum(arm_scmocha_call),
      n_af5 = sum(arm_scmocha_af5)
    ),
    by = sample
  ],
  id.vars = "sample",
  variable.name = "arm",
  value.name = "n_variants"
)
d_a[,
  arm := factor(
    data.table::fcase(
      arm == "n_mgatk" , ARM_MGATK   ,
      arm == "n_call"  , ARM_SCMOCHA ,
      default = ARM_SCMOCHA_AF5
    ),
    levels = ARM_LEVELS
  )
]

n_zero_mgatk <- d_a[arm == ARM_MGATK & n_variants == 0, .N]

d_a_fig <- fn_fig_subset(d_a)
n_zero_fig <- d_a_fig[arm == ARM_MGATK & n_variants == 0, .N]

# A zero is the headline here, so the axis has to be able to show it.
p_a <- d_a_fig |>
  ggplot(aes(x = sample, y = n_variants, fill = arm)) +
  geom_col(position = position_dodge2(width = 0.8, preserve = "single")) +
  geom_text(
    aes(label = scales::comma(n_variants)),
    position = position_dodge2(width = 0.8, preserve = "single"),
    vjust = -0.4,
    size = 2.6
  ) +
  scale_x_discrete(labels = sample_lab) +
  scale_y_continuous(
    labels = scales::comma,
    expand = expansion(mult = c(0, 0.18))
  ) +
  scale_fill_manual(values = color_arm) +
  fn_theme() +
  labs(
    title = "Variants retained by each arm, in every sample",
    subtitle = glue::glue(
      "original mgatk retains nothing in {n_zero_fig} of ",
      "{length(CROSS_SAMPLE_FIG_IDS)} samples shown: every S1 variant there ",
      "fails its strand r > {CUTOFF_STRAND_MGATK} floor \u00b7 linear axis, ",
      "so the zero bars are labelled rather than drawn \u00b7 the ",
      "{length(CROSS_SAMPLE_FIG_EXCLUDE)} 5' samples are left out of this ",
      "panel and kept in the tables"
    ),
    x = NULL,
    y = "Variants retained",
    fill = NULL
  )

# b: what the mgatk cell filter costs ---------------------------------------
d_b <- cell_inclusion[, .(
  sample,
  cells_total,
  cells_dropped_mgatk,
  detections_in_dropped_cells,
  variants_lost_by_cell_filter
)]
d_b[, frac_dropped := cells_dropped_mgatk / cells_total]

p_b <- fn_fig_subset(d_b) |>
  ggplot(aes(x = sample, y = frac_dropped)) +
  geom_col(aes(fill = sample), width = 0.6) +
  geom_text(
    aes(
      label = glue::glue(
        "{scales::comma(cells_dropped_mgatk)} /\n{scales::comma(cells_total)}"
      )
    ),
    vjust = -0.3,
    size = 2.8
  ) +
  scale_x_discrete(labels = sample_lab) +
  scale_y_continuous(
    labels = scales::percent,
    expand = expansion(mult = c(0, 0.2))
  ) +
  scale_fill_manual(values = color_sample) +
  fn_theme() +
  theme(legend.position = "none") +
  labs(
    title = "Cells discarded by the mgatk coverage filter",
    subtitle = glue::glue(
      "mean MT coverage > {CUTOFF_CELL_MEANCOV} \u00b7 scMOCHA applies no ",
      "cell filter \u00b7 labels are dropped over total cells"
    ),
    x = NULL,
    y = "Cells dropped"
  )

# c: where mgatk's own S1 variants sit against its strand floor --------------
# The gate that removes everything in four samples is the strand-correlation
# floor, so this panel shows the quantity that floor acts on. It is
# descriptive: the Spearman column tests whether strand correlation tracks
# coverage within a sample, and it does not in the deep sample, so no
# mechanism is claimed beyond what is drawn.
d_c <- variants[s1_mgatk == TRUE & !is.na(strand_mgatk)]

d_c_stat <- d_c[,
  {
    cv <- cov_mgatk
    keep <- !is.na(cv) & cv > 0
    rho <- NA_real_
    rho_p <- NA_real_
    if (sum(keep) >= 10L) {
      ct <- suppressWarnings(stats::cor.test(
        cv[keep],
        strand_mgatk[keep],
        method = "spearman"
      ))
      rho <- unname(ct$estimate)
      rho_p <- ct$p.value
    }
    .(
      n = .N,
      median_strand = stats::median(strand_mgatk),
      max_strand = max(strand_mgatk),
      n_above_floor = sum(strand_mgatk > CUTOFF_STRAND_MGATK),
      rho_strand_vs_coverage = rho,
      rho_p_value = rho_p
    )
  },
  by = sample
]
d_c_stat[, frac_above_floor := n_above_floor / n]

n_no_floor <- d_c_stat[n_above_floor == 0, .N]

log_info(
  "strand correlation vs mgatk floor: ",
  "{paste(d_c_stat$sample, d_c_stat$n_above_floor, sep = ' above=', collapse = '; ')}"
)

d_c_fig <- fn_fig_subset(d_c)
n_no_floor_fig <- d_c_fig[,
  sum(strand_mgatk > CUTOFF_STRAND_MGATK) == 0L,
  by = sample
][V1 == TRUE, .N]

p_c <- fn_or_empty(
  nrow(d_c_fig) > 0L,
  d_c_fig |>
    ggplot(aes(x = sample, y = strand_mgatk, color = sample)) +
    geom_jitter(width = 0.2, height = 0, size = 0.8, alpha = 0.55) +
    stat_summary(
      fun = stats::median,
      geom = "crossbar",
      width = 0.5,
      linewidth = 0.3,
      color = "grey25"
    ) +
    geom_hline(
      yintercept = CUTOFF_STRAND_MGATK,
      linetype = "dashed",
      color = "grey30"
    ) +
    scale_x_discrete(labels = sample_lab) +
    scale_color_manual(values = color_sample) +
    fn_theme() +
    theme(legend.position = "none") +
    labs(
      title = "Strand correlation of mgatk's own S1 variants",
      subtitle = glue::glue(
        "{nrow(d_c_fig)} variants with n_cells_conf_detected \u2265 ",
        "{CUTOFF_NCELLS_CONF} under mgatk's rule \u00b7 dashed line is its ",
        "strand r > {CUTOFF_STRAND_MGATK} floor \u00b7 in {n_no_floor_fig} ",
        "of {length(CROSS_SAMPLE_FIG_IDS)} samples shown not one variant ",
        "reaches it \u00b7 bar is the median"
      ),
      x = NULL,
      y = "Strand correlation (original mgatk)"
    ),
  "Strand correlation of mgatk's own S1 variants",
  "No mgatk S1 variant carries a strand correlation in any sample."
)

# d: which mgatk cutoff rejects each S1 variant ------------------------------
d_d <- variants[s1_mgatk == TRUE, .N, by = .(sample, excl_mgatk)]
d_d[, frac := N / sum(N), by = sample]

d_d_fig <- fn_fig_subset(d_d)

p_d <- fn_or_empty(
  nrow(d_d_fig) > 0L,
  d_d_fig |>
    ggplot(aes(x = sample, y = frac, fill = excl_mgatk)) +
    geom_col(width = 0.6) +
    scale_x_discrete(labels = sample_lab) +
    scale_y_continuous(labels = scales::percent, expand = c(0, 0)) +
    scale_fill_manual(values = color_exclusion) +
    fn_theme() +
    labs(
      title = "Which mgatk cutoff rejects each of its own S1 variants",
      subtitle = glue::glue(
        "every variant with n_cells_conf_detected \u2265 ",
        "{CUTOFF_NCELLS_CONF} under mgatk's rule, attributed to the first ",
        "criterion it fails \u00b7 {d_d_fig[, sum(N)]} variants over ",
        "{length(CROSS_SAMPLE_FIG_IDS)} samples shown"
      ),
      x = NULL,
      y = "Fraction of mgatk S1 variants",
      fill = NULL
    ),
  "Which mgatk cutoff rejects each of its own S1 variants",
  "No sample has an mgatk S1 variant."
)

# e: heteroplasmy of what the gate rejects, per sample -----------------------
d_e <- variants[gate_scmocha == TRUE & !is.na(af_carrier_median)]
d_e[,
  gate_call := factor(
    data.table::fifelse(
      gate_mgatk,
      "Passes mgatk VMR/strand gate",
      "Rejected by mgatk VMR/strand gate"
    ),
    levels = names(color_mgatk_gate)
  )
]

d_e_stat <- d_e[,
  {
    ok <- fn_testable(gate_call)
    p <- if (ok) {
      suppressWarnings(
        stats::wilcox.test(af_carrier_median ~ droplevels(gate_call))$p.value
      )
    } else {
      NA_real_
    }
    .(
      n_passed = sum(gate_mgatk),
      n_rejected = sum(!gate_mgatk),
      median_passed = stats::median(af_carrier_median[gate_mgatk]),
      median_rejected = stats::median(af_carrier_median[!gate_mgatk]),
      testable = ok,
      p_value = p
    )
  },
  by = sample
]

# A sample with no gate-passing variant produces no group, and a missing row is
# indistinguishable from a failed run. Every sample gets an explicit zero.
d_e_stat <- merge(
  data.table::data.table(sample = factor(SAMPLE_IDS, levels = SAMPLE_IDS)),
  d_e_stat,
  by = "sample",
  all.x = TRUE
)
d_e_stat[is.na(n_passed), `:=`(n_passed = 0L, n_rejected = 0L)]
d_e_stat[is.na(testable), testable := FALSE]

log_info(
  "per-sample gate test: ",
  "{paste(d_e_stat$sample, d_e_stat$testable, sep = ' testable=', collapse = '; ')}"
)

# Points with a median crossbar rather than a violin: several samples carry
# only a handful of variants, where a density estimate would be fiction.
d_e_fig <- fn_fig_subset(d_e)

p_e <- fn_or_empty(
  nrow(d_e_fig) > 0L,
  d_e_fig |>
    ggplot(aes(x = gate_call, y = af_carrier_median, color = gate_call)) +
    geom_jitter(width = 0.18, height = 0, size = 1.1, alpha = 0.7) +
    stat_summary(
      fun = stats::median,
      geom = "crossbar",
      width = 0.5,
      linewidth = 0.3,
      color = "grey25"
    ) +
    facet_wrap(~sample, nrow = 1, labeller = labeller(sample = sample_lab)) +
    scale_x_discrete(labels = c("passes", "rejected"), drop = FALSE) +
    scale_y_log10(labels = scales::label_number()) +
    scale_color_manual(values = color_mgatk_gate, drop = FALSE) +
    fn_theme() +
    theme(legend.position = "bottom") +
    labs(
      title = "Heteroplasmy of the variants mgatk's gate discards, per sample",
      subtitle = glue::glue(
        "variants passing the scMOCHA reliability gate, split by mgatk's ",
        "vmr > {CUTOFF_VMR_MGATK} and strand r > {CUTOFF_STRAND_MGATK} ",
        "\u00b7 the test is computable in ",
        "{d_e_stat[testable == TRUE & sample %in% CROSS_SAMPLE_FIG_IDS, .N]} ",
        "of {length(CROSS_SAMPLE_FIG_IDS)} samples shown; elsewhere mgatk's ",
        "gate passes nothing, so there is no comparison group \u00b7 bar is ",
        "the median"
      ),
      x = NULL,
      y = "Median AF across carrier cells (log scale)",
      color = NULL
    ),
  "Heteroplasmy of the variants mgatk's gate discards, per sample",
  "No variant passes the scMOCHA reliability gate in any sample."
)

# f, g: the call spectrum and its read support, pooled over the figure samples
# ---------------------------------------------------------------------------
# Built from the step 08 tables rather than recomputed, so the pooled panels
# and the per-sample ones are the same quantity by construction. Pooled over
# CROSS_SAMPLE_FIG_IDS only, matching every other panel in this step.
call_long_all <- data.table::rbindlist(
  lapply(CROSS_SAMPLE_FIG_IDS, fn_sample_tab, name = "08-call-af.tsv"),
  fill = TRUE
)
call_long_all[, measure := factor(measure, levels = CALL_AF_MEASURES)]
call_long_all[, af_bin := fn_call_af_bin(af)]

cells_all <- data.table::rbindlist(
  lapply(CROSS_SAMPLE_FIG_IDS, fn_sample_tab, name = "08-cell-af-depth.tsv"),
  fill = TRUE
)

pooled_note <- glue::glue(
  "{length(CROSS_SAMPLE_FIG_IDS)} samples pooled"
)

# Denominator for 07g: every cell in the figure samples, not just the ones
# holding a call.
n_cells_fig <- cell_inclusion[
  sample %in% CROSS_SAMPLE_FIG_IDS,
  sum(cells_total)
]

p_f <- fn_or_empty(
  nrow(call_long_all) > 0L,
  fn_plot_call_af_spectrum(call_long_all, pooled_note),
  "Allele frequency of every scMOCHA variant call",
  "No sample has a scMOCHA variant call."
)

p_g <- fn_or_empty(
  nrow(cells_all) > 0L,
  fn_plot_cell_af_depth(cells_all, pooled_note, n_cells_fig),
  "Read support behind each cell-level allele frequency",
  "No cell holds an alt read at a called variant position."
)

n_calls_pooled <- call_long_all[
  measure == levels(call_long_all$measure)[1],
  .N
]
log_info(
  "pooled call spectrum: {n_calls_pooled} calls over ",
  "{data.table::uniqueN(call_long_all$variant)} distinct variants; ",
  "{nrow(cells_all)} cell-variant pairs from ",
  "{nrow(unique(cells_all[, .(sample, barcode)]))} of {n_cells_fig} cells"
)

# overview ------------------------------------------------------------------
d_ov <- merge(
  SAMPLES,
  variants[,
    .(
      s0_scmocha = sum(candidate_scmocha),
      s0_mgatk = sum(candidate_mgatk),
      s1_scmocha = sum(s1_scmocha),
      s1_mgatk = sum(s1_mgatk),
      retained_mgatk = sum(arm_mgatk),
      retained_scmocha_call = sum(arm_scmocha_call),
      retained_scmocha_af5 = sum(arm_scmocha_af5)
    ),
    by = sample
  ],
  by.x = "sample_id",
  by.y = "sample"
)
d_ov <- merge(
  d_ov,
  cell_inclusion[, .(sample, cells_total, cells_dropped_mgatk)],
  by.x = "sample_id",
  by.y = "sample"
)
stopifnot(nrow(d_ov) == length(SAMPLE_IDS))

# The overview is keyed on sample_id rather than sample, so it is relabelled
# directly instead of through fn_relabel().
d_ov_out <- data.table::copy(d_ov)
d_ov_out[, chemistry := NULL]
d_ov_out[, sample := SAMPLE_DISPLAY[sample_id]]
data.table::setcolorder(d_ov_out, c("sample", "sample_id"))

# save ---------------------------------------------------------------------
export(
  d_ov_out,
  as.character(fs::path(paths$tabdir, "07-sample-overview.tsv"))
)
export(
  fn_relabel(d_a[, .(sample, arm, n_variants)]),
  as.character(fs::path(paths$tabdir, "07-arm-yield.tsv"))
)
export(
  fn_relabel(d_b),
  as.character(fs::path(paths$tabdir, "07-cell-filter.tsv"))
)
export(
  fn_relabel(d_c_stat),
  as.character(fs::path(paths$tabdir, "07-strand-support.tsv"))
)
export(
  fn_relabel(
    d_d[, .(sample, excluded_by = as.character(excl_mgatk), n_variants = N)]
  ),
  as.character(fs::path(paths$tabdir, "07-exclusion-reasons.tsv"))
)
export(
  fn_relabel(d_e_stat),
  as.character(fs::path(paths$tabdir, "07-gate-test.tsv"))
)

# The counts behind 07f, so "n of N calls below the cutoff" is auditable
# rather than read off the panel.
d_bins <- call_long_all[,
  .(n_calls = .N),
  by = .(measure = as.character(measure), af_bin)
]
d_bins[, fraction := n_calls / sum(n_calls), by = measure]
data.table::setorder(d_bins, measure, af_bin)
export(
  d_bins,
  as.character(fs::path(paths$tabdir, "07-call-af-bins.tsv"))
)

saveplot(
  as.character(fs::path(paths$figdir, "07a-arm-yield.pdf")),
  fn_wrap_labs(p_a, width = 95),
  width = 9,
  height = 5.5,
  device = cairo_pdf
)
saveplot(
  as.character(fs::path(paths$figdir, "07b-cell-filter.pdf")),
  fn_wrap_labs(p_b, width = 85),
  width = 8,
  height = 5,
  device = cairo_pdf
)
saveplot(
  as.character(fs::path(paths$figdir, "07c-strand-correlation.pdf")),
  fn_wrap_labs(p_c, width = 95),
  width = 9,
  height = 5.5,
  device = cairo_pdf
)
saveplot(
  as.character(fs::path(paths$figdir, "07d-exclusion-reasons.pdf")),
  fn_wrap_labs(p_d, width = 90),
  width = 8.5,
  height = 5.5,
  device = cairo_pdf
)
saveplot(
  as.character(fs::path(paths$figdir, "07e-gate-rejected-af.pdf")),
  fn_wrap_labs(p_e, width = 110),
  width = 13,
  height = 5,
  device = cairo_pdf
)
saveplot(
  as.character(fs::path(paths$figdir, "07f-call-af-spectrum.pdf")),
  fn_wrap_labs(p_f, width = 95),
  width = 8.5,
  height = 7,
  device = cairo_pdf
)
saveplot(
  as.character(fs::path(paths$figdir, "07g-cell-af-depth.pdf")),
  fn_wrap_labs(p_g, width = 95),
  width = 8.5,
  height = 5.5,
  device = cairo_pdf
)

log_info(
  "07: {length(SAMPLE_IDS)} samples in the tables, ",
  "{length(CROSS_SAMPLE_FIG_IDS)} in the figures; mgatk retains nothing in ",
  "{n_zero_mgatk}; no variant reaches the strand floor in {n_no_floor}; ",
  "gate test computable in {d_e_stat[testable == TRUE, .N]}. ",
  "7 figures and 7 tables written"
)
