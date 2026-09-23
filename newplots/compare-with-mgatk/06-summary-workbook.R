#!/usr/bin/env Rscript
# Metainfo ----------------------------------------------------------------
# @AUTHOR: Chun-Jie Liu
# @CONTACT: chunjie.sam.liu.at.gmail.com
# @DATE: 2026-09-09
# @DESCRIPTION: Assemble every table this stage produced, for every sample,
#               into one styled workbook for the response letter, plus a
#               rendered criterion comparison table. Runs after 07.
#               Usage: pixi run Rscript newplots/compare-with-mgatk/06-summary-workbook.R
# @VERSION: v0.2.0

# Reproducibility ----------------------------------------------------------
set.seed(9527)

# Library ------------------------------------------------------------------
suppressMessages({
  library(jutils)
  load_pkg(openxlsx2)
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
fs::dir_create(paths$tabdir)

# body ---------------------------------------------------------------------
criteria <- data.table::data.table(
  id = c("C1", "C2", "C3", "C4", "C5", "C6", "C7"),
  stage = c(
    "Cell inclusion",
    "Confident detection in one cell",
    "Variant retention",
    "Reliability gate",
    "Cell-by-variant AF matrix",
    "Reliability gate, scMOCHA side",
    "Downstream classification"
  ),
  original_mgatk = c(
    glue::glue("mean MT coverage > {CUTOFF_CELL_MEANCOV}"),
    glue::glue(
      "fwd alt >= {CUTOFF_ALT_STRAND} and rev alt >= {CUTOFF_ALT_STRAND}"
    ),
    glue::glue("n_cells_conf_detected >= {CUTOFF_NCELLS_CONF}"),
    glue::glue(
      "vmr > {CUTOFF_VMR_MGATK} and strand correlation > {CUTOFF_STRAND_MGATK}"
    ),
    "only variants passing C3 are written",
    "none beyond C4",
    "none; mgatk stops at the variant list"
  ),
  scmocha = c(
    "no cell filter; all barcodes retained",
    glue::glue(
      "fwd alt >= {CUTOFF_ALT_STRAND}, rev alt >= {CUTOFF_ALT_STRAND} and ",
      "fwd + rev alt >= {CUTOFF_ALT_READS}"
    ),
    glue::glue("n_cells_conf_detected >= {CUTOFF_NCELLS_CONF} (identical)"),
    "no VMR filter, no strand-correlation filter",
    "a raw matrix with every candidate variant is also written",
    glue::glue(
      "mis-alignment and RNA-editing position blacklist, then >= ",
      "{CUTOFF_NOTRELIABLE} cells at cell AF >= {CUTOFF_HETEROPLASMIC} with ",
      "cell depth >= {CUTOFF_MIN_READS}"
    ),
    glue::glue(
      "homoplasmic / heteroplasmic / somatic from cross-cell-type AF at ",
      "{CUTOFF_HETEROPLASMIC} and {CUTOFF_HOMOPLASMIC}, cluster reads >= ",
      "{CUTOFF_MIN_READS}"
    )
  ),
  direction = c(
    "scMOCHA more permissive",
    "scMOCHA more stringent",
    "same",
    "different basis",
    "scMOCHA more permissive",
    "different basis",
    "different question"
  )
)

definitions <- data.table::data.table(
  term = c(
    "S0 candidates",
    "S1",
    "S2",
    "detection",
    "population AF",
    "VMR",
    "strand correlation",
    "alt reads"
  ),
  meaning = c(
    "every non-reference allele observed on both strands across cells",
    glue::glue(
      "variants with n_cells_conf_detected >= {CUTOFF_NCELLS_CONF}; the rule ",
      "is identical, the underlying per-cell criterion is not"
    ),
    "variants surviving that caller's reliability gate",
    glue::glue(
      "one cell with AF >= {CUTOFF_HETEROPLASMIC} and depth >= ",
      "{CUTOFF_MIN_READS} at the variant position"
    ),
    paste(
      "total alt reads divided by total coverage across all cells;",
      "taken from the scMOCHA table for every variant so the two sets are",
      "compared on one measurement"
    ),
    "variance-to-mean ratio of per-cell AF across cells",
    "correlation between forward and reverse alt counts across cells",
    "reconstructed as AF x depth; the callers store AF and coverage"
  )
)

read_tab <- function(name) {
  f <- fs::path(paths$tabdir, name)
  if (!fs::file_exists(f)) {
    log_warn("missing cross-sample table, sheet skipped: {f}")
    return(NULL)
  }
  data.table::as.data.table(import(as.character(f)))
}

# One sheet per topic with a sample column, rather than one workbook per
# sample: five workbooks of the same 15 sheets would not be opened.
read_sample_tab <- function(name) {
  parts <- lapply(SAMPLE_IDS, function(s) {
    f <- fs::path(stage_paths(s)$tabdir, name)
    if (!fs::file_exists(f)) {
      log_warn("missing table for {s}, sample omitted from sheet: {name}")
      return(NULL)
    }
    d <- data.table::as.data.table(import(as.character(f)))
    # A header-only file reads back with every column typed logical, which
    # collides with the other samples on rbind. Those samples contribute no
    # rows anyway; the zero itself is recorded in the count tables.
    if (nrow(d) == 0L) {
      log_info("empty table for {s}, no rows contributed to sheet: {name}")
      return(NULL)
    }
    d[, sample := as.character(sample)]
    d[]
  })
  parts <- parts[!vapply(parts, is.null, logical(1))]
  if (length(parts) == 0L) {
    return(NULL)
  }
  # A column that is whole-numbered in a shallow sample and fractional in a
  # deep one comes back with a different class from each file.
  data.table::rbindlist(parts, fill = TRUE, ignore.attr = TRUE)
}

sheets <- list(
  "00_Samples" = SAMPLES,
  "01_Criteria" = criteria,
  "02_Definitions" = definitions,
  "03_Overview" = read_tab("07-sample-overview.tsv"),
  "04_Arm_yield" = read_tab("07-arm-yield.tsv"),
  "05_Cell_filter" = read_tab("07-cell-filter.tsv"),
  "06_Strand_correlation" = read_tab("07-strand-support.tsv"),
  "07_Exclusion_pooled" = read_tab("07-exclusion-reasons.tsv"),
  "08_Gate_test_per_sample" = read_tab("07-gate-test.tsv"),
  "09_Call_AF_bins_pooled" = read_tab("07-call-af-bins.tsv"),
  "10_Cell_inclusion" = read_sample_tab("02-cell-inclusion.tsv"),
  "11_Funnel_counts" = read_sample_tab("03-funnel-counts.tsv"),
  "12_Gate_crossapplied" = read_sample_tab("03-gate-crossapplied.tsv"),
  "13_Exclusion_reasons" = read_sample_tab("03-exclusion-reasons.tsv"),
  "14_Arm_combinations" = read_sample_tab("03-arm-combinations.tsv"),
  "15_AF_bins" = read_sample_tab("04-af-bins.tsv"),
  "16_AF_bins_own_definition" = read_sample_tab(
    "04-af-bins-own-definition.tsv"
  ),
  "17_Gate_test" = read_sample_tab("04-tests.tsv"),
  "18_Measure_sensitivity" = read_sample_tab("04-measure-sensitivity.tsv"),
  "19_VMR_strand_rejected" = read_sample_tab("05-vmr-strand-rejected.tsv"),
  "20_Arm_in_mgatk_plane" = read_sample_tab("05-arm-in-mgatk-plane.tsv"),
  "21_Arm_in_scmocha_plane" = read_sample_tab(
    "05-arm-in-scmocha-plane.tsv"
  ),
  "22_Read_support" = read_sample_tab("05-read-support.tsv"),
  # The cell-level companion 08-cell-af-depth.tsv is deliberately not a sheet:
  # it runs to several hundred thousand rows across the samples. It stays a
  # TSV next to the figure it backs.
  "23_Call_AF" = read_sample_tab("08-call-af.tsv"),
  "24_Variant_membership" = read_sample_tab("03-variant-membership.tsv")
)
sheets <- sheets[!vapply(sheets, is.null, logical(1))]

P_COLS <- c("p_value")
NUM_COLS <- c(
  "fraction",
  "frac",
  "frac_dropped",
  "rho",
  "background_alt_rate",
  "median_af_passed",
  "median_af_rejected",
  "median_rejected",
  "median_passed",
  "location_shift",
  "ci_low",
  "ci_high",
  "median_alt_reads"
)
INT_COLS <- c(
  "value",
  "n",
  "n_variants",
  "n_cells",
  "n_detections",
  "n_passed",
  "n_rejected",
  "cells_total",
  "cells_dropped_mgatk",
  "detections_in_dropped_cells",
  "variants_lost_by_cell_filter",
  "s0_scmocha",
  "s0_mgatk",
  "s1_scmocha",
  "s1_mgatk",
  "retained_mgatk",
  "retained_scmocha_call",
  "retained_scmocha_af5"
)

wb <- wb_workbook()

for (nm in names(sheets)) {
  dt <- sheets[[nm]]
  wb$add_worksheet(nm)$add_data(sheet = nm, x = dt, with_filter = TRUE)
  wb$set_col_widths(sheet = nm, cols = seq_len(ncol(dt)), widths = "auto")
  wb$freeze_pane(sheet = nm, first_active_row = 2)
  wb$add_fill(
    sheet = nm,
    dims = wb_dims(rows = 1, cols = seq_len(ncol(dt))),
    color = wb_color(hex = color_xlsx_hdr)
  )
  wb$add_font(
    sheet = nm,
    dims = wb_dims(rows = 1, cols = seq_len(ncol(dt))),
    bold = TRUE,
    color = wb_color(hex = color_xlsx_white)
  )

  rows <- seq_len(nrow(dt)) + 1L
  for (spec in list(
    list(cols = P_COLS, fmt = "0.00E+00"),
    list(cols = NUM_COLS, fmt = "0.0000"),
    list(cols = INT_COLS, fmt = "#,##0")
  )) {
    hit <- which(names(dt) %in% spec$cols)
    if (length(hit) > 0) {
      wb$add_numfmt(
        sheet = nm,
        dims = wb_dims(rows = rows, cols = hit),
        numfmt = spec$fmt
      )
    }
  }
}

# save ---------------------------------------------------------------------
outfile <- fs::path(paths$tabdir, "06-compare-with-mgatk.xlsx")
wb$save(as.character(outfile))

export(
  criteria,
  as.character(fs::path(
    paths$tabdir,
    "06-criteria-comparison.tsv"
  ))
)

log_info(
  "wrote {length(sheets)} sheets to {outfile}: ",
  "{paste(names(sheets), collapse = ', ')}"
)
