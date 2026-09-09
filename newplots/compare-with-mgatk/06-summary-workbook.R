#!/usr/bin/env Rscript
# Metainfo ----------------------------------------------------------------
# @AUTHOR: Chun-Jie Liu
# @CONTACT: chunjie.sam.liu.at.gmail.com
# @DATE: 2026-09-09
# @DESCRIPTION: Assemble every table this stage produced into one styled
#               workbook for the response letter, plus a rendered criterion
#               comparison table.
#               Usage: pixi run Rscript newplots/compare-with-mgatk/06-summary-workbook.R
# @VERSION: v0.1.0

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
    log_warn("missing table, sheet skipped: {f}")
    return(NULL)
  }
  data.table::as.data.table(import(as.character(f)))
}

sheets <- list(
  "00_Criteria" = criteria,
  "01_Definitions" = definitions,
  "02_Cell_inclusion" = read_tab("02-cell-inclusion.tsv"),
  "03_Funnel_counts" = read_tab("03-funnel-counts.tsv"),
  "04_Gate_crossapplied" = read_tab("03-gate-crossapplied.tsv"),
  "05_Exclusion_reasons" = read_tab("03-exclusion-reasons.tsv"),
  "06_Arm_combinations" = read_tab("03-arm-combinations.tsv"),
  "07_AF_bins" = read_tab("04-af-bins.tsv"),
  "08_AF_bins_own_definition" = read_tab("04-af-bins-own-definition.tsv"),
  "09_Gate_test" = read_tab("04-tests.tsv"),
  "10_Measure_sensitivity" = read_tab("04-measure-sensitivity.tsv"),
  "11_VMR_strand_rejected" = read_tab("05-vmr-strand-rejected.tsv"),
  "12_Arm_in_mgatk_plane" = read_tab("05-arm-in-mgatk-plane.tsv"),
  "13_Read_support" = read_tab("05-read-support.tsv"),
  "14_Variant_membership" = read_tab("03-variant-membership.tsv")
)
sheets <- sheets[!vapply(sheets, is.null, logical(1))]

P_COLS <- c("p_value")
NUM_COLS <- c(
  "fraction",
  "frac",
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
  "n_rejected"
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
