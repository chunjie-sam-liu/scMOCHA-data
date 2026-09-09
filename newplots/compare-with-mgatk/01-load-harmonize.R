#!/usr/bin/env Rscript
# Metainfo ----------------------------------------------------------------
# @AUTHOR: Chun-Jie Liu
# @CONTACT: chunjie.sam.liu.at.gmail.com
# @DATE: 2026-09-09
# @DESCRIPTION: Join the scMOCHA and original-mgatk variant_stats tables, tag
#               each variant with the stage it survives under each rule set,
#               and extract per-cell depth at candidate variant positions so
#               the scMOCHA reliability gate can be evaluated downstream.
#               Usage: pixi run Rscript newplots/compare-with-mgatk/01-load-harmonize.R
# @VERSION: v0.1.0

# Reproducibility ----------------------------------------------------------
set.seed(9527)

# Library ------------------------------------------------------------------
suppressMessages({
  library(jutils)
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
source(fs::path(
  fs::path_expand(Sys.getenv("REPODIR")),
  "newplots",
  "compare-with-mgatk",
  "config.R"
))

paths <- stage_paths()
inputs <- stage_inputs(paths)
fs::dir_create(paths$cachedir)

for (f in inputs) {
  if (!fs::file_exists(f)) {
    log_error("missing input: {f}")
    stop("missing input: ", f)
  }
}

# function ----------------------------------------------------------------
fn_read_stats <- function(file, caller) {
  data.table::fread(cmd = glue::glue("zcat {file}")) |>
    dplyr::mutate(caller = caller) |>
    data.table::as.data.table()
}

# The AF matrix header carries an empty first field, so header = TRUE would
# consume the first barcode as the column name and lose a cell.
fn_read_af_barcodes <- function(file) {
  data.table::fread(
    cmd = glue::glue("zcat {file} | cut -f1"),
    header = FALSE
  )[[1]]
}

# Streams the per-cell coverage file through awk so only the positions that
# carry a candidate variant ever reach R.
fn_read_coverage_at <- function(file, positions) {
  posfile <- fs::path(fs::path_expand("~/tmp/cmp-mgatk"), "positions.txt")
  fs::dir_create(fs::path_dir(posfile))
  data.table::fwrite(
    data.table::data.table(pos = sort(unique(positions))),
    posfile,
    col.names = FALSE
  )
  on.exit(fs::file_delete(posfile), add = TRUE)

  cmd <- glue::glue(
    "zcat {file} | ",
    "awk -F, 'NR==FNR{{keep[$1]=1; next}} ($1 in keep)' {posfile} -"
  )
  data.table::fread(
    cmd = cmd,
    header = FALSE,
    col.names = c(
      "position",
      "barcode",
      "depth"
    )
  )
}

# cut emits fields in file order regardless of the order requested, so the
# column names are taken from the header in that same order.
fn_read_af_subset <- function(file, variants) {
  con <- gzfile(file, "r")
  hdr <- strsplit(readLines(con, n = 1L), "\t")[[1]]
  close(con)

  idx <- sort(unique(match(variants, hdr)))
  stopifnot(!anyNA(idx))

  data.table::fread(
    cmd = glue::glue(
      "zcat {file} | cut -f{paste(c(1L, idx), collapse = ',')}"
    ),
    skip = 1L,
    header = FALSE,
    col.names = c("barcode", hdr[idx])
  )
}

# body ---------------------------------------------------------------------
stats_scmocha <- fn_read_stats(inputs$stats_scmocha, "scMOCHA")
stats_mgatk <- fn_read_stats(inputs$stats_mgatk, "Original mgatk")

log_info(
  "variant_stats rows: scMOCHA {nrow(stats_scmocha)}, ",
  "original mgatk {nrow(stats_mgatk)}"
)

# The two tables share position/nucleotide but every statistic differs,
# because the cell set and the confident-detection rule differ.
stats_joined <- merge(
  stats_scmocha[, .(
    variant,
    position,
    nucleotide,
    vmr_scmocha = vmr,
    mean_scmocha = mean,
    variance_scmocha = variance,
    nconf_scmocha = n_cells_conf_detected,
    n5_scmocha = n_cells_over_5,
    n10_scmocha = n_cells_over_10,
    n20_scmocha = n_cells_over_20,
    n95_scmocha = n_cells_over_95,
    maxhet_scmocha = max_heteroplasmy,
    strand_scmocha = strand_correlation,
    cov_scmocha = mean_coverage
  )],
  stats_mgatk[, .(
    variant,
    position,
    nucleotide,
    vmr_mgatk = vmr,
    mean_mgatk = mean,
    variance_mgatk = variance,
    nconf_mgatk = n_cells_conf_detected,
    n5_mgatk = n_cells_over_5,
    n10_mgatk = n_cells_over_10,
    n20_mgatk = n_cells_over_20,
    n95_mgatk = n_cells_over_95,
    maxhet_mgatk = max_heteroplasmy,
    strand_mgatk = strand_correlation,
    cov_mgatk = mean_coverage
  )],
  by = c("variant", "position", "nucleotide"),
  all = TRUE
)

stats_joined[, `:=`(
  candidate_scmocha = !is.na(nconf_scmocha),
  candidate_mgatk = !is.na(nconf_mgatk)
)]

# S1: the n_cells_conf_detected >= 3 gate, identical in both callers but
# computed from different confident-detection rules.
stats_joined[, `:=`(
  s1_scmocha = !is.na(nconf_scmocha) & nconf_scmocha >= CUTOFF_NCELLS_CONF,
  s1_mgatk = !is.na(nconf_mgatk) & nconf_mgatk >= CUTOFF_NCELLS_CONF
)]

# S2 for original mgatk only. scMOCHA's S2 needs per-cell depth and is built
# in step 03; it is not a function of these columns.
stats_joined[,
  s2_mgatk := s1_mgatk &
    !is.na(vmr_mgatk) &
    !is.na(strand_mgatk) &
    vmr_mgatk > CUTOFF_VMR_MGATK &
    strand_mgatk > CUTOFF_STRAND_MGATK
]

stats_joined[, blacklisted := fn_is_blacklisted(position)]

log_info(
  "candidates: scMOCHA {sum(stats_joined$candidate_scmocha)}, ",
  "mgatk {sum(stats_joined$candidate_mgatk)}"
)
log_info(
  "S1 (n_cells_conf >= {CUTOFF_NCELLS_CONF}): scMOCHA ",
  "{sum(stats_joined$s1_scmocha)}, mgatk {sum(stats_joined$s1_mgatk)}"
)
log_info("S2 mgatk (vmr + strand gate): {sum(stats_joined$s2_mgatk)}")

# Per-cell mean MT coverage, and the cell filter original mgatk applies.
cell_depth <- data.table::fread(
  inputs$depth,
  header = FALSE,
  col.names = c("barcode", "mean_coverage")
)
cell_depth[, kept_by_mgatk := mean_coverage > CUTOFF_CELL_MEANCOV]

barcodes_scmocha <- fn_read_af_barcodes(inputs$af_scmocha)
barcodes_mgatk <- fn_read_af_barcodes(inputs$af_mgatk)

# Barcode order differs between depthTable and the AF matrices, so these are
# set comparisons, not row-count comparisons.
stopifnot(
  setequal(cell_depth$barcode, barcodes_scmocha),
  setequal(cell_depth[kept_by_mgatk == TRUE, barcode], barcodes_mgatk)
)

log_info(
  "cells: {nrow(cell_depth)} total, {sum(cell_depth$kept_by_mgatk)} kept by ",
  "mgatk, {sum(!cell_depth$kept_by_mgatk)} dropped"
)

# Only S1 variants of either rule set are ever needed downstream.
variants_s1 <- stats_joined[s1_scmocha | s1_mgatk, variant]
positions_needed <- stats_joined[s1_scmocha | s1_mgatk, unique(position)]
log_info(
  "S1 union: {length(variants_s1)} variants at ",
  "{length(positions_needed)} positions"
)

cell_pos_coverage <- fn_read_coverage_at(inputs$coverage, positions_needed)
log_info("coverage rows retained: {nrow(cell_pos_coverage)}")

stopifnot(
  data.table::uniqueN(cell_pos_coverage$position) == length(positions_needed)
)

# The raw matrix is used, not the post-C3 one, so both rule sets' S1 variants
# are present for every cell.
cell_af <- fn_read_af_subset(inputs$af_raw, variants_s1)
stopifnot(
  setequal(cell_af$barcode, cell_depth$barcode),
  ncol(cell_af) == length(variants_s1) + 1L
)
log_info(
  "cell AF submatrix: {nrow(cell_af)} cells x {ncol(cell_af) - 1L} variants"
)

# save ---------------------------------------------------------------------
export(
  stats_joined,
  as.character(fs::path(
    paths$cachedir,
    "01-variant-joined.qs"
  ))
)
export(cell_depth, as.character(fs::path(paths$cachedir, "01-cell-depth.qs")))
export(
  cell_pos_coverage,
  as.character(fs::path(
    paths$cachedir,
    "01-cell-pos-coverage.qs"
  ))
)
export(cell_af, as.character(fs::path(paths$cachedir, "01-cell-af-s1.qs")))

log_info(
  "saved 01-variant-joined.qs ({nrow(stats_joined)} rows), ",
  "01-cell-depth.qs ({nrow(cell_depth)} rows), ",
  "01-cell-pos-coverage.qs ({nrow(cell_pos_coverage)} rows), ",
  "01-cell-af-s1.qs ({nrow(cell_af)} x {ncol(cell_af) - 1L}) to {paths$cachedir}"
)
