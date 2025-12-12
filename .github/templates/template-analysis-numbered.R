#!/usr/bin/env Rscript
# Metainfo ----------------------------------------------------------------
# @AUTHOR: Chun-Jie Liu
# @CONTACT: chunjie.sam.liu.at.gmail.com
# @DATE: 2024-12-12 10:00:00
# @DESCRIPTION: Template for numbered analysis script
# @VERSION: v0.0.1

# Library -----------------------------------------------------------------
load_pkg(jutils)
load_pkg(
  ggplot2,
  patchwork,
  prismatic,
  paletteer,
  data.table,
  dplyr,
  tidyr,
  glue,
  fs,
  GetoptLong,
  logger,
  qs,
  DBI,
  duckdb
)

# args --------------------------------------------------------------------
# s: string, i: integer, f: float, !: boolean
# @: array, %: hash
GetoptLong.options(help_style = "two-column")
VERSION = "v0.0.1"

verbose = FALSE
gseid = "GSE123456"
basedir = "/liulab/chunjie/data/scMOCHA"
ncores = 20

GetoptLong(
  "gseid=s",
  "GSE accession ID",
  "basedir=s",
  "base directory path",
  "ncores=i",
  "number of cores",
  "verbose!",
  "print detailed messages"
)

# src ---------------------------------------------------------------------
source("00-colors.R") # Load project color schemes

# header ------------------------------------------------------------------
log_threshold(TRACE)
log_layout(layout_glue_colors)

# log_info('Starting analysis...')
# log_debug('This is a debug message')
# log_success('Analysis complete!')
# log_warn('Warning message')
# log_error('Error message')

# function ----------------------------------------------------------------

# Custom function example
fn_process_data <- function(.data, threshold = 0.05) {
  # Function implementation
  .data |>
    dplyr::filter(value > threshold) |>
    dplyr::mutate(
      log_value = log10(value + 1)
    )
}

# load data ---------------------------------------------------------------

# Set paths
datadir <- fs::path(basedir, "data")
outdir <- fs::path(basedir, "analysis/zzz/MANUSCRIPTFIGURES")

# Ensure output directory exists
fs::dir_create(outdir, recurse = TRUE)

# Load input data (adjust filename as needed)
# input_data <- import(fs::path(datadir, "input.qs"))

# Alternative: DuckDB
# conn <- DBI::dbConnect(duckdb::duckdb(), fs::path(datadir, "database.duckdb"))
# input_data <- dplyr::tbl(conn, "table_name") |> as.data.table()
# DBI::dbDisconnect(conn, shutdown = TRUE)

# body --------------------------------------------------------------------

# Example analysis pipeline:
# processed_data <- input_data |>
#   dplyr::filter(!is.na(variant)) |>
#   dplyr::mutate(
#     vaf_pct = vaf * 100,
#     category = case_when(
#       vaf < 0.01 ~ "low",
#       vaf < 0.1 ~ "medium",
#       TRUE ~ "high"
#     )
#   ) |>
#   tidyr::nest(.by = c(gseid, celltype), .key = "nested") |>
#   dplyr::mutate(
#     stats = purrr::map(.x = nested, .f = \(.x) {
#       # Process each nested group
#       fn_process_data(.x)
#     })
#   ) |>
#   tidyr::unnest(stats)

# Example visualization:
# p <- ggplot(processed_data, aes(x = celltype, y = vaf_pct, fill = disease)) +
#   geom_boxplot() +
#   scale_fill_manual(values = color_disease) +
#   labs(
#     title = "Variant Allele Frequency by Cell Type and Disease",
#     x = "Cell Type",
#     y = "VAF (%)",
#     fill = "Disease"
#   ) +
#   theme_classic() +
#   theme(
#     axis.text.x = element_text(angle = 45, hjust = 1)
#   )

# Save plot
# ggsave(
#   filename = fs::path(outdir, "NN-analysis-name.pdf"),
#   plot = p,
#   width = 10,
#   height = 6
# )

# save --------------------------------------------------------------------

# Export processed data
# export(processed_data, fs::path(datadir, "processed_output.qs"))

# Alternative: Save to DuckDB
# conn <- DBI::dbConnect(duckdb::duckdb(), fs::path(datadir, "database.duckdb"))
# DBI::dbWriteTable(conn, "processed_table", processed_data, overwrite = TRUE)
# DBI::dbDisconnect(conn, shutdown = TRUE)

# log_success("Analysis complete! Output saved to {outdir}")
