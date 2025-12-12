#!/usr/bin/env Rscript
# Metainfo ----------------------------------------------------------------
# @AUTHOR: Chun-Jie Liu
# @CONTACT: chunjie.sam.liu.at.gmail.com
# @DATE: 2024-12-12 10:00:00
# @DESCRIPTION: Template for uppercase utility script
# @VERSION: v0.0.1

# Library -----------------------------------------------------------------
load_pkg(jutils)
load_pkg(
  data.table,
  dplyr,
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
input_file = ""
output_file = ""
basedir = "/liulab/chunjie/data/scMOCHA"

GetoptLong(
  "input_file=s",
  "input file path (required)",
  "output_file=s",
  "output file path (required)",
  "basedir=s",
  "base directory path",
  "verbose!",
  "print detailed messages"
)

# Validate required arguments
if (input_file == "" || output_file == "") {
  stop("Both input_file and output_file are required")
}

# header ------------------------------------------------------------------
log_threshold(TRACE)
log_layout(layout_glue_colors)

# log_info("Processing {input_file}...")

# function ----------------------------------------------------------------

# Utility function
fn_transform_data <- function(.data) {
  # Data transformation logic
  .data |>
    dplyr::mutate(
      # Add transformations here
      processed = TRUE
    )
}

# load data ---------------------------------------------------------------

# Check input file exists
if (!fs::file_exists(input_file)) {
  log_error("Input file not found: {input_file}")
  stop("Missing required input file")
}

# Load input data
input_data <- tryCatch(
  {
    import(input_file)
  },
  error = function(e) {
    log_error("Failed to read input file: {e$message}")
    stop(e)
  }
)

# log_info("Loaded {nrow(input_data)} rows from input file")

# body --------------------------------------------------------------------

# Process data
processed_data <- input_data |>
  fn_transform_data()

# log_success("Processing complete: {nrow(processed_data)} rows")

# save --------------------------------------------------------------------

# Ensure output directory exists
output_dir <- fs::path_dir(output_file)
fs::dir_create(output_dir, recurse = TRUE)

# Export processed data
tryCatch(
  {
    export(processed_data, output_file)
    # log_success("Output saved to {output_file}")
  },
  error = function(e) {
    log_error("Failed to write output file: {e$message}")
    stop(e)
  }
)
