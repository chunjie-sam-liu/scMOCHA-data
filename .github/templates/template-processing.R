#!/usr/bin/env Rscript
# Metainfo ----------------------------------------------------------------
# @AUTHOR: Chun-Jie Liu
# @CONTACT: chunjie.sam.liu.at.gmail.com
# @DATE: 2024-12-12 10:00:00
# @DESCRIPTION: Template for preprocessing/processing script
# @VERSION: v0.0.1

# Library -----------------------------------------------------------------
load_pkg(jutils)
load_pkg(
  data.table,
  dplyr,
  tidyr,
  glue,
  fs,
  parallel,
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
  "number of cores for parallel processing",
  "verbose!",
  "print detailed messages"
)

# header ------------------------------------------------------------------
log_threshold(TRACE)
log_layout(layout_glue_colors)

# log_info("Starting processing for {gseid}...")

# function ----------------------------------------------------------------

# Processing function for parallel execution
fn_process_sample <- function(.srrid, .gseid, .basedir) {
  # Sample-level processing logic

  # Construct paths
  input_path <- fs::path(.basedir, "data", .gseid, .srrid, "input.qs")
  output_path <- fs::path(.basedir, "data", .gseid, .srrid, "output.qs")

  # Check if input exists
  if (!fs::file_exists(input_path)) {
    # log_warn("Input file not found for {.srrid}")
    return(NULL)
  }

  # Process data
  result <- tryCatch(
    {
      data <- import(input_path)

      # Apply transformations
      processed <- data |>
        dplyr::filter(!is.na(value)) |>
        dplyr::mutate(
          processed_value = value * 2
        )

      # Save result
      export(processed, output_path)

      # log_success("Processed {.srrid}: {nrow(processed)} rows")
      return(list(srrid = .srrid, status = "success", nrow = nrow(processed)))
    },
    error = function(e) {
      # log_error("Failed to process {.srrid}: {e$message}")
      return(list(srrid = .srrid, status = "failed", error = e$message))
    }
  )

  return(result)
}

# load data ---------------------------------------------------------------

# Set paths
datadir <- fs::path(basedir, "data", gseid)
metafile <- fs::path(datadir, "metadata.qs")

# Load metadata
if (!fs::file_exists(metafile)) {
  log_error("Metadata file not found: {metafile}")
  stop("Missing required metadata file")
}

metadata <- import(metafile)
# log_info("Loaded metadata: {nrow(metadata)} samples")

# body --------------------------------------------------------------------

# Get list of samples to process
sample_list <- metadata$srrid

# Parallel processing
# log_info("Processing {length(sample_list)} samples using {ncores} cores...")

results <- parallel::mclapply(
  X = sample_list,
  FUN = function(.srrid) {
    fn_process_sample(.srrid, gseid, basedir)
  },
  mc.cores = ncores
)

# Summarize results
results_summary <- data.table::rbindlist(results, fill = TRUE)

# log_info("Processing summary:")
# log_info("  Success: {sum(results_summary$status == 'success', na.rm = TRUE)}")
# log_info("  Failed: {sum(results_summary$status == 'failed', na.rm = TRUE)}")

# save --------------------------------------------------------------------

# Save processing summary
summary_file <- fs::path(datadir, "processing_summary.qs")
export(results_summary, summary_file)

# log_success("Processing complete! Summary saved to {summary_file}")
