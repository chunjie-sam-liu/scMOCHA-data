#!/usr/bin/env Rscript
# =============================================================================
# Metainfo
# =============================================================================
# @AUTHOR: Chun-Jie Liu
# @DATE: 2026-01-01
# @DESCRIPTION:
#   Parallel cell-type-specific variant plotting with file locking
#   - JOY + HIST + CUMFRAC in one PDF
#   - DETAIL in a separate PDF
#   - Skip if PDF exists
#   - Parallel-safe (filelock)
# =============================================================================

# -----------------------------------------------------------------------------
# Library
# -----------------------------------------------------------------------------
load_pkg(jutils, future, furrr, filelock)


# -----------------------------------------------------------------------------
# Args
# -----------------------------------------------------------------------------
GetoptLong.options(help_style = "two-column")
VERSION <- "v0.3.0"

verbose <- TRUE
GetoptLong("verbose!", "print verbose log messages")

# -----------------------------------------------------------------------------
# Logger
# -----------------------------------------------------------------------------
logger::log_layout(logger::layout_simple)
logger::log_threshold(if (verbose) logger::TRACE else logger::INFO)

logger::log_info("Starting script (version: {VERSION})")

# -----------------------------------------------------------------------------
# Parallel plan (REAL parallelism)
# -----------------------------------------------------------------------------
dotenv(".env")
workers <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", 8))

# Use multicore (fork) on Linux for shared environment - much faster
# multisession spawns separate R processes that don't share globals
if (.Platform$OS.type == "unix") {
  future::plan(future::multicore, workers = workers)
} else {
  future::plan(future::multisession, workers = workers)
}

logger::log_info(
  "Parallel workers: {workers}, plan: {class(future::plan())[1]}"
)

# -----------------------------------------------------------------------------
# Paths
# -----------------------------------------------------------------------------
dotenv(".env")
outdir <- path(Sys.getenv("OUTDIR"))

outdirnotuse <- path(Sys.getenv("OUTDIRNOTUSE"))

plot_dir <- outdirnotuse / "celltype-specific-each"
unlink(plot_dir, recursive = TRUE)
fs::dir_create(plot_dir)

cli::cli_alert_info(
  "Plot output dir: {plot_dir}, file exists {file_exists(plot_dir)}"
)

# -----------------------------------------------------------------------------
# Load data
# -----------------------------------------------------------------------------
ALLVARIANTS_TEST <- import(
  outdir / "VARIANT-KRUSKAL-WALLIS-TEST.xlsx"
)

ALLVARIANTS_TEST_SIG <- ALLVARIANTS_TEST |>
  dplyr::filter(p.value < 0.05) |>
  dplyr::select(variant, gseid, srrid) |>
  dplyr::distinct()

logger::log_info("Total tasks: {nrow(ALLVARIANTS_TEST_SIG)}")

# -----------------------------------------------------------------------------
# Source plotting functions
# -----------------------------------------------------------------------------
source(
  path(
    Sys.getenv("HIGHRESDIR"),
    "plot_celltype_specific_variant.R"
  )
)
source(
  path(
    Sys.getenv("HIGHRESDIR"),
    "plot_individual_proportion.R"
  )
)

# =============================================================================
# Helper: safe PDF plotting with file lock
# =============================================================================
safe_pdf_with_lock <- function(file, width, height, plot_expr) {
  lockfile <- paste0(file, ".lock")

  lock <- tryCatch(
    filelock::lock(lockfile, timeout = 0),
    error = function(e) NULL
  )

  # Another worker is plotting
  if (is.null(lock)) {
    logger::log_trace("LOCKED (skip): {file}")
    return(invisible(NULL))
  }

  # Double-check existence after acquiring lock
  if (file.exists(file)) {
    logger::log_trace("EXISTS (skip): {file}")
    filelock::unlock(lock)
    unlink(lockfile)
    return(invisible(NULL))
  }

  tryCatch(
    {
      logger::log_info("Plotting: {file}")
      pdf(file = file, width = width, height = height)
      force(plot_expr)
    },
    error = function(e) {
      logger::log_error(
        "FAILED: {file} - {skip_formatter(conditionMessage(e))}"
      )
    },
    finally = {
      # Always close device if open
      if (grDevices::dev.cur() > 1) {
        dev.off()
      }

      # Always unlock
      try(filelock::unlock(lock), silent = TRUE)

      # Explicitly remove lock file
      if (file.exists(lockfile)) {
        unlink(lockfile)
      }
    }
  )
}

# =============================================================================
# Parallel plotting
# =============================================================================
furrr::future_pwalk(
  .l = list(
    ALLVARIANTS_TEST_SIG$variant,
    ALLVARIANTS_TEST_SIG$gseid,
    ALLVARIANTS_TEST_SIG$srrid
  ),
  .f = \(thevariant, thegseid, thesrrid) {
    # thevariant <- ALLVARIANTS_TEST_SIG$variant[[1]]
    # thegseid <- ALLVARIANTS_TEST_SIG$gseid[[1]]
    # thesrrid <- ALLVARIANTS_TEST_SIG$srrid[[1]]
    base <- glue::glue("{thevariant}-{thegseid}-{thesrrid}")

    main_pdf <- plot_dir / glue::glue("JOY-HIST-CUMFRAC-{base}.pdf")
    detail_pdf <- plot_dir / glue::glue("DETAIL-{base}.pdf")
    pseudo_bulk_pdf <- plot_dir / glue::glue("PSEUDO-BULK-{base}.pdf")

    # -------------------------------------------------------------------------
    # JOY + HIST + CUMFRAC (same PDF)
    # -------------------------------------------------------------------------
    safe_pdf_with_lock(
      file = main_pdf,
      width = 11,
      height = 6,
      {
        print(fn_plot_joy(thevariant, thegseid, thesrrid))
        print(fn_plot_hist(thevariant, thegseid, thesrrid))
        print(fn_plot_cumulative_fraction(thevariant, thegseid, thesrrid))
      }
    )

    # -------------------------------------------------------------------------
    # DETAIL plot (separate PDF)
    # -------------------------------------------------------------------------
    safe_pdf_with_lock(
      file = detail_pdf,
      width = 20,
      height = 12,
      {
        print(fn_plot_joy_celltype_detail(
          thevariant,
          thegseid,
          thesrrid
        ))
      }
    )

    # -------------------------------------------------------------------------
    # Pseudo bulk plot (separate PDF)
    # -------------------------------------------------------------------------

    safe_pdf_with_lock(
      file = pseudo_bulk_pdf,
      width = 15,
      height = 8,
      {
        print(fn_plot_variant_ratio(thevariant))
        print(fn_plot_hetero_pseudo_bulk(thevariant))
      }
    )
  },
  .options = furrr::furrr_options(seed = TRUE)
)

future::plan(future::sequential)
# -----------------------------------------------------------------------------
# Footer
# -----------------------------------------------------------------------------
logger::log_success("All parallel plotting tasks finished")
