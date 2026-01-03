#!/usr/bin/env Rscript
# =============================================================================
# Metainfo
# =============================================================================
# @AUTHOR: Chun-Jie Liu
# @DATE: 2026-01-01
# @DESCRIPTION:
#   Parallel cell-type-specific variant plotting
#   - JOY + HIST + CUMFRAC -> one PDF
#   - DETAIL -> separate PDF
#   - Skip if PDF exists
# =============================================================================

load_pkg(jutils)

library(future)
library(furrr)

# -----------------------------------------------------------------------------
# Args
# -----------------------------------------------------------------------------
GetoptLong.options(help_style = "two-column")
VERSION <- "v0.2.0"

verbose <- TRUE
GetoptLong("verbose!", "print verbose log messages")

# -----------------------------------------------------------------------------
# Logger
# -----------------------------------------------------------------------------
logger::log_layout(logger::layout_glue_colors)
logger::log_threshold(if (verbose) logger::TRACE else logger::INFO)

logger::log_info("Starting script (version: {VERSION})")

# -----------------------------------------------------------------------------
# Parallel plan (REAL parallelism)
# -----------------------------------------------------------------------------
workers <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", 8))
future::plan(multisession, workers = workers)

logger::log_info("Parallel workers: {workers}")

# -----------------------------------------------------------------------------
# Paths
# -----------------------------------------------------------------------------
outdir <- path("/home/liuc9/github/scMOCHA-data/analysis/zzz/MANUSCRIPTFIGURES")
outdirnotuse <- path(
  "/home/liuc9/github/scMOCHA-data/analysis/high-res/MANUSCRIPTFIGURES-notuse"
)

plot_dir <- outdirnotuse / "celltype-specific-each"
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

# -----------------------------------------------------------------------------
# Load data
# -----------------------------------------------------------------------------
ALLVARIANTS_TEST <- import(
  outdir / "VARIANT-KRUSKAL-WALLIS-TEST.xlsx"
)

ALLVARIANTS_TEST_SIG <- ALLVARIANTS_TEST |>
  dplyr::select(variant, gseid, srrid) |>
  dplyr::distinct()

logger::log_info("Total tasks: {nrow(ALLVARIANTS_TEST_SIG)}")

# -----------------------------------------------------------------------------
# Source plotting functions
# -----------------------------------------------------------------------------
source(
  "/home/liuc9/github/scMOCHA-data/analysis/high-res/plot_celltype_specific_variant.R"
)
source(
  "/home/liuc9/github/scMOCHA-data/analysis/high-res/plot_individual_proportion.R"
)

# -----------------------------------------------------------------------------
# Helper
# -----------------------------------------------------------------------------
safe_pdf <- function(file, width, height, expr) {
  tryCatch(
    {
      pdf(file = file, width = width, height = height)
      force(expr)
    },
    error = function(e) {
      logger::log_error("FAILED: {file}")
      logger::log_error(conditionMessage(e))
    },
    finally = {
      if (grDevices::dev.cur() > 1) dev.off()
    }
  )
}

# =============================================================================
# PARALLEL plotting
# =============================================================================

furrr::future_pwalk(
  .l = list(
    ALLVARIANTS_TEST_SIG$variant,
    ALLVARIANTS_TEST_SIG$gseid,
    ALLVARIANTS_TEST_SIG$srrid
  ),
  .f = \(thevariant, thegseid, thesrrid) {
    base <- glue::glue("{thevariant}-{thegseid}-{thesrrid}")

    main_pdf <- plot_dir / glue::glue("JOY-HIST-CUMFRAC-{base}.pdf")
    detail_pdf <- plot_dir / glue::glue("DETAIL-{base}.pdf")

    # JOY + HIST + CUMFRAC
    if (!file.exists(main_pdf)) {
      safe_pdf(
        main_pdf,
        11,
        6,
        {
          print(fn_plot_joy(thevariant, thegseid, thesrrid))
          print(fn_plot_hist(thevariant, thegseid, thesrrid))
          print(fn_plot_cumulative_fraction(thevariant, thegseid, thesrrid))
        }
      )
    }

    # DETAIL
    if (!file.exists(detail_pdf)) {
      safe_pdf(
        detail_pdf,
        20,
        12,
        {
          print(fn_plot_joy_celltype_detail(thevariant, thegseid, thesrrid))
        }
      )
    }
  },
  .options = furrr::furrr_options(seed = TRUE)
)

logger::log_success("Parallel plotting finished")
