#!/usr/bin/env Rscript
# Metainfo ----------------------------------------------------------------
# @AUTHOR: Chun-Jie Liu
# @CONTACT: chunjie.sam.liu.at.gmail.com
# @DATE: 2026-03-03 14:02:10
# @DESCRIPTION: this script is used for ...

# Reproducibility ----------------------------------------------------------
set.seed(1)
# Library -----------------------------------------------------------------

suppressMessages({
  load_pkg(jutils)
})

# Args --------------------------------------------------------------------

# s: string, i: integer, f: float, !: boolean, @: array, %: list
VERSION = "v0.0.1"

GetoptLong.options(help_style = "two-column")

# default: default value specified here.

nthread = 8
GetoptLong(
  "nthread=i",
  "Number of threads to use",
  "verbose",
  "Enable verbose logging"
)


# Logger ------------------------------------------------------------------

log_layout(layout_glue_colors)

if (isTRUE(verbose)) {
  log_threshold(TRACE)
  log_info("Verbose mode enabled")
} else {
  log_threshold(INFO)
}


# Source ---------------------------------------------------------------------

# Load data ---------------------------------------------------------------
load_pkg(jutils)
dotenv()
suppressMessages({
  conflicted::conflict_prefer("filter", "dplyr")
})
outdir <- path(Sys.getenv("OUTDIR"))
outdirnotuse <- path(Sys.getenv("OUTDIRNOTUSE"))

variants <- c(
  "13592C>T",
  "5031G>T",
  "8362T>G"
)


# Conn ---------------------------------------------------------------

# Function ----------------------------------------------------------------

# Main --------------------------------------------------------------------
library(clusterProfiler)

plotdir <- path(
  "/home/liuc9/github/scMOCHA-data/high-res-MANUSCRIPTFIGURES-notuse/AD/corr/go/plot"
)
fs::dir_create(plotdir)

variants |>
  purrr::walk(\(.variant) {
    log_info("Processing variant: {.variant}")
    .gofile <- outdirnotuse /
      "AD" /
      "corr" /
      "go" /
      glue::glue("ad-variant-{.variant}-go-enrichment.qs")

    if (!fs::file_exists(.gofile)) {
      log_warn("File not found: {.gofile}")
      return(invisible(NULL))
    }

    .go <- import(.gofile)
    .safe_variant <- gsub(">", "_", .variant)
    .outdir <- path(
      plotdir,
      glue::glue("variant-{.safe_variant}")
    )
    dir_create(.outdir)

    # .go is a tibble: celltype + pos_bp/cc/mf + neg_bp/cc/mf + *_plot columns
    .go |>
      dplyr::select(celltype, dplyr::ends_with("_plot")) |>
      tidyr::pivot_longer(
        cols = -celltype,
        names_to = "goname",
        values_to = "p"
      ) |>
      dplyr::filter(!purrr::map_lgl(p, is.null)) |>
      dplyr::mutate(
        a = purrr::pmap(
          list(celltype, goname, p),
          \(.celltype, .goname, .p) {
            .safe_celltype <- gsub(" ", "_", .celltype)
            .filename <- glue::glue(
              "{.safe_variant}_{.safe_celltype}_{.goname}.pdf"
            )
            log_info("Saving: {.filename}")
            ggplot2::ggsave(
              path = .outdir,
              filename = .filename,
              plot = .p,
              width = 10,
              height = 8,
              dpi = 300
            )
          }
        )
      )
  })

# Save  --------------------------------------------------------------

# Session info -------------------------------------------------------------
if (isTRUE(verbose)) {
  sessionInfo()
}
