#!/usr/bin/env Rscript
# Metainfo ----------------------------------------------------------------
# @AUTHOR: Chun-Jie Liu
# @CONTACT: chunjie.sam.liu.at.gmail.com
# @DATE: 2026-03-23 00:00:00
# @DESCRIPTION: Export COVID-19 GO enrichment plots to PDF files

# Reproducibility ----------------------------------------------------------
set.seed(1)

# Library -----------------------------------------------------------------

suppressMessages({
  load_pkg(jutils)
})

# Args --------------------------------------------------------------------

VERSION = "v0.0.1"

GetoptLong.options(help_style = "two-column")

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

# Load data ----------------------------------------------------------------

dotenv()

outdirnotuse <- path(Sys.getenv("OUTDIRNOTUSE"))
covid_outdir <- outdirnotuse / "COVID19"
plotdir <- covid_outdir / "corr" / "go" / "plot"
dir_create(plotdir)

variants <- import(covid_outdir / "COVID19-variant-top-ttest-cluster-variants.qs2")

# Main ---------------------------------------------------------------------

variants |>
  walk(
    \(.variant) {
      log_info("Exporting GO plots for {.variant}")
      gofile <- covid_outdir / "corr" / "go" / glue("covid-variant-{.variant}-go-enrichment.qs2")

      if (!file_exists(gofile)) {
        log_warn("GO file not found: {gofile}")
        return(invisible(NULL))
      }

      go_tbl <- import(gofile)
      safe_variant <- gsub(">", "_", .variant, fixed = TRUE)
      variant_dir <- plotdir / glue("variant-{safe_variant}")
      dir_create(variant_dir)

      go_tbl |>
        select(celltype, ends_with("_plot")) |>
        pivot_longer(
          cols = -celltype,
          names_to = "goname",
          values_to = "plot_obj"
        ) |>
        filter(!purrr::map_lgl(plot_obj, is.null)) |>
        mutate(
          saved = purrr::pmap(
            list(celltype, goname, plot_obj),
            function(celltype, goname, plot_obj) {
              safe_celltype <- gsub(" ", "_", celltype)
              filename <- glue("{safe_variant}_{safe_celltype}_{goname}.pdf")
              ggsave(
                path = variant_dir,
                filename = filename,
                plot = plot_obj,
                width = 10,
                height = 8,
                dpi = 300
              )
            }
          )
        )
    }
  )

if (isTRUE(verbose)) {
  sessionInfo()
}
