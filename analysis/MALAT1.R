#!/usr/bin/env Rscript
# Metainfo ----------------------------------------------------------------
# @AUTHOR: Chun-Jie Liu
# @CONTACT: chunjie.sam.liu.at.gmail.com
# @DATE: 2025-08-15 12:53:23
# @DESCRIPTION: filename
# @VERSION: v0.0.1

# Library -----------------------------------------------------------------

suppressPackageStartupMessages(library(magrittr))
library(ggplot2)
library(patchwork)
library(prismatic)
library(paletteer)
library(data.table)
#library(rlang)
library(glue)
library(parallel)
library(GetoptLong)
library(logger)

# args --------------------------------------------------------------------

# s: string, i: integer, f: float, !: boolean
# @: array
# %: hash
# default: default value specified here.
verbose <- FALSE
spec <- "
Usage: Rscript foorbar.R [options]
Options:

<verbose!> Print messages
"

GetoptLong.options(help_style = "two-column")
GetoptLong(spec, template_control = list(opt_width = 21))

# header ------------------------------------------------------------------

# future: :plan(future: :multisession, workers = 10)

# load data ---------------------------------------------------------------
gseid <- "GSE226602"
srrid <- "GSM7080049"
# thepath <- "/home/liuc9/github/scMOCHA-data/data/GSE226602/final/GSM7080049"
srrdir <- "/home/liuc9/github/scMOCHA-data/analysis/zzz/clean-data/gse_srrid_srrdir.csv" |>
  import()


library(Seurat)
sc <- readr::read_rds(
  file.path(
    thepath,
    "sc_azimuth.rds.gz"
  )
)
# load conn ---------------------------------------------------------------

# src ---------------------------------------------------------------------

# function ----------------------------------------------------------------
fn_plot <- function(gseid, srrid, sc) {
  sc$sct <- Seurat::SCTransform(
    object = sc$sc,
    vars.to.regress = c("percent.mt", "percent.ribo")
  )

  LayerData(sc$sct[["SCT"]], layer = "data")["MALAT1", ] |>
    tibble::enframe(
      name = "barcode",
      value = "MALAT1"
    ) -> malat1_expr

  tibble::tibble(
    barcode = colnames(sc$sc)
  ) |>
    dplyr::mutate(
      keep = barcode %in% colnames(sc$sc_filter$RNA),
    ) |>
    dplyr::mutate(
      kept = ifelse(keep, "kept", "filtered")
    ) -> cell_kept

  cell_kept |>
    dplyr::left_join(
      malat1_expr,
      by = "barcode"
    ) |>
    ggplot(aes(
      x = MALAT1
    )) +
    geom_density(
      aes(color = kept)
    ) +
    theme_minimal() +
    theme(
      legend.position = c(0.3, 0.8)
    ) +
    labs(
      x = "MALAT1 Expression",
      y = "Density",
      title = "{gseid}-{srrid}" |> glue::glue()
    )
}
# body --------------------------------------------------------------------
srrdir |>
  # head(20) |>
  dplyr::mutate(
    p = parallel::mcmapply(
      FUN = \(gseid, srrid, srrdir) {
        sc <- readr::read_rds(
          file.path(
            srrdir,
            "sc_azimuth.rds.gz"
          )
        )
        fn_plot(gseid, srrid, sc) -> p
        ggsave(
          plot = p,
          filename = file.path(
            srrdir,
            "MALAT1_expression_density.pdf"
          ),
          device = "pdf",
          width = 10,
          height = 6,
        )
        p
      },
      gseid = gseid,
      srrid = srrid,
      srrdir = srrdir,
      mc.cores = 20,
      SIMPLIFY = FALSE
    )
  ) -> srrdir_plot
# footer ------------------------------------------------------------------

# future: :plan(future: :sequential)

# save image --------------------------------------------------------------
