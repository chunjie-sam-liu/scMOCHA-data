#!/usr/bin/env Rscript
# Metainfo ----------------------------------------------------------------
# @AUTHOR: Chun-Jie Liu
# @CONTACT: chunjie.sam.liu.at.gmail.com
# @DATE: 2025-09-02 12:01:16
# @DESCRIPTION: filename
# @VERSION: v0.0.1

# Library -----------------------------------------------------------------

suppressPackageStartupMessages(library(magrittr))
library(ggplot2)
library(patchwork)
library(prismatic)
library(paletteer)
library(data.table)
# library(rlang)
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
gseid_srrid_variant_hetero_plot_ratio <- import(
  "/home/liuc9/github/scMOCHA-data/analysis/zzz/plot-real-somatic-variant/main-variants/gseid_srrid_variant_hetero_plot_ratio.qs"
)
# load conn ---------------------------------------------------------------

# src ---------------------------------------------------------------------

# function ----------------------------------------------------------------
fn_de <- function(
  thegseid,
  thesrrid,
  thevariant,
  forplot_
) {
  library(Seurat)
  .dir <- file.path(
    "/home/liuc9/github/scMOCHA-data/data/",
    thegseid,
    "final",
    thesrrid
  )
  sc <- import(
    file.path(
      .dir,
      "sc_azimuth.rds.gz"
    )
  )
  sc_azimuth <- sc$sc_azimuth
  rm(sc)
  gc()

  sc_azimuth@meta.data |>
    tibble::rownames_to_column("barcode") |>
    as.data.table() |>
    dplyr::left_join(
      forplot_ |>
        as.data.table() |>
        dplyr::mutate(barcode = as.character(barcode)),
    ) |>
    dplyr::mutate(
      barcode_new = glue::glue("{thegseid}-{thesrrid}-{barcode}")
    ) |>
    as.data.frame() -> d_merge

  new_names <- setNames(d_merge$barcode_new, d_merge$barcode)

  sc_azimuth <- RenameCells(
    sc_azimuth,
    new.names = new_names
  )

  sc_azimuth@meta.data <- d_merge |>
    tibble::column_to_rownames("barcode_new")

  sc_azimuth <- Seurat::SCTransform(
    sc_azimuth,
    assay = "RNA",
  )
  DefaultAssay(sc_azimuth) <- "SCT"

  .dir_de <- file.path(
    .dir,
    "de"
  )
  if (dir.exists(.dir_de) == FALSE) {
    dir.create(.dir_de)
  }

  export(
    sc_azimuth,
    file = file.path(.dir_de, "sc_azimuth.sct.qs")
  )

  markers <- Seurat::FindMarkers(
    object = sc_azimuth,
    ident.1 = "Heteroplasmy",
    ident.2 = "Sufficient reads",
    test.use = "wilcox",
    group.by = "cellvarianttype"
  )

  export(
    markers,
    file = file.path(
      .dir_de,
      "sc_azimuth.markers.hetero_vs_sufficient.qs"
    ),
  )

  rm(sc_azimuth)
  gc()

  markers
}

fn_de_plot <- function(
  markers,
  .cutoff_pval = 0.05,
  .cutoff_log2fc = 1
) {
  markers |>
    tibble::rownames_to_column("gene") |>
    dplyr::mutate(
      fdr = -log10(p_val_adj)
    ) |>
    dplyr::mutate(
      avg_log2FC = ifelse(
        abs(avg_log2FC) > 100,
        sign(avg_log2FC) * 100,
        avg_log2FC
      )
    ) |>
    dplyr::mutate(
      color = dplyr::case_when(
        p_val_adj < .cutoff_pval & avg_log2FC > .cutoff_log2fc ~ "red",
        p_val_adj < .cutoff_pval & avg_log2FC < -.cutoff_log2fc ~ "blue",
        TRUE ~ "grey"
      )
    ) -> forplot

  forplot |>
    ggplot(aes(
      x = avg_log2FC,
      y = fdr,
      color = color
    )) +
    geom_point(aes()) +
    ggrepel::geom_text_repel(
      data = forplot |>
        dplyr::filter(
          (p_val_adj < .cutoff_pval & avg_log2FC > .cutoff_log2fc) |
            (p_val_adj < .cutoff_pval & avg_log2FC < -.cutoff_log2fc)
        ) |>
        dplyr::slice(1:10),
      aes(label = gene),
      size = 3,
      max.overlaps = 20
    ) +
    scale_color_identity() +
    geom_vline(xintercept = c(-1, 1), linetype = "dashed") +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
    theme_classic() +
    labs(
      x = "Fold change(Heteroplasmy/Sufficient reads)",
      y = "FDR",
      # title = "{thegseid}-{thesrrid}-m.{thevariant}" |> glue::glue()
    ) +
    theme(
      plot.title = element_text(
        hjust = 0.5,
        face = "bold",
        color = "black",
        size = 16
      )
    )
}

# body --------------------------------------------------------------------

thevariant <- "3727T>C"
thesrrid <- ""

gseid_srrid_variant_hetero_plot_ratio |>
  # dplyr::filter(variant == thevariant) |>
  dplyr::filter(variant %in% c("3727T>C", "3728C>T")) |>
  # dplyr::select(tidyselect::contains("ratio"))
  dplyr::select(
    gseid,
    srrid,
    variant,
    forplot,
    # tidyselect::contains("ratio")
  ) -> filtered_data


thevariant <- "3728C>T"
thesrrid <- "GSM7080053"
thegseid <- "GSE226602"
forplot_ <- filtered_data$forplot[[1]]

# filtered_data |>
#   dplyr::filter(
#     variant == thevariant,
#     srrid == thesrrid
#   ) |>
#   dplyr::select(forplot) |>
#   tidyr::unnest(cols = c(forplot)) -> variant_cell_barcode

#
#
# ? plot --------------------------------------------------------------------
#
#

filtered_data |>
  head(10) |>
  dplyr::mutate(
    p = parallel::mcmapply(
      FUN = \(thegseid, thesrrid, thevariant, forplot_) {
        tryCatch(
          expr = {
            fn_de_plot(
              markers = fn_de(
                thegseid = thegseid,
                thesrrid = thesrrid,
                thevariant = thevariant,
                forplot_ = forplot_
              )
            ) +
              ggtitle(glue::glue("{thegseid}-{thesrrid}-m.{thevariant}"))
          },
          error = \(e) {
            message(glue::glue("{thegseid}-{thesrrid}-m.{thevariant} error"))
            return(NULL)
          }
        )
      },
      thegseid = gseid,
      thesrrid = srrid,
      thevariant = variant,
      forplot_ = forplot,
      SIMPLIFY = FALSE,
      USE.NAMES = FALSE,
      mc.cores = 8
    )
  ) |>
  dplyr::select(-forplot) -> filtered_data_plots

# footer------------------------------------------------------------------

# future: :plan(future: :sequential)

# save image --------------------------------------------------------------
