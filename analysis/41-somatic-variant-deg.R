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

fn_load_sc_and_sct <- function(.filepath, thegseid, thesrrid, forplot_) {
  .sc <- import(.filepath)
  sc_azimuth <- .sc$sc_azimuth
  rm(.sc)
  gc()
  sc_azimuth@meta.data |>
    tibble::rownames_to_column("barcode") |>
    as.data.table() |>
    dplyr::left_join(
      forplot_ |>
        as.data.table() |>
        dplyr::mutate(
          barcode = as.character(barcode)
        ),
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

  sc_azimuth[["SCT"]]@scale.data <- matrix()
  sc_azimuth
}


fn_sct <- function(
  thegseid,
  thesrrid,
  thevariant,
  forplot_
) {
  library(Seurat)
  .dir <- path(
    "/home/liuc9/github/scMOCHA-data/data/",
    thegseid,
    "final",
    thesrrid
  )
  .dir_de <- path(
    .dir,
    "de"
  )

  dir_create(.dir_de)

  .sct_filepath <- path(
    .dir_de,
    "sc_azimuth.sct.qs"
  )

  sc_azimuth <- if (file_exists(.sct_filepath)) {
    log_fatal("{.sct_filepath} exists, skip!" |> glue::glue())
    return(NULL)
  } else {
    sc_azimuth <- fn_load_sc_and_sct(
      .filepath = path(
        .dir,
        "sc_azimuth.rds.gz"
      ),
      thegseid = thegseid,
      thesrrid = thesrrid,
      forplot_ = forplot_
    )
    export(
      sc_azimuth,
      file = .sct_filepath
    )
    sc_azimuth
  }

  # markers <- Seurat::FindMarkers(
  #   object = sc_azimuth,
  #   ident.1 = "Heteroplasmy",
  #   ident.2 = "Sufficient reads",
  #   test.use = "wilcox",
  #   group.by = "cellvarianttype"
  # )

  # export(
  #   markers,
  #   file = file.path(
  #     .dir_de,
  #     "sc_azimuth.markers.hetero_vs_sufficient.{thevariant}.qs" |>
  #       glue::glue()
  #   ),
  # )

  # rm(sc_azimuth)
  # gc()

  # markers
}


fn_de_high_vs_low <- function(
  thegseid,
  thesrrid,
  thevariant,
  forplot_,
  .vs = c(0.5, 0.5)
) {
  library(Seurat)
  .dir <- file.path(
    "/home/liuc9/github/scMOCHA-data/data/",
    thegseid,
    "final",
    thesrrid
  )
  .dir_de <- file.path(
    .dir,
    "de"
  )

  sc_azimuth <- import(
    file.path(
      .dir_de,
      "sc_azimuth.sct.qs"
    )
  )

  sc_azimuth@meta.data |>
    as.data.table() |>
    dplyr::filter(cellvarianttype == "Heteroplasmy") |>
    dplyr::pull(af) |>
    quantile(probs = seq(0, 1, 0.05), na.rm = FALSE) -> .quant

  .high <- .quant[glue::glue("{.vs[1] * 100}%")]
  .low <- .quant[glue::glue("{.vs[2] * 100}%")]
  # scales::label_number(accuracy = 0.01)(median_af)
  .label_high <- glue::glue(
    "High={scales::label_number(accuracy = 1)(.vs[1] * 100)}% AF={scales::label_number(accuracy = 0.01)(.high)}"
  )
  .label_low <- glue::glue(
    "Low={scales::label_number(accuracy = 1)(.vs[2] * 100)}% AF={scales::label_number(accuracy = 0.01)(.low)}"
  )
  hetero_label <- glue::glue(
    "Heteroplasmy ({.label_high}) vs ({.label_low})"
  )

  sc_azimuth@meta.data |>
    dplyr::mutate(
      cellvarianttype2 = dplyr::case_when(
        cellvarianttype == "Heteroplasmy" &
          af >= .high ~
          glue::glue("{.label_high}"),
        cellvarianttype == "Heteroplasmy" &
          af < .low ~
          glue::glue("{.label_low}"),
        TRUE ~ as.character(cellvarianttype)
      )
    ) -> sc_azimuth@meta.data

  sc_azimuth@meta.data |>
    dplyr::count(cellvarianttype2) |>
    dplyr::arrange(cellvarianttype2) |>
    tibble::deframe() -> n_cellvarianttype2

  markers <- Seurat::FindMarkers(
    object = sc_azimuth,
    ident.1 = glue::glue("{.label_high}"),
    ident.2 = glue::glue("{.label_low}"),
    assay = "SCT",
    slot = "data",
    test.use = "wilcox",
    group.by = "cellvarianttype2"
  )

  export(
    markers,
    file = file.path(
      .dir_de,
      "sc_azimuth.markers.{hetero_label}.{thevariant}.qs" |>
        glue::glue()
    ),
  )

  rm(sc_azimuth)
  gc()

  .labs = labs(
    x = glue::glue(
      "Fold change Heteroplasmy ({.label_high} n={
          scales::label_comma()(n_cellvarianttype2[.label_high])
        }) vs Low ({.label_low} n={
          scales::label_comma()(n_cellvarianttype2[.label_low])
        })"
    ),
    y = "FDR",
    # title = "m.{thevariant}" |> glue::glue()
    title = "Markers: {hetero_label} (m.{thevariant})" |>
      glue::glue()
  )

  list(
    markers = markers,
    .labs = .labs,
    hetero_label = hetero_label
  ) -> m
  m
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
  # dplyr::filter(variant %in% c("3727T>C", "3728C>T")) |>
  # dplyr::select(tidyselect::contains("ratio"))
  dplyr::select(
    gseid,
    srrid,
    variant,
    forplot,
    # tidyselect::contains("ratio")
  ) -> filtered_data


# thevariant <- "3728C>T"
# thesrrid <- "GSM7080053"
# thegseid <- "GSE226602"
# forplot_ <- filtered_data$forplot[[1]]

# filtered_data |>
#   dplyr::filter(
#     variant == thevariant,
#     srrid == thesrrid
#   ) |>
#   dplyr::select(forplot) |>
#   tidyr::unnest(cols = c(forplot)) -> variant_cell_barcode

#
#
# ? save sct --------------------------------------------------------------------
#
#

filtered_data |>
  # head(6) |>
  dplyr::mutate(
    p = parallel::mcmapply(
      FUN = \(thegseid, thesrrid, thevariant, forplot_) {
        tryCatch(
          expr = {
            fn_sct(
              thegseid = thegseid,
              thesrrid = thesrrid,
              thevariant = thevariant,
              forplot_ = forplot_
            )
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
      mc.cores = 20
    )
  ) -> filtered_data_sct
#

#
# ? plot --------------------------------------------------------------------
#
#

filtered_data |>
  # head(6) |>
  dplyr::mutate(
    p = parallel::mcmapply(
      FUN = \(thegseid, thesrrid, thevariant, forplot_) {
        # thegseid <- filtered_data$gseid[[1]]
        # thesrrid <- filtered_data$srrid[[1]]
        # thevariant <- filtered_data$variant[[1]]
        # forplot_ <- filtered_data$forplot[[1]]
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
      mc.cores = 20
    )
  ) |>
  dplyr::select(-forplot) -> filtered_data_plots


filtered_data_plots |>
  dplyr::mutate(
    saveimage = parallel::mcmapply(
      FUN = \(p, gseid, srrid, variant) {
        if (is.null(p)) {
          return(NULL)
        }
        .outdir <- file.path(
          "/home/liuc9/github/scMOCHA-data/analysis/zzz/plot-real-somatic-variant/main-variants",
          variant,
          "deg"
        )
        dir.create(
          .outdir,
          recursive = TRUE,
          showWarnings = FALSE
        )
        ggsave(
          filename = file.path(
            .outdir,
            "{gseid}-{srrid}-m.{variant}.deg.hetero_vs_sufficient.pdf" |>
              glue::glue()
          ),
          plot = p,
          width = 6,
          height = 5
        )
      },
      p = p,
      gseid = gseid,
      srrid = srrid,
      variant = variant,
      mc.cores = 8,
      SIMPLIFY = FALSE
    )
  )

#
#
# ? plot hetero high vs low --------------------------------------------------------------------
#
#

vss <- list(
  c(0.5, 0.5),
  c(0.6, 0.4),
  c(0.7, 0.3),
  c(0.8, 0.2),
  c(0.9, 0.1)
)


vss |>
  purrr::map(
    .f = \(vs) {
      filtered_data |>
        dplyr::mutate(
          p = parallel::mcmapply(
            FUN = \(thegseid, thesrrid, thevariant, forplot_, vs) {
              tryCatch(
                expr = {
                  m <- fn_de_high_vs_low(
                    thegseid = thegseid,
                    thesrrid = thesrrid,
                    thevariant = thevariant,
                    forplot_ = forplot_,
                    .vs = vs
                  )
                  p <- fn_de_plot(
                    markers = m$markers
                  ) +
                    m$.labs
                  if (is.null(p)) {
                    return(NULL)
                  }
                  .outdir <- file.path(
                    "/home/liuc9/github/scMOCHA-data/analysis/zzz/plot-real-somatic-variant/main-variants",
                    thevariant,
                    "deg"
                  )
                  dir.create(
                    .outdir,
                    recursive = TRUE,
                    showWarnings = FALSE
                  )
                  sanitize_filename <- function(x) {
                    x %>%
                      gsub(">", "GT", ., fixed = TRUE) %>%
                      gsub("<", "LT", ., fixed = TRUE) %>%
                      gsub("%", "pct", ., fixed = TRUE) %>%
                      gsub("=", "-", ., fixed = TRUE) %>%
                      gsub("[()]", "", .) %>%
                      gsub("[[:space:]]+", "_", .) %>%
                      trimws()
                  }
                  .outfilename <- "{thegseid}-{thesrrid}-m.{thevariant}.deg.{p$labels$title}.pdf" |>
                    glue::glue() |>
                    fs::path_sanitize() |>
                    sanitize_filename()
                  ggsave(
                    path = .outdir,
                    filename = .outfilename,
                    plot = p,
                    width = 10,
                    height = 10
                  )
                },
                error = \(e) {
                  message(glue::glue(
                    "{thegseid}-{thesrrid}-m.{thevariant} error"
                  ))
                  return(NULL)
                }
              )
            },
            thegseid = gseid,
            thesrrid = srrid,
            thevariant = variant,
            forplot_ = forplot,
            vs = list(vs),
            SIMPLIFY = FALSE,
            USE.NAMES = FALSE,
            mc.cores = 20
          )
        ) |>
        dplyr::select(-forplot) -> filtered_data_plots
    }
  ) -> m

# footer------------------------------------------------------------------

# future: :plan(future: :sequential)

# save image --------------------------------------------------------------
