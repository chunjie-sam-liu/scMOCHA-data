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
#' Very important function
#' @example fn_plot_cell_af_depth_forplot("10398A>G", "GSM5494107")
fn_plot_cell_af_depth_forplot <- function(thevariant, thesrrid) {
  conn <- DBI::dbConnect(
    duckdb::duckdb(),
    "/home/liuc9/github/scMOCHA-data/analysis/zzz/clean-data/all_hetero_af.cell.duckdb.1.2.1",
    read_only = TRUE
  )
  source("analysis/00-colors.R")

  colorcode <- setNames(names(color_variantcell), color_variantcell)

  dplyr::tbl(
    conn,
    "allvariants_cell"
  ) |>
    dplyr::filter(
      srrid == thesrrid,
      variant == thevariant
    ) |>
    dplyr::collect() |>
    dplyr::mutate(
      variant_type = dplyr::case_match(
        variant_type,
        "colorful" ~ "red",
        "black" ~ "darkblue",
        "white" ~ "white",
        "grey" ~ "gray",
        NA ~ "white"
      )
    ) |>
    dplyr::mutate(
      variant_type = factor(
        variant_type,
        levels = color_variantcell
      )
    ) |>
    dplyr::arrange(
      variant_type,
      -af
    ) -> forplot_

  forplot_ |>
    dplyr::mutate(
      barcode = factor(
        barcode,
        levels = forplot_$barcode
      )
    ) |>
    dplyr::mutate(
      celltype = gsub(
        "_",
        " ",
        celltype
      )
    ) |>
    dplyr::mutate(
      celltype = factor(
        celltype,
        names(color_celltype)
      )
    ) |>
    dplyr::mutate(
      af = ifelse(
        af < 0.01,
        NA_real_,
        af
      )
    ) |>
    # dplyr::mutate(
    #   variant_type = as.character(variant_type),
    # ) |>
    dplyr::mutate(
      depth = log2(depth + 1) # log2 transform to reduce skewness
    ) |>
    dplyr::mutate(
      cellvarianttype = colorcode[variant_type]
    ) |>
    dplyr::mutate(
      cellvarianttype = factor(
        cellvarianttype,
        levels = colorcode
      )
    ) -> forplot
  DBI::dbDisconnect(conn, shutdown = TRUE)
  forplot
}

fn_load_sc_and_sct <- function(.filepath, thegseid, thesrrid, forplot_) {
  .sc <- import(.filepath)
  sc_azimuth <- .sc$sc_azimuth
  rm(.sc)
  gc()
  sc_azimuth@meta.data |>
    tibble::rownames_to_column("barcode") |>
    as.data.table() |>
    # dplyr::left_join(
    #   forplot_ |>
    #     as.data.table() |>
    #     dplyr::mutate(
    #       barcode = as.character(barcode)
    #     ),
    # ) |>
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
library(Seurat)
fn_de <- function(
  thegseid,
  thesrrid,
  thevariant,
  .vs = c("Heteroplasmy", "Sufficient reads"),
  .celltype = NA_character_
) {
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
  .dir_markers <- path(
    .dir_de,
    "markers"
  )
  dir_create(.dir_markers)
  forplot_ <- fn_plot_cell_af_depth_forplot(
    thevariant = thevariant,
    thesrrid = thesrrid
  )

  .sct_filepath <- path(
    .dir_de,
    "sc_azimuth.sct.qs"
  )
  sc_azimuth <- if (is.na(.celltype)) {
    import(.sct_filepath)
  } else {
    import(.sct_filepath) |>
      subset(
        subset = predicted.celltype.l1 == .celltype
      )
  }
  sc_azimuth@meta.data |>
    as.data.table() |>
    dplyr::left_join(
      forplot_ |>
        as.data.table() |>
        dplyr::mutate(
          barcode = as.character(barcode)
        ),
    ) |>
    as.data.frame() -> .d_merge
  sc_azimuth@meta.data <- .d_merge

  sc_azimuth@meta.data |>
    dplyr::count(cellvarianttype) |>
    dplyr::arrange(cellvarianttype) |>
    tibble::deframe() -> n_cellvarianttype

  if (all(.vs == c("Heteroplasmy", "Sufficient reads"))) {
    markers <- Seurat::FindMarkers(
      object = sc_azimuth,
      ident.1 = "Heteroplasmy",
      ident.2 = "Sufficient reads",
      test.use = "wilcox",
      group.by = "cellvarianttype"
    )
    .prefixout <- path(
      "markers.hetero_vs_sufficient.{thegseid}.{thesrrid}.{thevariant}.{ifelse(is.na(.celltype), 'all', .celltype)}_" |>
        glue::glue()
    )
    log_fatal("{.dir_markers}/{.prefixout}  run!" |> glue::glue())
    export(
      markers,
      file = path(
        .dir_markers,
        .prefixout
      ),
      format = "qs"
    )
    return(
      list(
        markers = markers,
        labs = labs(
          x = "Fold change Heteroplasmy (n={
          scales::label_comma()(n_cellvarianttype['Heteroplasmy'])
        }) vs Sufficient reads (n={
          scales::label_comma()(n_cellvarianttype['Sufficient reads'])
        })" |>
            glue::glue(),
          y = "FDR",
          title = "Markers: Heteroplasmy vs Sufficient Reads (m.{thevariant}) {ifelse(is.na(.celltype), '', .celltype)}" |>
            glue::glue()
        ),
        hetero_label = "Heteroplasmy vs Sufficient reads",
        celltype = .celltype,
        prefixout = .prefixout
      )
    )
  }

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

  .prefixout <- path(
    "markers.{hetero_label}.{thegseid}.{thesrrid}.{thevariant}.{ifelse(is.na(.celltype), 'all', .celltype)}_" |>
      glue::glue()
  )
  export(
    markers,
    file = path(
      .dir_markers,
      .prefixout
    ),
    format = "qs"
  )
  return(
    list(
      markers = markers,
      labs = labs(
        x = glue::glue(
          "Fold change Heteroplasmy ({.label_high} n={
          scales::label_comma()(n_cellvarianttype2[.label_high])
        }) vs Low ({.label_low} n={
          scales::label_comma()(n_cellvarianttype2[.label_low])
        })"
        ),
        y = "FDR",
        # title = "m.{thevariant}" |> glue::glue()
        title = "Markers: {hetero_label} (m.{thevariant}) {ifelse(is.na(.celltype), '', .celltype)}" |>
          glue::glue()
      ),
      hetero_label = hetero_label,
      celltype = .celltype,
      prefixout = .prefixout
    )
  )
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

  forplot |> dplyr::count(color) |> tibble::deframe() -> n_color

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
      ),
      plot.subtitle = element_text(
        hjust = 0.5,
        face = "bold",
        color = "black",
        size = 14
      ),
    ) +
    labs(
      subtitle = "Up={ifelse(is.na(n_color['red']), 0, n_color['red'])}, Down={ifelse(is.na(n_color['blue']), 0, n_color['blue'])}" |>
        glue::glue()
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

filtered_data |> dplyr::count(variant) |> dplyr::arrange(variant)

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
            # thegseid <- filtered_data$gseid[[1]]
            # thesrrid <- filtered_data$srrid[[1]]
            # thevariant <- filtered_data$variant[[1]]
            # forplot_ <- filtered_data$forplot[[1]]
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

vss <- list(
  c("Heteroplasmy", "Sufficient reads"),
  c(0.5, 0.5),
  c(0.6, 0.4),
  c(0.7, 0.3),
  c(0.8, 0.2),
  c(0.9, 0.1)
)
celltypes <- c(
  NA_character_,
  "CD4 T",
  "CD8 T",
  "B",
  "NK",
  "other T",
  "Mono",
  "DC",
  "other"
)
celltypes |>
  purrr::map(
    .f = \(.celltype) {
      filtered_data |>
        dplyr::select(-forplot) |>
        dplyr::mutate(celltype = .celltype) -> .d
      vss |>
        purrr::map(
          .f = \(.vs) {
            .d |>
              dplyr::mutate(vs = list(.vs))
          }
        ) |>
        dplyr::bind_rows()
    }
  ) |>
  dplyr::bind_rows() -> filtered_data_new


filtered_data_new |>
  # head(6) |>
  dplyr::mutate(
    p = parallel::mcmapply(
      FUN = \(thegseid, thesrrid, thevariant, celltype, vs) {
        # thegseid <- filtered_data_new$gseid[[1]]
        # thesrrid <- filtered_data_new$srrid[[1]]
        # thevariant <- filtered_data_new$variant[[1]]
        # .celltype <- filtered_data_new$celltype[[1]]
        # .vs <- filtered_data_new$vs[[1]]

        tryCatch(
          expr = {
            m <- fn_de(
              thegseid = thegseid,
              thesrrid = thesrrid,
              thevariant = thevariant,
              .celltype = celltype,
              .vs = vs
            )
            p <- fn_de_plot(
              markers = m$markers
            ) +
              m$labs
            if (is.null(p)) {
              return(NULL)
            }
            m$p <- p
            .outdir <- path(
              "/home/liuc9/github/scMOCHA-data/analysis/zzz/plot-real-somatic-variant/main-variants",
              thevariant,
              "deg_new"
            )
            dir_create(.outdir)
            ggsave(
              path = .outdir,
              filename = "{m$prefixout}.pdf" |> glue::glue(),
              plot = p,
              width = 13,
              height = 7
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
      celltype = celltype,
      vs = vs,
      SIMPLIFY = FALSE,
      USE.NAMES = FALSE,
      mc.cores = 20
    )
  ) -> filtered_data_plots

#
#
# ? plot hetero high vs low --------------------------------------------------------------------
#
#

# vss <- list(
#   c(0.5, 0.5),
#   c(0.6, 0.4),
#   c(0.7, 0.3),
#   c(0.8, 0.2),
#   c(0.9, 0.1)
# )

# vss |>
#   purrr::map(
#     .f = \(vs) {
#       filtered_data |>
#         dplyr::mutate(
#           p = parallel::mcmapply(
#             FUN = \(thegseid, thesrrid, thevariant, forplot_, vs) {
#               tryCatch(
#                 expr = {
#                   m <- fn_de_high_vs_low(
#                     thegseid = thegseid,
#                     thesrrid = thesrrid,
#                     thevariant = thevariant,
#                     forplot_ = forplot_,
#                     .vs = vs
#                   )
#                   p <- fn_de_plot(
#                     markers = m$markers
#                   ) +
#                     m$.labs
#                   if (is.null(p)) {
#                     return(NULL)
#                   }
#                   .outdir <- file.path(
#                     "/home/liuc9/github/scMOCHA-data/analysis/zzz/plot-real-somatic-variant/main-variants",
#                     thevariant,
#                     "deg"
#                   )
#                   dir.create(
#                     .outdir,
#                     recursive = TRUE,
#                     showWarnings = FALSE
#                   )
#                   sanitize_filename <- function(x) {
#                     x %>%
#                       gsub(">", "GT", ., fixed = TRUE) %>%
#                       gsub("<", "LT", ., fixed = TRUE) %>%
#                       gsub("%", "pct", ., fixed = TRUE) %>%
#                       gsub("=", "-", ., fixed = TRUE) %>%
#                       gsub("[()]", "", .) %>%
#                       gsub("[[:space:]]+", "_", .) %>%
#                       trimws()
#                   }
#                   .outfilename <- "{thegseid}-{thesrrid}-m.{thevariant}.deg.{p$labels$title}.pdf" |>
#                     glue::glue() |>
#                     fs::path_sanitize() |>
#                     sanitize_filename()
#                   ggsave(
#                     path = .outdir,
#                     filename = .outfilename,
#                     plot = p,
#                     width = 10,
#                     height = 10
#                   )
#                 },
#                 error = \(e) {
#                   message(glue::glue(
#                     "{thegseid}-{thesrrid}-m.{thevariant} error"
#                   ))
#                   return(NULL)
#                 }
#               )
#             },
#             thegseid = gseid,
#             thesrrid = srrid,
#             thevariant = variant,
#             forplot_ = forplot,
#             vs = list(vs),
#             SIMPLIFY = FALSE,
#             USE.NAMES = FALSE,
#             mc.cores = 20
#           )
#         ) |>
#         dplyr::select(-forplot) -> filtered_data_plots
#     }
#   ) -> m

# footer------------------------------------------------------------------

# future: :plan(future: :sequential)

# save image --------------------------------------------------------------
