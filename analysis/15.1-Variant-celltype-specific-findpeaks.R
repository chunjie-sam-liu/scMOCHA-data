#!/usr/bin/env Rscript
# Metainfo ----------------------------------------------------------------
# @AUTHOR: Chun-Jie Liu
# @CONTACT: chunjie.sam.liu.at.gmail.com
# @DATE: 2025-06-20 19:17:51
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
cleandatadir <- "/home/liuc9/github/scMOCHA-data/data/zzz/clean-data"
dbdir <- "/home/liuc9/github/scMOCHA-data/analysis/zzz/db"
ks_test_dir <- file.path(dbdir, "all_hetero_af.cell.ks_test")
plotdir <- "/home/liuc9/github/scMOCHA-data/analysis/zzz/plot-celltype-specific-variant"

# gseid_srrid_ks_load <- import(
#   file.path(
#     ks_test_dir,
#     "a_gseid_srrid_ks_load.nocellaf.qs"
#   )
# )



# ALLVARIANTS <- import(file.path(
#   "/home/liuc9/github/scMOCHA-data/analysis/zzz/clean-data/", "all_variant.qs"
# ))

# META <- import("/home/liuc9/github/scMOCHA-data/analysis/zzz/clean-data/gse_dataset_metadata_full.sex_pred.qs") |>
#   dplyr::select(
#     gseid, srrid, Age_new, Age_group,
#     Haplogroup,
#     disease, Chemistry, sex_pred
#   ) |>
#   dplyr::mutate(
#     Haplogroup = purrr::map_chr(
#       .x = Haplogroup,
#       .f = \(.x) {
#         # if (stringr::str_starts(.x, "L")) {
#         #   gsub("L", "L0", .x)
#         # }
#         gsub("\\d+.*", "", .x)
#       }
#     )
#   )


# load conn --------------------------------------------------------------------

conn <- DBI::dbConnect(
  duckdb::duckdb(),
  "/home/liuc9/github/scMOCHA-data/analysis/zzz/clean-data/all_hetero_af.cell.duckdb.1.2.1"
)
DBI::dbListTables(conn)
DBI::dbDisconnect(conn, shutdown = TRUE)

tbl_all_hetero_af_cell <- dplyr::tbl(
  conn,
  "all_hetero_af_cell"
)
tbl_all_hetero_altdepth_cell <- dplyr::tbl(
  conn,
  "all_hetero_altdepth_cell"
)
tbl_all_hetero_sumdepth_cell <- dplyr::tbl(
  conn,
  "all_hetero_sumdepth_cell"
)
tbl_barcode <- dplyr::tbl(
  conn,
  "barcode"
)
tbl_allvariants <- dplyr::tbl(
  conn,
  "allvariants"
)
tbl_meta <- dplyr::tbl(conn, "meta") |>
  dplyr::select(
    gseid, srrid, Age_new, Age_group,
    Haplogroup,
    disease, Chemistry, sex_pred
  )

# DBI::dbDisconnect(conn, shutdown = TRUE)

# src ---------------------------------------------------------------------
source("./analysis/00-colors.R")

# function ----------------------------------------------------------------

fn_plot_ggdist <- function(
    thevariant,
    thegseid,
    thesrrid) {
  library(ggdist)
  tbl_all_hetero_af_cell |>
    dplyr::filter(
      # gseid == thegseid,
      srrid == thesrrid,
      variant == thevariant,
      af > 0
    ) |>
    dplyr::collect() ->
  .d
  # thevariant <- "7833T>C"
  .variant <- .d$variant[1]
  .gseid <- .d$gseid[1]
  .srrid <- .d$srrid[1]

  .d |>
    dplyr::filter(af > 0) |>
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
        levels = names(color_celltype) |> rev()
      )
    ) ->
  forplot_

  forplot_ |>
    ggplot(aes(
      x = af,
      y = celltype,
      fill = celltype
    )) +
    ggdist::stat_halfeye(
      scale = 1
    ) +
    ggdist::stat_interval(
      show.legend = FALSE,
    ) +
    stat_summary(geom = "point", fun = median, show.legend = FALSE) +
    scale_fill_manual(
      values = color_celltype,
      na.value = "grey50"
    ) +
    scale_color_manual(values = MetBrewer::met.brewer("VanGogh3")) +
    # scale_color_brewer() +
    guides(col = "none") +
    ggridges::theme_ridges() +
    # ggridges::stat_density_ridges(
    #   quantile_lines = TRUE, quantiles = 2
    # ) +
    theme(
      legend.position = "none",
      plot.title = element_text(
        hjust = 0.5,
        size = 16
      ),
    ) +
    labs(
      title = paste0(.variant, "\n(", .gseid, "-", .srrid, ")"),
      x = "Allele Frequency",
      y = "Cell Type"
    )
}

fn_plot_joy <- function(
    thevariant,
    thegseid,
    thesrrid) {
  tbl_all_hetero_af_cell |>
    dplyr::filter(
      # gseid == thegseid,
      srrid == thesrrid,
      variant == thevariant,
      af > 0
    ) |>
    dplyr::collect() ->
  .d
  # thevariant <- "7833T>C"
  .variant <- .d$variant[1]
  .gseid <- .d$gseid[1]
  .srrid <- .d$srrid[1]

  .d |>
    dplyr::filter(af > 0) |>
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
        levels = names(color_celltype) |> rev()
      )
    ) ->
  forplot_

  forplot_ |>
    ggplot(aes(
      x = af,
      y = celltype,
      fill = celltype
    )) +
    ggridges::geom_density_ridges(
      # scale = 3,
      # alpha = 0.8,
      rel_min_height = 0.01,
      size = 0.1
    ) +
    scale_fill_manual(
      values = color_celltype,
      na.value = "grey50"
    ) +
    ggridges::theme_ridges() +
    # ggridges::stat_density_ridges(
    #   quantile_lines = TRUE, quantiles = 2
    # ) +
    theme(
      legend.position = "none",
      plot.title = element_text(
        hjust = 0.5,
        size = 16
      ),
    ) +
    labs(
      title = paste0(.variant, "\n(", .gseid, "-", .srrid, ")"),
      x = "Allele Frequency",
      y = "Cell Type"
    )
}

fn_plot_joy_celltype_level2_level3 <- function(
    thevariant,
    thegseid,
    thesrrid,
    thecelltype,
    thecelltype_prefix,
    thecelltype_level) {
  tbl_barcode |>
    dplyr::filter(
      srrid == thesrrid
    ) |>
    dplyr::collect() ->
  celltypedetail

  tbl_all_hetero_af_cell |>
    dplyr::filter(
      gseid == thegseid,
      srrid == thesrrid,
      variant == thevariant,
      af > 0
    ) |>
    dplyr::select(-celltype) |>
    dplyr::collect() |>
    dplyr::left_join(
      celltypedetail,
      by = c("barcode")
    ) |>
    dplyr::rename(
      plotcelltype = "celltype_{thecelltype_level}" |> glue::glue(),
    ) ->
  thevariant_data


  thevariant_data |>
    dplyr::filter(
      celltype == thecelltype,
      af > 0
    ) |>
    dplyr::filter(grepl(thecelltype_prefix, plotcelltype)) |>
    dplyr::mutate(
      plotcelltype = factor(
        plotcelltype
      )
    ) ->
  forplot

  levels(forplot$plotcelltype)


  color_celltype_detail <- log(seq(1, exp(1), length.out = length(levels(forplot$plotcelltype)))) |>
    purrr::map_chr(
      ~ prismatic::clr_lighten(
        color_celltype[thecelltype],
        .x
      )
    )
  names(color_celltype_detail) <- levels(forplot$plotcelltype)

  forplot |>
    ggplot(aes(
      x = af,
      y = plotcelltype,
      fill = plotcelltype
    )) +
    ggridges::geom_density_ridges(
      # scale = 3,
      # alpha = 0.8,
      rel_min_height = 0.01,
      size = 0.1
    ) +
    scale_fill_manual(
      values = color_celltype_detail,
      na.value = "grey50"
    ) +
    ggridges::theme_ridges() +
    theme(
      legend.position = "none",
      plot.title = element_text(
        hjust = 0.5,
        # size = 16
      ),
    ) +
    labs(
      title = "{thecelltype}-{thecelltype_level}-{thevariant}\n({thegseid}-{thesrrid})" |> glue::glue(),
      x = "Allele Frequency",
      y = "Cell Type"
    )
}

fn_plot_joy_celltype_detail <- function(
    thevariant,
    thegseid,
    thesrrid) {
  tibble::tibble(
    thevariant = thevariant,
    thegseid = thegseid,
    thesrrid = thesrrid,
    thecelltype = c(
      c("B", "CD4 T", "CD8 T", "other T") |> rep(times = 2),
      c("NK", "DC", "Mono", "other") |> rep(times = 2)
    ),
    thecelltype_prefix = c(
      c("B", "CD4", "CD8", "") |> rep(times = 2),
      c("NK", "DC", "Mono", "") |> rep(times = 2)
    ),
    thecelltype_level = c("l2", "l3") |> rep(each = 4) |> rep(times = 2)
  ) ->
  thevariant_celltype_df

  thevariant_celltype_df |>
    dplyr::mutate(
      p = parallel::mcmapply(
        thevariant = thevariant,
        thegseid = thegseid,
        thesrrid = thesrrid,
        thecelltype = thecelltype,
        thecelltype_prefix = thecelltype_prefix,
        thecelltype_level = thecelltype_level,
        FUN = fn_plot_joy_celltype_level2_level3,
        mc.cores = 5,
        SIMPLIFY = FALSE
      )
    ) ->
  plot_thevariant_celltype_list
}

find_density_peaks <- function(x, bw = "nrd0", min_height = 0.01, min_prominence = 0.1, min_distance = 0.05, ...) {
  # Calculate density
  dx <- density(x, bw = bw, ...)

  # Find local maxima
  n <- length(dx$y)
  peaks <- c()

  # Check each point (excluding endpoints)
  for (i in 2:(n - 1)) {
    if (dx$y[i] > dx$y[i - 1] && dx$y[i] > dx$y[i + 1] && dx$y[i] >= min_height) {
      peaks <- c(peaks, i)
    }
  }

  # Filter peaks by prominence and distance
  if (length(peaks) > 0) {
    peaks_info <- data.frame(
      idx = peaks,
      x = dx$x[peaks],
      y = dx$y[peaks]
    )

    # Calculate prominence for each peak
    peaks_info$prominence <- 0
    for (i in 1:nrow(peaks_info)) {
      peak_idx <- peaks_info$idx[i]
      peak_height <- peaks_info$y[i]

      # Find minimum heights on both sides
      left_min <- min(dx$y[1:peak_idx])
      right_min <- min(dx$y[peak_idx:n])
      baseline <- max(left_min, right_min)

      peaks_info$prominence[i] <- peak_height - baseline
    }

    # Filter by prominence
    peaks_info <- peaks_info[peaks_info$prominence >= min_prominence, ]

    # Filter by minimum distance between peaks
    if (nrow(peaks_info) > 1) {
      # Sort by height (keep strongest peaks when too close)
      peaks_info <- peaks_info[order(peaks_info$y, decreasing = TRUE), ]

      keep <- rep(TRUE, nrow(peaks_info))
      for (i in 1:(nrow(peaks_info) - 1)) {
        if (keep[i]) {
          for (j in (i + 1):nrow(peaks_info)) {
            if (abs(peaks_info$x[i] - peaks_info$x[j]) < min_distance) {
              keep[j] <- FALSE
            }
          }
        }
      }
      peaks_info <- peaks_info[keep, ]
    }

    # Sort by density value (highest first) and remove helper columns
    peaks_info <- peaks_info[order(peaks_info$y, decreasing = TRUE), c("x", "y", "prominence")]
  } else {
    peaks_info <- data.frame(x = numeric(0), y = numeric(0), prominence = numeric(0))
  }

  return(list(
    peaks = peaks_info,
    density_obj = dx
  ))
}

find_cell_density_peak <- function(x, pheight = 0.1, pprominence = 0.2) {
  dx <- density(x)
  # Find multiple peaks with adaptive filtering for bimodal distributions
  max_density <- max(dx$y)
  adaptive_min_height <- max_density * pheight # 10% of maximum density
  adaptive_min_prominence <- max_density * pprominence # 20% of maximum density

  peaks_result <- find_density_peaks(x, min_height = adaptive_min_height, min_prominence = adaptive_min_prominence, min_distance = 0.2)

  data.table(
    npeaks = nrow(peaks_result$peaks),
    peaks = list(peaks_result$peaks)
  )
}


# body --------------------------------------------------------------------

thevariant <- "3727T>C"
tbl_all_hetero_af_cell |>
  dplyr::filter(
    variant == thevariant,
    af > 0,
  ) |>
  dplyr::left_join(
    tbl_all_hetero_altdepth_cell,
    by = c("gseid", "srrid", "variant", "barcode")
  ) |>
  dplyr::left_join(
    tbl_all_hetero_sumdepth_cell,
    by = c("gseid", "srrid", "variant", "barcode")
  ) |>
  dplyr::left_join(
    tbl_barcode,
    by = c("gseid", "srrid", "barcode", "celltype")
  ) |>
  dplyr::select(
    variant, srrid, barcode, celltype, af, altdepth, sumdepth,
    celltype_l2, celltype_l3
  ) |>
  dplyr::collect() ->
thevariant_data

thevariant_data |>
  tidyr::nest(
    .by = c(variant, srrid, celltype),
    .key = "data"
  ) |>
  dplyr::mutate(
    peaks = parallel::mclapply(
      X = data,
      FUN = \(.x) {
        tryCatch(
          find_cell_density_peak(.x$af),
          error = function(e) {
            message(glue::glue("Error processing {unique(.x$variant)} in {unique(.x$srrid)}: {e$message}"))
            data.table(npeaks = 0, peaks = list(data.frame(x = numeric(0), y = numeric(0), prominence = numeric(0))))
          }
        )
      },
      mc.cores = 5
    )
  ) ->
thevariant_data_peaks

thevariant_data_peaks |>
  tidyr::unnest(cols = peaks) |>
  dplyr::select(
    srrid, celltype, npeaks
  ) |>
  tidyr::spread(
    key = celltype,
    value = npeaks
  ) |>
  dplyr::mutate(
    total_peaks = rowSums(dplyr::across(where(is.numeric)), na.rm = TRUE)
  ) |>
  dplyr::arrange(-total_peaks) ->
thevariant_data_peaks_summary

thevariant_data_peaks_summary |>
  dplyr::filter(total_peaks < 10)

fn_plot_joy(
  thevariant = thevariant,
  thesrrid = "GSM4697614"
)
fn_plot_joy_celltype_detail(
  thevariant = thevariant,
  thesrrid = "GSM4697614"
) -> a


# footer ------------------------------------------------------------------

# future: :plan(future: :sequential)

# save image --------------------------------------------------------------
DBI::dbDisconnect(conn, shutdown = TRUE)
