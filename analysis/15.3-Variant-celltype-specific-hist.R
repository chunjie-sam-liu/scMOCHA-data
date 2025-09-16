#!/usr/bin/env Rscript
# Metainfo ----------------------------------------------------------------
# @AUTHOR: Chun-Jie Liu
# @CONTACT: chunjie.sam.liu.at.gmail.com
# @DATE: 2025-09-15 13:22:36
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
library(scales)

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

# load conn ---------------------------------------------------------------

conn <- DBI::dbConnect(
  duckdb::duckdb(),
  "/home/liuc9/github/scMOCHA-data/analysis/zzz/clean-data/all_hetero_af.cell.duckdb.1.2.1"
)
DBI::dbListTables(conn)
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
tbl_allvariants |> dplyr::glimpse()

tbl_gseid_srrid_srrdir <- dplyr::tbl(
  conn,
  "gseid_srrid_srrdir"
)

tbl_meta <- dplyr::tbl(
  conn,
  "meta"
)

tbl_gseid_srrid_variant <- dplyr::tbl(
  conn,
  "gseid_srrid_variant"
) |>
  dplyr::collect() |>
  dplyr::mutate(
    variant = purrr::map(
      .x = variant_alltype,
      ~ {
        .x |>
          jsonlite::fromJSON() -> .xx
        .hetero <- .xx$heteroplasmic_variant
        if (length(.hetero) == 0) {
          return(NULL)
        } else {
          return(list(.hetero))
        }
      }
    )
  ) |>
  dplyr::select(gseid, srrid, variant) |>
  tidyr::unnest(cols = variant) |>
  tidyr::unnest(cols = variant)

tbl_gseid_srrid_variant_celltype_ks_test <- dplyr::tbl(
  conn,
  "gseid_srrid_variant_celltype_ks_test"
) |>
  dplyr::collect() |>
  dplyr::semi_join(
    tbl_gseid_srrid_variant,
    by = c("gseid", "srrid", "variant")
  )

# src ---------------------------------------------------------------------

source("./analysis/00-colors.R")

# function ----------------------------------------------------------------

fn_plot_hist <- function(
  thevariant,
  thegseid,
  thesrrid
) {
  tbl_all_hetero_af_cell |>
    dplyr::filter(
      # gseid == thegseid,
      srrid == thesrrid,
      variant == thevariant,
      af > 0
    ) |>
    dplyr::collect() -> .d
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
        levels = names(color_celltype)
      )
    ) -> forplot_

  forplot_ |>
    dplyr::count(celltype) |>
    dplyr::mutate(
      label = glue::glue(
        "n={scales::label_comma()(n)}"
      )
    ) -> .forlabel

  fn_xy_breaks_limits(forplot_$af, step = 0.2) -> .xbl
  # fn_xy_breaks_limits(forplot_$celltype, step = 1) -> .ybl

  forplot_ |>
    ggplot(aes(
      x = af,
      # y = celltype,
      fill = celltype
    )) +
    geom_histogram(binwidth = 0.1) +
    scale_fill_manual(
      values = color_celltype,
    ) +
    geom_text(
      data = .forlabel,
      aes(
        x = 0.5,
        y = Inf,
        label = label
      ),
      vjust = 1.5
    ) +
    scale_y_continuous(
      labels = scales::label_number(accuracy = 1),
      expand = expansion(mult = c(0.01, 0.01)),
    ) +
    scale_x_continuous(
      # limits = .xbl$limits,
      breaks = seq(0, 1, 0.2),
      labels = scales::label_number(accuracy = 0.1),
      expand = expansion(add = c(0.01, 0.01)),
    ) +
    theme(
      # text = element_text(family = "Times New Roman"),
      legend.position = "none",
      # axis.text.y = element_blank(),
      axis.line = element_line(color = "black"),
      strip.background = element_rect(
        fill = "white",
        color = "black"
      ),
      plot.title = element_text(
        hjust = 0.5,
        # size = 16
      ),
      panel.background = element_blank(),
    ) +
    ggh4x::facet_wrap2(
      ~celltype,
      nrow = 1,
      # ncol = 8,
      strip.position = "top",
      strip = ggh4x::strip_themed(
        background_x = ggh4x::elem_list_rect(
          fill = color_celltype
        ),
        text_x = ggh4x::elem_list_text(
          colour = "white",
          face = c("bold")
        ),
        by_layer_y = TRUE,
      ),
      scales = "free_y",
    ) +
    labs(
      title = glue::glue("m.{.variant}({.gseid}-{.srrid})"),
      x = glue::glue("m.{.variant} Heteroplasmy Level"),
      y = "Cell Count"
    )

  # forplot_ |>
  #   ggplot(aes(
  #     x = af,
  #     color = celltype
  #   )) +
  #   stat_ecdf(
  #     geom = "step",
  #     size = 1
  #   ) +
  #   scale_color_manual(
  #     values = color_celltype,
  #   )
}

fn_plot_joy_celltype_level2_level3 <- function(
  thevariant,
  thegseid,
  thesrrid,
  thecelltype,
  thecelltype_prefix,
  thecelltype_level
) {
  celltypedetail <- import(
    "/mnt/isilon/u01_project/large-scale/liuc9/raw/{thegseid}/final/{thesrrid}/sc_azimuth_celltype.csv" |>
      glue::glue()
  )

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
    ) -> thevariant_data

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
    ) -> forplot

  levels(forplot$plotcelltype)

  color_celltype_detail <- log(seq(
    1,
    exp(1),
    length.out = length(levels(forplot$plotcelltype))
  )) |>
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
      title = "{thecelltype}-{thecelltype_level}-{thevariant}\n({thegseid}-{thesrrid})" |>
        glue::glue(),
      x = "Allele Frequency",
      y = "Cell Type"
    )
}

fn_plot_joy_celltype_detail <- function(
  thevariant,
  thegseid,
  thesrrid
) {
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
  ) -> thevariant_celltype_df

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
    ) -> plot_thevariant_celltype_list
}

# body --------------------------------------------------------------------
pheno <- import(
  "/mnt/isilon/u01_project/large-scale/ting/raw/GSE235050/GSE235050.pheno.csv"
) |>
  dplyr::mutate(
    srrid = geo_accession,
    status = `status:ch1`,
    sex = `Sex:ch1`,
    age = `age:ch1`,
    treatment = `treatment:ch1`
  ) |>
  dplyr::select(srrid, status, sex, age, treatment)
thevariant <- "3727T>C"
thegseid <- "GSE235050"
tbl_gseid_srrid_variant_celltype_ks_test |>
  dplyr::filter(
    variant == thevariant,
    gseid == thegseid,
    # p.value < 0.05
  ) |>
  dplyr::left_join(
    pheno,
    by = c("srrid")
  )


tbl_gseid_srrid_variant_celltype_ks_test |>
  dplyr::filter(
    gseid == thegseid
  ) |>
  dplyr::count(variant) |>
  dplyr::arrange(-n, variant) |>
  head(20) |>
  print(n = Inf)
# footer ------------------------------------------------------------------

# future: :plan(future: :sequential)

# save image --------------------------------------------------------------
