#!/usr/bin/env Rscript
# Metainfo ----------------------------------------------------------------
# @AUTHOR: Chun-Jie Liu
# @CONTACT: chunjie.sam.liu.at.gmail.com
# @DATE: 2025-08-04 17:14:17
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
library(future.apply)

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

# load conn ---------------------------------------------------------------

conn_all_variant_cell <- DBI::dbConnect(
  duckdb::duckdb(),
  dbdir = "/home/liuc9/github/scMOCHA-data/analysis/zzz/clean-data/all_variant_cell.duckdb.1.2.1"
)
DBI::dbListTables(conn_all_variant_cell)

conn_all_hetero_af <- DBI::dbConnect(
  duckdb::duckdb(),
  "/home/liuc9/github/scMOCHA-data/analysis/zzz/clean-data/all_hetero_af.cell.duckdb.1.2.1"
)
DBI::dbListTables(conn_all_hetero_af)

conn_all_hetero_af |>
  dplyr::tbl("gseid_srrid_variant")

tbl_allvariants <- conn_all_hetero_af |>
  dplyr::tbl("allvariants")


tbl_gseid_srrid_variant <- conn_all_hetero_af |>
  dplyr::tbl("gseid_srrid_variant")

tbl_gseid_srrid_variant |>
  dplyr::collect() |>
  dplyr::mutate(
    a = purrr::map(
      .x = variant_alltype,
      ~ {
        jsonlite::fromJSON(.x) |>
          purrr::pluck("heteroplasmic_variant") -> .v
        if (length(.v) == 0) {
          return(NULL)
        } else {
          return(tibble::tibble(variant = .v))
        }
      }
    )
  ) |>
  dplyr::select(-variant_alltype) |>
  tidyr::unnest(cols = c(a)) -> gseid_srrid_variant_hetero

# src ---------------------------------------------------------------------

# function ----------------------------------------------------------------

fn_plot_somatic_variant <- function(thevariant, thesrrid) {
  dplyr::tbl(conn_all_variant_cell, "all_variant_cell") |>
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
        levels = c("red", "darkblue", "gray", "white")
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
    ) -> forplot
  source("analysis/00-colors.R")

  thetheme <- theme(
    panel.background = element_blank(),
    panel.grid = element_blank(),
    axis.ticks = element_blank(),
    axis.text = element_blank(),
    axis.title.x = element_blank(),
    plot.margin = margin(t = 0, b = 0, unit = "cm"),
  )

  forplot |>
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
    ggplot(aes(
      x = barcode,
      y = 1,
      fill = celltype
    )) +
    geom_col() +
    scale_fill_manual(
      name = "Cell Type",
      values = color_celltype,
    ) +
    scale_y_continuous(expand = expansion(mult = c(0, 0)), ) +
    thetheme +
    # theme(panel.background = element_rect(color = "red")) +
    labs(
      y = "Cell Type",
    ) -> p1_celltype

  viridis::viridis_pal(option = "D")(10) |> color()

  forplot |>
    dplyr::mutate(
      af = ifelse(
        af < 0.01,
        NA_real_,
        af
      )
    ) |>
    ggplot(aes(
      x = barcode,
      y = af,
      fill = af
    )) +
    geom_col() +
    scale_fill_gradient2(
      name = "Allele Frequency",
      high = "#FDE725FF",
      mid = "#21908CFF",
      low = "#440154FF"
    ) +
    scale_y_continuous(
      expand = expansion(mult = c(0, 0)),
    ) +
    theme(
      panel.background = element_blank(),
      panel.grid = element_blank(),
      axis.line = element_line(colour = "black"),
      axis.ticks.x = element_blank(),
      axis.text.x = element_blank(),
      axis.title.x = element_blank(),
      plot.margin = margin(t = 0, b = 0, unit = "cm"),
    ) +
    labs(
      y = "Allele Frequency",
    ) -> p2_af

  forplot |>
    dplyr::mutate(
      variant_type = as.character(variant_type),
    ) |>
    ggplot(aes(
      x = barcode,
      y = 1,
      fill = variant_type
    )) +
    geom_col() +
    scale_y_continuous(expand = expansion(mult = c(0, 0)), ) +
    scale_fill_identity(
      guide = "legend",
      name = "Variant cell",
      breaks = c("red", "darkblue", "gray", "white"),
      labels = c(
        "Heteroplasmy",
        "Suficcient reads",
        "No sufficient reads",
        "No reads"
      )
    ) +
    thetheme +
    labs(
      y = "Variant cells",
    ) -> p3_variant_cells

  forplot |>
    dplyr::mutate(
      depth = log2(depth + 1) # log2 transform to reduce skewness
    ) |>
    ggplot(aes(
      x = barcode,
      y = depth,
      fill = depth
    )) +
    geom_col() +
    scale_fill_gradient(
      name = "log2(depth + 1)",
      high = "gold",
      low = "white"
    ) +
    scale_y_continuous(
      expand = expansion(mult = c(0, 0)),
    ) +
    theme(
      panel.background = element_blank(),
      panel.grid = element_blank(),
      axis.line = element_line(colour = "black"),
      axis.ticks.x = element_blank(),
      axis.text.x = element_blank(),
      axis.title.x = element_blank(),
      plot.margin = margin(t = 0, b = 0, unit = "cm"),
    ) +
    labs(
      y = "Log2(Depth + 1)",
    ) -> p4_depth

  wrap_plots(
    p2_af,
    plot_spacer(),
    p4_depth,
    plot_spacer(),
    p3_variant_cells,
    plot_spacer(),
    p1_celltype,
    ncol = 1,
    heights = c(15, -1.05, 15, -1.05, 10, -1.05, 10),
    guides = "collect"
  ) +
    plot_annotation(
      title = glue::glue(
        "Variant {thevariant} in {thesrrid}"
      ),
      theme = theme(
        plot.title = element_text(hjust = 0.5, size = 16, face = "bold")
      )
    ) -> p_all

  p_all
}

# body --------------------------------------------------------------------
thevariant <- "3173G>A"
thevariants <- c(
  "3173G>A",
  "3176A>T",
  "3178T>A",
  "3727T>C",
  "3728C>T"
)

gseid_srrid_variant_hetero |>
  dplyr::filter(variant %in% thevariants) |>
  dplyr::mutate(
    p = parallel::mcmapply(
      FUN = fn_plot_somatic_variant,
      thevariant = variant,
      thesrrid = srrid,
      SIMPLIFY = FALSE,
      mc.cores = 20
    )
  ) -> gseid_srrid_variant_hetero_plot


gseid_srrid_variant_hetero_plot |>
  dplyr::mutate(
    # save image to /home/liuc9/github/scMOCHA-data/analysis/zzz/plot-real-somatic-variant/main-variants
    a = parallel::mcmapply(
      FUN = \(.gseid, .srrid, .variant, .p) {
        .filename = "{.gseid}_{.srrid}_{.variant}.pdf" |> glue::glue()
        ggsave(
          plot = .p,
          filename = .filename,
          path = "/home/liuc9/github/scMOCHA-data/analysis/zzz/plot-real-somatic-variant/main-variants",
          width = 13,
          height = 8
        )
      },
      .gseid = gseid,
      .srrid = srrid,
      .variant = variant,
      .p = p,
      SIMPLIFY = FALSE,
      mc.cores = 20
    )
  )

# footer ------------------------------------------------------------------

# future: :plan(future: :sequential)

# save image --------------------------------------------------------------
