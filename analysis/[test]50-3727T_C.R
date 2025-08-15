#!/usr/bin/env Rscript
# Metainfo ----------------------------------------------------------------
# @AUTHOR: Chun-Jie Liu
# @CONTACT: chunjie.sam.liu.at.gmail.com
# @DATE: 2025-08-05 16:24:15
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

# load conn ---------------------------------------------------------------
conn_all_hetero_af <- DBI::dbConnect(
  duckdb::duckdb(),
  "/home/liuc9/github/scMOCHA-data/analysis/zzz/clean-data/all_hetero_af.cell.duckdb.1.2.1"
)
DBI::dbListTables(conn_all_hetero_af)

thevariant <- "3727T>C"

tbl_all_hetero_af_bulk <- dplyr::tbl(
  conn_all_hetero_af,
  "all_hetero_af_bulk"
)

tbl_all_hetero_af_cluster <- dplyr::tbl(
  conn_all_hetero_af,
  "all_hetero_af_cluster"
)

tbl_gseid_srrid_variant <- dplyr::tbl(
  conn_all_hetero_af,
  "gseid_srrid_variant"
)

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

# body --------------------------------------------------------------------
thevariant <- "3727T>C"
gseid_srrid_variant_hetero |>
  dplyr::filter(variant == thevariant) -> individuals_with_variant

dplyr::bind_rows(
  tbl_all_hetero_af_bulk |>
    dplyr::filter(variant == thevariant) |>
    dplyr::filter(srrid %in% individuals_with_variant$srrid) |>
    dplyr::collect(),
  tbl_all_hetero_af_cluster |>
    dplyr::filter(variant == thevariant) |>
    dplyr::filter(srrid %in% individuals_with_variant$srrid) |>
    dplyr::collect()
) -> all_hetero_af_thevariant

source("/home/liuc9/github/scMOCHA-data/analysis/00-colors.R")
color_celltype_bulk <- c(
  "Pseudo-bulk" = "red",
  color_celltype
)

all_hetero_af_thevariant |>
  dplyr::mutate(
    barcode = gsub(barcode, pattern = "_", replacement = " "),
    barcode = ifelse(barcode == "bulk", "Pseudo-bulk", barcode),
  ) |>
  dplyr::mutate(
    barcode = factor(
      barcode,
      levels = names(color_celltype_bulk)
    ),
  ) -> forplot

forplot |>
  dplyr::filter(barcode == "Pseudo-bulk") |>
  dplyr::arrange(-af) -> rank_pseudo_bulk


rank_pseudo_bulk |> dplyr::filter(af >= 0.05)

forplot |>
  dplyr::mutate(
    srrid = factor(srrid, levels = rank_pseudo_bulk$srrid),
  ) |>
  ggplot(aes(x = srrid, y = af, fill = barcode)) +
  geom_col() +
  geom_hline(
    aes(yintercept = 0.05),
    linetype = 20,
    color = "red"
  ) +
  geom_hline(
    aes(yintercept = 0.1),
    linetype = 21,
    color = "black"
  ) +
  scale_fill_manual(
    name = "Cell type",
    values = color_celltype_bulk
  ) +
  scale_y_continuous(
    name = "Heteroplasmy frequency",
    limits = c(0, 0.4),
    breaks = seq(0, 0.4, by = 0.05),
    labels = scales::percent_format(accuracy = 1),
    expand = expansion(mult = c(0.001, 0.01))
  ) +
  ggh4x::facet_grid2(
    ~barcode,
    strip = ggh4x::strip_themed(
      background_x = ggh4x::elem_list_rect(
        fill = color_celltype_bulk,
        color = NA
      ),
      text_x = ggh4x::elem_list_text(
        colour = "white",
        face = c("bold")
      )
    ),
    switch = "x",
  ) +
  theme(
    plot.margin = margin(t = 0.2, b = 0.1, l = 0.1, r = 0.2, unit = "cm"),
    panel.grid = element_blank(),
    panel.background = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.line = element_line(color = "black"),
    legend.position = c(0.2, 0.6),
    strip.placement = "outside"
  ) +
  labs(
    x = "Individuals",
    title = glue::glue(
      "Heteroplasmy frequency of {thevariant} in 34 individuals with the variant"
    ),
    subtitle = "Pseudo-bulk and cell type specific"
  ) -> p_variant


ggsave(
  plot = p_variant,
  filename = "3727T_C_celltype_AF.pdf",
  path = "/home/liuc9/github/scMOCHA-data/analysis/zzz/plot-real-somatic-variant/main-variants/3727T_C",
  width = 10,
  height = 6
)

# footer ------------------------------------------------------------------

# future: :plan(future: :sequential)

# save image --------------------------------------------------------------
