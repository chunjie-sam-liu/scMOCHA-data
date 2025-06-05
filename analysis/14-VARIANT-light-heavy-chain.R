#!/usr/bin/env Rscript
# Metainfo ----------------------------------------------------------------
# @AUTHOR: Chun-Jie Liu
# @CONTACT: chunjie.sam.liu.at.gmail.com
# @DATE: 2025-06-05 10:53:21
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

# src ---------------------------------------------------------------------

# header ------------------------------------------------------------------


# future: :plan(future: :multisession, workers = 10)


# load data ---------------------------------------------------------------

cleandatadir <- "/home/liuc9/github/scMOCHA-data/data/zzz/clean-data"
all_variant <- import(
  file.path(cleandatadir, "all_variant.qs")
)
all_variant |>
  dplyr::filter(
    issomatic == "heteroplasmic"
  ) ->
v_hete

all_variant |>
  dplyr::filter(
    issomatic == "homoplasmic"
  ) ->
v_homo

all_variant |>
  dplyr::filter(
    issomatic %in% c("homoplasmic", "heteroplasmic")
  ) ->
v_homo_hete

# function ----------------------------------------------------------------

fn_variant_L_H_strand <- function(.v) {
  # .v <- v_hete
  Ori <- c(
    "pos1" = 210,
    "pos2" = 16172
  )
  Other <- c(211, 16171)
  L_ref <- c("C", "A")
  L_variant <- c(
    "C>A",
    "C>G",
    "C>T",
    "A>C",
    "A>G",
    "A>T"
  )
  H_ref <- c("G", "T")
  H_variant <- c(
    "G>A",
    "G>C",
    "G>T",
    "T>A",
    "T>C",
    "T>G"
  )
  variant_reverse <- c(
    "G>A" = "C>T",
    "G>C" = "C>G",
    "G>T" = "C>A",
    "T>A" = "A>T",
    "T>C" = "A>G",
    "T>G" = "A>C",
    "C>A" = "G>T",
    "C>G" = "G>C",
    "C>T" = "G>A",
    "A>C" = "T>G",
    "A>G" = "T>C",
    "A>T" = "T>A"
  )
  spontaneous_deamination <- c(
    "A>G",
    "C>T",
    # complementary
    "T>C",
    "G>A"
  )
  ros_damage <- c(
    "G>T",
    "G>C",
    # complementary
    "C>A",
    "C>G"
  )
  other_damage <- c(
    "T>G",
    "T>A",
    # complementary
    "A>C",
    "A>T"
  )
  .v |>
    dplyr::mutate(
      variant_short = gsub(
        "[0-9]*",
        "",
        variant
      )
    ) |>
    tidyr::separate(
      variant_short,
      into = c("ref", "alt"),
      remove = FALSE,
    ) |>
    dplyr::mutate(
      variant_location = ifelse(
        dplyr::between(
          Position, Other[1], Other[2]
        ),
        "Other region",
        "Ori region"
      )
    ) |>
    dplyr::mutate(
      L_H_strand = ifelse(
        variant_short %in% L_variant,
        "L",
        "H"
      )
    ) |>
    dplyr::mutate(
      deamination_ros = dplyr::case_when(
        variant_short %in% spontaneous_deamination ~ "Spontaneous deamination",
        variant_short %in% ros_damage ~ "ROS damage",
        variant_short %in% other_damage ~ "Other damage",
        TRUE ~ "Unknown"
      )
    ) |>
    dplyr::mutate(
      variant_six = purrr::map2_chr(
        .x = variant_short,
        .y = L_H_strand,
        .f = function(x, y) {
          if (y == "L") {
            return(x)
          } else {
            return(variant_reverse[x])
          }
        }
      )
    )
}

fn_plot_pie <- function(.d, .colors = NULL) {
  .d |>
    dplyr::select(group = 1, n) |>
    dplyr::arrange(-n) |>
    dplyr::mutate(csum = rev(cumsum(rev(n)))) %>%
    dplyr::mutate(pos = n / 2 + dplyr::lead(csum, 1)) %>%
    dplyr::mutate(pos = dplyr::if_else(is.na(pos), n / 2, pos)) %>%
    dplyr::mutate(percentage = n / sum(n)) |>
    dplyr::mutate(group = factor(group, levels = group)) ->
  .dd

  .scalefill <- if (is.null(.colors)) {
    ggsci::scale_fill_aaas(
      name = NULL
    )
  } else {
    scale_fill_manual(
      name = NULL,
      values = .colors
    )
  }
  .scalecolor <- if (is.null(.colors)) {
    ggsci::scale_color_aaas(
      name = NULL
    )
  } else {
    scale_color_manual(
      name = NULL,
      values = .colors
    )
  }

  .dd |>
    ggplot(aes(
      x = "",
      y = n,
    )) +
    geom_bar(
      aes(fill = group),
      stat = "identity",
      width = 1,
      color = "white",
      show.legend = FALSE
    ) +
    .scalefill +
    ggrepel::geom_label_repel(
      aes(
        y = pos,
        label = glue::glue("{group}\n{n} ({scales::percent(percentage)})"),
        color = group,
      ),
      size = 6,
      nudge_x = 1,
      nudge_y = 0,
      show.legend = FALSE,
      max.overlaps = Inf,
    ) +
    .scalecolor +
    coord_polar(theta = "y", start = 0) +
    theme_void() +
    theme(
      plot.title = element_text(
        hjust = 0.5,
        size = 22,
      ),
      # legend.position = "none"
    )
}

# body --------------------------------------------------------------------


v_hete_L_H_strand <- fn_variant_L_H_strand(v_hete)
v_homo_hete_L_H_strand <- fn_variant_L_H_strand(v_homo_hete)

v_hete_L_H_strand |>
  dplyr::count(variant_location)

v_hete_L_H_strand |>
  ggplot(aes(
    x = variant_six,
    fill = L_H_strand
  )) +
  geom_bar()

v_hete_L_H_strand |>
  dplyr::count(deamination_ros) |>
  fn_plot_pie()

# footer ------------------------------------------------------------------

# future: :plan(future: :sequential)

# save image --------------------------------------------------------------
