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
    "C>T",
    "A>G",
    # complementary
    "G>A",
    "T>C"
  )
  ros_damage <- c(
    "G>T",
    "G>C",
    # complementary
    "C>A",
    "T>A"
  )
  other_damage <- c(
    "T>G",
    "T>A",
    # complementary
    "A>C",
    "A>G"
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

# body --------------------------------------------------------------------


v_hete_L_H_strand <- fn_variant_L_H_strand(v_hete)

# footer ------------------------------------------------------------------

# future: :plan(future: :sequential)

# save image --------------------------------------------------------------
