#!/usr/bin/env Rscript
# Metainfo ----------------------------------------------------------------
# @AUTHOR: Chun-Jie Liu
# @CONTACT: chunjie.sam.liu.at.gmail.com
# @DATE: 2025-09-12 14:00:19
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
d <- import(
  "/mnt/isilon/u01_project/large-scale/ting/raw/GSE279945/final/GSM8583894/variant_info_from_heatmap.qs"
)
# load conn ---------------------------------------------------------------

# src ---------------------------------------------------------------------

# function ----------------------------------------------------------------

# body --------------------------------------------------------------------

d |>
  # head(20) |>
  tidyr::replace_na(list(
    AF = 0,
    AR = 0,
    CF = 0,
    CR = 0,
    GF = 0,
    GR = 0,
    TF = 0,
    TR = 0
  )) |>
  dplyr::mutate(
    strandbias = parallel::mcmapply(
      variant = variant,
      variant_type = variant_type,
      AF = AF,
      AR = AR,
      CF = CF,
      CR = CR,
      GF = GF,
      GR = GR,
      TF = TF,
      TR = TR,
      FUN = function(variant, variant_type, AF, AR, CF, CR, GF, GR, TF, TR) {
        if (variant_type != "colorful") {
          return(
            tibble::tibble(
              pvalue = NA_real_,
              strand_ratio = NA_real_
            )
          )
        }
        tryCatch(
          expr = {
            ref <- gsub("\\d*|>.*", "", variant)
            alt <- gsub(".*>", "", variant)

            # rf <- ifelse(is.na(get(paste0(ref, "F"))), 0, get(paste0(ref, "F")))
            # rr <- ifelse(is.na(get(paste0(ref, "R"))), 0, get(paste0(ref, "R")))
            # af <- ifelse(is.na(get(paste0(alt, "F"))), 0, get(paste0(alt, "F")))
            # ar <- ifelse(is.na(get(paste0(alt, "R"))), 0, get(paste0(alt, "R")))
            rf <- get(paste0(ref, "F"))
            rr <- get(paste0(ref, "R"))
            af <- get(paste0(alt, "F"))
            ar <- get(paste0(alt, "R"))

            table <- matrix(c(rf, rr, af, ar), nrow = 2, byrow = T)
            colnames(table) <- c("Forward", "Reverse")
            rownames(table) <- c("Ref", "Alt")
            result <- fisher.test(table)
            strand_ratio <- max(af, ar) / (af + ar)
            return(
              tibble::tibble(
                pvalue = result$p.value,
                strand_ratio = strand_ratio
              )
            )
          },
          error = \(e) {
            return(
              tibble::tibble(
                pvalue = NA_real_,
                strand_ratio = NA_real_
              )
            )
          }
        )
      },
      mc.cores = 20,
      SIMPLIFY = FALSE
    )
  ) |>
  tidyr::unnest(cols = strandbias) -> d_strandbias

d_strandbias |> dplyr::filter(!is.na(pvalue)) |> dplyr::filter(pvalue < 0.05)

# footer ------------------------------------------------------------------

# future: :plan(future: :sequential)

# save image --------------------------------------------------------------
