#!/usr/bin/env Rscript
# Metainfo ----------------------------------------------------------------
# @AUTHOR: Chun-Jie Liu
# @CONTACT: chunjie.sam.liu.at.gmail.com
# @DATE: 2025-02-03 17:08:06
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
log_threshold(TRACE)
log_layout(layout_glue_colors)

# future: :plan(future: :multisession, workers = 10)

# function ----------------------------------------------------------------


# load data ---------------------------------------------------------------


basedir <- "/home/liuc9/github/scMOCHA-data/data"
outdir <- "/home/liuc9/github/scMOCHA-data/data/out_variant_check"

thepaths <- c(
  "/mnt/isilon/u01_project/large-scale/liuc9/raw/GSE226602/final/GSM7080044",
  "/mnt/isilon/u01_project/large-scale/liuc9/raw/GSE163668/final/GSM4995425",
  "/mnt/isilon/u01_project/large-scale/liuc9/raw/GSE163668/final/GSM4995448",
  "/mnt/isilon/u01_project/large-scale/liuc9/raw/GSE279945/final/GSM8583898",
  "/home/liuc9/github/scMOCHA-data/data/scfoundation/GSE162117/final/GSM4933442"
)

sc5p_pe_variant <- c(
  "1255T>C", "1314C>T", "1315G>A", "1380G>T", "1382A>T", "1397T>A", "1670A>G", "2191A>C", "2285T>C", "2289G>T", "2442T>C", "3173G>A", "3176A>T", "3178T>A", "3727T>C",
  "3728C>T", "3734A>G", "7428G>A", "11560A>G", "13752T>G", "13954C>A", "14082C>G", "15666T>C"
)
sc5p_r2_variant <- c(
  "153A>G", "195T>C", "827A>G", "1002C>T", "2352T>C", "3547A>G", "3766T>C", "4820G>A", "4977T>C", "6164C>T", "6473C>T", "8362T>G", "8598T>C", "8730A>G", "9196G>A",
  "9497T>C", "10604T>A", "10819A>G", "11177C>T", "14212T>C", "14905G>A", "15047G>A", "15535C>T", "15747T>C"
)
sc3pv2_variant <- c(
  "3010G>A", "6260G>A", "8251G>A", "9055G>A", "9150A>G", "9698T>C", "9950T>C", "9974C>T", "10211C>T", "11840C>T", "12016C>A", "15662A>G"
)

sample_chem <- tibble::tibble(
  gseid = c("GSE226602", "GSE163668", "GSE163668", "GSE279945", "GSE162117"),
  gsmid = c("GSM7080044", "GSM4995425", "GSM4995448", "GSM8583898", "GSM4933442"),
  thepath = thepaths,
  thevariant = list(sc5p_pe_variant, sc5p_r2_variant, NULL, NULL, sc3pv2_variant),
  chemistry = c("SC5P-PE", "SC5P-R2", "SC5P-R2", "SC3Pv3", "SC3Pv2"),
  color = c("red", "blue", "black", "black", "green")
) |>
  dplyr::mutate(
    gsmid_label = glue::glue("{gsmid} ({chemistry})")
  ) |>
  dplyr::mutate(
    gsmid_label = factor(gsmid_label, levels = gsmid_label),
    gsmid = factor(gsmid, levels = gsmid),
    chemistry = factor(chemistry, levels = chemistry |> unique())
  )

pcc <- readr::read_tsv(file = "https://raw.githubusercontent.com/chunjie-sam-liu/chunjie-sam-liu.life/master/public/data/pcc.tsv") |>
  dplyr::arrange(cancer_types)

stats |>
  dplyr::filter(
    strand_correlation > 0.65
  ) |>
  dplyr::filter(
    vmr > 0.01
  ) |>
  dplyr::filter(mean_coverage > 10)
# body --------------------------------------------------------------------

stats <- data.table::fread("/mnt/isilon/u01_project/large-scale/liuc9/raw/GSE226602/final/GSM7080044/cell.variant_stats.tsv.gz")
all_variants <- readxl::read_excel("/home/liuc9/github/scMOCHA-data/data/GSE226602/final/GSM7080044/cell_variant_annotation.xlsx")

stats |>
  # dplyr::filter(variant %in% all_variants$variant) |>
  dplyr::mutate(
    vmr_log = log10(vmr)
  ) |>
  dplyr::mutate(
    color = ifelse(variant %in% sc5p_pe_variant, "red", "black")
  ) |>
  ggplot(aes(x = strand_correlation, y = vmr_log)) +
  geom_point(aes(color = color)) +
  scale_color_identity() +
  geom_vline(
    xintercept = 0.65,
    linetype = 20,
    color = "red"
  ) +
  geom_hline(
    yintercept = log10(0.01),
    linetype = 20,
    color = "red"
  )


# footer ------------------------------------------------------------------

# future: :plan(future: :sequential)

# save image --------------------------------------------------------------
