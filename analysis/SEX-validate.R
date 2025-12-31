#!/usr/bin/env Rscript
# Metainfo ----------------------------------------------------------------
# @AUTHOR: Chun-Jie Liu
# @CONTACT: chunjie.sam.liu.at.gmail.com
# @DATE: 2025-05-15 22:01:45
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

# load data ---------------------------------------------------------------
sex_pred <- import(
  "/home/liuc9/github/scMOCHA-data/analysis/zzz/clean-data/gse_srrid_srrdir_sex.qs"
) |>
  dplyr::select(
    gseid,
    srrid,
    srrdir,
    sex_pred = sex
  )
sex_real <- import(
  "/home/liuc9/github/scMOCHA-data/analysis/zzz/clean-data/gse_dataset_metadata_full.qs"
) |>
  dplyr::select(gseid, srrid, sex_real = Gender)

# body --------------------------------------------------------------------

sex_pred |>
  dplyr::mutate(
    sex_pred = factor(sex_pred, levels = c("Male", "Female"))
  ) |>
  dplyr::left_join(sex_real, by = c("gseid", "srrid")) -> sex_pred_real


sex_pred_real |>
  dplyr::count(sex_pred) |>
  dplyr::mutate(sex_pred_str = glue::glue("{sex_pred}\n(n={n})")) %>%
  dplyr::select(-n) |>
  dplyr::mutate(
    sex_pred_str = factor(sex_pred_str, levels = sex_pred_str)
  ) -> sex_pred_str

sex_pred_real |>
  dplyr::count(sex_real) |>
  dplyr::mutate(sex_real_str = glue::glue("{sex_real}\n(n={n})")) %>%
  dplyr::select(-n) |>
  dplyr::mutate(
    sex_real_str = factor(sex_real_str, levels = sex_real_str)
  ) -> sex_real_str

sex_pred_real |>
  dplyr::select(sex_pred, sex_real) |>
  dplyr::count(sex_pred, sex_real) |>
  dplyr::left_join(sex_pred_str, by = "sex_pred") |>
  dplyr::left_join(sex_real_str, by = "sex_real") -> for_sankey_plot

source("/home/liuc9/github/scMOCHA-data/analysis/00-colors.R")
library(ggalluvial)
for_sankey_plot |>
  dplyr::mutate(sex_pred = factor(sex_pred, levels = c("Female", "Male"))) |>
  ggplot(aes(
    axis1 = sex_real_str,
    axis2 = sex_pred_str,
    y = n
  )) +
  ggalluvial::geom_alluvium(
    aes(fill = sex_real),
    width = 1 / 12,
    alpha = 0.8,
  ) +
  ggalluvial::geom_stratum(
    # aes(fill = Chemistry),
    width = 0.3
  ) +
  scale_fill_manual(
    name = "Sex",
    values = color_gender
  ) +
  geom_text(
    stat = "stratum",
    aes(label = after_stat(stratum))
  ) +
  scale_x_discrete(
    limits = c("sex_real_str", "sex_pred_str"),
    labels = c("Sex (Real)", "Sex (Predicted)")
  ) +
  scale_y_continuous(
    expand = expand_scale(mult = c(0.01, 0.02))
  ) +
  # ggsci::scale_fill_aaas() +
  theme(
    axis.text.y = element_blank(),
    axis.text.x = element_text(color = "black", size = 12),
    axis.ticks = element_blank(),
    axis.title = element_blank(),
    panel.background = element_rect(fill = NA, colour = NA),
    legend.position = "top"
  ) +
  guides(fill = guide_legend(title = "Sex", nrow = 1)) -> sex_pred_real_plot

ggsave(
  filename = "SEX-validate.pdf",
  plot = sex_pred_real_plot,
  device = "pdf",
  path = "/home/liuc9/github/scMOCHA-data/analysis/zzz/",
  height = 4,
  width = 6,
)


sex_pred_real |>
  dplyr::filter(sex_real != "Unknown") |>
  dplyr::mutate(
    sex_pred = as.character(sex_pred),
    sex_real = as.character(sex_real)
  ) |>
  dplyr::filter(sex_pred != sex_real) |>
  dplyr::arrange(sex_real)

# footer ------------------------------------------------------------------

# future: :plan(future: :sequential)

# save image --------------------------------------------------------------
