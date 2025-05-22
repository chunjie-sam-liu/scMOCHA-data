#!/usr/bin/env Rscript
# Metainfo ----------------------------------------------------------------
# @AUTHOR: Chun-Jie Liu
# @CONTACT: chunjie.sam.liu.at.gmail.com
# @DATE: 2025-05-22 17:46:41
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

theme_cor <- function() {
  theme( # plot.margin = unit(c(0.5,0.5,0.5,0.5), "cm"),
    plot.title = element_text(size = rel(1.3), vjust = 2, hjust = 0.5, lineheight = 0.8),

    # axis
    axis.title.x = element_text(face = "bold", size = 16),
    axis.title.y = element_text(face = "bold", size = 16, angle = 90),
    axis.text = element_text(size = rel(1.1)),
    axis.text.x = element_text(hjust = 0.5, vjust = 0, size = 14),
    axis.text.y = element_text(vjust = 0.5, hjust = 0, size = 14),
    axis.line = element_line(colour = "black"),

    # ticks
    axis.ticks = element_line(colour = "black"),

    # legend
    legend.title = element_text(size = rel(1.1), face = "bold"),
    legend.text = element_text(size = rel(1.1), face = "bold"),
    # legend.position = "bottom",
    # legend.position = "none",
    # legend.direction = "horizontal",
    legend.background = element_blank(),
    legend.key = element_rect(fill = NA, colour = NA),

    # strip
    strip.text = element_text(size = rel(1.3)),

    # panel
    panel.background = element_blank(),
    # aspect.ratio = 1,
    complete = T
  )
}

# load data ---------------------------------------------------------------

outdir <- "/home/liuc9/github/scMOCHA-data/stats/stats/zzz/disease/go"
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
library(clusterProfiler)
variant_go_all <- import(file.path(outdir, "variant_go_all.qs"))

expr <- import("/mnt/isilon/u01_project/large-scale/liuc9/raw/zzz/db/EXPR/gse_srrid_celltype_gene_expr.qs") |>
  dplyr::filter(celltype == "Mono")

v_3173G_A <- import("/home/liuc9/github/scMOCHA-data/stats/stats/zzz/ad/ad-celltype-variant-af-3173G>A.qs") |>
  dplyr::filter(celltype == "Mono")

import(
  "/home/liuc9/github/scMOCHA-data/stats/stats/zzz/clean-data/gse_dataset_metadata_full.qs"
) |>
  dplyr::select(srrid, disease, Chemistry) |>
  dplyr::filter(
    Chemistry == "SC5P-PE"
  ) |>
  dplyr::filter(
    disease %in% c(
      "Alzheimer's Disease",
      "Healthy"
    )
  ) ->
metadata

# body --------------------------------------------------------------------

fn_load_corr <- function(.variant) {
  import(
    "/home/liuc9/github/scMOCHA-data/stats/stats/zzz/ad/ad-celltype-variant-af-{.variant}-corr.fst" |> glue::glue()
  ) |>
    dplyr::filter(celltype == "Mono") |>
    dplyr::filter(pval < 0.05) |>
    dplyr::arrange(desc(corr)) |>
    dplyr::filter(abs(corr) > 0.3)
}

corr_3173G_A <- fn_load_corr("3173G>A") |>
  dplyr::arrange(desc(corr))


expr |>
  dplyr::left_join(
    v_3173G_A,
    by = c("celltype")
  ) |>
  dplyr::inner_join(
    corr_3173G_A,
    by = c("genename", "celltype")
  ) |>
  dplyr::arrange(corr) |>
  dplyr::filter(!grepl("ENSG0", genename)) ->
expr_v_3173G_A


source("/home/liuc9/github/scMOCHA-data/stats/stats/00-colors.R")
expr_v_3173G_A |>
  dplyr::mutate(
    corr_plot = purrr::pmap(
      list(
        .expr = expr,
        .af = af,
        .genename = genename,
        .corr = corr,
        .pval = pval
      ),
      .f = \(.expr, .af, .genename, .corr, .pval) {
        # .expr <- expr_v_3173G_A$expr[[1]]
        # .af <- expr_v_3173G_A$af[[1]]
        # .genename <- expr_v_3173G_A$genename[[1]]
        # .corr <- expr_v_3173G_A$corr[[1]]
        # .pval <- expr_v_3173G_A$pval[[1]]

        .expr |>
          dplyr::left_join(
            .af,
            by = c("gseid", "srrid")
          ) |>
          dplyr::inner_join(
            metadata,
            by = c("srrid")
          ) ->
        .expr_af
        cor_test <- cor.test(~ expr + af, data = .expr_af, method = "pearson")

        tryCatch(
          {
            .expr_af |>
              ggplot(aes(
                x = af,
                y = expr
              )) +
              geom_point(
                aes(
                  color = disease,
                ),
                position = position_jitter(width = 0.1, height = 0.1),
              ) +
              scale_color_manual(
                values = color_disease,
                name = "Disease",
              ) +
              geom_smooth(
                color = "red",
                method = "lm"
              ) +
              scale_x_continuous(
                # limits = c(0, 96),
                # breaks = seq(0, 96, by = 10),
                # labels = seq(0, 96, by = 10),
                expand = expansion(mult = c(0.01, 0.01))
              ) +
              scale_y_continuous(
                # limits = c(0, 90),
                # breaks = seq(0, 90, by = 10),
                # labels = seq(0, 90, by = 10),
                expand = expansion(mult = c(0.01, 0))
              ) +
              theme_cor() +
              theme(
                # plot.title = element_blank(),
                legend.position = "bottom"
              ) +
              labs(
                title = human_read_latex_pval(
                  .x = human_read(.pval),
                  .s = glue::glue("R={round(.corr,3)}")
                ),
                x = "3173G>A Allele frequency",
                y = "{.genename} normalized gene expression" |> glue::glue(),
              )
          },
          error = function(e) {
            # message("Error in cor.test: ", e)
            return(NULL)
          }
        )
      }
    )
  ) ->
expr_v_3173G_A_plot

# footer ------------------------------------------------------------------

# future: :plan(future: :sequential)

# save image --------------------------------------------------------------
