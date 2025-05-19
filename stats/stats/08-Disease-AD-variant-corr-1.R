#!/usr/bin/env Rscript
# Metainfo ----------------------------------------------------------------
# @AUTHOR: Chun-Jie Liu
# @CONTACT: chunjie.sam.liu.at.gmail.com
# @DATE: 2025-05-19 13:56:34
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

fn_go_enrich <- function(cancer_sgene, ont = c("BP", "CC", "MF")) {
  .ont <- match.arg(ont)

  gobp <- clusterProfiler::enrichGO(
    gene = cancer_sgene |> unique(),
    OrgDb = "org.Hs.eg.db",
    keyType = "SYMBOL",
    # keyType = "ENSEMBL",
    ont = .ont
  )

  gobp
}

fn_plot_go <- function(.go, .topn = Inf, .ont = c("BP", "CC", "MF")) {
  .ont <- match.arg(.ont)


  base_fill <- c("BP" = "#AE1700", "CC" = "#DF8F44FF", "MF" = "#00A1D5FF")
  ont_fullname <- c(
    "BP" = "Biological Process",
    "CC" = "Cellular Component",
    "MF" = "Molecular Function"
  )

  .ont_fill <- base_fill[.ont]
  x_label <- ont_fullname[.ont]

  .go %>%
    tibble::as_tibble() %>%
    dplyr::mutate(
      Description = stringr::str_wrap(
        stringr::str_to_sentence(string = Description),
        width = 60
      )
    ) %>%
    dplyr::mutate(adjp = -log10(p.adjust)) %>%
    dplyr::select(ID, Description, adjp, Count, geneID) %>%
    dplyr::arrange(adjp, Count) %>%
    dplyr::mutate(
      Description = factor(Description, levels = Description)
    ) ->
  .go_bp_for_plot

  if (!is.infinite(.topn)) {
    .go_bp_for_plot |>
      tail(.topn) ->
    .go_bp_for_plot
  }


  .go_bp_for_plot |>
    ggplot(aes(x = Description, y = adjp)) +
    geom_col(fill = .ont_fill, color = NA, width = 0.7) +
    geom_text(aes(label = Count), hjust = 4, color = "white", size = 5) +
    labs(y = "-log10(Adj. P value)", x = x_label) +
    scale_y_continuous(expand = c(0, 0.02)) +
    coord_flip() +
    theme(
      panel.background = element_rect(fill = NA),
      panel.grid = element_blank(),
      axis.line.x = element_line(color = "black"),
      axis.line.y = element_line(color = "black"),
      # axis.title.y = element_blank(),
      axis.text.y = element_text(color = "black", size = 13, hjust = 1),
      axis.ticks.length.y = unit(3, units = "mm"),
      axis.text.x = element_text(color = "black", size = 12),
      axis.title = element_text(colour = "black", size = 16, face = "bold")
    )
}



# load data ---------------------------------------------------------------
expr_v_3173G_A_corr <- import("/home/liuc9/github/scMOCHA-data/stats/stats/zzz/ad/ad-celltype-variant-af-3173G>A-corr.fst")
expr_v_1670A_G_corr <- import(
  "/home/liuc9/github/scMOCHA-data/stats/stats/zzz/ad/ad-celltype-variant-af-1670A>G-corr.fst"
)

# body --------------------------------------------------------------------


expr_v_1670A_G_corr |>
  dplyr::filter(celltype == "Mono") |>
  dplyr::filter(pval < 0.05) ->
mono_corr

mono_corr |>
  dplyr::filter(corr > 0.3) ->
mono_corr_03_pos



mono_corr_03_pos_bp <- fn_go_enrich(cancer_sgene = unique(mono_corr_03_pos$genename), "BP")
mono_corr_03_pos_cc <- fn_go_enrich(cancer_sgene = unique(mono_corr_03_pos$genename), "CC")
mono_corr_03_pos_mf <- fn_go_enrich(cancer_sgene = unique(mono_corr_03_pos$genename), "MF")

mono_corr |>
  dplyr::filter(corr < -0.3) ->
mono_corr_03_neg

mono_corr_03_neg_bp <- fn_go_enrich(cancer_sgene = unique(mono_corr_03_neg$genename), "BP")
mono_corr_03_neg_cc <- fn_go_enrich(cancer_sgene = unique(mono_corr_03_neg$genename), "CC")
mono_corr_03_neg_mf <- fn_go_enrich(cancer_sgene = unique(mono_corr_03_neg$genename), "MF")

fn_plot_go(mono_corr_03_neg_bp, 10, "BP")
fn_plot_go(mono_corr_03_neg_cc, 10, "CC")
fn_plot_go(mono_corr_03_neg_mf, 10, "MF")

# footer ------------------------------------------------------------------

# future: :plan(future: :sequential)

# save image --------------------------------------------------------------
