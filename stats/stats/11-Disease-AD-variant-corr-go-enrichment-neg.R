#!/usr/bin/env Rscript
# Metainfo ----------------------------------------------------------------
# @AUTHOR: Chun-Jie Liu
# @CONTACT: chunjie.sam.liu.at.gmail.com
# @DATE: 2025-05-21 17:58:36
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
    geom_text(aes(label = Count), hjust = 1, color = "white", size = 5) +
    labs(y = "-log10(Adj. P value)", x = x_label) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
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

fn_load_corr <- function(.variant) {
  import(
    "/home/liuc9/github/scMOCHA-data/stats/stats/zzz/ad/ad-celltype-variant-af-{.variant}-corr.fst" |> glue::glue()
  ) |>
    dplyr::filter(celltype == "Mono") |>
    dplyr::filter(pval < 0.05) |>
    dplyr::arrange(desc(corr)) |>
    dplyr::filter(abs(corr) > 0.3)
}

fn_variant_go <- function(.variant) {
  # .variant <- variants[[1]]
  .corr <- fn_load_corr(.variant)
  .corr_pos <- .corr |> dplyr::filter(corr > 0.3)
  .corr_neg <- .corr |> dplyr::filter(corr < -0.3)

  .pos_bp <- fn_go_enrich(cancer_sgene = unique(.corr_pos$genename), "BP")
  .pos_cc <- fn_go_enrich(cancer_sgene = unique(.corr_pos$genename), "CC")
  .pos_mf <- fn_go_enrich(cancer_sgene = unique(.corr_pos$genename), "MF")

  .neg_bp <- fn_go_enrich(cancer_sgene = unique(.corr_neg$genename), "BP")
  .neg_cc <- fn_go_enrich(cancer_sgene = unique(.corr_neg$genename), "CC")
  .neg_mf <- fn_go_enrich(cancer_sgene = unique(.corr_neg$genename), "MF")

  tibble::tibble(
    variant = .variant,
    pos_bp = list(.pos_bp),
    pos_cc = list(.pos_cc),
    pos_mf = list(.pos_mf),
    neg_bp = list(.neg_bp),
    neg_cc = list(.neg_cc),
    neg_mf = list(.neg_mf),
    pos_bp_plot = list(fn_plot_go(.pos_bp, 20, "BP")),
    pos_cc_plot = list(fn_plot_go(.pos_cc, 20, "CC")),
    pos_mf_plot = list(fn_plot_go(.pos_mf, 20, "MF")),
    neg_bp_plot = list(fn_plot_go(.neg_bp, 20, "BP")),
    neg_cc_plot = list(fn_plot_go(.neg_cc, 20, "CC")),
    neg_mf_plot = list(fn_plot_go(.neg_mf, 20, "MF"))
  )
}

# load data ---------------------------------------------------------------


# body --------------------------------------------------------------------
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
    # legend.negition = "bottom",
    # legend.negition = "none",
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

outdir <- "/home/liuc9/github/scMOCHA-data/stats/stats/zzz/disease"
variants <- c(
  "1670A>G",
  "1397T>A",
  "3173G>A",
  "3176A>T",
  "3178T>A"
)
library(clusterProfiler)
variant_go_all <- import(
  file = file.path(outdir, "variant_go_all.qs")
)

variant_go_all |>
  dplyr::mutate(
    a = purrr::map2(
      .x = variant,
      .y = neg_bp,
      .f = ~ {
        fn_plot_go(.y, 20, "BP") +
          labs(title = .x) +
          theme(
            plot.title = element_text(size = 20)
          )
      }
    )
  ) |>
  dplyr::pull(a) ->
neg_bp_plot_all

names(neg_bp_plot_all) <- variant_go_all$variant

neg_bp_plot_all |>
  wrap_plots() +
  plot_annotation(
    title = "negitive correlation with expression",
    theme = theme(plot.title = element_text(size = 20))
  ) ->
neg_bp_plot_all_patch

ggsave(
  path = outdir,
  filename = file.path("neg_bp_plot_all.pdf"),
  plot = neg_bp_plot_all_patch,
  width = 20,
  height = 10,
  dpi = 300
)

variants |>
  purrr::map(
    ~ {
      ggsave(
        path = outdir,
        filename = file.path("neg_bp_{.x}.pdf" |> glue::glue()),
        plot = neg_bp_plot_all[[.x]],
        width = 10,
        height = 7,
        dpi = 300
      )
    }
  )


# footer ------------------------------------------------------------------

# future: :plan(future: :sequential)

# save image --------------------------------------------------------------
