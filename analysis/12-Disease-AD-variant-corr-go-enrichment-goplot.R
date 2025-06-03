#!/usr/bin/env Rscript
# Metainfo ----------------------------------------------------------------
# @AUTHOR: Chun-Jie Liu
# @CONTACT: chunjie.sam.liu.at.gmail.com
# @DATE: 2025-05-22 14:55:39
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
outdir <- "/home/liuc9/github/scMOCHA-data/analysis/zzz/disease/go"
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

variants <- c(
  "1670A>G",
  "1397T>A",
  "3173G>A",
  "3176A>T",
  "3178T>A"
)
library(clusterProfiler)
variant_go_all <- import(file.path(outdir, "variant_go_all.qs"))


# body --------------------------------------------------------------------


variant_go_all |>
  dplyr::select(
    variant,
    dplyr::contains("plot")
  ) |>
  tidyr::pivot_longer(
    cols = -variant,
    names_to = "goname",
    values_to = "p"
  ) |>
  tidyr::separate(
    col = goname,
    into = c("posneg", "gotype", ".s"),
    sep = "_",
    remove = FALSE
  ) ->
forsaveplots


# ! save go for each plot --------------------------------------------------------------------

forsaveplots |>
  dplyr::mutate(
    a = parallel::mcmapply(
      FUN = \(variant, goname, p) {
        .filename <- "{variant}_{goname}.pdf" |> glue::glue()
        ggsave(
          path = outdir,
          filename = .filename,
          plot = p,
          width = 10,
          height = 6,
          dpi = 300
        )
        1
      },
      variant = variant,
      goname = goname,
      p = p,
      mc.cores = 10,
      SIMPLIFY = FALSE
    )
  )


# save plots and make plots together ---------------------------------
forsaveplots |>
  tidyr::nest(
    .by = gotype,
    .key = "go"
  ) |>
  dplyr::mutate(
    a = purrr::map2(
      .x = gotype,
      .y = go,
      .f = \(.gotype, .go) {
        # pos

        .go |>
          dplyr::filter(posneg == "pos") |>
          dplyr::pull(p) |>
          wrap_plots(ncol = 3) +
          plot_annotation(
            title = "Positive correlation with expression",
            theme = theme(plot.title = element_text(size = 20))
          ) ->
        .pos_plot
        .posout_filename <- "pos_{.gotype}_all.pdf" |> glue::glue()
        ggsave(
          path = outdir,
          filename = .posout_filename,
          plot = .pos_plot,
          width = 20,
          height = 10,
          dpi = 300
        )
        message(.gotype)
        .go |>
          dplyr::filter(posneg == "neg") |>
          dplyr::pull(p) |>
          wrap_plots(ncol = 3) +
          plot_annotation(
            title = "Negative correlation with expression",
            theme = theme(plot.title = element_text(size = 20))
          ) ->
        .neg_plot
        .posneg_filename <- "neg_{.gotype}_all.pdf" |> glue::glue()
        ggsave(
          path = outdir,
          filename = .posneg_filename,
          plot = .neg_plot,
          width = 20,
          height = 10,
          dpi = 300
        )
      }
    )
  )



# ! don't run below --------------------------------------------------------------------


#
#
#
#
#
#
#
#
#
# ! 3173G>A--------------------------------------------------------------------

variant_go_all |>
  dplyr::filter(
    variant == "3173G>A"
  ) |>
  dplyr::pull(neg_cc) |>
  as.data.frame() |>
  dplyr::filter(
    grepl("mitochondrial", Description, ignore.case = TRUE)
  ) |>
  dplyr::pull(geneID) |>
  stringr::str_split("/") |>
  unlist() |>
  sort() |>
  unique() ->
genes_pathway

expr <- import("/mnt/isilon/u01_project/large-scale/liuc9/raw/zzz/db/EXPR/gse_srrid_celltype_gene_expr.qs") |>
  dplyr::filter(celltype == "Mono")

v_3173G_A <- import("/home/liuc9/github/scMOCHA-data/analysis/zzz/ad/ad-celltype-variant-af-3173G>A.qs") |>
  dplyr::filter(celltype == "Mono")

source("/home/liuc9/github/scMOCHA-data/analysis/00-colors.R")

import(
  "/home/liuc9/github/scMOCHA-data/analysis/zzz/clean-data/gse_dataset_metadata_full.qs"
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

fn_load_corr <- function(.variant) {
  import(
    "/home/liuc9/github/scMOCHA-data/analysis/zzz/ad/ad-celltype-variant-af-{.variant}-corr.fst" |> glue::glue()
  ) |>
    dplyr::filter(celltype == "Mono") |>
    dplyr::filter(pval < 0.05) |>
    dplyr::arrange(desc(corr)) |>
    dplyr::filter(abs(corr) > 0.3)
}

import(
  "/home/liuc9/github/scMOCHA-data/analysis/zzz/ad/ad-celltype-variant-af-3173G>A-corr.fst" |> glue::glue()
) |>
  dplyr::arrange(corr) |>
  dplyr::filter(pval < 0.05) |>
  dplyr::filter(abs(corr) > 0.3) |>
  dplyr::filter(!celltype %in% c("other_T", "other", "DC")) |>
  dplyr::filter(!grepl("ENSG", genename))

corr_3173G_A <- fn_load_corr("3173G>A") |>
  dplyr::filter(genename %in% c(genes_pathway)) |>
  dplyr::arrange(desc(corr))

expr |>
  dplyr::left_join(
    v_3173G_A,
    by = c("celltype")
  ) |>
  dplyr::left_join(
    corr_3173G_A,
    by = c("genename", "celltype")
  ) |>
  dplyr::arrange(desc(corr)) |>
  dplyr::filter(genename %in% genes_pathway) ->
expr_v_3173G_A

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

expr_v_3173G_A |>
  dplyr::mutate(
    corr_plot = purrr::pmap(
      list(
        .expr = expr,
        .af = af,
        .genename = genename
      ),
      .f = \(.expr, .af, .genename) {
        # .expr <- expr_v_3173G_A$expr[[1]]
        # .af <- expr_v_3173G_A$af[[1]]
        # .genename <- expr_v_3173G_A$genename[[1]]
        # cor_test <- cor(.expr, .af, method = "pearson")

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
                  .x = human_read(cor_test$p.value),
                  .s = glue::glue("R={round(cor_test$estimate,3)}")
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

expr_v_3173G_A_plot |>
  dplyr::mutate(
    a = purrr::map2(
      .x = genename,
      .y = corr_plot,
      .f = ~ {
        .filename <- glue::glue("{.x}.pdf")
        ggsave(
          filename = .filename,
          path = "/home/liuc9/github/scMOCHA-data/analysis/zzz/disease/corr/3173G_A",
          plot = .y,
          width = 7,
          height = 6,
          dpi = 300
        )
      }
    )
  )

expr_v_3173G_A |>
  dplyr::arrange(corr) |>
  dplyr::select(genename, corr, pval) |>
  export(
    file = "/home/liuc9/github/scMOCHA-data/analysis/zzz/disease/corr/3173G_A/3173G_A_corr.csv"
  )

expr_v_3173G_A |>
  dplyr::arrange(corr) |>
  dplyr::select(genename, corr, pval) |>
  dplyr::slice(1:25) |>
  # dplyr::filter(corr < -0.35) |>
  dplyr::mutate(
    variant = "3173G>A",
  ) ->
expr_v_3173G_A_top10

variants |>
  purrr::map(
    ~ {
      import(
        "/home/liuc9/github/scMOCHA-data/analysis/zzz/ad/ad-celltype-variant-af-{.x}-corr.fst" |> glue::glue()
      ) |>
        dplyr::filter(genename %in% expr_v_3173G_A_top10$genename) |>
        dplyr::mutate(variant = .x)
    }
  ) |>
  dplyr::bind_rows() ->
variant_corr_top10

variant_corr_top10 |>
  dplyr::filter(variant == "3173G>A") |>
  dplyr::mutate(
    celltype = gsub("_", " ", celltype),
  ) |>
  dplyr::mutate(
    celltype = factor(celltype, levels = names(color_celltype) |> rev())
  ) |>
  dplyr::mutate(
    genename = factor(genename, levels = unique(expr_v_3173G_A_top10$genename))
  ) |>
  dplyr::mutate(
    mark = dplyr::case_when(
      pval < 0.001 ~ "***",
      pval < 0.01 ~ "**",
      pval < 0.05 ~ "*",
      TRUE ~ ""
    )
  ) ->
variant_corr_top10_3173G_A_forplot

variant_corr_top10_3173G_A_forplot |>
  ggplot(aes(
    x = genename,
    y = celltype
  )) +
  geom_tile(aes(fill = corr)) +
  geom_text(
    aes(label = mark),
    size = 5,
    color = "black"
  ) +
  scale_fill_gradient2(
    breaks = round(seq(-0.35, 0.4, length.out = 5), digits = 2),
    labels = format(seq(-0.4, 0.4, length.out = 5), digits = 2),
    low = "#00fefe",
    mid = "white",
    high = "#fe0000"
  ) +
  theme_cor() +
  theme(
    axis.text.y = element_text(hjust = 1, size = 14, face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, face = "bold"),
    axis.line = element_blank(),
    axis.ticks = element_blank(),
    legend.position = "right",
    axis.title = element_blank()
  ) +
  guides(
    fill = guide_legend(
      # legend title
      title = "Pearsons' correlation (r)",
      title.position = "right",
      title.theme = element_text(angle = -90, size = 14, family = "Times"),
      title.vjust = -0.5,
      title.hjust = 0.5,

      # legend label
      label = TRUE,
      label.position = "left",
      label.theme = element_text(size = 14, family = "Times"),
      label.hjust = 0.5,
      label.vjust = 0.5,

      # legend key
      keywidth = 1,
      keyheight = 1.8,
      reverse = TRUE
    )
  ) +
  coord_fixed(ratio = 1) +
  labs(
    x = "Gene",
    y = "Cell type"
  ) ->
variant_corr_top10_3173G_A_plot

ggsave(
  path = "/home/liuc9/github/scMOCHA-data/analysis/zzz/disease/corr/",
  filename = file.path("variant_corr_top10_3173G_A_celltype.pdf"),
  plot = variant_corr_top10_3173G_A_plot,
  width = 12,
  height = 5,
  dpi = 300
)


variant_corr_top10 |>
  # dplyr::filter(variant == "3173G>A") |>
  dplyr::filter(celltype == "Mono") |>
  dplyr::mutate(
    celltype = gsub("_", " ", celltype),
  ) |>
  dplyr::mutate(
    celltype = factor(celltype, levels = names(color_celltype) |> rev())
  ) |>
  dplyr::mutate(
    genename = factor(genename, levels = unique(expr_v_3173G_A_top10$genename))
  ) |>
  dplyr::mutate(
    mark = dplyr::case_when(
      pval < 0.001 ~ "***",
      pval < 0.01 ~ "**",
      pval < 0.05 ~ "*",
      TRUE ~ ""
    )
  ) |>
  ggplot(aes(
    x = genename,
    y = variant
  )) +
  geom_tile(aes(fill = corr)) +
  geom_text(
    aes(label = mark),
    size = 5,
    color = "black"
  ) +
  scale_fill_gradient2(
    breaks = round(seq(-0.35, 0.4, length.out = 5), digits = 2),
    labels = format(seq(-0.4, 0.4, length.out = 5), digits = 2),
    low = "#00fefe",
    mid = "white",
    high = "#fe0000"
  ) +
  theme_cor() +
  theme(
    axis.text.y = element_text(hjust = 1, size = 14, face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, face = "bold"),
    axis.line = element_blank(),
    axis.ticks = element_blank(),
    legend.position = "right",
    axis.title = element_blank()
  ) +
  guides(
    fill = guide_legend(
      # legend title
      title = "Pearsons' correlation (r)",
      title.position = "right",
      title.theme = element_text(angle = -90, size = 14, family = "Times"),
      title.vjust = -0.5,
      title.hjust = 0.5,

      # legend label
      label = TRUE,
      label.position = "left",
      label.theme = element_text(size = 14, family = "Times"),
      label.hjust = 0.5,
      label.vjust = 0.5,

      # legend key
      keywidth = 1,
      keyheight = 1.8,
      reverse = TRUE
    )
  ) +
  coord_fixed(ratio = 1) +
  labs(
    x = "Gene",
    y = "Cell type"
  ) ->
variant_corr_top10_5variant_plot

ggsave(
  path = "/home/liuc9/github/scMOCHA-data/analysis/zzz/disease/corr/",
  filename = file.path("variant_corr_top10_5_variant.pdf"),
  plot = variant_corr_top10_5variant_plot,
  width = 12,
  height = 5,
  dpi = 300
)

# footer ------------------------------------------------------------------

# future: :plan(future: :sequential)

# save image --------------------------------------------------------------
