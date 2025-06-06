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

outdir <- "/home/liuc9/github/scMOCHA-data/analysis/zzz/plot-ad/go"
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

library(clusterProfiler)
variant_go_all <- import(file.path(outdir, "variant_go_all.qs"))

expr <- import("/mnt/isilon/u01_project/large-scale/liuc9/raw/zzz/db/EXPR/gse_srrid_celltype_gene_expr.qs") |>
  dplyr::filter(celltype == "Mono")

v_3173G_A <- import("/home/liuc9/github/scMOCHA-data/analysis/zzz/ad/ad-celltype-variant-af-3173G>A.qs") |>
  dplyr::filter(celltype == "Mono")

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

# body --------------------------------------------------------------------

fn_load_corr <- function(.variant) {
  import(
    "/home/liuc9/github/scMOCHA-data/analysis/zzz/ad/ad-celltype-variant-af-{.variant}-corr.fst" |> glue::glue()
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



# ? plot scatter and correlation --------------------------------------------------------------------

source("/home/liuc9/github/scMOCHA-data/analysis/00-colors.R")

expr_v_3173G_A |>
  dplyr::mutate(
    corr_plot = parallel::mcmapply(
      FUN = \(.expr, .af, .genename, .corr, .pval) {
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
      },
      .expr = expr,
      .af = af,
      .genename = genename,
      .corr = corr,
      .pval = pval,
      mc.cores = 10,
      SIMPLIFY = FALSE
    )
  ) ->
expr_v_3173G_A_plot



# ? save plots --------------------------------------------------------------------

expr_v_3173G_A_plot |>
  dplyr::mutate(
    a = parallel::mcmapply(
      FUN = \(.x, .y) {
        .filename <- glue::glue("{.x}.pdf")
        ggsave(
          filename = .filename,
          path = "/home/liuc9/github/scMOCHA-data/analysis/zzz/plot-ad/corr/3173G_A",
          plot = .y,
          width = 7,
          height = 6,
          dpi = 300
        )
      },
      .x = genename,
      .y = corr_plot,
      mc.cores = 10,
      SIMPLIFY = FALSE
    )
  )

expr_v_3173G_A_plot |>
  export(
    file = "/home/liuc9/github/scMOCHA-data/analysis/zzz/plot-ad/corr/3173G_A/3173G_A_corr.qs"
  )




# ? all correlation with 3173G>A --------------------------------------------------------------------


genebiotype <- import("/home/liuc9/github/scMOCHA-data/config/Homo_sapiens.GRCh38.107.gtf.id_name_length_genetype.fst")

expr_v_3173G_A_plot |>
  dplyr::filter(AD >= 20, Healthy >= 20) |>
  dplyr::left_join(
    genebiotype,
    by = c("genename" = "gene_name")
  ) ->
expr_v_3173G_A_plot_filtered


expr_v_3173G_A_plot_filtered |>
  tibble::rowid_to_column() |>
  dplyr::select(-c(expr, af)) |>
  dplyr::mutate(
    label = ifelse(
      abs(corr) > 0.53 & Gene_type == "protein_coding",
      genename,
      NA_character_
    )
  ) |>
  ggplot(aes(
    x = rowid,
    y = corr
  )) +
  geom_col() +
  ggrepel::geom_text_repel(
    aes(label = label)
  ) +
  scale_x_continuous(
    limits = c(0, 1400),
    breaks = seq(0, 1400, 200),
    expand = expansion(mult = c(0.01, 0.02))
  ) +
  scale_y_continuous(
    limits = c(-0.7, 0.7),
    breaks = seq(-0.7, 0.7, 0.2),
    expand = expansion(mult = c(0.02, 0.01))
  ) +
  theme_cor() +
  theme(
    # axis.title.x = element_blank()
  ) +
  labs(
    x = "Gene rank",
    y = "Corr. between gene expr and 3173G>A AF in Mono"
  ) ->
plot_corr_3173G_A

ggsave(
  filename = "corr_gene_3173G_A_in_mono.pdf",
  path = "/home/liuc9/github/scMOCHA-data/analysis/zzz/plot-ad/corr/",
  plot = plot_corr_3173G_A,
  width = 10,
  height = 6,
  dpi = 300
)




# ? all variant and celltype --------------------------------------------------------------------

variants <- c(
  "1670A>G",
  "1397T>A",
  "3173G>A",
  "3176A>T",
  "3178T>A"
)

variants |>
  purrr::map(
    ~ {
      import(
        "/home/liuc9/github/scMOCHA-data/analysis/zzz/ad/ad-celltype-variant-af-{.x}-corr.fst" |> glue::glue()
      ) |>
        dplyr::mutate(variant = .x)
    }
  ) |>
  dplyr::bind_rows() ->
variant_corr_celltype

expr_v_3173G_A_plot_filtered |>
  tibble::rowid_to_column() |>
  dplyr::select(-c(expr, af)) |>
  dplyr::mutate(
    label = ifelse(
      abs(corr) > 0.5 & Gene_type == "protein_coding",
      genename,
      NA_character_
    )
  ) |>
  dplyr::filter(!is.na(label)) ->
topcorrgenes

topcorrgenes |> dplyr::pull(genename)

variant_corr_celltype |>
  dplyr::filter(genename %in% topcorrgenes$genename) |>
  dplyr::filter(variant == "3173G>A") |>
  dplyr::mutate(
    celltype = gsub("_", " ", celltype),
  ) |>
  dplyr::mutate(
    celltype = factor(celltype, levels = names(color_celltype) |> rev())
  ) |>
  dplyr::mutate(
    genename = factor(genename, levels = topcorrgenes$genename)
  ) |>
  dplyr::mutate(
    mark = dplyr::case_when(
      pval < 0.001 ~ "***",
      pval < 0.01 ~ "**",
      pval < 0.05 ~ "*",
      TRUE ~ ""
    )
  ) ->
variant_topcorrgenes_3173G_A_forplot

variant_topcorrgenes_3173G_A_forplot |>
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
    breaks = round(seq(-0.6, 0.55, length.out = 5), digits = 2),
    labels = format(seq(-0.6, 0.6, length.out = 5), digits = 2),
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
variant_topcorrgenes_3173G_A_tileplot

ggsave(
  path = "/home/liuc9/github/scMOCHA-data/analysis/zzz/plot-ad/corr/",
  filename = file.path("variant_topcorrgenes_3173G_A_tileplot.pdf"),
  plot = variant_topcorrgenes_3173G_A_tileplot,
  width = 12,
  height = 5,
  dpi = 300
)


# ? mito genes --------------------------------------------------------------------


# variant_go_all |>
#   dplyr::filter(
#     variant == "3173G>A"
#   ) |>
#   dplyr::pull(neg_cc)


# footer ------------------------------------------------------------------

# future: :plan(future: :sequential)

# save image --------------------------------------------------------------
