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
outdir <- "/home/liuc9/github/scMOCHA-data/stats/stats/zzz/disease"

variants <- c(
  "1670A>G",
  "1397T>A",
  "3173G>A",
  "3176A>T",
  "3178T>A"
)




# ! body --------------------------------------------------------------------
# fn_variant_go("3173G>A") -> a
variants |>
  purrr::map_dfr(
    fn_variant_go,
    .id = "variant"
  ) ->
variant_go

variant_go |>
  dplyr::bind_rows() |>
  dplyr::mutate(
    variant = variants
  ) ->
variant_go_all


# export(
#   variant_go_all,
#   file = file.path(outdir, "variant_go_all.qs")
# )

variant_go_all |>
  dplyr::mutate(
    a = purrr::map2(
      .x = variant,
      .y = pos_bp,
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
pos_bp_plot_all

names(pos_bp_plot_all) <- variant_go_all$variant

pos_bp_plot_all |>
  wrap_plots() +
  plot_annotation(
    title = "Positive correlation with expression",
    theme = theme(plot.title = element_text(size = 20))
  ) ->
pos_bp_plot_all_patch

ggsave(
  path = outdir,
  filename = file.path("pos_bp_plot_all.pdf"),
  plot = pos_bp_plot_all_patch,
  width = 20,
  height = 10,
  dpi = 300
)

variants |>
  purrr::map(
    ~ {
      ggsave(
        path = outdir,
        filename = file.path("pos_bp_{.x}.pdf" |> glue::glue()),
        plot = pos_bp_plot_all[[.x]],
        width = 10,
        height = 7,
        dpi = 300
      )
    }
  )



# ! 3173G>A --------------------------------------------------------------------
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

mt_dna_repair_genes <- c(
  "OGG1",
  "MUTYH",
  "NEIL1", "NEIL2", "NTHL1",
  "APEX1",
  "POLG",
  "LIG3",
  "TWINKLE",
  "TFAM",
  "SIRT3", "SIRT4", "SIRT5",
  "PRDX3", "PRDX5", "GPX1",
  "MPV17",
  "TWNK",
  "TYMP"
)


variant_go_all |>
  dplyr::filter(variant == "3173G>A") |>
  dplyr::pull(pos_bp) |>
  as.data.frame() |>
  dplyr::filter(grepl("dna repair", Description, ignore.case = TRUE)) |>
  dplyr::pull(geneID) |>
  stringr::str_split("/") |>
  unlist() |>
  sort() |>
  unique() ->
dna_repair_genes_pathway
dna_repair_genes <- c(dna_repair_genes_pathway, mt_dna_repair_genes)


expr <- import("/mnt/isilon/u01_project/large-scale/liuc9/raw/zzz/db/EXPR/gse_srrid_celltype_gene_expr.qs") |>
  dplyr::filter(celltype == "Mono")

v_3173G_A <- import("/home/liuc9/github/scMOCHA-data/stats/stats/zzz/ad/ad-celltype-variant-af-3173G>A.qs") |>
  dplyr::filter(celltype == "Mono")

source("/home/liuc9/github/scMOCHA-data/stats/stats/00-colors.R")

color_disease

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


corr_3173G_A <- fn_load_corr("3173G>A") |>
  dplyr::filter(genename %in% c(dna_repair_genes)) |>
  dplyr::arrange(desc(corr))

corr_3173G_A |>
  dplyr::filter(genename %in% mt_dna_repair_genes)

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
  dplyr::filter(genename %in% dna_repair_genes) ->
expr_v_3173G_A

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
          path = .path,
          plot = .y,
          width = 7,
          height = 6,
          dpi = 300
        )
      }
    )
  )

expr_v_3173G_A |>
  dplyr::select(genename, corr, pval) |>
  export(
    file = "/home/liuc9/github/scMOCHA-data/stats/stats/zzz/disease/corr/3173G_A/3173G_A_corr.csv"
  )


expr_v_3173G_A |>
  dplyr::select(genename, corr, pval) |>
  dplyr::slice(1:10) |>
  dplyr::mutate(
    variant = "3173G>A",
  ) ->
expr_v_3173G_A_top10

variants |>
  purrr::map(
    ~ {
      import(
        "/home/liuc9/github/scMOCHA-data/stats/stats/zzz/ad/ad-celltype-variant-af-{.x}-corr.fst" |> glue::glue()
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
  path = outdir,
  filename = file.path("variant_corr_top10_3173G_A_celltype.pdf"),
  plot = variant_corr_top10_3173G_A_plot,
  width = 8,
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
variant_corr_top10_5_variant_plot

ggsave(
  path = outdir,
  filename = file.path("variant_corr_top10_5_variant.pdf"),
  plot = variant_corr_top10_5_variant_plot,
  width = 8,
  height = 5,
  dpi = 300
)
# footer ------------------------------------------------------------------

# future: :plan(future: :sequential)

# save image --------------------------------------------------------------
