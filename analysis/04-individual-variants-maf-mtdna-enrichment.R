#!/usr/bin/env Rscript
# Metainfo ----------------------------------------------------------------
# @AUTHOR: Chun-Jie Liu
# @CONTACT: chunjie.sam.liu.at.gmail.com
# @DATE: 2025-04-22 14:48:15
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

fn_plot_mtdna <- function() {
  # mt_exons_df <- "/home/liuc9/github/scMOCHA/fasta/mt_exons.df.rds.gz"

  LENGTH <- 16569
  # rCRS <- Biostrings::readDNAStringSet("/home/liuc9/github/scMOCHA-data/config/rCRS.MT.fasta")
  gtf_gene_df <- readr::read_rds("/home/liuc9/github/scMOCHA-data/config/mtdna_genes_dloop.rds.gz")


  library(gggenes)
  ggplot(
    gtf_gene_df,
    aes(
      xmin = start,
      xmax = end,
      y = seqnames,
    )
  ) +
    # geom_gene_arrow() +
    geom_gene_arrow(
      aes(
        fill = COLOR
      ),
      arrowhead_height = unit(3, "mm"), arrowhead_width = unit(1, "mm"),
    ) +
    scale_fill_identity(
      name = "Gene type",
      guide = "legend",
      labels = c("MT rRNA", "Protein coding", "MT tRNA", "MT OLR", "D-Loop")
    ) +
    # scale_fill_brewer(
    #   palette = "Set1",
    #   name = "Gene type",
    #   labels = c("D-Loop", "MT rRNA", "MT tRNA", "Protein coding")
    # ) +
    ggrepel::geom_text_repel(
      aes(
        x = (start + end) / 2,
        label = gsub(
          pattern = "MT-",
          replacement = "",
          x = gene_name
        ),
      ),
      color = "black",
      # fill = "white",
      # nudge_x =1,
      # nudge_y =0.001,
      size = 3,
      show.legend = F,
      max.overlaps = Inf,
    ) +
    # scale_color_brewer(palette = "Set1") +
    scale_x_continuous(
      limits = c(0, LENGTH),
      breaks = c(seq(0, LENGTH, 1000), LENGTH),
      labels = c(seq(0, LENGTH, 1000), LENGTH),
      expand = expansion(mult = c(0, 0.01)),
    ) +
    scale_y_discrete(
      expand = expansion(mult = c(0, 0), add = c(0, 0))
    ) +
    # theme_genes() +
    theme(
      legend.position = "bottom",
      axis.title = element_blank(),
      axis.text.y = element_blank(),
      # axis.text.x = element_text(size = 14),
      # legend.text = element_text(size = 14),
      panel.background = element_blank(),
      panel.grid = element_blank(),
      axis.ticks.y = element_blank(),
      axis.ticks.x = element_line(color = "black"),
      axis.line.x = element_line(color = "black"),
      axis.text.x = element_text(
        vjust = -1,
      ),
    )
}

# load data ---------------------------------------------------------------

basedir <- "/home/liuc9/github/scMOCHA-data/data"
outdir <- "/home/liuc9/github/scMOCHA-data/analysis/zzz/plot-basic"

gse_dataset_metadata_full <- import(
  "/home/liuc9/github/scMOCHA-data/analysis/zzz/clean-data/gse_dataset_metadata_full.qs"
)

gse_data <- import(
  "/home/liuc9/github/scMOCHA-data/analysis/zzz/clean-data/gse_data.qs"
)


pcc <- import(file = "https://raw.githubusercontent.com/chunjie-sam-liu/chunjie-sam-liu.life/master/public/data/pcc.tsv") |>
  dplyr::arrange(cancer_types)


# thegseid <- "GSE168453"


# body --------------------------------------------------------------------


gse_data |>
  dplyr::select(gseid, srrid, chemistry, anno, hetero, haplo_variant, haplo_violin, somatic_variant, celltype_ratio) |>
  dplyr::left_join(
    gse_dataset_metadata_full |> dplyr::select(-gseid),
    by = c("srrid" = "srrid")
  ) |>
  dplyr::arrange(disease, Chemistry) ->
gse_data_haplo_variant



# ! all variants --------------------------------------------------------------------
regions_missalignment_error <- c(66:71, 300:316, 513:525, 3106:3107, 12418:12425, 16182:16194)
regions_rare_heteroplasmic_variants <- c(499, 538, 545, 10953, 12684)
variants_tobe_excluded <- c(
  regions_missalignment_error,
  regions_rare_heteroplasmic_variants
)

gse_data_haplo_variant |>
  dplyr::mutate(
    heteroplasmic = purrr::map(
      somatic_variant,
      .f = \(.x) {
        tibble::tibble(
          variant = .x$somatic
        ) |>
          dplyr::mutate(
            pos = stringr::str_extract(variant, "\\d+") |> as.integer(),
          ) |>
          dplyr::filter(
            !pos %in% variants_tobe_excluded,
          ) ->
        .xx
        .xx$variant -> heteroplasmic_variant
        c(.x$high_af, .x$haplo) |> unique() -> homoplasmic_variant
        .x$heteroplasmic_variant <- heteroplasmic_variant
        .x$homoplasmic_variant <- homoplasmic_variant

        .x
      }
    )
  ) |>
  dplyr::mutate(
    n_heteroplasmic = purrr::map(
      heteroplasmic,
      .f = \(.x) {
        tibble::tibble(
          n_heteroplasmic = length(.x$heteroplasmic_variant),
          n_homoplasmic = length(.x$homoplasmic_variant),
        )
      }
    )
  ) |>
  tidyr::unnest(cols = n_heteroplasmic) ->
gse_data_variant_heteroplasmic



gse_data_variant_heteroplasmic$heteroplasmic |>
  purrr::map(
    .f = \(.x) {
      .x$heteroplasmic_variant
    }
  ) |>
  purrr::reduce(
    union
  ) ->
heteroplasmic_variant

gse_data_variant_heteroplasmic$heteroplasmic |>
  purrr::map(
    .f = \(.x) {
      .x$homoplasmic_variant
    }
  ) |>
  purrr::reduce(
    union
  ) ->
homoplasmic_variant

gse_data_variant_heteroplasmic |>
  dplyr::select(srrid, hetero) |>
  tidyr::unnest(cols = hetero) |>
  dplyr::group_by(srrid, variant) |>
  dplyr::summarise(
    af = mean(af, na.rm = TRUE),
  ) |>
  dplyr::ungroup() |>
  dplyr::group_by(variant) |>
  dplyr::summarise(
    af = mean(af, na.rm = TRUE),
  ) ->
variant_mean_af

gse_data_haplo_variant |>
  dplyr::select(gseid, srrid, chemistry, haplo_variant) |>
  tidyr::unnest(cols = haplo_variant) ->
all_variants

all_variants |>
  dplyr::filter(Position %in% c(3243, 8344, 8993, 13513, 14709, 10191, 14459, 3460))

all_variants |>
  dplyr::filter(Position %in% c(1555, 1494, 961, 2336, 3090))

# all_variants |>
#   dplyr::filter(Disease != "") |>
#   View()


all_variants |>
  dplyr::select(Position, variant, aachange, Disease, `Gnomad Frequency`) |>
  dplyr::mutate(
    Disease = ifelse(is.na(Disease), "", Disease),
  ) |>
  dplyr::distinct() |>
  dplyr::arrange(Position) ->
variant_type


all_variants |>
  dplyr::count(variant) |>
  dplyr::left_join(
    variant_type,
    by = "variant"
  ) |>
  dplyr::mutate(
    issomatic = ifelse(variant %in% heteroplasmic_variant, "heteroplasmic", "other"),
  ) |>
  dplyr::mutate(
    issomatic = ifelse(variant %in% homoplasmic_variant, "homoplasmic", issomatic),
  ) |>
  dplyr::arrange(
    desc(n)
  ) |>
  dplyr::group_by(Position) |>
  dplyr::mutate(
    issomatic = ifelse(
      dplyr::n() > 1,
      "multiple",
      issomatic
    )
  ) |>
  dplyr::ungroup() |>
  dplyr::left_join(
    variant_mean_af,
    by = "variant"
  ) ->
variant_count

\(){
  export(
    variant_count,
    file = file.path(
      "/home/liuc9/github/scMOCHA-data/analysis/zzz/clean-data/", "all_variant.csv"
    ),
    format = "both"
  )
  # readr::write_rds(
  #   variant_count,
  #   file = file.path(
  #     "/home/liuc9/github/scMOCHA-data/analysis/zzz/clean-data/", "all_variant.rds"
  #   )
  # )
}

# variant_count |> dplyr::filter(Position == 3933)
variant_count <- import(
  file.path(
    "/home/liuc9/github/scMOCHA-data/analysis/zzz/clean-data/", "all_variant.qs"
  )
)

# ! homoplasmic variant --------------------------------------------------------------------



variant_count |>
  dplyr::filter(
    issomatic == "homoplasmic"
  ) |>
  dplyr::group_by(Position) |>
  dplyr::filter(dplyr::n() == 1) |>
  dplyr::ungroup() |>
  dplyr::mutate(
    freq = n / nrow(gse_data),
    `Gnomad Frequency` = `Gnomad Frequency` / 100
  ) ->
variant_count_homoplasmic

variant_count_homoplasmic |>
  ggpubr::ggscatter(
    x = "freq",
    y = "Gnomad Frequency",
    # conf.int = TRUE,
    cor.coef = TRUE, # Add correlation coefficient. see ?stat_cor
    cor.coeff.args = list(
      method = "pearson",
      label.x = 0.1,
      label.sep = "\n",
      size = 5
    )
  ) +
  # geom_point(alpha = 0.7) +
  ggpointdensity::geom_pointdensity(
    adjust = 0.01,
    show.legend = FALSE,
  ) +
  viridis::scale_color_viridis() +
  scale_x_continuous(limits = c(0, 1), expand = expansion(add = c(0.01, 0))) +
  scale_y_continuous(limits = c(0, 1), expand = expansion(add = c(0.01, 0))) +
  labs(
    x = "Germline variant popultation frequency",
    y = "Gnomad Frequency"
  ) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "darkgray") +
  theme(
    axis.title = element_text(size = 16),
  ) ->
p_homoplasmic_variant_correlates_with_gnomad
ggsave(
  filename = file.path(
    "/home/liuc9/github/scMOCHA-data/analysis/zzz/plot-heteroplasmic", "homoplasmic_variant_correlates_with_gnomad.pdf"
  ),
  plot = p_homoplasmic_variant_correlates_with_gnomad,
  width = 6,
  height = 5,
  dpi = 300
)

variant_count_homoplasmic |>
  ggplot(aes(
    x = Position,
    y = n
  )) +
  geom_bar(stat = "identity", color = "red") +
  geom_text(
    data = variant_count_homoplasmic |>
      head(5) |>
      dplyr::mutate(
        label = glue::glue("{variant}, {aachange} \n FQ={n}/{nrow(gse_data)}, GF={`Gnomad Frequency`}\n{Disease}")
      ),
    aes(
      label = label,
    ),
    color = "black",
    size = 3,
    vjust = -0.5,
    hjust = 0.5,
    show.legend = FALSE
  ) +
  scale_x_continuous(
    expand = expansion(mult = c(0, 0.01)),
    limits = c(1, 16569),
    breaks = c(seq(0, 17000, 1000), 16569),
    labels = c(seq(0, 17000, 1000), 16569),
  ) +
  scale_y_continuous(
    expand = expansion(mult = c(0, 0.2)),
    label = scales::label_number()
  ) +
  # scale_fill_identity(
  #   name = "Sample"
  # ) +
  theme(
    plot.margin = margin(t = 0, b = 0, unit = "cm"),
    panel.background = element_blank(),
    panel.grid = element_blank(),
    axis.line.y.left = element_line(color = "black"),
    # axis.line.x.bottom = element_line(color = "black"),
    axis.ticks.x = element_blank(),
    axis.text.x = element_blank(),
    axis.line.x = element_blank(),
    axis.title.x = element_blank(),
    # legend.position = c(0.8, 0.5),
    legend.position = "none",
    legend.key = element_blank(),
    axis.title.y = element_text(size = 16, color = "black"),
    axis.text.y = element_text(color = "black"),
    legend.text = element_text(
      size = 14,
      color = "black"
    ),
    legend.title = element_text(
      size = 16,
      colour = "black"
    ),
    strip.background = element_blank(),
    strip.text = element_text(
      size = 8,
      color = "black",
      face = "bold"
    )
  ) +
  labs(
    y = "# of Samples",
  ) ->
p_variant_count_homoplasmic

wrap_plots(
  p_variant_count_homoplasmic,
  fn_plot_mtdna(),
  ncol = 1,
  heights = c(15, 1)
)

ggsave(
  filename = file.path(
    "/home/liuc9/github/scMOCHA-data/analysis/zzz/plot-heteroplasmic", "homoplasmic_variant_distribution__samples.pdf"
  ),
  plot = wrap_plots(
    p_variant_count_homoplasmic,
    fn_plot_mtdna(),
    ncol = 1,
    heights = c(15, 1)
  ),
  width = 20,
  height = 10,
  dpi = 300
)


variant_count_homoplasmic |>
  ggplot(aes(
    x = Position,
    y = af,
    color = af
  )) +
  geom_bar(stat = "identity") +
  viridis::scale_color_viridis(direction = 1) +
  scale_x_continuous(
    expand = expansion(mult = c(0, 0.01)),
    limits = c(1, 16569),
    breaks = c(seq(0, 17000, 1000), 16569),
    labels = c(seq(0, 17000, 1000), 16569),
  ) +
  scale_y_continuous(
    expand = expansion(add = c(0, 0.01)),
    limits = c(0, 1),
    label = scales::label_number()
  ) +
  theme(
    plot.margin = margin(t = 0, b = 0, unit = "cm"),
    panel.background = element_blank(),
    panel.grid = element_blank(),
    axis.line.y.left = element_line(color = "black"),
    # axis.line.x.bottom = element_line(color = "black"),
    axis.ticks.x = element_blank(),
    axis.text.x = element_blank(),
    axis.line.x = element_blank(),
    axis.title.x = element_blank(),
    # legend.position = c(0.8, 0.5),
    legend.position = "none",
    legend.key = element_blank(),
    axis.title.y = element_text(size = 16, color = "black"),
    axis.text.y = element_text(color = "black"),
    legend.text = element_text(
      size = 14,
      color = "black"
    ),
    legend.title = element_text(
      size = 16,
      colour = "black"
    ),
    strip.background = element_blank(),
    strip.text = element_text(
      size = 8,
      color = "black",
      face = "bold"
    )
  ) +
  labs(
    y = "Mean heteroplasmic frequency",
  ) ->
p_variant_af_homoplasmic


wrap_plots(
  p_variant_af_homoplasmic,
  fn_plot_mtdna(),
  ncol = 1,
  heights = c(15, 1)
)

ggsave(
  filename = file.path(
    "/home/liuc9/github/scMOCHA-data/analysis/zzz/plot-heteroplasmic", "homoplasmic_variant_distribution_af_samples.pdf"
  ),
  plot = wrap_plots(
    p_variant_af_homoplasmic,
    fn_plot_mtdna(),
    ncol = 1,
    heights = c(15, 1)
  ),
  width = 20,
  height = 10,
  dpi = 300
)

# ! heteroplasmic --------------------------------------------------------------------


variant_count |>
  dplyr::filter(
    issomatic == "heteroplasmic"
  ) |>
  dplyr::group_by(Position) |>
  dplyr::filter(dplyr::n() == 1) |>
  dplyr::ungroup() |>
  dplyr::mutate(
    freq = n / nrow(gse_data),
    `Gnomad Frequency` = `Gnomad Frequency` / 100
  ) ->
variant_count_heteroplasmic

variant_count_heteroplasmic |>
  # dplyr::filter(aachange == "rRNA") |>
  ggpubr::ggscatter(
    x = "freq",
    y = "Gnomad Frequency",
    # conf.int = TRUE,
    cor.coef = TRUE, # Add correlation coefficient. see ?stat_cor
    cor.coeff.args = list(
      method = "pearson",
      label.x = 0.1,
      label.y = 0.1,
      label.sep = "\n",
      size = 5
    )
  ) +
  # geom_point(alpha = 0.7) +
  ggpointdensity::geom_pointdensity(
    adjust = 0.01,
    show.legend = FALSE,
  ) +
  viridis::scale_color_viridis() +
  scale_x_continuous(limits = c(0, 1), expand = expansion(add = c(0.01, 0))) +
  scale_y_continuous(limits = c(0, 1), expand = expansion(add = c(0.01, 0))) +
  labs(
    x = "Somatic variant popultation frequency",
    y = "Gnomad Frequency"
  ) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "darkgray") +
  theme(
    axis.title = element_text(size = 16),
  ) +
  coord_fixed(
    xlim = c(0, 0.2),
    ylim = c(0, 0.2)
  ) ->
p_heteroplasmic_variant_correlates_with_gnomad
ggsave(
  filename = file.path(
    "/home/liuc9/github/scMOCHA-data/analysis/zzz/plot-heteroplasmic", "heteroplasmic_variant_correlates_with_gnomad.pdf"
  ),
  plot = p_heteroplasmic_variant_correlates_with_gnomad,
  width = 6,
  height = 5,
  dpi = 300
)

variant_count_heteroplasmic |>
  ggplot(aes(
    x = Position,
    y = n
  )) +
  geom_bar(stat = "identity", color = "red") +
  geom_text(
    data = variant_count_heteroplasmic |>
      head(5) |>
      dplyr::mutate(
        label = glue::glue("{variant}, {aachange} \n FQ={n}/{nrow(gse_data)}, GF={`Gnomad Frequency`}\n{Disease}")
      ),
    aes(
      label = label,
    ),
    color = "black",
    size = 3,
    vjust = -0.5,
    hjust = 0.5,
    show.legend = FALSE
  ) +
  scale_x_continuous(
    expand = expansion(mult = c(0, 0.01)),
    limits = c(1, 16569),
    breaks = c(seq(0, 17000, 1000), 16569),
    labels = c(seq(0, 17000, 1000), 16569),
  ) +
  scale_y_continuous(
    expand = expansion(mult = c(0, 0.2)),
    label = scales::label_number()
  ) +
  # scale_fill_identity(
  #   name = "Sample"
  # ) +
  theme(
    plot.margin = margin(t = 0, b = 0, unit = "cm"),
    panel.background = element_blank(),
    panel.grid = element_blank(),
    axis.line.y.left = element_line(color = "black"),
    # axis.line.x.bottom = element_line(color = "black"),
    axis.ticks.x = element_blank(),
    axis.text.x = element_blank(),
    axis.line.x = element_blank(),
    axis.title.x = element_blank(),
    # legend.position = c(0.8, 0.5),
    legend.position = "none",
    legend.key = element_blank(),
    axis.title.y = element_text(size = 16, color = "black"),
    axis.text.y = element_text(color = "black"),
    legend.text = element_text(
      size = 14,
      color = "black"
    ),
    legend.title = element_text(
      size = 16,
      colour = "black"
    ),
    strip.background = element_blank(),
    strip.text = element_text(
      size = 8,
      color = "black",
      face = "bold"
    )
  ) +
  labs(
    y = "# of Samples",
  ) ->
p_variant_count_heteroplasmic

wrap_plots(
  p_variant_count_heteroplasmic,
  fn_plot_mtdna(),
  ncol = 1,
  heights = c(15, 1)
)

ggsave(
  filename = file.path(
    "/home/liuc9/github/scMOCHA-data/analysis/zzz/plot-heteroplasmic", "heteroplasmic_variant_distribution__samples.pdf"
  ),
  plot = wrap_plots(
    p_variant_count_heteroplasmic,
    fn_plot_mtdna(),
    ncol = 1,
    heights = c(15, 1)
  ),
  width = 20,
  height = 10,
  dpi = 300
)

variant_count_heteroplasmic |>
  ggplot(aes(
    x = Position,
    y = af,
    color = af
  )) +
  geom_bar(stat = "identity") +
  geom_text(
    data = variant_count_heteroplasmic |>
      dplyr::arrange(desc(af)) |>
      head(5) |>
      dplyr::mutate(
        label = glue::glue("{variant}, {aachange} \n FQ={n}/{nrow(gse_data)}, GF={`Gnomad Frequency`}\n{Disease}")
      ),
    aes(
      label = label,
    ),
    color = "black",
    size = 3,
    vjust = -0.5,
    hjust = 0.5,
    show.legend = FALSE
  ) +
  viridis::scale_color_viridis(direction = 1) +
  scale_x_continuous(
    expand = expansion(mult = c(0, 0.01)),
    limits = c(1, 16569),
    breaks = c(seq(0, 17000, 1000), 16569),
    labels = c(seq(0, 17000, 1000), 16569),
  ) +
  scale_y_continuous(
    expand = expansion(add = c(0, 0.1)),
    limits = c(0, 1),
    label = scales::label_number()
  ) +
  theme(
    plot.margin = margin(t = 0, b = 0, unit = "cm"),
    panel.background = element_blank(),
    panel.grid = element_blank(),
    axis.line.y.left = element_line(color = "black"),
    # axis.line.x.bottom = element_line(color = "black"),
    axis.ticks.x = element_blank(),
    axis.text.x = element_blank(),
    axis.line.x = element_blank(),
    axis.title.x = element_blank(),
    # legend.position = c(0.8, 0.5),
    legend.position = "none",
    legend.key = element_blank(),
    axis.title.y = element_text(size = 16, color = "black"),
    axis.text.y = element_text(color = "black"),
    legend.text = element_text(
      size = 14,
      color = "black"
    ),
    legend.title = element_text(
      size = 16,
      colour = "black"
    ),
    strip.background = element_blank(),
    strip.text = element_text(
      size = 8,
      color = "black",
      face = "bold"
    )
  ) +
  labs(
    y = "Mean heteroplasmic frequency",
  ) ->
p_variant_af_heteroplasmic


wrap_plots(
  p_variant_af_heteroplasmic,
  fn_plot_mtdna(),
  ncol = 1,
  heights = c(15, 1)
)

ggsave(
  filename = file.path(
    "/home/liuc9/github/scMOCHA-data/analysis/zzz/plot-heteroplasmic", "heteroplasmic_variant_distribution_af_samples.pdf"
  ),
  plot = wrap_plots(
    p_variant_af_heteroplasmic,
    fn_plot_mtdna(),
    ncol = 1,
    heights = c(15, 1)
  ),
  width = 20,
  height = 10,
  dpi = 300
)


# body --------------------------------------------------------------------


# footer ------------------------------------------------------------------

# future: :plan(future: :sequential)

# save image --------------------------------------------------------------
