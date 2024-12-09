#!/usr/bin/env Rscript --vanilla
# Metainfo ----------------------------------------------------------------
# @AUTHOR: Chun-Jie Liu
# @CONTACT: chunjie.sam.liu.at.gmail.com
# @DATE: `r date()`
# @DESCRIPTION: filename

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
datadir <- "/home/liuc9/github/scMOCHA-data/data/GSE226602"
outdir <- file.path(datadir, "out")

# body --------------------------------------------------------------------

anno <- readr::read_rds(file.path(outdir, "GSE226602.scmocha.out.rds.gz"))
meta <- data.table::fread(file.path(outdir, "GSE226602.pheno.select.csv"))

anno |>
  dplyr::select(-c(srrdir, dir_exists)) |>
  dplyr::inner_join(meta, by = "srrid") ->
anno_meta

anno_meta |>
  dplyr::mutate(
    nmut = purrr::map_int(
      .x = anno,
      .f = \(.x) {
        if (is.null(.x)) {
          return(NA_integer_)
        }
        nrow(.x)
      }
    )
  ) |>
  dplyr::mutate(
    nmut_somatic = purrr::map_int(
      .x = somatic_variant,
      .f = \(.x) {
        if (is.null(.x$somatic)) {
          return(NA_integer_)
        }
        length(.x$somatic)
      }
    )
  ) |>
  dplyr::mutate(
    haplogroup = purrr::map2(
      .x = anno,
      .y = srrid,
      .f = \(.x, .y) {
        message(.y)
        if (is.null(.x)) {
          return(
            tibble::tibble(
              Haplogroup = NA_character_,
              Verbose_haplogroup = NA_character_
            )
          )
        }
        .x |>
          dplyr::select(Haplogroup, Verbose_haplogroup) |>
          dplyr::filter(!is.na(Haplogroup)) |>
          dplyr::filter(Haplogroup != "") |>
          dplyr::distinct() |>
          dplyr::mutate_all(.funs = as.character) ->
        .xx

        if (nrow(.xx) == 0) {
          tibble::tibble(
            Haplogroup = NA_character_,
            Verbose_haplogroup = NA_character_
          )
        } else {
          .xx
        }
      }
    )
  ) |>
  tidyr::unnest(cols = haplogroup) |>
  tidyr::unnest(
    cols = cell_stats
  ) |>
  dplyr::mutate(
    ratio = round(`number of cells after filtering` / `estimated number of cells`, 2)
  ) ->
anno_meta_info

anno_meta_info |>
  dplyr::arrange(disease, age, gender) |>
  dplyr::select(
    Sample = srrid,
    Age = age,
    `# of variants` = nmut,
    Haplogroup = Haplogroup,
    `# of somatic variants` = nmut_somatic,
    Gender = gender,
    Disease = disease,
    `# cells after filter` = `number of cells after filtering`,
  ) ->
metadata_clean
# save to xslx and tsv
data.table::fwrite(
  x = metadata_clean,
  file = file.path(outdir, "GSE226602.age.somatic.csv")
)
writexl::write_xlsx(
  x = metadata_clean,
  path = file.path(outdir, "GSE226602.age.somatic.xlsx")
)

# ggstats correlation plot ------------------------------------------------



ggstatsplot::ggscatterstats(
  data = metadata_clean,
  x = Age,
  y = `# of somatic variants`,
  title = "Age",
  xlab = "Age",
  ylab = "Number of Somatic Variants",
  ggtheme = ggplot2::theme_minimal()
) -> p
p



ggstatsplot::ggscatterstats(
  data = metadata_clean |> dplyr::filter(Disease == "Healthy Control"),
  x = Age,
  y = `# of somatic variants`,
  title = "Age",
  xlab = "Age",
  ylab = "Number of Somatic Variants",
  ggtheme = ggplot2::theme_minimal()
) -> p
p

# Age
# # of cells after filtering
anno_meta_info |>
  dplyr::glimpse()

anno_meta_info |>
  dplyr::mutate(
    avg_depth = purrr::map_dbl(
      .x = depth,
      .f = \(.d) {
        mean(.d$depth)
      }
    )
  ) |>
  dplyr::mutate(
    somatic_variant = purrr::map(
      .x = somatic_variant,
      .f = \(.x) {
        .x$somatic
      }
    )
  ) |>
  dplyr::select(
    srrid,
    ncells = `number of cells after filtering`,
    numi = `median UMI counts per cell`,
    avg_depth,
    nmut_somatic,
    age,
    gender,
    disease,
    haplo_violin,
    somatic_variant
  ) ->
anno_meta_info_clean


anno_meta_info_clean |>
  ggstatsplot::ggscatterstats(
    x = ncells,
    y = nmut_somatic,
    title = "Number of Cells",
    xlab = "",
    ylab = "Number of Somatic Variants",
  ) -> p_ncells
p_ncells

anno_meta_info_clean |>
  ggstatsplot::ggscatterstats(
    x = numi,
    y = nmut_somatic,
    title = "Number of UMI",
    xlab = "",
    ylab = "Number of Somatic Variants",
  ) -> p_numi
p_numi

anno_meta_info_clean |>
  ggstatsplot::ggscatterstats(
    x = avg_depth,
    y = nmut_somatic,
    title = "Average Depth",
    xlab = "",
    ylab = "Number of Somatic Variants",
  ) -> p_avg_depth
p_avg_depth

anno_meta_info_clean |>
  ggstatsplot::ggscatterstats(
    x = age,
    y = nmut_somatic,
    title = "Age",
    xlab = "",
    ylab = "Number of Somatic Variants",
  ) -> p_age
p_age

anno_meta_info_clean |> dplyr::glimpse()

anno_meta_info_clean |>
  ggstatsplot::ggbetweenstats(
    x = gender,
    y = nmut_somatic,
    xlab = "",
    ylab = "Number of Somatic Variants",
    title = "Gender",
  ) -> p_gender
p_gender

anno_meta_info_clean |>
  ggstatsplot::ggbetweenstats(
    x = disease,
    y = nmut_somatic,
    xlab = "",
    ylab = "Number of Somatic Variants",
    title = "Disease"
  ) -> p_disease
p_disease

wrap_plots(list(p_ncells, p_numi, p_avg_depth, p_age, p_gender, p_disease), ncol = 3) -> p_combined
p_combined

outdir_plot <- file.path(outdir, "plot")
dir.create(outdir_plot, showWarnings = FALSE, recursive = TRUE)
ggsave(
  path = outdir_plot,
  filename = "correlation_somatic_variants.pdf",
  plot = p_combined,
  width = 20,
  height = 10
)


anno_meta_info_clean |>
  ggstatsplot::ggscatterstats(
    x = avg_depth,
    y = ncells
  )


anno_meta_info_clean |>
  dplyr::mutate(
    somatic_variant = purrr::map(
      .x = somatic_variant,
      .f = \(.x) {
        as.character(.x)
      }
    )
  ) |>
  dplyr::mutate(
    cell_variant = purrr::map2(
      .x = haplo_violin,
      .y = somatic_variant,
      .f = \(.x, .y) {
        .x |>
          dplyr::filter(
            variant %in% .y
          )
      }
    )
  ) ->
anno_meta_info_clean_cell_variant

anno_meta_info_clean_cell_variant |>
  dplyr::mutate(
    somatic_variant = purrr::map(
      .x = somatic_variant,
      .f = \(.x) {
        tibble::tibble(
          variant = .x
        )
      }
    )
  ) |>
  tidyr::unnest(cols = somatic_variant) |>
  dplyr::group_by(variant) |>
  tidyr::nest() |>
  dplyr::ungroup() |>
  dplyr::mutate(
    srrid = purrr::map(
      .x = data,
      .f = function(.x) {
        .x |> dplyr::pull(srrid)
      }
    )
  ) |>
  dplyr::select(-data) ->
forupset

library(ggupset)
forupset |>
  ggplot(aes(x = srrid)) +
  geom_bar(width = 0.6) +
  geom_text(
    stat = "count",
    aes(label = after_stat(count)),
    vjust = -0.5,
    color = "black",
    size = 3,
    fontface = "bold"
  ) +
  scale_x_upset(order_by = "degree") +
  scale_y_continuous(
    expand = expansion(mult = c(0, 0.1), add = 0)
  ) +
  theme_combmatrix(
    combmatrix.label.make_space = TRUE,
    # combmatrix.panel.point.color.fill = d$color[1],
    combmatrix.panel.line.size = 0,
    combmatrix.label.text = element_text(
      # size = 12,
      color = "black",
      face = "bold"
    ),
    combmatrix.label.extra_spacing = 1,
    combmatrix.panel.striped_background.color.one = "white",
    combmatrix.panel.striped_background.color.two = "grey",
  ) +
  labs(
    y = "# of Variants",
    x = "",
    # title = .x
  ) +
  theme(
    panel.background = element_blank(),
    panel.grid = element_blank(),
    axis.line = element_line(size = 0.5, color = "black"),
    axis.title.y = element_text(
      size = 16,
      color = "black",
      face = "bold"
    ),
    axis.text.y = element_text(
      size = 14,
      color = "black"
    ),
    plot.title = element_text(
      hjust = 0.5,
      color = "black",
      size = 16,
      face = "bold"
    )
  ) ->
p_upset

# save upset plot
ggsave(
  path = outdir_plot,
  filename = "upset_somatic_variants.pdf",
  plot = p_upset,
  width = 27,
  height = 13
)

anno_meta_info_clean_cell_variant |>
  dplyr::select(-haplo_violin, -somatic_variant) |>
  tidyr::unnest(cols = cell_variant) ->
anno_meta_info_clean_cell_variant_unnest

anno_meta_info_clean_cell_variant_unnest |> dplyr::glimpse()

fn_plot_mtdna <- function() {
  mt_exons_df <- "/home/liuc9/github/scMOCHA/fasta/mt_exons.df.rds.gz"


  gtf_gene_df <-
    readr::read_rds(
      file = mt_exons_df
    )
  library(gggenes)
  ggplot(gtf_gene_df, aes(xmin = start, xmax = end, y = seqnames)) +
    # geom_gene_arrow() +
    geom_gene_arrow(
      aes(
        fill = gene_biotype
      ),
      arrowhead_height = unit(3, "mm"), arrowhead_width = unit(1, "mm")
    ) +
    scale_fill_brewer(
      palette = "Set1",
      name = "Gene type",
      labels = c("MT rRNA", "MT tRNA", "Protein coding")
    ) +
    ggrepel::geom_text_repel(
      aes(x = (start + end) / 2, label = gene_name, color = gene_biotype),
      # fill = "white",
      # nudge_x =1,
      # nudge_y = -0.1,
      size = 3,
      show.legend = F,
      max.overlaps = Inf,
    ) +
    scale_color_brewer(palette = "Set1") +
    scale_x_continuous(
      limits = c(0, 17000),
      breaks = seq(0, 17000, 1000),
      expand = expansion(mult = c(0, 0.03)),
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
      axis.line.x = element_line(color = "black"),
      axis.text.x = element_text(
        vjust = -1,
      ),
    ) ->
  pg
  pg
}

forupset |>
  dplyr::mutate(n = purrr::map_int(srrid, length)) |>
  dplyr::arrange(desc(n)) |>
  dplyr::filter(n >= 43) ->
sel_variants



anno_meta_info_clean_cell_variant_unnest |>
  dplyr::filter(variant %in% sel_variants$variant) |>
  dplyr::mutate(
    af = ifelse(af == 0, NA_real_, af)
  ) |>
  dplyr::mutate(
    af = ifelse(depth < log2(10), NA_real_, af)
  ) |>
  dplyr::filter(af > 0) ->
theforplot
pcc <- readr::read_tsv(file = "https://raw.githubusercontent.com/chunjie-sam-liu/chunjie-sam-liu.life/master/public/data/pcc.tsv") |>
  dplyr::arrange(cancer_types)
library(ggh4x)


theforplot |>
  ggplot() +
  ggh4x::facet_wrap2(
    ~cluster,
    ncol = 1,
    strip.position = "right",
    strip = ggh4x::strip_themed(
      background_y = elem_list_rect(
        fill = pcc$color
      ),
      text_y = elem_list_text(
        colour = "white",
        face = c("bold")
      ),
      by_layer_y = FALSE,
    )
  ) +
  ggbeeswarm::geom_quasirandom(
    aes(
      x = pos,
      y = af,
      color = af
    ),
    size = 1,
    dodge.width = .75,
    alpha = .5,
  ) +
  scale_color_gradient2(
    name = "AF",
    low = "white",
    mid = "red",
    high = "#3B0049",
    midpoint = 0.5,
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
    legend.position = "right",
    legend.key = element_blank(),
    axis.title.y = element_text(color = "black"),
    # axis.title.y = element_blank(),
    axis.text.y = element_text(color = "black"),
    # legend.text = element_text(
    #   size = 14,
    #   color = "black"
    # ),
    # legend.title = element_text(
    #   size = 16,
    #   colour = "black"
    # ),
    # strip.background = element_blank(),
    # strip.text = element_text(
    #   # size = 8,
    #   color = "black",
    #   face = "bold"
    # )
  ) +
  labs(y = "AF") ->
p_af_cell

# p_af_cell
ggsave(
  path = outdir_plot,
  filename = "af_allcell.pdf",
  plot = wrap_plots(
    p_af_cell,
    fn_plot_mtdna(),
    ncol = 1,
    heights = c(1.3, 0.1)
  ),
  width = 25,
  height = 12
)

# footer ------------------------------------------------------------------

# future: :plan(future: :sequential)

# save image --------------------------------------------------------------
