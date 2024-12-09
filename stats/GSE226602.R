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
    Haplogroup = Haplogroup,
    Age = age,
    Gender = gender,
    Disease = disease,
    `# cells after filter` = `number of cells after filtering`,
    `# of variants` = nmut,
    `# of somatic variants` = nmut_somatic,
  ) ->
metadata_clean

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
  dplyr::select(-haplo_violin, -somatic_variant) |>
  tidyr::unnest(cols = cell_variant) ->
anno_meta_info_clean_cell_variant_unnest

anno_meta_info_clean_cell_variant_unnest |> dplyr::glimpse()


anno_meta_info_clean_cell_variant_unnest |>
  dplyr::filter(variant == "3176A>T", cluster == "B") |>
  dplyr::mutate(
    depth = exp(depth)
  ) |>
  dplyr::mutate(
    depth = ifelse(depth > 500, 500, depth)
  ) |>
  ggstatsplot::ggscatterstats(
    x = depth,
    y = af
  )





# footer ------------------------------------------------------------------

# future: :plan(future: :sequential)

# save image --------------------------------------------------------------
