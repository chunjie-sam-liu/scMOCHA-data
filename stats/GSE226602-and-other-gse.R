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
basedir <- "/home/liuc9/github/scMOCHA-data/data"

outdir <- file.path(basedir, "out")

gseids <- c("GSE149689", "GSE155223", "GSE155673", "GSE157344", "GSE163668", "GSE166992", "GSE171555", "GSE181279", "GSE226602")

gseids_meta <- tibble::tibble(
  GSE_ID = c("GSE163668", "GSE149689", "GSE155223", "GSE155673", "GSE157344", "GSE166992", "GSE171555", "GSE226602", "GSE181279"),
  samples = c(38, 20, 18, 12, 33, 9, 48, 50, 5),
  Disease = c("COVID-19", "COVID-19", "COVID-19", "COVID-19", "COVID-19", "COVID-19", "COVID-19", "AD", "AD"),
  Source = c("PBMC", "PBMC", "PBMC", "PBMC", "PBMC", "PBMC", "PBMC", "PBMC", "PBMC"),
  Chemistry = c("SC5P-R2", "SC3Pv3", "SC5P-R2", "SC3Pv3", "SC3Pv3", "SC5P-PE", "SC5P-R2", "SC5P-PE", "SC5P-PE"),
  Publication = c("Nature, 2021", "Exp Mol Med, 2022", "Cell Rep, 2023", "Science, 2020", "Nat Commun, 2021", "Cell Rep, 2021", "Med, 2021", "Neuron, 2024", "Front Immunol., 2021")
)

# body --------------------------------------------------------------------

# load gse cell ratio and variant data

tibble::tibble(
  gseid = gseids
) |>
  dplyr::mutate(
    cell_ratio_variant = purrr::map(
      .x = gseid,
      .f = \(.gseid) {
        data.table::fread(
          file.path(basedir, .gseid, "out", glue::glue("{.gseid}.cell_ratio_and_variant_clean.csv"))
        )
      }
    )
  ) |>
  dplyr::left_join(
    gseids_meta,
    by = c("gseid" = "GSE_ID")
  ) ->
gse_cell_ratio_variant_meta

gse_cell_ratio_variant_meta |>
  dplyr::mutate(
    `Avg # of somatic variants` = purrr::map_dbl(
      cell_ratio_variant,
      ~ mean(.x$`# of somatic variants`)
    )
  ) ->
gse_cell_ratio_variant_meta_xlsx

gse_cell_ratio_variant_meta_xlsx |>
  dplyr::select(-cell_ratio_variant) |>
  writexl::write_xlsx(
    path = file.path(outdir, "gses_cell_ratio_variant_meta.xlsx")
  )

gse_cell_ratio_variant_meta_xlsx |>
  dplyr::arrange(dplyr::desc(`Avg # of somatic variants`)) |>
  dplyr::mutate(
    label = glue::glue("{gseid} ({round(`Avg # of somatic variants`, 2)})")
  ) |>
  dplyr::mutate(
    label = factor(label, levels = label)
  ) ->
gseid_ranked

gse_cell_ratio_variant_meta |>
  dplyr::left_join(
    gseid_ranked |> dplyr::select(gseid, label),
    by = "gseid"
  ) |>
  tidyr::unnest(cols = cell_ratio_variant) |>
  dplyr::mutate(
    Chemistry = factor(Chemistry, levels = c("SC3Pv3", "SC5P-R2", "SC5P-PE") |> rev()),
  ) ->
forplot

cor.test(~ `Depth mean` + `# of somatic variants`, data = forplot) |> broom::tidy() -> cor_test_all
cor.test(~ `Depth mean` + `# of somatic variants`, data = forplot, subset = `Depth mean` < 250000) |> broom::tidy() -> cor_test_250k

forplot |>
  ggplot(aes(
    x = `Depth mean`,
    y = `# of somatic variants`,
  )) +
  geom_point(aes(
    shape = Chemistry,
    color = label
  )) +
  geom_smooth(
    method = "loess", se = FALSE, color = "black",
    linetype = 21,
  ) +
  geom_hline(
    yintercept = 10,
    linetype = 21,
    color = "red"
  ) +
  scale_color_brewer(
    name = "GSE ID",
    palette = "Set1",
    direction = -1
  ) +
  scale_x_continuous(
    labels = scales::label_number(),
    limits = c(0, 350000),
    breaks = seq(0, 350000, 50000),
  ) +
  scale_y_continuous(
    labels = scales::label_number(),
    limits = c(0, 80),
    breaks = seq(0, 80, 10),
  ) +
  theme(
    panel.background = element_blank(),
    panel.grid = element_blank(),
    axis.line = element_line(size = 0.5, color = "black"),
    axis.title = element_text(
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
  ) +
  labs(
    x = "Average Depth",
    y = "Number of Somatic Variants",
    title = glue::glue("All samples test, Pearson's r = {round(cor_test_all$estimate, 2)}, p-value = {scales::pvalue(cor_test_all$p.value)}"),
    subtitle = glue::glue("Exclude 3 samples, Pearson's r = {round(cor_test_250k$estimate, 2)}, p-value = {scales::pvalue(cor_test_250k$p.value)}")
  ) ->
p_somatic_variant

ggsave(
  filename = file.path(outdir, "somatic_variant.pdf"),
  plot = p_somatic_variant,
  width = 9,
  height = 6,
  dpi = 300
)

# footer ------------------------------------------------------------------

# future: :plan(future: :sequential)

# save image --------------------------------------------------------------
