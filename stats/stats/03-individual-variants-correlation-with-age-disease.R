#!/usr/bin/env Rscript
# Metainfo ----------------------------------------------------------------
# @AUTHOR: Chun-Jie Liu
# @CONTACT: chunjie.sam.liu.at.gmail.com
# @DATE: 2025-03-06 13:32:57
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


# load data ---------------------------------------------------------------


basedir <- "/home/liuc9/github/scMOCHA-data/data"


# body --------------------------------------------------------------------

chem_levels <- c("SC3Pv2", "SC3Pv3", "SC5P-R2", "SC5P-PE") |> rev()
ggsci::pal_aaas()(3) |> prismatic::color()
chem_colors <- viridis::viridis_pal(option = "D")(4) |>
  prismatic::color()

anno_meta_info_clean <- readr::read_rds(
  "/home/liuc9/github/scMOCHA-data/stats/stats/zzz/clean-data/gse_dataset_metadata_full.rds"
) |>
  dplyr::mutate(
    Chemistry = factor(
      Chemistry,
      levels = chem_levels
    ),
  )

fn_eda_ggpubr <- function(anno_meta_info_clean) {
  p_ncells <- tryCatch(
    expr = {
      anno_meta_info_clean |>
        ggpubr::ggscatter(
          x = "# cells after filter",
          y = "# of somatic variants",
          title = "Number of Cells",
          xlab = "",
          ylab = "Number of Somatic Variants",
          color = "Chemistry",
          add = "reg.line",
          palette = chem_colors,
          # conf.int = TRUE,
          cor.coef = TRUE, # Add correlation coefficient. see ?stat_cor
          cor.coeff.args = list(
            method = "pearson",
            label.x = 3,
            label.sep = "\n"
          )
        )
    },
    error = function(e) {
      log_fatal("Error in plotting p_ncells: ", e)
      ggplot()
    }
  )

  p_numi <- tryCatch(
    expr = {
      anno_meta_info_clean |>
        dplyr::mutate(
          `Median UMI cell` = as.numeric(`Median UMI/cell`)
        ) |>
        ggpubr::ggscatter(
          x = "Median UMI cell",
          y = "# of somatic variants",
          title = "Number of UMI",
          xlab = "",
          ylab = "Number of Somatic Variants",
          color = "Chemistry",
          add = "reg.line",
          palette = chem_colors,
          # conf.int = TRUE,
          cor.coef = TRUE, # Add correlation coefficient. see ?stat_cor
          cor.coeff.args = list(
            method = "pearson",
            label.x = 3,
            label.sep = "\n"
          )
        )
    },
    error = function(e) {
      log_fatal("Error in plotting p_numi: ", e)
      ggplot()
    }
  )

  p_avg_depth <- tryCatch(
    expr = {
      anno_meta_info_clean |>
        ggpubr::ggscatter(
          x = "Depth read mean",
          y = "# of somatic variants",
          title = "Average Depth",
          xlab = "",
          ylab = "Number of Somatic Variants",
          color = "Chemistry",
          add = "reg.line",
          palette = chem_colors,
          # conf.int = TRUE,
          cor.coef = TRUE, # Add correlation coefficient. see ?stat_cor
          cor.coeff.args = list(
            method = "pearson",
            label.x = 3,
            label.sep = "\n"
          )
        )
    },
    error = function(e) {
      log_fatal("Error in plotting p_avg_depth: ", e)
      ggplot()
    }
  )

  p_age <- tryCatch(
    expr = {
      anno_meta_info_clean |>
        ggpubr::ggscatter(
          x = "Age_new",
          y = "# of somatic variants",
          title = "Age",
          xlab = "",
          ylab = "Number of Somatic Variants",
          color = "Chemistry",
          add = "reg.line",
          palette = chem_colors,
          # conf.int = TRUE,
          cor.coef = TRUE, # Add correlation coefficient. see ?stat_cor
          cor.coeff.args = list(
            method = "pearson",
            label.x = 3,
            label.sep = "\n"
          )
        )
    },
    error = function(e) {
      log_fatal("Error in plotting p_age: ", e)
      ggplot()
    }
  )

  p_gender <- tryCatch(
    expr = {
      anno_meta_info_clean |>
        ggpubr::ggboxplot(
          x = "Gender",
          y = "# of somatic variants",
          xlab = "",
          ylab = "Number of Somatic Variants",
          title = "Gender",
          color = "Chemistry",
          palette = chem_colors,
          add = "jitter",
        )
    },
    error = function(e) {
      log_fatal("Error in plotting p_gender: ", e)
      ggplot()
    }
  )
  p_disease <- tryCatch(
    expr = {
      anno_meta_info_clean |>
        dplyr::mutate(
          disease = dplyr::case_when(
            disease %in% c("Alzheimer's Disease", "Healthy", "COVID-19", "Unknown") ~ disease,
            TRUE ~ "Other"
          )
        ) |>
        dplyr::mutate(
          disease = factor(
            disease,
            levels = c(
              "Healthy",
              "Alzheimer's Disease",
              "COVID-19",
              "Unknown",
              "Other"
            )
          )
        ) |>
        ggpubr::ggboxplot(
          x = "disease",
          y = "# of somatic variants",
          xlab = "",
          ylab = "Number of Somatic Variants",
          title = "Disease",
          color = "Chemistry",
          palette = chem_colors,
          add = "jitter"
        )
    },
    error = function(e) {
      log_fatal("Error in plotting p_disease: ", e)
      ggplot()
    }
  )

  wrap_plots(list(p_ncells, p_numi, p_avg_depth, p_age, p_gender, p_disease), ncol = 3) -> p_combined
  p_combined
}

zzz_out <- "/home/liuc9/github/scMOCHA-data/stats/stats/zzz"
fn_eda_ggpubr(anno_meta_info_clean)
ggsave(
  path = zzz_out,
  filename = "ggpubr_correlation_somatic_variants.pdf",
  plot = fn_eda_ggpubr(anno_meta_info_clean),
  width = 20,
  height = 10
)

ggsave(
  path = zzz_out,
  filename = "ggpubr_correlation_somatic_variants_cutoff2.pdf",
  plot = fn_eda_ggpubr(anno_meta_info_clean |> dplyr::filter(`# of somatic variants` >= 2)),
  width = 20,
  height = 10
)

ggsave(
  path = zzz_out,
  filename = "ggpubr_correlation_somatic_variants_cutoff1.pdf",
  plot = fn_eda_ggpubr(anno_meta_info_clean |> dplyr::filter(`# of somatic variants` >= 1)),
  width = 20,
  height = 10
)


# footer ------------------------------------------------------------------

# future: :plan(future: :sequential)

# save image --------------------------------------------------------------
