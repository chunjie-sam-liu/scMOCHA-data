#!/usr/bin/env Rscript
# Metainfo ----------------------------------------------------------------
# @AUTHOR: Chun-Jie Liu
# @CONTACT: chunjie.sam.liu.at.gmail.com
# @DATE: 2025-12-18 12:02:16
# @DESCRIPTION: this script is used for ...

# Library -----------------------------------------------------------------

load_pkg(jutils)

# args --------------------------------------------------------------------

# s: string, i: integer, f: float, !: boolean, @: array, %: list
GetoptLong.options(help_style = "two-column")
VERSION = "v0.0.1"

# default: default value specified here.

verbose = TRUE

GetoptLong("verbose!", "print messages")


logger::log_threshold(logger::TRACE)
logger::log_layout(logger::layout_glue_colors)

# header ------------------------------------------------------------------

# load data ---------------------------------------------------------------
outdir <- path("/home/liuc9/github/scMOCHA-data/analysis/zzz/MANUSCRIPTFIGURES")
METAFULL <- import(outdir / "SAMPLES-METADATA-FULL.xlsx")
# source color
source(path(
  "/home/liuc9/github/scMOCHA-data/analysis/high-res/src/00-colors.R"
))

# load conn ---------------------------------------------------------------

# src ---------------------------------------------------------------------

# function ----------------------------------------------------------------

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
          palette = color_chemistry,
          # conf.int = TRUE,
          cor.coef = TRUE, # Add correlation coefficient. see ?stat_cor
          cor.coeff.args = list(
            method = "pearson",
            label.x = 3,
            label.sep = "\n"
          )
        ) +
        theme(
          legend.position = "none",
        )
    },
    error = function(e) {
      log_fatal("Error in plotting p_ncells: {e$message}")
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
          palette = color_chemistry,
          # conf.int = TRUE,
          cor.coef = TRUE, # Add correlation coefficient. see ?stat_cor
          cor.coeff.args = list(
            method = "pearson",
            label.x = 3,
            label.sep = "\n"
          )
        ) +
        theme(
          legend.position = "none",
        )
    },
    error = function(e) {
      log_fatal("Error in plotting p_numi: {e$message}")
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
          palette = color_chemistry,
          # conf.int = TRUE,
          cor.coef = TRUE, # Add correlation coefficient. see ?stat_cor
          cor.coeff.args = list(
            method = "pearson",
            label.x = 3,
            label.sep = "\n"
          )
        ) +
        theme(
          legend.position = "none",
        )
    },
    error = function(e) {
      log_fatal("Error in plotting p_avg_depth: {e$message}")
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
          palette = color_chemistry,
          add = "reg.line",
          # conf.int = TRUE,
          cor.coef = TRUE, # Add correlation coefficient. see ?stat_cor
          cor.coeff.args = list(
            method = "pearson",
            label.x = 3,
            label.sep = "\n"
          )
        ) +
        theme(
          legend.position = "right",
        )
    },
    error = function(e) {
      log_fatal("Error in plotting p_age: {e$message}")
      ggplot()
    }
  )

  p_gender <- tryCatch(
    expr = {
      anno_meta_info_clean |>
        ggpubr::ggboxplot(
          x = "sex_pred",
          y = "# of somatic variants",
          xlab = "",
          ylab = "Number of Somatic Variants",
          title = "Gender",
          color = "Chemistry",
          palette = color_chemistry,
          add = "jitter",
        ) +
        theme(
          legend.position = "none",
        )
    },
    error = function(e) {
      log_fatal("Error in plotting p_gender: {e$message}")
      ggplot()
    }
  )
  p_disease <- tryCatch(
    expr = {
      anno_meta_info_clean |>
        dplyr::mutate(
          disease = dplyr::case_when(
            disease %in%
              c(
                "Alzheimer's Disease",
                "Healthy",
                "COVID-19",
                "Unknown"
              ) ~ disease,
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
          palette = color_chemistry,
          add = "jitter"
        ) +
        theme(
          legend.position = "none",
        )
    },
    error = function(e) {
      log_fatal("Error in plotting p_disease: {e$message}")
      ggplot()
    }
  )

  wrap_plots(
    list(p_ncells, p_numi, p_avg_depth, p_age, p_gender, p_disease),
    ncol = 3,
    guides = "collect",
  ) -> p_combined
  p_combined
}

fn_eda_ggpubr(METAFULL)

# body --------------------------------------------------------------------

# footer ------------------------------------------------------------------

# save image --------------------------------------------------------------
