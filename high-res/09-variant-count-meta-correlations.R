#!/usr/bin/env Rscript
# Metainfo ----------------------------------------------------------------
# @AUTHOR: Chun-Jie Liu
# @CONTACT: chunjie.sam.liu.at.gmail.com
# @DATE: 2025-12-18 12:02:16
# @DESCRIPTION: this script is used for ...

# Library -----------------------------------------------------------------

load_pkg(jutils)

dotenv(".env")

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
  Sys.getenv("HIGHRESDIR"),
  "00-colors.R"
))

# load conn ---------------------------------------------------------------

# src ---------------------------------------------------------------------

# function ----------------------------------------------------------------

fn_eda_ggpubr <- function(
  anno_meta_info_clean,
  varianttype = "# of somatic variants"
) {
  ylab_varianttype <- switch(
    varianttype,
    "# of somatic variants" = "Number of Somatic Variants",
    "# of homoplasmic variants" = "Number of Homoplasmic Variants",
    "# of heteroplasmic variants" = "Number of Heteroplasmic Variants",
    "# of variants" = "Number of Total Variants",
    stop("Invalid varianttype")
  )

  p_ncells <- tryCatch(
    expr = {
      anno_meta_info_clean |>
        ggpubr::ggscatter(
          x = "# cells after filter",
          y = varianttype,
          title = "Number of Cells",
          xlab = "",
          ylab = ylab_varianttype,
          # color = "Chemistry",
          add = "reg.line",
          add.params = list(color = "blue"),
          # palette = color_chemistry,
          # conf.int = TRUE,
          cor.coef = TRUE, # Add correlation coefficient. see ?stat_cor
          cor.coeff.args = list(
            method = "pearson",
            label.x = 3,
            label.sep = "\n",
            p.accuracy = 0.001
          )
        ) +
        scale_x_continuous(
          labels = scales::label_comma()
        ) +
        theme(
          plot.title = element_text(hjust = 0.5),
          legend.position = "none",
          aspect.ratio = 1
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
          y = varianttype,
          title = "Number of UMI",
          xlab = "",
          ylab = ylab_varianttype,
          # color = "Chemistry",
          add = "reg.line",
          add.params = list(color = "blue"),
          # palette = color_chemistry,
          # conf.int = TRUE,
          cor.coef = TRUE, # Add correlation coefficient. see ?stat_cor
          cor.coeff.args = list(
            method = "pearson",
            label.x = 3,
            label.sep = "\n",
            p.accuracy = 0.001
          )
        ) +
        scale_x_continuous(
          labels = scales::label_comma()
        ) +
        theme(
          plot.title = element_text(hjust = 0.5),
          legend.position = "none",
          aspect.ratio = 1
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
          y = varianttype,
          title = "Average Depth",
          xlab = "",
          ylab = ylab_varianttype,
          # color = "Chemistry",
          add = "reg.line",
          add.params = list(color = "blue"),
          # palette = color_chemistry,
          # conf.int = TRUE,
          cor.coef = TRUE, # Add correlation coefficient. see ?stat_cor
          cor.coeff.args = list(
            method = "pearson",
            label.x = 3,
            label.sep = "\n",
            p.accuracy = 0.001
          )
        ) +
        scale_x_continuous(
          labels = scales::label_comma()
        ) +
        theme(
          plot.title = element_text(hjust = 0.5),
          legend.position = "none",
          aspect.ratio = 1
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
          y = varianttype,
          title = "Age",
          xlab = "",
          ylab = ylab_varianttype,
          # color = "Chemistry",
          # palette = color_chemistry,
          add = "reg.line",
          add.params = list(color = "blue"),
          # conf.int = TRUE,
          cor.coef = TRUE, # Add correlation coefficient. see ?stat_cor
          cor.coeff.args = list(
            method = "pearson",
            label.x = 3,
            label.sep = "\n",
            p.accuracy = 0.001
          )
        ) +
        scale_x_continuous(
          labels = scales::label_comma()
        ) +
        theme(
          plot.title = element_text(hjust = 0.5),
          legend.position = "none",
          aspect.ratio = 1
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
          x = "SEXPRED",
          y = varianttype,
          xlab = "",
          ylab = ylab_varianttype,
          title = "Sex",
          color = "SEXPRED",
          # palette = color_chemistry,
          add = "jitter",
        ) +
        scale_color_manual(
          values = color_gender
        ) +
        ggpubr::stat_compare_means(
          method = "t.test",
          label = "p.format"
        ) +
        theme(
          plot.title = element_text(hjust = 0.5),
          legend.position = "none",
          aspect.ratio = 1
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
          y = varianttype,
          xlab = "",
          ylab = ylab_varianttype,
          title = "Disease",
          color = "disease",
          # palette = color_chemistry,
          add = "jitter"
        ) +
        scale_color_manual(
          values = color_disease
        ) +
        ggpubr::stat_compare_means(
          method = "anova",
          label.y = max(
            anno_meta_info_clean$`# of somatic variants`,
            na.rm = TRUE
          ) *
            1.1
        ) +
        theme(
          plot.title = element_text(hjust = 0.5),
          legend.position = "none",
          aspect.ratio = 1,
          axis.text.x = element_text(angle = 15, hjust = 1)
        )
    },
    error = function(e) {
      log_fatal("Error in plotting p_disease: {e$message}")
      ggplot()
    }
  )

  {
    pdf(
      file = outdir /
        glue::glue(
          "VARIANT-COUNT-vs-metadata-{ylab_varianttype}.pdf"
        ),
      width = 6,
      height = 4
    )
    print(p_ncells)
    print(p_numi)
    print(p_avg_depth)
    print(p_age)
    print(p_gender)
    print(p_disease)
    dev.off()
  }

  wrap_plots(
    list(p_ncells, p_numi, p_avg_depth, p_age, p_gender, p_disease),
    ncol = 3,
    guides = "collect",
  ) -> p_combined
  p_combined
}

fn_eda_ggpubr(METAFULL, "# of somatic variants")
fn_eda_ggpubr(METAFULL, "# of homoplasmic variants")
fn_eda_ggpubr(METAFULL, "# of heteroplasmic variants")
fn_eda_ggpubr(METAFULL, "# of variants")

# body --------------------------------------------------------------------

# footer ------------------------------------------------------------------

# save image --------------------------------------------------------------
