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

# chem_levels <- c("SC3Pv2", "SC3Pv3", "SC5P-R2", "SC5P-PE") |> rev()
# ggsci::pal_aaas()(3) |> prismatic::color()
# color_chemistry <- viridis::viridis_pal(option = "D")(4) |>
#   prismatic::color()

source("/home/liuc9/github/scMOCHA-data/stats/stats/00-colors.R")

sex_pred <- import("/home/liuc9/github/scMOCHA-data/stats/stats/zzz/clean-data/gse_srrid_srrdir_sex.qs") |>
  dplyr::select(
    srrid,
    sex_pred = sex
  ) |>
  dplyr::mutate(
    sex_pred = factor(sex_pred, names(color_gender))
  )

anno_meta_info_clean <- import(
  "/home/liuc9/github/scMOCHA-data/stats/stats/zzz/clean-data/gse_dataset_metadata_full.fst"
) |>
  dplyr::mutate(
    Chemistry = factor(
      Chemistry,
      levels = names(color_chemistry)
    ),
  ) |>
  dplyr::left_join(sex_pred, by = "srrid")

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
      log_fatal("Error in plotting p_age: ", e)
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
          palette = color_chemistry,
          add = "jitter"
        ) +
        theme(
          legend.position = "none",
        )
    },
    error = function(e) {
      log_fatal("Error in plotting p_disease: ", e)
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




# ! age --------------------------------------------------------------------


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

cor_test <- cor.test(
  x = anno_meta_info_clean$`Age_new`,
  y = anno_meta_info_clean$`# of somatic variants`,
  method = "pearson"
)

anno_meta_info_clean |>
  ggplot(aes(
    x = `Age_new`,
    y = `# of somatic variants`
  )) +
  geom_point(
    aes(
      color = Chemistry,
    ),
    position = position_jitter(width = 0.1, height = 0.1),
  ) +
  scale_color_manual(
    values = color_chemistry,
    name = "Chemistry",
    labels = names(color_chemistry)
  ) +
  geom_smooth(
    color = "red",
    method = "lm"
  ) +
  scale_x_continuous(
    limits = c(0, 96),
    breaks = seq(0, 96, by = 10),
    labels = seq(0, 96, by = 10),
    expand = expansion(mult = c(0.01, 0.01))
  ) +
  scale_y_continuous(
    limits = c(0, 90),
    breaks = seq(0, 90, by = 10),
    labels = seq(0, 90, by = 10),
    expand = expansion(mult = c(0.01, 0))
  ) +
  annotate(
    geom = "text",
    x = 20,
    y = 80,
    label = human_read_latex_pval(
      .x = human_read(cor_test$p.value),
      .s = glue::glue("R={round(cor_test$estimate,3)}")
    ),
    size = 6
  ) +
  theme_cor() +
  theme(
    plot.title = element_blank(),
    legend.position = "top"
  ) +
  labs(
    title = "Number of Somatic Variants correlates with Age",
    x = "Age",
    y = "# of somatic variants"
  ) ->
p_age_cor
ggsave(
  path = zzz_out,
  filename = "ggpubr_correlation_age_somatic_variants.pdf",
  plot = p_age_cor,
  width = 8,
  height = 6
)



# ! age correlation with ggstatsplot --------------------------------------

library(ggstatsplot)

# Create age correlation plot with ggscatterstats
ggscatterstats(
  data = anno_meta_info_clean,
  x = Age_new,
  y = `# of somatic variants`,
  xlab = "Age",
  ylab = "Number of Somatic Variants",
  title = "Relationship between Age and Number of Somatic Variants",
  xfill = "#CC79A7",
  yfill = "#009E73",
  marginal = TRUE,
  point.args = list(alpha = 0.7, size = 3),
  smooth.line.args = list(color = "red", linewidth = 1, method = "lm"),
  bf.message = FALSE,
  caption = "Each point represents a single dataset"
) +
  scale_x_continuous(
    limits = c(0, 96),
    breaks = seq(0, 96, by = 10),
    labels = seq(0, 96, by = 10),
    expand = expansion(mult = c(0.01, 0.01))
  ) +
  scale_y_continuous(
    limits = c(0, 90),
    breaks = seq(0, 90, by = 10),
    labels = seq(0, 90, by = 10),
    expand = expansion(mult = c(0.01, 0))
  ) +
  annotate(
    geom = "text",
    x = 20,
    y = 80,
    label = human_read_latex_pval(
      .x = human_read(cor_test$p.value),
      .s = glue::glue("R={round(cor_test$estimate,3)}")
    ),
    size = 6
  ) +
  # theme_cor() +
  theme(
    plot.title = element_blank(),
    panel.grid = element_blank(),
    panel.background = element_blank(),
    axis.line = element_line(colour = "black"),
  ) +
  labs(
    title = "Number of Somatic Variants correlates with Age",
    x = "Age",
    y = "# of somatic variants"
  ) ->
p_age_cor_ggstatsplot

# Save the plot
ggsave(
  path = zzz_out,
  filename = "ggpubr_ggstatsplot_correlation_age_somatic_variants.pdf",
  plot = p_age_cor_ggstatsplot,
  width = 12,
  height = 8
)


# ! sex --------------------------------------------------------------------
library(ggnewscale)
anno_meta_info_clean |>
  ggplot(aes(
    x = sex_pred,
    y = `# of somatic variants`
  )) +
  geom_boxplot(
    aes(color = sex_pred),
    width = 0.5,
    show.legend = FALSE,
    outlier.shape = NA,
  ) +
  scale_color_manual(
    values = color_gender
  ) +
  new_scale_color() +
  geom_point(
    aes(color = Chemistry),
    position = position_jitter(width = 0.1, height = 0.1),
  ) +
  scale_color_manual(
    values = color_chemistry,
    name = "Chemistry",
    labels = names(color_chemistry)
  ) +
  ggsignif::geom_signif(
    aes(
      y = `# of somatic variants`,
    ),
    comparisons = list(
      c("Female", "Male")
    ),
    y_position = 80
  ) +
  theme_cor() +
  theme(legend.position = "top", ) +
  labs(
    x = "",
    y = "# of somatic variants",
  ) ->
p1

anno_meta_info_clean |>
  dplyr::filter(
    Chemistry == "SC5P-PE",
    `# of somatic variants` >= 1
  ) |>
  ggplot(aes(
    x = sex_pred,
    y = `# of somatic variants`
  )) +
  geom_boxplot(
    width = 0.5,
    aes(color = sex_pred),
    show.legend = FALSE,
    outlier.shape = NA,
  ) +
  scale_color_manual(
    values = color_gender
  ) +
  new_scale_color() +
  geom_point(
    aes(color = Chemistry),
    position = position_jitter(width = 0.1, height = 0.1),
    show.legend = FALSE
  ) +
  scale_color_manual(
    values = color_chemistry,
    name = "Chemistry",
    labels = names(color_chemistry)
  ) +
  ggsignif::geom_signif(
    aes(
      y = `# of somatic variants`,
    ),
    comparisons = list(
      c("Female", "Male")
    ),
    y_position = 80
  ) +
  theme_cor() +
  labs(
    x = "",
    y = "# of somatic variants",
  ) ->
p2

wrap_plots(
  p1,
  p2,
  ncol = 2
) +
  # guide_area() +
  plot_layout(guides = "collect") &
  theme(
    legend.direction = "horizontal",
    legend.box = "horizontal",
    legend.position = "top",
    # plot.background = element_rect(color = "black"),
  ) ->
p_collect

# Save the plot
ggsave(
  path = zzz_out,
  filename = "ggpubr_ggstatsplot_correlation_sex_somatic_variants.pdf",
  plot = p_collect,
  width = 8,
  height = 5
)


# footer ------------------------------------------------------------------

# future: :plan(future: :sequential)

# save image --------------------------------------------------------------
