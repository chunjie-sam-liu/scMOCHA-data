#!/usr/bin/env Rscript
# Metainfo ----------------------------------------------------------------
# @AUTHOR: Chun-Jie Liu
# @CONTACT: chunjie.sam.liu.at.gmail.com
# @DATE: 2025-03-13 11:29:32
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
basedir <- "/home/liuc9/github/scMOCHA-data/data/scfoundation2/PBMC"


gseids = (gseids <- c(
  "GSE143353", # done
  "GSE147794", # finished, some are not success, need rerun, done
  "GSE148215", # under run, done
  "GSE153421", # under run, done
  "GSE163314", # under run, done
  "GSE163633", # under run, some are not success, need rerun, under run, done
  "GSE164690", # under run, some are not success, need rerun, done
  "GSE167825", # under run, some are not success, need rerun, under run, done
  "GSE168453", # under run, some are not success, need rerun, under run, done
  "GSE174125", # under run, done
  "GSE184703" # under run, some are not success, need rerun, under run, done
))


gseids_meta_scfoundation <- tibble::tibble(
  GSE_ID = c(
    # scfoundation2
    "GSE143353", # done
    "GSE147794", # finished, some are not success, need rerun, done
    "GSE148215", # under run, done
    "GSE153421", # under run, done
    "GSE163314", # under run, done
    "GSE163633", # under run, some are not success, need rerun, under run, done
    "GSE164690", # under run, some are not success, need rerun, done
    "GSE167825", # under run, some are not success, need rerun, under run, done
    "GSE168453", # under run, some are not success, need rerun, under run, done
    "GSE174125", # under run, done
    "GSE184703" # under run, some are not success, need rerun, under run, done
  ),
) |>
  dplyr::mutate(
    samples = purrr::map(
      GSE_ID,
      .f = \(.x) {
        basedir <- "/home/liuc9/github/scMOCHA-data/data/scfoundation2/PBMC"
        data.table::fread(
          file.path(
            basedir,
            .x,
            "out",
            glue::glue("{.x}.cell_ratio_and_variant_clean.csv")
          )
        ) -> .d
        tibble::tibble(
          samples = nrow(.d),
          Disease = "-",
          Source = "PBMC",
          Chemistry = unique(.d$Chemistry)[[1]],
          Publication = "-"
        )
      }
    ),
  ) |>
  tidyr::unnest(cols = samples)

# body --------------------------------------------------------------------
gseids_meta <- dplyr::bind_rows(
  gseids_meta_scfoundation
)

tibble::tibble(
  gseid = gseids
) |>
  dplyr::mutate(
    cell_ratio_variant = purrr::map(
      .x = gseid,
      .f = \(.gseid) {
        basedir <- "/home/liuc9/github/scMOCHA-data/data/scfoundation2/PBMC"
        data.table::fread(
          file.path(
            basedir,
            .gseid,
            "out",
            glue::glue("{.gseid}.cell_ratio_and_variant_clean.csv")
          )
        ) |>
          dplyr::select(-Chemistry)
      }
    )
  ) |>
  dplyr::mutate(
    anno = purrr::map(
      .x = gseid,
      .f = \(.gseid) {
        readr::read_rds(
          file.path(
            basedir,
            .gseid,
            "out",
            glue::glue("{.gseid}.scmocha.out.rds.gz")
          )
        )
      }
    )
  ) -> gse_data_loaded

outdir <- file.path(basedir, "out")
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

gse_data_loaded |>
  dplyr::select(-anno) |>
  tidyr::unnest(cols = cell_ratio_variant) |>
  dplyr::group_by(gseid) |>
  dplyr::summarise(
    `Avg. mutation` = mean(`# of variants`, na.rm = TRUE),
    `Avg. somatic mutation` = mean(`# of somatic variants`, na.rm = TRUE),
    `Avg. total reads` = mean(`Total reads`, na.rm = TRUE),
    `Avg. mapped reads` = mean(`Depth read mean`, na.rm = TRUE),
    `Avg. call depth` = mean(`Depth mean`, na.rm = TRUE),
  ) |>
  dplyr::left_join(
    # gseids_meta,
    gseids_meta,
    by = c("gseid" = "GSE_ID")
  ) |>
  dplyr::slice(match(gseids, gseid)) |>
  writexl::write_xlsx(
    path = file.path(outdir, "gses_meta_read.xlsx")
  )

gse_data_loaded |>
  dplyr::left_join(
    gseids_meta,
    by = c("gseid" = "GSE_ID")
  ) -> gse_cell_ratio_variant_meta

# save gse cell ratio and variant data ------------------------------------
gse_cell_ratio_variant_meta |>
  tidyr::unnest(cols = cell_ratio_variant) -> gse_cell_ratio_variant_meta_xlsx

gse_cell_ratio_variant_meta_xlsx |>
  dplyr::select(-anno) |>
  dplyr::arrange(`# of somatic variants`) |>
  writexl::write_xlsx(
    path = file.path(outdir, "gses_cell_ratio_variant_meta.xlsx")
  )

gse_cell_ratio_variant_meta_xlsx |>
  dplyr::group_by(gseid, Chemistry) |>
  dplyr::summarise(
    `Avg # of somatic variants` = mean(`# of somatic variants`, na.rm = TRUE),
    `# of samples` = dplyr::n()
  ) |>
  dplyr::ungroup() |>
  dplyr::arrange(dplyr::desc(`Avg # of somatic variants`)) |>
  dplyr::mutate(
    label = glue::glue(
      "{gseid} ({Chemistry}, {`# of samples`}, {round(`Avg # of somatic variants`, 2)})"
    )
  ) |>
  dplyr::mutate(
    label = factor(label, levels = label)
  ) -> gseid_ranked


# plot average depth and number of somatic variants ------------------------

chem_levels <- c("SC3Pv2", "SC3Pv3", "SC5P-R2", "SC5P-PE") |> rev()
chem_colors <- viridis::viridis_pal(option = "D")(4) |>
  prismatic::color()

gse_cell_ratio_variant_meta |>
  dplyr::left_join(
    gseid_ranked |> dplyr::select(gseid, label),
    by = "gseid"
  ) |>
  tidyr::unnest(cols = cell_ratio_variant) |>
  dplyr::mutate(
    Chemistry = factor(Chemistry, levels = chem_levels),
  ) -> forplot

cor.test(~ `Depth mean` + `# of somatic variants`, data = forplot) |>
  broom::tidy() -> cor_test_all

pcc <- readr::read_tsv(
  file = "https://raw.githubusercontent.com/chunjie-sam-liu/chunjie-sam-liu.life/master/public/data/pcc.tsv"
) |>
  dplyr::arrange(cancer_types)

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
    method = "loess",
    se = FALSE,
    color = "black",
    linetype = 21,
  ) +
  geom_hline(
    yintercept = 10,
    linetype = 21,
    color = "red"
  ) +
  # ggsci::scale_color_aaas(
  #   name = "GSE ID",
  # ) +
  scale_color_manual(
    name = "GSE ID",
    values = pcc$color
  ) +
  scale_x_continuous(
    labels = scales::label_number(),
    limits = c(0, 10000),
    breaks = seq(0, 10000, 1000),
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
    ),
    legend.box = "horizontal",
  ) +
  labs(
    x = "Average Depth",
    y = "Number of Somatic Variants",
    title = glue::glue(
      "All samples test, Pearson's r = {round(cor_test_all$estimate, 2)}, p-value = {scales::pvalue(cor_test_all$p.value)}"
    )
  ) -> p_somatic_variant

ggsave(
  filename = file.path(outdir, "somatic_variant_with_avg_depth.pdf"),
  plot = p_somatic_variant,
  width = 17,
  height = 7,
  dpi = 300
)

# footer ------------------------------------------------------------------

# future: :plan(future: :sequential)

# save image --------------------------------------------------------------
