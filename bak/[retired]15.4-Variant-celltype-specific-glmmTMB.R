#!/usr/bin/env Rscript
# Metainfo ----------------------------------------------------------------
# @AUTHOR: Chun-Jie Liu
# @CONTACT: chunjie.sam.liu.at.gmail.com
# @DATE: 2025-06-16 14:32:45
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
library(glue)
library(parallel)
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

# load data ---------------------------------------------------------------

cleandatadir <- "/home/liuc9/github/scMOCHA-data/data/zzz/clean-data"


dbdir <- "/home/liuc9/github/scMOCHA-data/analysis/zzz/db"
ks_test_dir <- file.path(dbdir, "all_hetero_af.cell.ks_test")
plotdir <- "/home/liuc9/github/scMOCHA-data/analysis/zzz/plot-celltype-specific-variant"


META <- import(
  "/home/liuc9/github/scMOCHA-data/analysis/zzz/clean-data/gse_dataset_metadata_full.sex_pred.qs"
) |>
  dplyr::select(
    gseid,
    srrid,
    Age_new,
    Age_group,
    Haplogroup,
    disease,
    Chemistry,
    sex_pred
  ) |>
  dplyr::mutate(
    Haplogroup = purrr::map_chr(
      .x = Haplogroup,
      .f = \(.x) {
        # if (stringr::str_starts(.x, "L")) {
        #   gsub("L", "L0", .x)
        # }
        gsub("\\d+.*", "", .x)
      }
    )
  )


conn <- DBI::dbConnect(
  duckdb::duckdb(),
  dbdir = "/home/liuc9/github/scMOCHA-data/analysis/zzz/clean-data/all_hetero_af.cell.duckdb.1.2.1" |>
    glue::glue()
)

all_hetero_af_cell_tbl <- dplyr::tbl(conn, "all_hetero_af_cell")

all_hetero_af_cell_tbl |>
  dplyr::select(variant) |>
  dplyr::distinct() |>
  dplyr::collect() -> all_hetero_af_cell_variants
all_hetero_af_cell_tbl |>
  dplyr::select(celltype) |>
  dplyr::distinct() |>
  dplyr::collect() -> all_hetero_af_cell_celltypes


# function ----------------------------------------------------------------

# body --------------------------------------------------------------------

META |>
  dplyr::select(-Age_group) |>
  dplyr::mutate(
    disease = as.character(disease),
  ) |>
  dplyr::mutate(
    disease = ifelse(disease == "Unknown", NA_character_, disease),
  ) |>
  dplyr::mutate(
    Chemistry = as.character(Chemistry),
  ) -> META_age

library(glmmTMB)
library(DHARMa)

all_hetero_af_cell_variants |>
  dplyr::mutate(
    model = parallel::mclapply(
      X = variant,
      FUN = function(.thevariant) {
        # .thevariant <- all_hetero_af_cell_variants$variant[[1]]
        all_hetero_af_cell_tbl |>
          dplyr::filter(
            variant == .thevariant,
            # celltype == thecelltype,
            af > 0
          ) |>
          as.data.table() -> .d
        .d |>
          dplyr::inner_join(
            META_age,
            by = c("gseid", "srrid")
          ) |>
          dplyr::mutate(
            af = ifelse(af == 1, 0.9999, af),
          ) -> .dd

        model <- glmmTMB(
          af ~ celltype +
            Age_new +
            disease +
            Chemistry +
            Haplogroup +
            sex_pred +
            (1 | srrid),
          family = beta_family(),
          data = .dd
        )
        model
      },
      mc.cores = 50
    )
  ) -> all_hetero_af_cell_variants_models

export(
  all_hetero_af_cell_variants_models,
  file = file.path(
    "/home/liuc9/github/scMOCHA-data/analysis/zzz/db/all_hetero_af.cell.glmmTMB",
    "all_hetero_af.cell.variants.models.qs"
  )
)


all_hetero_af_cell_variants_models |>
  # head(2) |>
  dplyr::mutate(
    params = parallel::mclapply(
      X = model,
      FUN = \(.model) {
        # summary(.model)
        if (performance::check_singularity(.model)) {
          logger::log_warn("Model is singular, skipping parameter extraction.")
          return(NULL)
        }
        parameters::model_parameters(.model) |> as.data.table()
      },
      mc.cores = 20
    )
  ) -> all_hetero_af_cell_variants_models_params

all_hetero_af_cell_variants_models_params |>
  tidyr::unnest(cols = c(params)) |>
  dplyr::filter(
    p < 0.05
  ) -> all_hetero_af_cell_variants_models_params_filtered

lymphoid_cells <- c(
  "B",
  "CD4_T",
  "CD8_T",
  "other_T",
  "NK"
)
myeloid_cells <- c(
  "DC",
  "Mono"
)

all_hetero_af_cell_variants_models_params_filtered |>
  dplyr::filter(
    grepl("celltype", Parameter) |
      grepl("Age_new", Parameter)
  ) |>
  dplyr::select(
    variant,
    Parameter,
    Coefficient,
    p
  ) |>
  tidyr::nest(
    .by = variant,
    .key = "params"
  ) |>
  dplyr::mutate(
    params_new = parallel::mclapply(
      X = params,
      FUN = function(.x) {
        cell_age <- gsub("celltype", "", .x$Parameter)
        tibble::tibble(
          lymphoid_cell = intersect(lymphoid_cells, cell_age) |> length(),
          myeloid_cell = intersect(myeloid_cells, cell_age) |> length(),
          cell_age = intersect("Age_new", cell_age) |> length(),
        )
      },
      mc.cores = 20
    )
  ) |>
  dplyr::select(-params) |>
  tidyr::unnest(
    cols = c(params_new)
  ) -> all_hetero_af_cell_variants_models_params_filtered_sorted

all_hetero_af_cell_variants_models_params_filtered_sorted |>
  dplyr::filter(cell_age > 0)


source("analysis/00-colors.R")

.thevariant <- "3727T>C"
.thevariant <- "12367A>G"
.thevariant <- "4886C>T"
.thevariant <- "1082A>G"
.thevariant <- "8852G>A"
.thevariant <- "2091A>T"

all_hetero_af_cell_tbl |>
  dplyr::filter(
    variant == .thevariant,
    # celltype == thecelltype,
    af > 0
  ) |>
  as.data.table() |>
  dplyr::inner_join(
    META_age,
    by = c("gseid", "srrid")
  ) |>
  dplyr::filter(af > 0) |>
  dplyr::mutate(
    celltype = gsub(
      "_",
      " ",
      celltype
    )
  ) |>
  dplyr::mutate(
    celltype = factor(
      celltype,
      levels = names(color_celltype) |> rev()
    )
  ) |>
  ggplot(aes(
    x = af,
    y = celltype,
    fill = celltype
  )) +
  ggjoy::geom_joy(
    # scale = 0.9,
    # alpha = 0.8,
    rel_min_height = 0.01,
    size = 0.1
  ) +
  scale_fill_manual(
    values = color_celltype,
    na.value = "grey50"
  ) +
  ggjoy::theme_joy() +
  theme(
    legend.position = "none",
    plot.title = element_text(
      hjust = 0.5,
      size = 16
    ),
  ) +
  labs(
    x = "Allele Frequency",
    y = "Cell Type"
  )

# footer ------------------------------------------------------------------

# future: :plan(future: :sequential)

# save image --------------------------------------------------------------
