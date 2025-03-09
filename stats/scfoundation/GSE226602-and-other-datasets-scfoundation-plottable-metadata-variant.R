#!/usr/bin/env Rscript
# Metainfo ----------------------------------------------------------------
# @AUTHOR: Chun-Jie Liu
# @CONTACT: chunjie.sam.liu.at.gmail.com
# @DATE: 2025-03-07 16:21:52
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

fn_load_hetero <- function(.filename) {
  # .filename <- file.path(
  #   sc_dir,
  #   "mgatk_out/final",
  #   "sc.cell_heteroplasmic_df.tsv.gz"
  # )

  data.table::fread(input = .filename) -> .d

  data.table::setnames(.d, "V1", "barcode")
  .d <- data.table::melt(.d, id.vars = "barcode", variable.name = "variant", value.name = "af")
  .d[, pos := gsub(pattern = ">|[AGCT]", replacement = "", x = variant)]
  .d[, pos := as.integer(pos)]
  .d
}

fn_load_cluster <- function(.filename) {
  data.table::fread(
    input = .filename,
    sep = "\t",
    col.names = c("barcode", "tag", "celltype")
  ) |>
    dplyr::arrange(celltype) |>
    dplyr::mutate(celltype = factor(celltype)) |>
    dplyr::select(-tag)
}

fn_load_coverage <- function(.filename) {
  data.table::fread(
    input = .filename,
    sep = ",",
    col.names = c("pos", "barcode", "depth")
  )
}

fn_load_meta <- function(.filename) {
  data.table::fread(
    input = .filename,
    sep = "\t"
  ) |>
    dplyr::rename(
      barcode = cellbarcode
    ) |>
    dplyr::select(-orig.ident)
}

fn_forplot <- function(.af, .coverage, .meta) {
  # print(.meta)
  .af |>
    dplyr::select(barcode, cluster, dplyr::contains(">")) |>
    tidyr::pivot_longer(
      cols = -c(barcode, cluster),
      names_to = "variant",
      values_to = "af"
    ) |>
    dplyr::group_by(barcode, cluster) |>
    dplyr::summarise(s_af = sum(af, na.rm = T)) |>
    dplyr::ungroup() |>
    dplyr::arrange(cluster, -s_af) ->
  .rank

  .af |>
    dplyr::select(barcode, dplyr::contains(">")) |>
    tidyr::pivot_longer(
      cols = -barcode,
      names_to = "variant",
      values_to = "af"
    ) |>
    dplyr::mutate(
      pos = gsub(pattern = "([[:digit:]]*).*", "\\1", variant) |>
        as.numeric()
    ) |>
    dplyr::left_join(
      .coverage,
      by = c("barcode", "pos")
    ) |>
    tidyr::replace_na(
      replace = list(
        af = 0
      )
    ) |>
    dplyr::mutate(af = ifelse(is.na(depth), NA, af)) |>
    # dplyr::mutate(af = ifelse(depth < log2(10), -0.1, af)) |>
    dplyr::mutate(af = ifelse(depth < 10, -0.1, af)) |>
    dplyr::arrange(pos) ->
  .forplot

  .coverage |>
    dplyr::group_by(barcode) |>
    dplyr::summarise(sum_depth = sum(depth, na.rm = TRUE)) ->
  .coverage_cell

  list(
    rank = .rank,
    forplot = .forplot,
    meta = .meta,
    coverage_cell = .coverage_cell
  )
}

# load data ---------------------------------------------------------------

basedir <- "/home/liuc9/github/scMOCHA-data/data"
foundation_out <- file.path(basedir, "scfoundation/out")

gse_dataset_metadata_full <- data.table::fread(
  file.path(foundation_out, "gse_dataset_metadata_full.csv")
)
gseids <- c(
  "GSE155673",
  "GSE157344",
  "GSE149689",
  "GSE171555",
  "GSE155223",
  "GSE163668",
  "GSE175524",
  "GSE206283",
  "GSE226598",
  "GSE261140",
  "GSE279945",
  "GSE214865",
  "GSE220189",
  "GSE233844",
  "GSE175499",
  "GSE149313",
  "GSE154386",
  "GSE159117",
  "GSE188632",
  "GSE166992",
  "GSE162117",
  "GSE226602",
  "GSE161354",
  "GSE235050",
  "GSE181279"
)

# body --------------------------------------------------------------------
tibble::tibble(
  gseid = gseids
) |>
  dplyr::mutate(
    anno = purrr::map(
      .x = gseid,
      .f = \(.gseid) {
        .anno <- readr::read_rds(
          file.path(basedir, .gseid, "out", glue::glue("{.gseid}.scmocha.out.rds.gz"))
        )
      }
    )
  ) ->
gse_data_loaded


gse_data_loaded |>
  tidyr::unnest(cols = anno) ->
gse_data

gse_data |>
  dplyr::select(somatic_variant) |>
  dplyr::mutate(
    somatic_variant_new = purrr::map(
      .x = somatic_variant,
      .f = \(.x) {
        .x$somatic
      }
    )
  ) |>
  dplyr::pull(somatic_variant_new) |>
  purrr::reduce(union) ->
somatic_variants

gse_data |>
  dplyr::select(gseid, srrid, srrdir) |>
  dplyr::mutate(
    raw = parallel::mclapply(
      X = srrdir,
      FUN = function(.srrdir) {
        barcode_cluster_file <- "barcode_cluster.tsv"
        cluster_umap <- fn_load_cluster(
          .filename = file.path(.srrdir, barcode_cluster_file)
        )

        cell_hetero_raw_file <- "cell.cell_heteroplasmic_df_raw.tsv.gz"
        cell_hetero_raw <- fn_load_hetero(
          .filename = file.path(.srrdir, cell_hetero_raw_file)
        ) |>
          dplyr::filter(variant %in% somatic_variants)

        cell_coverage_file <- "cell.coverage.txt.gz"
        cell_coverage <- fn_load_coverage(
          .filename = file.path(.srrdir, cell_coverage_file)
        )

        cell_meta_data_file <- "cell_meta_data.tsv"
        metadata <- fn_load_meta(
          .filename = file.path(.srrdir, cell_meta_data_file)
        )


        cell_raw_cluster_af <- cluster_umap |>
          dplyr::left_join(cell_hetero_raw, by = "barcode") |>
          dplyr::rename(cluster = celltype) |>
          tidyr::pivot_wider(
            names_from = variant,
            values_from = af
          )

        cell_raw_cluster_forplot <- fn_forplot(
          .af = cell_raw_cluster_af,
          .coverage = cell_coverage,
          .meta = metadata
        )

        list(
          cluster_umap = cluster_umap,
          cell_hetero_raw = cell_hetero_raw,
          cell_coverage = cell_coverage,
          metadata = metadata,
          cell_raw_cluster_af = cell_raw_cluster_af,
          cell_raw_cluster_forplot = cell_raw_cluster_forplot
        )
      },
      mc.cores = 100
    )
  ) ->
gse_data_af

readr::write_rds(
  gse_data_af,
  file.path(
    foundation_out,
    "GSE226602-and-other-datasets-scfoundation-plottable-metadata-variant.R.gse_data_af.rds.gz"
  )
)

# footer ------------------------------------------------------------------

# future: :plan(future: :sequential)

# save image --------------------------------------------------------------
save.image(
  file = file.path(
    foundation_out,
    "/home/liuc9/github/scMOCHA-data/stats/scfoundation/GSE226602-and-other-datasets-scfoundation-plottable-metadata-variant.R.RData"
  )
)
