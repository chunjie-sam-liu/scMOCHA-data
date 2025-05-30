#!/usr/bin/env Rscript
# Metainfo ----------------------------------------------------------------
# @AUTHOR: Chun-Jie Liu
# @CONTACT: chunjie.sam.liu.at.gmail.com
# @DATE: 2025-02-24 16:22:48
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
project_filename <- "/home/liuc9/github/scMOCHA-data/data/scfoundation/out/project_source.csv"
project_source <- data.table::fread(project_filename)

# body --------------------------------------------------------------------

sqlite_file <- "/mnt/isilon/xing_lab/liuc9/refdata/sradb/SRAmetadb.sqlite"

sra_con <- DBI::dbConnect(
  RSQLite::SQLite(),
  sqlite_file
)
# DBI::dbDisconnect(sra_con)

DBI::dbListTables(sra_con)

sra_table <- dplyr::tbl(sra_con, "sra")
study_table <- dplyr::tbl(sra_con, "study")
sample_table <- dplyr::tbl(sra_con, "sample")

col_inter <- setdiff(intersect(colnames(sra_table), colnames(study_table)), "study_accession")


proj_IDs <- project_source$proj_ID

sra_table |>
  dplyr::filter(study_name %in% proj_IDs) |>
  dplyr::inner_join(
    study_table |>
      dplyr::select(-col_inter),
    by = "study_accession"
  ) |>
  as.data.table() ->
sra_df

sra_df |>
  dplyr::pull(sample_accession) |>
  unique() ->
sample_accessions


sample_df <- sample_table |>
  dplyr::filter(
    sample_accession %in% sample_accessions
  ) |>
  as.data.table()

DBI::dbDisconnect(sra_con)

cleaned_sample_df <- sample_df %>%
  dplyr::select(-which(apply(is.na(.), 2, all))) |>
  dplyr::select(sample_accession, sample_attribute)


cleaned_sample_df |>
  dplyr::mutate(
    sa = parallel::mclapply(
      X = sample_attribute,
      FUN = function(.x) {
        .x |>
          stringr::str_split(" \\|\\| ") |>
          _[[1]] ->
        .xx

        tibble::tibble(
          xx = .xx
        ) |>
          tidyr::separate(
            col = xx,
            into = c("key", "value"),
            sep = ": ",
            remove = TRUE
          ) ->
        .d

        .d |> dplyr::filter(key == "source_name") -> .dd

        tibble::tibble(
          source_name = .dd$value,
          sample_attribute_new = list(.d)
        ) |>
          dplyr::mutate(
            source_name = stringr::str_trim(source_name)
          )
      },
      mc.cores = 20
    )
  ) |>
  tidyr::unnest(cols = sa) ->
cleaned_sample_source_name

# study_name
# experiment_name

sra_df |>
  dplyr::select(-sample_attribute) |>
  dplyr::inner_join(
    cleaned_sample_source_name,
    by = c("sample_accession" = "sample_accession")
  ) |>
  dplyr::mutate(
    proj_ID = study_name,
    samp_ID = experiment_name
  ) |>
  dplyr::relocate(
    proj_ID, samp_ID,
    .before = 1
  ) ->
cleaned_sample_df_sra

cleaned_sample_df_sra |> dplyr::glimpse()



project_source |>
  dplyr::left_join(
    cleaned_sample_df_sra,
    by = c("proj_ID" = "proj_ID", "samp_ID" = "samp_ID")
  ) ->
project_source_sra

project_source_sra |> dplyr::glimpse()

readr::write_rds(
  project_source_sra,
  "/home/liuc9/github/scMOCHA-data/data/scfoundation/out/project_source_sra.rds.gz"
)

# footer ------------------------------------------------------------------

# future: :plan(future: :sequential)

# save image --------------------------------------------------------------
