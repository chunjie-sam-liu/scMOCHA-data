#!/usr/bin/env Rscript
# Metainfo ----------------------------------------------------------------
# @AUTHOR: Chun-Jie Liu
# @CONTACT: chunjie.sam.liu.at.gmail.com
# @DATE: 2025-01-02 12:49:40
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

basedir <- "/home/liuc9/github/scMOCHA-data/data/scfoundation"

readr::read_lines("/home/liuc9/github/scMOCHA-data/data/scfoundation/cmds.sh") |> gsub(
  pattern = "Rscript /home/liuc9/github/scMOCHA-data/src/01-sra-metadata.R -g | -b /home/liuc9/github/scMOCHA-data/data/scfoundation",
  replacement = "",
  x = _
) -> gseids

# body --------------------------------------------------------------------
sqlite_file <- "/mnt/isilon/xing_lab/liuc9/refdata/sradb/SRAmetadb.sqlite"

sra_con <- DBI::dbConnect(
  RSQLite::SQLite(),
  sqlite_file
)
DBI::dbListTables(sra_con)

sra_table <- dplyr::tbl(
  sra_con, "sra"
)
study_table <- dplyr::tbl(
  sra_con, "study"
)


sra_table |>
  dplyr::filter(study_name %in% gseids) |>
  dplyr::inner_join(study_table, by = "study_accession") |>
  as.data.table() ->
sra_df

sra_df |>
  dplyr::pull(sample_accession) |>
  unique() ->
sample_accessions

sample_table <- dplyr::tbl(
  sra_con, "sample"
)

sample_df <- sample_table |>
  dplyr::filter(
    sample_accession %in% sample_accessions
  ) |>
  as.data.table()

cleaned_sample_df <- sample_df %>%
  dplyr::select(-which(apply(is.na(.), 2, all))) |>
  dplyr::select(sample_accession, sample_attribute)

cleaned_sample_df$sample_attribute[[1]] |>
  stringr::str_split(pattern = "\\|\\|", simplify = T) |>
  stringr::str_split(pattern = ": ") |>
  purrr::map(.f = \(.x) {
    .x[[1]]
  }) |>
  purrr::reduce(.f = c) ->
new_columns

cleaned_sample_df |>
  tidyr::separate(
    col = sample_attribute,
    into = new_columns,
    sep = " \\|\\|"
  ) %>%
  dplyr::mutate_at(
    new_columns,
    ~ stringr::str_remove(., ".*:")
  ) ->
cleaned_sample_df
DBI::dbDisconnect(sra_con)

sra_df |>
  dplyr::inner_join(
    cleaned_sample_df,
    by = "sample_accession"
  ) |>
  dplyr::relocate(
    run_accession,
    .before = 1
  ) ->
cleaned_sample_df_sra
#
# footer ------------------------------------------------------------------

# future: :plan(future: :sequential)

# save image --------------------------------------------------------------
