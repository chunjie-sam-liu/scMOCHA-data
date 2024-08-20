#!/usr/bin/env Rscript
# Metainfo ----------------------------------------------------------------

# @AUTHOR: Chun-Jie Liu
# @CONTACT: chunjie.sam.liu.at.gmail.com
# @DATE: Fri Aug 16 17:45:09 2024
# @DESCRIPTION: filename

# Library -----------------------------------------------------------------

suppressPackageStartupMessages(library(magrittr))
library(ggplot2)
library(patchwork)
library(prismatic)
library(paletteer)
library(data.table)
#library(rlang)
library(GetoptLong)
library(logger)

# args --------------------------------------------------------------------

# s: string, i: integer, f: float, !: boolean
# @: array
# %: hash
# default: default value specified here.

gseid <- "GSE163668"

verbose <- FALSE
spec <- "
Usage: Rscript foorbar.R [options]

Options:
<gseid=s> gseid
<verbose!> Print messages
"

GetoptLong.options(help_style = "two-column")
GetoptLong(spec, template_control = list(opt_width = 21))

# src ---------------------------------------------------------------------


# header ------------------------------------------------------------------
log_threshold(TRACE)
log_layout(layout_glue_colors)

# log_info('Starting the script...')
# log_debug('This is the second log line')
# log_trace('Note that the 2nd line is being placed right after the 1st one.')
# log_success('Doing pretty well so far!')
# log_warn('But beware, as some errors might come :/')
# log_error('This is a problem')
# log_debug('Note that getting an error is usually bad')
# log_error('This is another problem')
# log_fatal('The last problem')

# future::plan(future::multisession, workers = 10)

# function ----------------------------------------------------------------


fn_get_html <- function(theid, prefix, datadir) {
  # .url <- glue::glue("https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc={theid}") # geo
  # .url <- glue::glue("https://www.ncbi.nlm.nih.gov/Traces/study/?acc={theid}") # sraruntable
  .url <- glue::glue("{prefix}={theid}")


  # .html <- rvest::read_html(.url)
  httr::GET(.url)

  .of <- file.path(
    datadir,
    "{theid}.html" |> glue::glue()
  )

  xml2::write_html(
    x = .html,
    file = .of
  )
  .of
}

fn_parse_bioprojectid <- function(.of) {
  rvest::read_html(.of) |>
    rvest::html_elements("td") %>%
    .[rvest::html_text2(.) == "BioProject"] ->
    .ele

  .ele |>
    rvest::html_element(xpath = "./parent::tr") |>
    rvest::html_elements("td") %>%
    .[[2]] ->
    .the_ele

  .projid <- rvest::html_text2(.the_ele)
  .projid


}

fn_parse_sratable <- function(.of) {
  rvest::read_html(.of) |>
    rvest::html_table()

}


# load data ---------------------------------------------------------------

basedir <- "/home/liuc9/github/scMOCHA-data/data"

# body --------------------------------------------------------------------



datadir <- file.path(
  basedir, gseid
)

dir.create(
  path = datadir,
  showWarnings = F,
  recursive = T
)

# get gsm -----------------------------------------------------------------


gse <- GEOquery::getGEO(
  GEO = gseid,
  destdir = datadir,
  GSEMatrix = TRUE
)


gse[[1]] |>
  Biobase::phenoData() |>
  Biobase::pData() ->
  gse_pheno

data.table::fwrite(
  x = gse_pheno,
  file = file.path(
    datadir,
    "{gseid}.pheno.csv" |> glue::glue()
  )
)


# strsplit(
#   x = gse@header$relation[[2]],
#   split = "term=",
# )[[1]][[2]] ->
#   srpid

# gse[[1]] |>
#   Biobase::


# sradb -------------------------------------------------------------------


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
  dplyr::filter(study_name == gseid) |>
  dplyr::inner_join(study_table, by = "study_accession") |>
  as.data.table() ->
  sra_df
#
sample_accessions <- sra_df |>
  dplyr::pull(sample_accession)

sample_table <- dplyr::tbl(
  sra_con, "sample"
)

sample_df <- sample_table |>
  dplyr::filter(
    sample_accession %in% sample_accessions
  ) |>
  as.data.table()

cleaned_sample_df <- sample_df  %>%
  dplyr::select(-which(apply(is.na(.), 2, all))) |>
  dplyr::select(sample_accession, sample_attribute)

cleaned_sample_df$sample_attribute[[1]] |>
  stringr::str_split(pattern = "\\|\\|", simplify = T) |>
  stringr::str_split(pattern = ": ") |>
  purrr::map(.f = \(.x) {.x[[1]]}) |>
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
    run_accession, .before = 1
  ) ->
  cleaned_sample_df_sra



data.table::fwrite(
  x = cleaned_sample_df_sra,
  file = file.path(
    datadir,
    "{gseid}.metadata.csv" |> glue::glue()
  )
)

# cleaned_sample_df_sra <- data.table::fread("/scr1/users/liuc9/mitochondrial/realdata/06-bigdata/GSE163668/GSE163668.metadata.csv")
cleaned_sample_df_sra |>
  dplyr::select(run_accession, experiment_name, experiment_accession) |>
  data.table::fwrite(
    file = file.path(
      datadir,
      "{gseid}.metadata.gsm.csv" |> glue::glue()
    )
  )





# footer ------------------------------------------------------------------

# future::plan(future::sequential)

# save image --------------------------------------------------------------
