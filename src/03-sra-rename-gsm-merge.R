#!/usr/bin/env Rscript
# Metainfo ----------------------------------------------------------------

# @AUTHOR: Chun-Jie Liu
# @CONTACT: chunjie.sam.liu.at.gmail.com
# @DATE: Sun Aug 18 21:18:32 2024
# @DESCRIPTION: filename

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


# gseid <- "GSE226602"

# s: string, i: integer, f: float, !: boolean
# @: array
# %: hash
# default: default value specified here.
basedir <- "/mnt/isilon/u01_project/large-scale/liuc9/raw"
verbose <- FALSE
spec <- "
Usage: Rscript foorbar.R [options]

Options:
<gseid=s> gseid
<basedir=s> basedir
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

# rename ------------------------------------------------------------------

# https://www.10xgenomics.com/support/single-cell-gene-expression/documentation/steps/sequencing/sequencing-requirements-for-single-cell-3

# v4 r1=28, i7=10, i5=10, r2=90
# v3/v3.1 r1=28, i7=10, i5=10,r2=90
# v3/v3.1 r1=28, i7=8, i5=0,r2=91
# v2, r1=26, i7=8, i5=0, r2=98

rename_code <- c(
  "26" = "R1",
  "28" = "R1",
  "30" = "R1",
  "8" = "I1",
  "10" = "I1",
  "90" = "R2",
  "91" = "R2",
  "98" = "R2",
  "150" = "R2"
)

rename_code2 <- c()
# function ----------------------------------------------------------------

fn_get_fastq_read_length <- function(.fastq) {
  # .fastq <- .fastqs$fastq[[1]]

  .line <- readr::read_lines(file = .fastq, n_max = 1)
  stringr::str_split(.line, "length=", simplify = T)[[2]]
}

fn_rename <- function(.srrdir) {
  .fastqs <-
    tibble::tibble(
      from = list.files(
        path = .srrdir,
        pattern = "fastq$",
        full.names = T
      )
    )


  .fastqs |>
    dplyr::mutate(
      rl = purrr::map_chr(
        .x = from,
        .f = fn_get_fastq_read_length
      )
    ) |>
    dplyr::arrange(rl) ->
  .fastqs_rl

  # .fastqs_rl |>
  #   dplyr::mutate(
  #     read_type = plyr::revalue(
  #       x = rl,
  #       replace = rename_code
  #     )
  #   ) ->
  #   .fastqs_rl_rt

  .fastqs_rl |>
    dplyr::mutate(
      rl = as.integer(rl)
    ) |>
    dplyr::mutate(
      rt = dplyr::case_when(
        rl < 20 ~ "I1",
        rl <= 30 ~ "R1",
        TRUE ~ "R2"
      )
    ) |>
    dplyr::group_by(rt) |>
    dplyr::mutate(n = dplyr::n()) |>
    dplyr::mutate(idx = 1:dplyr::n()) |>
    dplyr::ungroup() |>
    tidyr::separate(
      col = rt,
      into = c("ir", "irn"),
      sep = -1,
      remove = F
    ) |>
    dplyr::mutate(read_type = ifelse(
      n == 1,
      glue::glue("{ir}{irn}"),
      glue::glue("{ir}{idx}")
    )) ->
  .fastqs_rl_rt

  .fastqs_rl_rt
}

# load data ---------------------------------------------------------------

log_warn(gseid)
# basedir <- "/home/liuc9/github/scMOCHA-data/data"
# basedir <- "/mnt/isilon/u01_project/large-scale/liuc9/raw"
datadir <- file.path(
  basedir, gseid
)

dir.create(
  path = datadir,
  showWarnings = F,
  recursive = T
)

runfile <- data.table::fread(
  file.path(
    datadir,
    "{gseid}.runfile.csv" |> glue::glue()
  )
)

# body --------------------------------------------------------------------
gsm <- data.table::fread(
  file.path(
    datadir,
    "{gseid}.SraRunTable.GSM" |> glue::glue()
  )
)

gsm |>
  dplyr::rename(
    srrid = run_accession
  ) |>
  dplyr::inner_join(
    runfile,
    by = "srrid"
  ) |>
  dplyr::group_by(experiment_name) |>
  tidyr::nest() |>
  dplyr::ungroup() ->
gsm_nest


gsm_nest |>
  dplyr::mutate(
    rename = purrr::map2(
      .x = experiment_name,
      .y = data,
      .f = \(.x, .y) {
        # .x <- gsm_nest$experiment_name[[1]]
        # .y <- gsm_nest$data[[1]]


        gsmdir <- file.path(
          datadir,
          "gsm",
          .x
        )

        dir.create(
          gsmdir,
          showWarnings = F,
          recursive = T
        )


        .y |>
          tibble::rowid_to_column() ->
        .y_idx

        .y_idx |>
          dplyr::mutate(
            rename = purrr::map2(
              .x = rowid,
              .y = srrdir,
              .f = \(.rowid, .srrdir) {
                # .rowid <- .y_idx$rowid[[1]]
                # .srrdir <- .y_idx$srrdir[[1]]

                .rt <- fn_rename(.srrdir)

                .rt |>
                  dplyr::mutate(
                    targetname = "{.x}_S1_L00{.rowid}_{read_type}_001.fastq" |> glue::glue()
                  ) |>
                  dplyr::mutate(
                    to = file.path(gsmdir, targetname)
                  ) ->
                .rt_from_to

                .rt_from_to |>
                  dplyr::mutate(
                    a = purrr::map2(
                      .x = from,
                      .y = to,
                      .f = \(.from, .to) {
                        if (file.exists(.to)) {
                          file.remove(.to)
                        }
                        log_success(.from)
                        log_success(.to)
                        file.symlink(
                          from = .from,
                          to = .to
                        )
                      }
                    )
                  )
              }
            )
          )
      }
    )
  )


# footer ------------------------------------------------------------------

# future::plan(future::sequential)

# save image --------------------------------------------------------------
