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
#library(rlang)
library(GetoptLong)
library(logger)

# args --------------------------------------------------------------------


gseid <- "GSE163668"

# s: string, i: integer, f: float, !: boolean
# @: array
# %: hash
# default: default value specified here.
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


# load data ---------------------------------------------------------------

log_warn(gseid)
basedir <- "/home/liuc9/github/scMOCHA-data/data"
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
    "{gseid}.metadata.gsm.csv" |> glue::glue()
  )
)

gsm |>
  dplyr::rename(
    srrid = run_accession
  ) |>
  dplyr::inner_join(
    runfile, by = "srrid"
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
          .x
        )

        dir.create(
          gsmdir,
          showWarnings = F,
          recursive = T
        )


        .y |>
          tibble::rowid_to_column() |>
          dplyr::mutate(
            rename = purrr::map2(
              .x = rowid,
              .y = srrdir,
              .f = \(.rowid, .srrdir) {
                # .rowid <- a$rowid[[1]]
                # .srrdir <- a$srrdir[[1]]

                .srrid <- basename(.srrdir)
                .from_R1 <- file.path(
                  .srrdir,
                  "{.srrid}_1.fastq" |> glue::glue()
                )
                .to_R1 <- file.path(
                  gsmdir,
                  "{.x}_S1_L00{.rowid}_R1_001.fastq" |> glue::glue()
                )
                if(file.exists(.to_R1)) {
                  file.remove(.to_R1)
                }
                file.symlink(
                  from = .from_R1,
                  to = .to_R1
                )
                .from_R2 <- file.path(
                  .srrdir,
                  "{.srrid}_2.fastq" |> glue::glue()
                )
                .to_R2 <- file.path(
                  gsmdir,
                  "{.x}_S1_L00{.rowid}_R2_001.fastq" |> glue::glue()
                )
                if(file.exists(.to_R2)) {
                  file.remove(.to_R2)
                }
                file.symlink(
                  from = .from_R2,
                  to = .to_R2
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
