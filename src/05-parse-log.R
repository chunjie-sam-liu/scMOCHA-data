#!/usr/bin/env Rscript
# Metainfo ----------------------------------------------------------------

# @AUTHOR: Chun-Jie Liu
# @CONTACT: chunjie.sam.liu.at.gmail.com
# @DATE: Wed Sep  4 14:07:48 2024
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

fn_parse_log <- function(logfile) {
  output_lines <- character()

  con <- file(logfile, "r")
  on.exit(close(con))

  while (length(line <- readLines(con, n = 1, warn = FALSE)) > 0) {
    if (grepl("scMOCHA.output_dir_tar_gz", line)) {
      cleaned_line <- stringr::str_remove_all(line, "\"| |,|scMOCHA.output_dir_tar_gz|\\:")
      output_lines <- c(output_lines, cleaned_line)
    }
  }

  output_lines
}

# load data ---------------------------------------------------------------
# basedir <- "/home/liuc9/github/scMOCHA-data/data"
basedir <- "/mnt/isilon/u01_project/large-scale/liuc9/raw"
datadir <- file.path(
  basedir, gseid
)

dir.create(
  path = datadir,
  showWarnings = F,
  recursive = T
)

targzdir <- file.path(
  datadir, "targz"
)
dir.create(
  path = targzdir,
  showWarnings = F,
  recursive = T
)

# body --------------------------------------------------------------------
# parse log file ----------------------------------------------------------


logfile <- file.path(
  datadir,
  "04.{gseid}.batch.log" |> glue::glue()
)

log_warn("load the logfile ", logfile)

targz <- fn_parse_log(logfile = logfile)
# log_warn("Found the tar.gz files", targz)

readr::write_lines(
  x = targz,
  file = file.path(
    datadir,
    "{gseid}.scmocha.targz.txt" |> glue::glue()
  )
)


# cp tar gz file into targz for backup ------------------------------------

cmd_cp_targz <- glue::glue("cp {targz} {targzdir} &")

readr::write_lines(
  c(cmd_cp_targz),
  file = file.path(
    datadir,
    "05.{gseid}.scmocha.cptargz.sh" |> glue::glue()
  )
)

# uncompress --------------------------------------------------------------
untargzdir <- file.path(
  datadir, "final"
)
dir.create(
  path = untargzdir,
  showWarnings = F,
  recursive = T
)

cped_targzs <- glue::glue("{targzdir}/{basename(targz)}")
cmd_untar <- glue::glue("tar -zxvf {cped_targzs} -C {untargzdir} &")
readr::write_lines(
  c(cmd_untar),
  file = file.path(
    datadir,
    "07.{gseid}.scmocha.untargz.sh" |> glue::glue()
  )
)

# rm fastq file -----------------------------------------------------------
gsms <- basename(targz) |>
  gsub(".tar.gz", "", x = _)

srarun <- data.table::fread(
  file.path(
    datadir,
    "{gseid}.SraRunTable.GSM" |> glue::glue()
  )
)

srarun |>
  dplyr::filter(
    experiment_name %in% gsms
  ) ->
toberemoved

runfile <- data.table::fread(
  file.path(
    datadir,
    "{gseid}.runfile.csv" |> glue::glue()
  )
)

runfile |>
  # dplyr::filter(srafile_exist) |>
  dplyr::filter(srrid %in% toberemoved$run_accession) |>
  dplyr::select(srrdir) |>
  dplyr::mutate(
    rm_cmd = purrr::map_chr(
      .x = srrdir,
      .f = \(.srrdir) {
        .srrid <- basename(.srrdir)
        cmd <- "rm {.srrdir}/{.srrid}*.fastq" |> glue::glue()
        cmd
      }
    )
  ) |>
  dplyr::select(rm_cmd) ->
rm_cmds


cmd_rm_fastq <- glue::glue("{rm_cmds$rm_cmd} &")

readr::write_lines(
  c(cmd_rm_fastq),
  file = file.path(
    datadir,
    "06.{gseid}.scmocha.clear.sh" |> glue::glue()
  )
)


# footer ------------------------------------------------------------------

# future::plan(future::sequential)

# save image --------------------------------------------------------------
