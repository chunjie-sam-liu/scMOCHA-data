#!/usr/bin/env Rscript
# Metainfo ----------------------------------------------------------------

# @AUTHOR: Chun-Jie Liu
# @CONTACT: chunjie.sam.liu.at.gmail.com
# @DATE: Sun Aug 18 20:56:24 2024
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
# body --------------------------------------------------------------------

sratable <- data.table::fread(
  file.path(
    datadir,
    "{gseid}.metadata.csv" |> glue::glue()
  )
)

sratable |>
  dplyr::select(run_accession = Run, experiment_name = `Sample Name`, experiment_accession = Experiment) |>
  data.table::fwrite(
    file = file.path(
      datadir,
      "{gseid}.metadata.gsm.csv" |> glue::glue()
    )
  )

sratable |>
  dplyr::select(srrid = Run) |>
  dplyr::mutate(
    srrdir = file.path(
      datadir, srrid
    )
  ) |>
  dplyr::mutate(
    srrdir_exists = file.exists(srrdir)
  ) |>
  dplyr::mutate(
    prefetch = "prefetch --max-size 50G {srrid} --output-directory {datadir}" |> glue::glue()
  ) |>
  dplyr::mutate(
    srafile = file.path(
      srrdir,
      glue::glue("{srrid}.sralite")
    )
  ) |>
  dplyr::mutate(
    srafile_exist = file.exists(srafile)
  ) ->
  srafiles

readr::write_lines(
  glue::glue("{ srafiles$prefetch} &"),
  file = file.path(
    datadir,
    "00.{gseid}.prefetch.sh" |> glue::glue()
  )
)

cmd_slrm <- c(
  "#!/usr/bin/env bash",
  "# @AUTHOR: Chun-Jie Liu",
  "# @CONTACT: chunjie.sam.liu.at.gmail.com",
  "# @DATE: {lubridate::now()}" |> glue::glue(),
  "",
  "#SBATCH --signal=USR2",
  "#SBATCH --ntasks=1",
  "#SBATCH --cpus-per-task=10",
  "#SBATCH --mem=50G",
  "#SBATCH --time=720:00:00",
  "#SBATCH --output={datadir}/01.{gseid}.dump.job.%j" |> glue::glue(),
  "#module load R/4.1.0"
)


srafiles |>
  dplyr::mutate(
    dump_cmd = purrr::map2_chr(
      .x = srrdir,
      .y = srafile,
      .f = \(.x, .y) {
        .srrid <- basename(.x)

        # dir.create(
        #   path = .x,
        #   showWarnings = F,
        #   recursive = T
        # )

        cmd_dump <- c(
          "fasterq-dump {.y} --temp /scr1/users/liuc9/tmp/fasterq_dump  --include-technical --mem 50G --threads 10 --split-files --outdir {.x}" |> glue::glue()
        )

        cmd <- c(
          # cmd_slrm,
          cmd_dump
        )

        # dump_slrm_file <- file.path(
        #   .x,
        #   "dump_{.srrid}.slrm" |> glue::glue()
        # )
        # readr::write_lines(
        #   cmd,
        #   file = dump_slrm_file
        # )

        cmd
      }
    )
  ) ->
  srafile_dump

readr::write_lines(
  c(
    cmd_slrm,
    srafile_dump$dump_cmd
  ),
  file = file.path(
    datadir,
    "01.{gseid}.dump.slrm" |> glue::glue()
  )
)

readr::write_lines(
  glue::glue("{ srafile_dump$dump_cmd} &"),
  file = file.path(
    datadir,
    "01.{gseid}.dump.sh" |> glue::glue()
  )
)

data.table::fwrite(
  x = srafile_dump,
  file = file.path(
    datadir,
    "{gseid}.runfile.csv" |>glue::glue()
  )
)






# footer ------------------------------------------------------------------

# future::plan(future::sequential)

# save image --------------------------------------------------------------
