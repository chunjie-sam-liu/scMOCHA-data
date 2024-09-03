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

gseid <- "GSE226602"

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
  ) ->
  sratable_prefetch

readr::write_lines(
  glue::glue("{ sratable_prefetch$prefetch} &"),
  file = file.path(
    datadir,
    "00.{gseid}.prefetch.sh" |> glue::glue()
  )
)

sratable_prefetch |>
  dplyr::mutate(
    srafile = file.path(
      srrdir,
      glue::glue("{srrid}.sra")
    )
  ) |>
  dplyr::mutate(
    srafile_exist = file.exists(srafile)
  ) ->
  srafiles

readr::write_lines(
  glue::glue("file {srafiles$srafile}"),
  file = file.path(
    datadir,
    "00.{gseid}.prefetch.check.sh" |> glue::glue()
  )
)

srafiles |>
  dplyr::mutate(
    dump_cmd = purrr::map2_chr(
      .x = srrdir,
      .y = srafile,
      .f = \(.x, .y) {
        .srrid <- basename(.x)
        cmd_dump <- c(
          "fasterq-dump {.y} --temp /scr1/users/liuc9/tmp/fasterq_dump  --include-technical --mem 50G --threads 10 --split-files --outdir {.x}" |> glue::glue()
        )

        cmd <- c(
          cmd_dump
        )
        cmd
      }
    )
  ) ->
  srafile_dump

# save to runfile ---------------------------------------------------------


data.table::fwrite(
  x = srafile_dump,
  file = file.path(
    datadir,
    "{gseid}.runfile.csv" |>glue::glue()
  )
)


# dump sh -----------------------------------------------------------------



readr::write_lines(
  glue::glue("{ srafile_dump$dump_cmd} &"),
  file = file.path(
    datadir,
    "01.{gseid}.dump.sh" |> glue::glue()
  )
)

# dump slrm ---------------------------------------------------------------


dir.create(
  file.path(
    datadir,
    "errout"
  ),
  showWarnings = F,
  recursive = T
)

slrm_header <- c(
  "#!/usr/bin/env bash",
  "# @AUTHOR: Chun-Jie Liu",
  "# @CONTACT: chunjie.sam.liu.at.gmail.com",
  "# @DATE: {lubridate::now()}" |> glue::glue(),
  "",
  "#SBATCH --job-name=01.{gseid}.dump" |> glue::glue(),
  "#SBATCH --output={datadir}/errout/01.{gseid}.dump._%A-%a.out" |> glue::glue(),
  "#SBATCH --error={datadir}/errout/01.{gseid}.dump._%A-%a.err" |> glue::glue(),
  "#SBATCH --cpus-per-task=10",
  "#SBATCH --mem=50G",
  "#SBATCH --array=1-{length(srafile_dump$dump_cmd)}" |> glue::glue(),
  "#SBATCH --time=720:00:00",
  "",
  ""
)

slrm_array <- c(
  "declare -a cmds",
  "while IFS= read -r line; do",
  "  line=$(echo ${line}|sed 's/ *&$//')",
  "  cmds+=(\"${line}\")",
  "done < \"{datadir}/01.{gseid}.dump.sh\"" |> glue::glue(),
  "",
  "",
  "index=$((SLURM_ARRAY_TASK_ID - 1))",
  'cmd="${cmds[$index]}"',
  "",
  "",
  "srun ${cmd}"
)

readr::write_lines(
  c(slrm_header, slrm_array),
  file = file.path(
    datadir,
    "01.{gseid}.dump.slrm" |> glue::glue()
  )
)










# footer ------------------------------------------------------------------

# future::plan(future::sequential)

# save image --------------------------------------------------------------
