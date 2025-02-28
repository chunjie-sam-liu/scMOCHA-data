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

# function ----------------------------------------------------------------


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

erroutdir <- file.path(
  datadir,
  "errout"
)
dir.create(erroutdir, showWarnings = F, recursive = T)

# body --------------------------------------------------------------------

sratable <- data.table::fread(
  file.path(
    datadir,
    "{gseid}.SraRunTable" |> glue::glue()
  ),
  sep = ","
)

# GSM samples -------------------------------------------------------------

sample_name_ <- ifelse("Sample Name" %in% colnames(sratable), "Sample Name", "SampleName")

sratable |>
  dplyr::select(
    run_accession = Run,
    experiment_name = sample_name_,
    experiment_accession = Experiment
  ) |>
  data.table::fwrite(
    file = file.path(
      datadir,
      "{gseid}.SraRunTable.GSM" |> glue::glue()
    )
  )

# data.table::fread(
#   file.path(
#     datadir,
#     "Metadata_for_{gseid}.txt" |> glue::glue()
#   )
# ) |>
#   dplyr::select(run_accession = srr_ids, experiment_name = gsm_id, experiment_accession = srx_id) |>
#   data.table::fwrite(
#     file = file.path(
#       datadir,
#       "{gseid}.SraRunTable.GSM" |> glue::glue()
#     )
#   )

# Prefetch ----------------------------------------------------------------
sradir <- file.path(
  datadir,
  "sra"
)
dir.create(sradir, showWarnings = F, recursive = T)

sratable |>
  dplyr::select(srrid = Run) |>
  dplyr::mutate(
    srrdir = file.path(
      sradir, srrid
    )
  ) |>
  dplyr::mutate(
    srrdir_exists = file.exists(srrdir)
  ) |>
  dplyr::mutate(
    prefetch = "prefetch -p --max-size 100G {srrid} --output-directory {sradir} 1>{erroutdir}/prefetch.{srrid}.err 2>{erroutdir}/prefetch.{srrid}.err " |> glue::glue()
  ) ->
sratable_prefetch



readr::write_lines(
  glue::glue("{ sratable_prefetch$prefetch} &"),
  file = file.path(
    datadir,
    "00.{gseid}.prefetch.sh" |> glue::glue()
  )
)


# Check sra file ----------------------------------------------------------


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
  glue::glue("vdb-validate {srafiles$srafile} 1> {srafiles$srafile}.validate 2>&1 &"),
  file = file.path(
    datadir,
    "01.{gseid}.prefetch.check.sh" |> glue::glue()
  )
)

# Generate dump scripts ---------------------------------------------------


srafiles |>
  dplyr::mutate(
    dump_cmd = purrr::map2_chr(
      .x = srrdir,
      .y = srafile,
      .f = \(.x, .y) {
        .srrid <- basename(.x)

        cmd_dump <- c(
          "fasterq-dump {.y} --temp /mnt/isilon/u01_project/large-scale/liuc9/tmp  --include-technical --mem 50G --threads 10 --split-files --outdir {.x} 1>{erroutdir}/fasterq_dump.{.srrid}.err 2>{erroutdir}/fasterq_dump.{.srrid}.err" |> glue::glue()
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
    "{gseid}.runfile.csv" |> glue::glue()
  )
)


# dump sh -----------------------------------------------------------------



readr::write_lines(
  glue::glue("{ srafile_dump$dump_cmd} &"),
  file = file.path(
    datadir,
    "02.{gseid}.dump.sh" |> glue::glue()
  )
)

# dump slrm ---------------------------------------------------------------

slrm_header <- c(
  "#!/usr/bin/env bash",
  "# @AUTHOR: Chun-Jie Liu",
  "# @CONTACT: chunjie.sam.liu.at.gmail.com",
  "# @DATE: {lubridate::now()}" |> glue::glue(),
  "",
  "#SBATCH --job-name=02.{gseid}.dump" |> glue::glue(),
  "#SBATCH --output={datadir}/errout/02.{gseid}.dump._%A-%a.err" |> glue::glue(),
  "#SBATCH --error={datadir}/errout/02.{gseid}.dump._%A-%a.err" |> glue::glue(),
  "#SBATCH --cpus-per-task=10",
  "#SBATCH --mem=80G",
  "#SBATCH --array=1-{length(srafile_dump$dump_cmd)}" |> glue::glue(),
  "#SBATCH --time=720:00:00",
  "",
  ""
)

slrm_array <- c(
  "declare -a cmds",
  "while IFS= read -r line; do",
  "  line=${line%%1>*}",
  "  cmds+=(\"${line}\")",
  "done < \"{datadir}/02.{gseid}.dump.sh\"" |> glue::glue(),
  "",
  "",
  "index=$((SLURM_ARRAY_TASK_ID - 1))",
  'cmd="${cmds[$index]}"',
  'echo "Executing command ${cmd}"',
  "",
  "",
  "srun ${cmd}"
)

readr::write_lines(
  c(slrm_header, slrm_array),
  file = file.path(
    datadir,
    "02.{gseid}.dump.slrm" |> glue::glue()
  )
)






# footer ------------------------------------------------------------------

# future::plan(future::sequential)

# save image --------------------------------------------------------------
