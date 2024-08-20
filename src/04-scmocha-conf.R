#!/usr/bin/env Rscript
# Metainfo ----------------------------------------------------------------

# @AUTHOR: Chun-Jie Liu
# @CONTACT: chunjie.sam.liu.at.gmail.com
# @DATE: Sun Aug 18 23:50:47 2024
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
  dplyr::select(srrid = experiment_name) |>
  dplyr::distinct() |>
  dplyr::mutate(
    scmocha = purrr::map(
      .x = srrid,
      .f = \(.x) {
        .ydir <- file.path(
          datadir,
          .x
        )
        dir.create(.ydir, showWarnings = F, recursive = T)
        .srrid <- .x

        conf <- list(
          "scMOCHA.output_id" = "{.srrid}" |> glue::glue(),
          "scMOCHA.fastqs" = "{.ydir}" |> glue::glue(),
          "scMOCHA.sample_id" = "{.srrid}" |> glue::glue(),
          "scMOCHA.transcriptome" = "/home/liuc9/data/refdata/mgatk_index/Human",
          "scMOCHA.rCRS" = "/home/liuc9/github/scMOCHA/fasta/rCRS.MT.fasta",
          "scMOCHA.mt_exons_df" = "/home/liuc9/github/scMOCHA/fasta/mt_exons.df.rds.gz",
          "scMOCHA.mt_features_gmoviz" = "/home/liuc9/github/scMOCHA/fasta/mt_features.grange.gmoviz.rds.gz",
          "scMOCHA.output_dir" = "{.srrid}" |> glue::glue(),
          "scMOCHA.chrM" = "MT",
          "scMOCHA.low_coverage_threshold" = 10,
          "scMOCHA.npcs" = 10,
          "scMOCHA.reso" = 0.1,
          "scMOCHA.cellrefname" = "pbmcref",
          "scMOCHA.celllevel" = "celltype.l1",
          "scMOCHA.memory" = "50 GB",
          "scMOCHA.boot_disk_size_gb" = "12",
          "scMOCHA.disk_space" = "50",
          "scMOCHA.cpu" = "10",
          "scMOCHA.scmocha_version" = "latest",
          "scMOCHA.docker" = "chunjiesamliu/scmocha",
          "scMOCHA.partition" = "defq",
          "scMOCHA.account" = "liuc9",
          "scMOCHA.IMAGE" = "/scr1/users/liuc9/sif/scmocha_latest.sif",
          "scMOCHA.perlscript" = "/home/liuc9/github/scMOCHA/bin/get_variants_info.pl",
          "scMOCHA.jar_path" = "/scr1/users/liuc9/tools/haplogrep3",
          "scMOCHA.sqlite_path" = "/mnt/isilon/xing_lab/liuc9/refdata/mitomaster/mitomap_sqlite_20230525.sqlite3",
          "scMOCHA.nFeature_RNA_min" =  500,
          "scMOCHA.nFeature_RNA_max" =  8000,
          "scMOCHA.x10_version" = "v3"
        )

        conf_file <- file.path(
          .ydir,
          "{.srrid}.json" |> glue::glue()
        )

        jsonlite::write_json(
          x = conf,
          path = conf_file,
          auto_unbox = TRUE
        )

        .jsonfile <- file.path(
          .ydir, "{.srrid}.json" |> glue::glue()
        )
        .errfile <- file.path(
          .ydir, "{.srrid}.err" |> glue::glue()
        )
        .logfile <- file.path(
          .ydir, "{.srrid}.log" |> glue::glue()
        )
        runwdl_sh_file <- file.path(
          .ydir, "runwdl_{.srrid}.sh" |> glue::glue()
        )

        runwdl_cmd <- c(
          "#!/usr/bin/env bash",
          "# @AUTHOR: Chun-Jie Liu",
          "# @CONTACT: chunjie.sam.liu.at.gmail.com",
          "# @DATE: {lubridate::now()}" |> glue::glue(),
          "",
          "module load Java/15.0.1",
          "nohup java -Dconfig.file=/home/liuc9/github/scMOCHA/config/slurm.conf \\",
          "-jar /home/liuc9/tools/cromwell-78.jar \\",
          "run /home/liuc9/github/scMOCHA/scMOCHA.wdl \\",
          "-i {.jsonfile} 1>{.logfile} 2>{.errfile} " |> glue::glue()
        )

        readr::write_lines(
          x = runwdl_cmd,
          file = runwdl_sh_file
        )

        tibble::tibble(
          # srrid = .x,
          srrdir = .ydir,
          scmocha_sh = runwdl_sh_file
        )

      }
    )
  ) |>
  tidyr::unnest(cols = scmocha) ->
  conf_scmocha

data.table::fwrite(
  x = conf_scmocha,
  file = file.path(
    datadir,
    "{gseid}.scmocha.csv" |> glue::glue()
  )
)

readr::write_lines(
  glue::glue("bash {conf_scmocha$scmocha_sh} &"),
  file = file.path(
    datadir,
    "02.{gseid}.runwdl.sh" |> glue::glue()
  )
)

# create job array --------------------------------------------------------

dir.create(
  path = file.path(
    datadir, "errout"
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
  "#SBATCH --job-name=02.{gseid}.runwdl" |> glue::glue(),
  "#SBATCH --output={datadir}/errout/02.{gseid}.runwdl.vscode-terminal:/565d84ffd6ee234e22ec6d0b37934647/7_%A-%a.out" |> glue::glue(),
  "#SBATCH --error={datadir}/errout/02.{gseid}.runwdl._%A-%a.err" |> glue::glue(),
  "#SBATCH --cpus-per-task=10",
  "#SBATCH --mem=50G",
  "#SBATCH --array=1-{length(conf_scmocha$scmocha_sh)}" |> glue::glue(),
  "#SBATCH --time=720:00:00",
  "",
  ""
)

slrm_array <- c(
  "#input_files=({paste0(conf_scmocha$scmocha_sh, collapse = ' ')})" |> glue::glue(),
  "input_files=($(sed 's/bash \\(.*\\) &/\\1/' {datadir}/02.{gseid}.runwdl.sh))" |> glue::glue(),
  "",
  "",
  "index=$((SLURM_ARRAY_TASK_ID - 1))",
  'file="${input_files[$index]}"',
  "",
  "",
  'srun bash "${file}"'
)

readr::write_lines(
  c(slrm_header, slrm_array),
  file = file.path(
    datadir,
    "02.1.{gseid}.runwdl.slrm" |> glue::glue()
  )
)



# for batch scMOCHA solving database conflicts ----------------------------

# input srr id
readr::write_lines(
  x = conf_scmocha$srrid,
  file = file.path(
    datadir,
    "{gseid}.srrid.list" |> glue::glue()
  )
)

# input srr dir
readr::write_lines(
  x = conf_scmocha$srrdir,
  file = file.path(
    datadir,
    "{gseid}.srrdir.list" |> glue::glue()
  )
)


inputs_json <- list(
  "scMOCHABatch.output_id_list" = file.path(
    datadir,
    "{gseid}.srrid.list" |> glue::glue()
  ),
  "scMOCHABatch.fastqs_list" = file.path(
    datadir,
    "{gseid}.srrdir.list" |> glue::glue()
  ),
  "scMOCHABatch.sample_id_list" = file.path(
    datadir,
    "{gseid}.srrid.list" |> glue::glue()
  ),

  "scMOCHABatch.transcriptome" = "/home/liuc9/data/refdata/mgatk_index/Human",
  "scMOCHABatch.rCRS" = "/home/liuc9/github/scMOCHA/fasta/rCRS.MT.fasta",
  "scMOCHABatch.mt_exons_df" = "/home/liuc9/github/scMOCHA/fasta/mt_exons.df.rds.gz",
  "scMOCHABatch.mt_features_gmoviz" = "/home/liuc9/github/scMOCHA/fasta/mt_features.grange.gmoviz.rds.gz",

  "scMOCHABatch.output_dir_list" = file.path(
    datadir,
    "{gseid}.srrid.list" |> glue::glue()
  ),

  "scMOCHABatch.chrM" = "MT",
  "scMOCHABatch.low_coverage_threshold" = 10,
  "scMOCHABatch.npcs" = 10,
  "scMOCHABatch.reso" = 0.1,
  "scMOCHABatch.cellrefname" = "pbmcref",
  "scMOCHABatch.celllevel" = "celltype.l1",
  "scMOCHABatch.memory" = "50GB",
  "scMOCHABatch.boot_disk_size_gb" = "12",
  "scMOCHABatch.disk_space" = "50",
  "scMOCHABatch.cpu" = "10",
  "scMOCHABatch.scmocha_version" = "latest",
  "scMOCHABatch.docker" = "chunjiesamliu/scmocha",
  "scMOCHABatch.partition" = "defq",
  "scMOCHABatch.account" = "liuc9",
  "scMOCHABatch.IMAGE" = "/scr1/users/liuc9/sif/scmocha_latest.sif",
  "scMOCHABatch.perlscript" = "/home/liuc9/github/scMOCHA/bin/get_variants_info.pl",
  "scMOCHABatch.jar_path" = "/scr1/users/liuc9/tools/haplogrep3",
  "scMOCHABatch.sqlite_path" = "/mnt/isilon/xing_lab/liuc9/refdata/mitomaster/mitomap_sqlite_20230525.sqlite3",
  "scMOCHABatch.nFeature_RNA_min" =  500,
  "scMOCHABatch.nFeature_RNA_max" =  8000,
  "scMOCHABatch.x10_version" = "v3"
)


inputs_json_file <- file.path(
  datadir,
  "02.2.{gseid}.batch.json" |> glue::glue()
)

errfile <- file.path(
  datadir,
  "02.2.{gseid}.batch.err" |> glue::glue()
)
logfile <- file.path(
  datadir,
  "02.2.{gseid}.batch.log" |> glue::glue()
)

runwdl_batch_cmd <- c(
  "#!/usr/bin/env bash",
  "# @AUTHOR: Chun-Jie Liu",
  "# @CONTACT: chunjie.sam.liu.at.gmail.com",
  "# @DATE: {lubridate::now()}" |> glue::glue(),
  "",
  "module load Java/15.0.1",
  "nohup java -Dconfig.file=/home/liuc9/github/scMOCHA/config/slurm.conf \\",
  "-jar /home/liuc9/tools/cromwell-78.jar \\",
  "run /home/liuc9/github/scMOCHA-data/scMOCHA.batch.wdl \\",
  "-i {inputs_json_file} 1>{logfile} 2>{errfile} &" |> glue::glue()
)

jsonlite::write_json(
  x = inputs_json,
  path = inputs_json_file,
  auto_unbox = TRUE
)


readr::write_lines(
  x = runwdl_batch_cmd,
  file = file.path(
    datadir,
    "02.2.{gseid}.batch.sh" |> glue::glue()
  )
)


# 249 * 0.7# footer ------------------------------------------------------------------

# future::plan(future::sequential)

# save image --------------------------------------------------------------
