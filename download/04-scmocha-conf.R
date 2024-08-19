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
          "-i {.jsonfile} 1>{.logfile} 2>{.errfile} &" |> glue::glue()
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

# footer ------------------------------------------------------------------

# future::plan(future::sequential)

# save image --------------------------------------------------------------
