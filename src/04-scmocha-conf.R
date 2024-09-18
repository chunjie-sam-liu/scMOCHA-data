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
# library(rlang)
library(GetoptLong)
library(logger)

# args --------------------------------------------------------------------

# gseid <- "GSE163668"


# s: string, i: integer, f: float, !: boolean
# @: array
# %: hash
# default: default value specified here.



verbose <- FALSE

# cell ranger chemistry auto detector
# Assay configuration. NOTE: by default the assay configuration is detected
# automatically, which is the recommended mode. You usually will not need to specify a
# chemistry. Options are: 'auto' for autodetection, 'threeprime' for Single Cell 3',
#           'fiveprime' for  Single Cell 5', 'SC3Pv1' or 'SC3Pv2' or 'SC3Pv3' or 'SC3Pv4' for
# Single Cell 3' v1/v2/v3/v4, 'SC3Pv3LT' for Single Cell 3' v3 LT, 'SC3Pv3HT' for
# Single Cell 3' v3 HT, 'SC5P-PE' or 'SC5P-PE-v3' or 'SC5P-R2' or 'SC5P-R2-v3', for
#           Single Cell 5', paired-end/R2-only, 'SC-FB' for Single Cell Antibody-only 3' v2 or
#           5'. To analyze the GEX portion of multiome data, chemistry must be set to 'ARC-v1'
# [default: auto]
chemistry <- "auto"

spec <- "
Usage: Rscript foorbar.R [options]

Options:
<gseid=s> gseid
<chemistry=s> Assay configuration. NOTE: by default the assay configuration is detected
          automatically, which is the recommended mode. You usually will not need to specify a
          chemistry. Options are: 'auto' for autodetection, 'threeprime' for Single Cell 3',
          'fiveprime' for  Single Cell 5', 'SC3Pv1' or 'SC3Pv2' or 'SC3Pv3' or 'SC3Pv4' for
          Single Cell 3' v1/v2/v3/v4, 'SC3Pv3LT' for Single Cell 3' v3 LT, 'SC3Pv3HT' for
          Single Cell 3' v3 HT, 'SC5P-PE' or 'SC5P-PE-v3' or 'SC5P-R2' or 'SC5P-R2-v3', for
          Single Cell 5', paired-end/R2-only, 'SC-FB' for Single Cell Antibody-only 3' v2 or
          5'. To analyze the GEX portion of multiome data, chemistry must be set to 'ARC-v1'
          [default: auto]
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


# basedir <- "/mnt/isilon/u01_project/large-scale/liuc9/scMOCHA-data/data"
basedir <- "/mnt/isilon/u01_project/large-scale/liuc9/raw"
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
          "scMOCHA.chemistry" = "{chemistry}" |> glue::glue(),
          "scMOCHA.transcriptome" = "/mnt/isilon/u01_project/large-scale/liuc9/refdata/mgatk_index/Human",
          "scMOCHA.rCRS" = "/mnt/isilon/u01_project/large-scale/liuc9/scMOCHA/fasta/rCRS.MT.fasta",
          "scMOCHA.mt_exons_df" = "/mnt/isilon/u01_project/large-scale/liuc9/scMOCHA/fasta/mt_exons.df.rds.gz",
          "scMOCHA.mt_features_gmoviz" = "/mnt/isilon/u01_project/large-scale/liuc9/scMOCHA/fasta/mt_features.grange.gmoviz.rds.gz",
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
          "scMOCHA.perlscript" = "/mnt/isilon/u01_project/large-scale/liuc9/scMOCHA/bin/get_variants_info.pl",
          "scMOCHA.jar_path" = "/mnt/isilon/u01_project/large-scale/liuc9/tools/haplogrep3",
          "scMOCHA.sqlite_path" = "/mnt/isilon/u01_project/large-scale/liuc9/refdata/mitomaster/mitomap_sqlite_20230525.sqlite3",
          "scMOCHA.nFeature_RNA_min" = 500,
          "scMOCHA.nFeature_RNA_max" = 8000
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
          "nohup java -Dconfig.file=/mnt/isilon/u01_project/large-scale/liuc9/scMOCHA/config/slurm.conf \\",
          "-jar /mnt/isilon/u01_project/large-scale/liuc9/tools/cromwell-78.jar \\",
          "run /mnt/isilon/u01_project/large-scale/liuc9/scMOCHA/scMOCHA.wdl \\",
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
    "03.{gseid}.runwdl.sh_notrun" |> glue::glue()
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
  "#SBATCH --output={datadir}/errout/02.{gseid}.runwdl._%A-%a.err" |> glue::glue(),
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
  "input_files=($(sed 's/bash \\(.*\\) &/\\1/' {datadir}/03.{gseid}.runwdl.sh_notrun))" |> glue::glue(),
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
    "03.{gseid}.runwdl.slrm_notrun" |> glue::glue()
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
  "scMOCHABatch.chemistry" = "{chemistry}" |> glue::glue(),
  "scMOCHABatch.transcriptome" = "/mnt/isilon/u01_project/large-scale/liuc9/refdata/mgatk_index/Human",
  "scMOCHABatch.rCRS" = "/mnt/isilon/u01_project/large-scale/liuc9/scMOCHA/fasta/rCRS.MT.fasta",
  "scMOCHABatch.mt_exons_df" = "/mnt/isilon/u01_project/large-scale/liuc9/scMOCHA/fasta/mt_exons.df.rds.gz",
  "scMOCHABatch.mt_features_gmoviz" = "/mnt/isilon/u01_project/large-scale/liuc9/scMOCHA/fasta/mt_features.grange.gmoviz.rds.gz",
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
  "scMOCHABatch.perlscript" = "/mnt/isilon/u01_project/large-scale/liuc9/scMOCHA/bin/get_variants_info.pl",
  "scMOCHABatch.jar_path" = "/mnt/isilon/u01_project/large-scale/liuc9/tools/haplogrep3",
  "scMOCHABatch.sqlite_path" = "/mnt/isilon/u01_project/large-scale/liuc9/refdata/mitomaster/mitomap_sqlite_20230525.sqlite3",
  "scMOCHABatch.nFeature_RNA_min" = 500,
  "scMOCHABatch.nFeature_RNA_max" = 8000
)


inputs_json_file <- file.path(
  datadir,
  "04.{gseid}.batch.json" |> glue::glue()
)

errfile <- file.path(
  datadir,
  "04.{gseid}.batch.err" |> glue::glue()
)
logfile <- file.path(
  datadir,
  "04.{gseid}.batch.log" |> glue::glue()
)

runwdl_batch_cmd <- c(
  "#!/usr/bin/env bash",
  "# @AUTHOR: Chun-Jie Liu",
  "# @CONTACT: chunjie.sam.liu.at.gmail.com",
  "# @DATE: {lubridate::now()}" |> glue::glue(),
  "",
  "module load Java/15.0.1",
  "nohup java -Dconfig.file=/mnt/isilon/u01_project/large-scale/liuc9/scMOCHA/config/slurm.conf \\",
  "-jar /mnt/isilon/u01_project/large-scale/liuc9/tools/cromwell-78.jar \\",
  "run /mnt/isilon/u01_project/large-scale/liuc9/scMOCHA-data/scMOCHA.batch.wdl \\",
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
    "04.{gseid}.batch.sh" |> glue::glue()
  )
)


# 249 * 0.7# footer ------------------------------------------------------------------

# future::plan(future::sequential)

# save image --------------------------------------------------------------
