#!/usr/bin/env Rscript
# Metainfo ----------------------------------------------------------------

# @AUTHOR: Chun-Jie Liu
# @CONTACT: chunjie.sam.liu.at.gmail.com
# @DATE: Thu Aug 22 12:47:40 2024
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

# s: string, i: integer, f: float, !: boolean
# @: array
# %: hash
# default: default value specified here.
verbose <- FALSE
spec <- "
Usage: Rscript foorbar.R [options]

Options:
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
  logs <- readr::read_lines(
    file = logfile,
    skip_empty_rows = T
  )

  output_index <- stringr::str_which(logs, '"outputs": \\{')
  valid_json <- gsub("^\"|\"$", "", logs[output_index + 1])
  targz <- jsonlite::fromJSON(paste0("{", valid_json, "}"))
  targz$scMOCHABatch.output_dir_tar_gzs
}

# load data ---------------------------------------------------------------

basedir <- "/home/liuc9/github/scMOCHA-data/data"
gseid <- "GSE226602"
datadir <- file.path(
  basedir, gseid
)
# body --------------------------------------------------------------------



# parse log files ---------------------------------------------------------


logfile <- file.path(
  datadir,
  "02.2.{gseid}.batch.log" |> glue::glue()
)
targz <- fn_parse_log(logfile = logfile)

tibble::tibble(
  targz = targz
) |>
  dplyr::mutate(
    srrdir = gsub(".tar.gz", "", x = targz)
  ) |>
  dplyr::mutate(
    srrid = basename(srrdir)
  ) |>
  dplyr::select(srrdir, srrid) ->
srr_out


outdir <- file.path(
  datadir,
  "output"
)
dir.create(outdir, showWarnings = F, recursive = T)



readr::write_lines(
  x = targz,
  file = file.path(
    outdir,
    "{gseid}.scmocha.out.targz.txt" |> glue::glue()
  )
)

srr_out |>
  dplyr::mutate(
    cell_stats = purrr::map(
      .x = srrdir,
      .f = \(.srrdir) {
        .cs <- readxl::read_xlsx(
          file.path(.srrdir, "qc_cell_stats.xlsx")
        )
        .depth <- data.table::fread(
          file.path(.srrdir, "possorted_genome_bam.MT.depth"),
          col.names = c("chrom", "pos", "depth")
        ) |>
          dplyr::select(-chrom)

        .celltype_ratio <- data.table::fread(
          file.path(.srrdir, "celltype_ratio.tsv")
        )
        .cva <- data.table::fread(
          file.path(.srrdir, "cell_variant_annotation.tsv")
        )
        tibble::tibble(
          cell_stats = list(.cs),
          depth = list(.depth),
          celltype_ratio = list(.celltype_ratio),
          anno = list(.cva)
        )
      }
    )
  ) |>
  tidyr::unnest(cols = cell_stats) ->
srr_out_cell_stats


readr::write_rds(
  srr_out_cell_stats,
  file.path(
    outdir,
    "scmocha_suc_out.rds"
  )
)
# footer ------------------------------------------------------------------

# future::plan(future::sequential)

# save image --------------------------------------------------------------
