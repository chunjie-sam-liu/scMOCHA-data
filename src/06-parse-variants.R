#!/usr/bin/env Rscript --vanilla
# Metainfo ----------------------------------------------------------------

# @AUTHOR: Chun-Jie Liu
# @CONTACT: chunjie.sam.liu.at.gmail.com
# @DATE: Fri Sep 27 10:41:12 2024
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
# gseid <- "GSE149689"

spec <- "
Usage: Rscript foorbar.R [options]
Options:

<gseid=s> gseid, required
<verbose!> Print messages
"

GetoptLong.options(help_style = "two-column")
GetoptLong(spec, template_control = list(opt_width = 21))

# src ---------------------------------------------------------------------


# header ------------------------------------------------------------------
log_threshold(TRACE)
log_layout(layout_glue_colors)

# future::plan(future::multisession, workers = 10)

# function ----------------------------------------------------------------
fn_somatic_variant <- function(.haplo_variant, .haplo_violin, .n_cells = 10) {
  # .haplo_variant <- srr_out_cell_stats$haplo_variant[[1]]
  # .haplo_violin <- srr_out_cell_stats$haplo_violin[[1]]

  # 1. filter by haplogrep marker variant
  .haplo_variant |>
    dplyr::filter(fill != "white") |>
    dplyr::pull(variant) ->
  .v_haplo

  # 2. filter by n_cells
  .n_cells <- 10
  .haplo_violin |>
    dplyr::count(variant) |>
    dplyr::filter(n < .n_cells) |>
    dplyr::pull(variant) ->
  .v_n_cells

  # 3. tRNA p9 and RNA editing position
  .editing_pos <- c(
    585, 1610, 3238, 4271, 5520, 7526, 8303, # tRNA p9
    9999, 10413, 12146, 12274, 14734, 15896, # tRNA p9
    295, 2617, 13710 # RNA editing
  )

  .haplo_variant |>
    dplyr::filter(Position %in% .editing_pos) |>
    dplyr::pull(variant) ->
  .v_editing

  .haplo_variant |>
    dplyr::filter(!variant %in% c(.v_haplo, .v_n_cells, .v_editing)) |>
    dplyr::pull(variant) ->
  .v_somatic

  list(
    haplo = .v_haplo,
    n_cells = .v_n_cells,
    editing = .v_editing,
    somatic = .v_somatic
  )
}

# load data ---------------------------------------------------------------

basedir <- "/mnt/isilon/u01_project/large-scale/liuc9/raw"
datadir <- file.path(
  basedir, gseid
)

targzdir <- file.path(
  datadir, "targz"
)

outdir <- file.path(
  datadir, "out"
)
dir.create(outdir, showWarnings = F, recursive = T)

srrid_list <- readr::read_lines(
  file = file.path(
    datadir,
    "{gseid}.srrid.list" |> glue::glue()
  )
)
# body --------------------------------------------------------------------

tibble::tibble(
  srrid = srrid_list
) |>
  dplyr::mutate(
    srrdir = file.path(targzdir, srrid)
  ) |>
  dplyr::mutate(
    dir_exists = dir.exists(srrdir)
  ) ->
srr_out

srr_out |>
  dplyr::mutate(
    cell_stats = parallel::mclapply(
      X = srrdir,
      FUN = purrr::safely(\(.srrdir) {
        if (!dir.exists(.srrdir)) {
          return(NULL)
        }
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

        .cva |>
          dplyr::mutate(
            v = glue::glue("{Position}{Ref}>{Alt}")
          ) |>
          dplyr::pull(v) ->
        .v

        .cva$Position -> .pos

        .hetero <- data.table::fread(
          file.path(.srrdir, "cluster.cell_heteroplasmic_df.tsv.gz")
        ) |>
          dplyr::rename(celltype = V1) |>
          tidyr::gather(-celltype, key = variant, value = af) |>
          dplyr::filter(variant %in% .v)

        .cov <- data.table::fread(
          file.path(.srrdir, "cluster.coverage.txt.gz"),
          col.names = c("pos", "celltype", "count")
        ) |>
          dplyr::filter(pos %in% .pos)

        .haplo_variant <- data.table::fread(
          file.path(.srrdir, "cluster_cell_violin_haplo_variant.csv")
        )

        .haplo_violin <- data.table::fread(
          file.path(.srrdir, "cluster_cell_violin_haplo_forplot.csv")
        )

        .somatic <- readr::read_rds(
          file.path(.srrdir, "somatic_variant.rds")
        )

        tibble::tibble(
          cell_stats = list(.cs),
          depth = list(.depth),
          celltype_ratio = list(.celltype_ratio),
          anno = list(.cva),
          hetero = list(.hetero),
          coverage = list(.cov),
          haplo_variant = list(.haplo_variant),
          haplo_violin = list(.haplo_violin),
          somatic_variant = list(.somatic)
        )
      }),
      mc.cores = 20
    )
  ) |>
  dplyr::mutate(
    cell_stats = purrr::map(cell_stats, "result")
  ) |>
  tidyr::unnest(cols = cell_stats) ->
srr_out_cell_stats

log_success("{gseid} save to {outdir}/{gseid}.scmocha.out.rds" |> glue::glue())
readr::write_rds(
  srr_out_cell_stats,
  file.path(
    outdir,
    "{gseid}.scmocha.out.rds" |> glue::glue()
  )
)

# variants ----------------------------------------------------------------

# srr_out_cell_stats -> variant
srr_out_cell_stats |>
  dplyr::mutate(
    nmut = purrr::map(
      .x = haplo_variant,
      .f = \(.x) {
        if (is.null(.x)) {
          return(NA_integer_)
        }
        tibble::tibble(
          nmut = nrow(.x),
          nmut_nohaplo = nrow(.x |> dplyr::filter(color != "white"))
        )
      }
    )
  ) |>
  tidyr::unnest(cols = nmut) |>
  dplyr::mutate(
    haplogroup = purrr::map2(
      .x = anno,
      .y = srrid,
      .f = \(.x, .y) {
        log_info(gseid, " ", .y)
        if (is.null(.x)) {
          return(
            tibble::tibble(
              Haplogroup = NA_character_,
              Verbose_haplogroup = NA_character_
            )
          )
        }
        .x |>
          dplyr::select(Haplogroup, Verbose_haplogroup) |>
          dplyr::filter(!is.na(Haplogroup)) |>
          dplyr::filter(Haplogroup != "") |>
          dplyr::distinct() |>
          dplyr::mutate_all(.funs = as.character) ->
        .xx

        if (nrow(.xx) == 0) {
          tibble::tibble(
            Haplogroup = NA_character_,
            Verbose_haplogroup = NA_character_
          )
        } else {
          .xx
        }
      }
    )
  ) |>
  tidyr::unnest(cols = haplogroup) |>
  # dplyr::inner_join(
  #   pheno, by = "srrid"
  # ) |>
  tidyr::unnest(
    cols = cell_stats
  ) |>
  dplyr::mutate(
    ratio = round(`number of cells after filtering` / `estimated number of cells`, 2)
  ) ->
metadata_anno



metadata_anno |>
  dplyr::select(
    srrid,
    `Median UMI/cell` = `median UMI counts per cell`,
    `Median genes/cell` = `median genes per cell`,
    `# of cells` = `estimated number of cells`,
    `# cells after filter` = `number of cells after filtering`,
    `Cell ratio` = ratio,
    `# of variants` = nmut,
    Haplogroup = Haplogroup,
    `# of somatic variants` = nmut_nohaplo
  ) ->
metadata_clean

metadata_clean |>
  writexl::write_xlsx(
    path = file.path(
      outdir,
      "{gseid}.cell_ratio_and_variant_clean.xlsx" |> glue::glue()
    )
  )
log_success("save metadata to {outdir}/{gseid}.cell_ratio_and_variant_clean.xlsx")

data.table::fwrite(
  x = metadata_clean,
  file = file.path(
    outdir,
    "{gseid}.cell_ratio_and_variant_clean.csv" |> glue::glue()
  )
)

log_success("save metadata to {outdir}/{gseid}.cell_ratio_and_variant_clean.csv")
# plot cell estimates -----------------------------------------------------



# footer ------------------------------------------------------------------

# future::plan(future::sequential)

# save image --------------------------------------------------------------
