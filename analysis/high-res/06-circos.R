#!/usr/bin/env Rscript
# Metainfo ----------------------------------------------------------------
# @AUTHOR: Chun-Jie Liu
# @CONTACT: chunjie.sam.liu.at.gmail.com
# @DATE: 2025-12-15 23:33:43
# @DESCRIPTION: this script is used for ...

# Library -----------------------------------------------------------------

load_pkg(jutils)

# args --------------------------------------------------------------------

# s: string, i: integer, f: float, !: boolean, @: array, %: list
GetoptLong.options(help_style = "two-column")
VERSION = "v0.0.1"

# default: default value specified here.

verbose = TRUE

GetoptLong("verbose!", "print messages")


logger::log_threshold(logger::TRACE)
logger::log_layout(logger::layout_glue_colors)

# header ------------------------------------------------------------------

# load data ---------------------------------------------------------------
allvariants <- import(
  "/home/liuc9/github/scMOCHA-data/analysis/high-res/MANUSCRIPTFIGURES/SAMPLE-VARIANT-CLASSIFICATION-CLUSTER-BULK-AF.xlsx"
) |>
  dplyr::mutate(
    coord = parallel::mclapply(
      X = variant,
      FUN = \(.v) {
        # .v <- gse_data_variant_classification_clusteraf_bulkaf$variant[[1]]
        pos <- gsub(
          pattern = "(\\d+)([A-Z])>([A-Z])",
          replacement = "\\1",
          x = .v
        ) |>
          as.integer()
        ref <- gsub(
          pattern = "(\\d+)([A-Z])>([A-Z])",
          replacement = "\\2",
          x = .v
        )
        alt <- gsub(
          pattern = "(\\d+)([A-Z])>([A-Z])",
          replacement = "\\3",
          x = .v
        )
        data.table(
          seqnames = "MT",
          start = pos,
          end = pos,
          ref = ref,
          alt = alt
        )
      },
      mc.cores = 10
    )
  ) |>
  tidyr::unnest(
    cols = coord
  )
# load conn ---------------------------------------------------------------

# src ---------------------------------------------------------------------

# function ----------------------------------------------------------------

source("/home/liuc9/github/scMOCHA-data/analysis/00-colors.R")
fn_plot_circos <- function(
  color_circos_track = c(
    "phastCons100way" = "#FFD700",
    "gnomad" = "#0000FF",
    "homo_paf" = "#03DC62",
    "hete_paf" = "#7A0202",
    "somatic_paf" = "#000000",
    "homo_af" = "#A6FCCB",
    "hete_af" = "#FF0000",
    "somatic_af" = "#040404",
    "gene_name_bg" = "#EAF7FF",
    "coverage" = "#3FB5FF"
  )
) {
  start.degree = 90
  canvas.xlim = c(-1, 1)
  canvas.ylim = c(-1, 1)
  gap.degree = 1
  LENGTH <- 16569

  library(circlize)
  # ! init --------------------------------------------------------------------

  circos.clear()
  circos.par(
    start.degree = start.degree,
    canvas.xlim = canvas.xlim,
    canvas.ylim = canvas.ylim,
    gap.degree = gap.degree
  )

  # ! axis --------------------------------------------------------------------
  coordinate_df <- data.table(
    seqnames = rep("MT", LENGTH),
    start = 1:LENGTH,
    end = 1:LENGTH
  )

  circos.genomicInitialize(
    data = coordinate_df,
    plotType = "axis",
    axis.labels.cex = 0.8 * par("cex"),
    major.by = 1000
  )
  # ! track: phastCons100way --------------------------------------------------------------------
  phastCons100way <- import(
    "/home/liuc9/github/scMOCHA-data/config/chrM.phastCons100way.wigFix.qs"
  )

  circos.genomicTrack(
    phastCons100way,
    track.height = 0.03,
    ylim = c(0, 1),
    track.margin = c(0, 0.005),
    cell.padding = c(0, 0, 0, 0),
    bg.border = NA,
    bg.col = NA,
    panel.fun = function(region, value, ...) {
      circos.genomicLines(
        region = region,
        value = value,
        col = color_circos_track["phastCons100way"],
        lwd = 0.5
      )
    }
  )

  # ! gnomad --------------------------------------------------------------------

  gnomad <- import(
    "/home/liuc9/github/scMOCHA-data/analysis/zzz/db/gnomad.qs"
  ) |>
    dplyr::filter(filters == "PASS") |>
    dplyr::select(position, af = af_hom, ac = ac_hom) |>
    dplyr::arrange(position) |>
    dplyr::mutate(
      seqnames = "MT",
      start = position,
      end = position
    ) |>
    dplyr::select(
      seqnames,
      start,
      end,
      af,
      ac
    ) |>
    dplyr::filter(af > 0.01)

  circos.genomicTrack(
    gnomad,
    track.height = 0.1,
    ylim = c(0, 1),
    track.margin = c(0, 0.01),
    cell.padding = c(0, 0, 0, 0),
    bg.border = NA,
    bg.col = NA,
    panel.fun = function(region, value, ...) {
      pos = region$start
      val = value$af
      circos.barplot(
        value = val,
        pos = pos,
        col = color_circos_track["gnomad"],
        border = color_circos_track["gnomad"]
      )
    }
  )

  # ! homoplasmic paf--------------------------------------------------------------------
  N_INDIVIDUALS <- 571
  allvariants |>
    dplyr::filter(variant_type == "homo") |>
    tidyr::nest(.by = c(seqnames, start, end, ref, alt, variant)) |>
    dplyr::arrange(start) |>
    dplyr::mutate(
      n = purrr::map_int(data, nrow),
      af = purrr::map_dbl(data, \(.df) {
        mean(.df$Bulk, na.rm = TRUE)
      })
    ) |>
    dplyr::select(-data) |>
    dplyr::mutate(paf = n / N_INDIVIDUALS) |>
    dplyr::mutate(af = pmax(0, pmin(1, af, na.rm = TRUE))) -> forplot_homo

  circos.genomicTrack(
    forplot_homo,
    track.height = 0.1,
    ylim = c(0, 1),
    track.margin = c(0, 0.01),
    cell.padding = c(0, 0, 0, 0),
    bg.border = NA,
    bg.col = NA,
    panel.fun = function(region, value, ...) {
      pos = region$start
      val = value$paf
      circos.barplot(
        value = val,
        pos = pos,
        col = color_circos_track["homo_paf"],
        border = color_circos_track["homo_paf"]
      )
    }
  )

  # ! heteroplasmic paf--------------------------------------------------------------------
  N_INDIVIDUALS <- 571
  allvariants |>
    dplyr::filter(variant_type == "hete") |>
    tidyr::nest(.by = c(seqnames, start, end, ref, alt, variant)) |>
    dplyr::arrange(start) |>
    dplyr::mutate(
      n = purrr::map_int(data, nrow),
      af = purrr::map_dbl(data, \(.df) {
        mean(.df$Bulk, na.rm = TRUE)
      })
    ) |>
    dplyr::select(-data) |>
    dplyr::mutate(paf = n / N_INDIVIDUALS) |>
    dplyr::mutate(af = pmax(0, pmin(1, af, na.rm = TRUE))) -> forplot_hete

  circos.genomicTrack(
    forplot_hete,
    track.height = 0.1,
    ylim = c(0, 1),
    track.margin = c(0, 0.01),
    cell.padding = c(0, 0, 0, 0),
    bg.border = NA,
    bg.col = NA,
    panel.fun = function(region, value, ...) {
      pos = region$start
      val = value$paf
      circos.barplot(
        value = val,
        pos = pos,
        col = color_circos_track["hete_paf"],
        border = color_circos_track["hete_paf"]
      )
    }
  )

  # ! somatic paf--------------------------------------------------------------------
  N_INDIVIDUALS <- 571
  allvariants |>
    dplyr::filter(variant_type == "somatic") |>
    tidyr::nest(.by = c(seqnames, start, end, ref, alt, variant)) |>
    dplyr::arrange(start) |>
    dplyr::mutate(
      n = purrr::map_int(data, nrow),
      af = purrr::map_dbl(data, \(.df) {
        mean(.df$Bulk, na.rm = TRUE)
      })
    ) |>
    dplyr::select(-data) |>
    dplyr::mutate(paf = n / N_INDIVIDUALS) |>
    dplyr::mutate(af = pmax(0, pmin(1, af, na.rm = TRUE))) -> forplot_somatic

  circos.genomicTrack(
    forplot_somatic,
    track.height = 0.1,
    ylim = c(0, 1),
    track.margin = c(0, 0.01),
    cell.padding = c(0, 0, 0, 0),
    bg.border = NA,
    bg.col = NA,
    panel.fun = function(region, value, ...) {
      pos = region$start
      val = value$paf
      circos.barplot(
        value = val,
        pos = pos,
        col = color_circos_track["somatic_paf"],
        border = color_circos_track["somatic_paf"]
      )
    }
  )

  # ! gene name--------------------------------------------------------------------
  # ! gene region --------------------------------------------------------------------

  gtf_gene_df <- import(
    "/home/liuc9/github/scMOCHA-data/config/mtdna_genes_dloop.qs"
  )

  circos.genomicTrack(
    gtf_gene_df,
    ylim = c(0, 1),
    track.height = 0.1,
    track.margin = c(0, 0.02),
    cell.padding = c(0, 0, 0, 0),
    bg.border = color_circos_track["gene_name_bg"],
    bg.col = color_circos_track["gene_name_bg"],
    panel.fun = function(region, value, ...) {
      genetypes <- sort(unique(value$TYPE))
      for (genotype in genetypes) {
        r <- region[value$TYPE == genotype, ]
        v <- value[value$TYPE == genotype, ]
        v_gene_name <- gsub(
          pattern = "MT-|-",
          replacement = "",
          x = v$gene_name
        )

        vcol <- value[value$TYPE == genotype, ]$COLOR
        circos.genomicRect(
          region = r,
          value = v,
          col = vcol,
          border = "white",
          lty = 1
        )
        circos.genomicText(
          region = r,
          value = v,
          y = 0,
          adj = c(0, 0.5),
          labels = v_gene_name,
          facing = "clockwise",
          niceFacing = TRUE,
          cex = 0.6,
          col = "black",
        )
      }
    },
  )

  # ! homoplasmic af--------------------------------------------------------------------

  circos.genomicTrack(
    forplot_homo,
    track.height = 0.1,
    ylim = c(0, 1),
    track.margin = c(0, 0.01),
    cell.padding = c(0, 0, 0, 0),
    bg.border = NA,
    bg.col = NA,
    panel.fun = function(region, value, ...) {
      pos = region$start
      val = value$af
      circos.genomicPoints(
        region = region,
        value = value,
        pch = 2,
        cex = 0.5,
        col = color_circos_track["homo_af"],
        bg.border = color_circos_track["homo_af"]
      )
    }
  )

  # ! heteroplasmic af--------------------------------------------------------------------

  circos.genomicTrack(
    forplot_hete,
    track.height = 0.1,
    ylim = c(0, 1),
    track.margin = c(0, 0.01),
    cell.padding = c(0, 0, 0, 0),
    bg.border = NA,
    bg.col = NA,
    panel.fun = function(region, value, ...) {
      pos = region$start
      val = value$af
      # circos.barplot(
      #   value = val,
      #   pos = pos,
      #   col = "red",
      #   border = "red"
      # )
      circos.genomicPoints(
        region = region,
        value = value,
        pch = 3,
        cex = 0.5,
        col = color_circos_track["hete_af"],
        bg.border = color_circos_track["hete_af"]
      )
    }
  )

  # ! somatic af--------------------------------------------------------------------

  circos.genomicTrack(
    forplot_somatic,
    track.height = 0.1,
    ylim = c(0, 1),
    track.margin = c(0, 0.01),
    cell.padding = c(0, 0, 0, 0),
    bg.border = NA,
    bg.col = NA,
    panel.fun = function(region, value, ...) {
      pos = region$start
      val = value$af
      # circos.barplot(
      #   value = val,
      #   pos = pos,
      #   col = "red",
      #   border = "red"
      # )
      circos.genomicPoints(
        region = region,
        value = value,
        pch = 3,
        cex = 0.5,
        col = color_circos_track["somatic_af"],
        bg.border = color_circos_track["somatic_af"]
      )
    }
  )

  circos.clear()
}

# body --------------------------------------------------------------------

# footer ------------------------------------------------------------------

# save image --------------------------------------------------------------
