#!/usr/bin/env Rscript
# Metainfo ----------------------------------------------------------------
# @AUTHOR: Chun-Jie Liu
# @CONTACT: chunjie.sam.liu.at.gmail.com
# @DATE: 2025-12-15 23:33:43
# @DESCRIPTION: this script is used for ...

# Library -----------------------------------------------------------------

load_pkg(jutils)

dotenv(".env")

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
load_pkg(circlize)

conflicted::conflict_prefer("filter", "dplyr")

allvariants <- import(
  path(
    Sys.getenv("OUTDIR"),
    "SAMPLE-VARIANT-CLASSIFICATION-CLUSTER-BULK-AF.xlsx"
  )
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

allvariants |>
  dplyr::select(variant_type, variant) |>
  dplyr::distinct() |>
  dplyr::count(variant_type)
# load conn ---------------------------------------------------------------

# src ---------------------------------------------------------------------

# function ----------------------------------------------------------------

source(path(
  Sys.getenv("HIGHRESDIR"),
  "00-colors.R"
))
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
    path(Sys.getenv("REPODIR"), "config/chrM.phastCons100way.wigFix.qs")
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
    path(Sys.getenv("ZZZDIR"), "db/gnomad.qs")
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

  # circos.genomicTrack(
  #   forplot_somatic,
  #   track.height = 0.1,
  #   ylim = c(0, 1),
  #   track.margin = c(0, 0.01),
  #   cell.padding = c(0, 0, 0, 0),
  #   bg.border = NA,
  #   bg.col = NA,
  #   panel.fun = function(region, value, ...) {
  #     pos = region$start
  #     val = value$paf
  #     circos.barplot(
  #       value = val,
  #       pos = pos,
  #       col = color_circos_track["somatic_paf"],
  #       border = color_circos_track["somatic_paf"]
  #     )
  #   }
  # )

  # ! gene name--------------------------------------------------------------------
  # ! gene region --------------------------------------------------------------------

  gtf_gene_df <- import(
    path(Sys.getenv("REPODIR"), "config/mtdna_genes_dloop.qs")
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
          cex = 1,
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
    track.margin = c(0, 0.02),
    cell.padding = c(0, 0, 0, 0),
    bg.border = NA,
    bg.col = NA,
    panel.fun = function(region, value, ...) {
      pos = region$start
      val = value$af
      circos.genomicPoints(
        region = region,
        value = val,
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
    bg.border = color_circos_track["hete_af"],
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
        value = val,
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
        value = val,
        pch = 3,
        cex = 0.5,
        col = color_circos_track["somatic_af"],
        bg.border = color_circos_track["somatic_af"]
      )
    }
  )

  circos.clear()

  # ! legend ------------------------------------------------------------------
  legend(
    x = "bottomleft",
    inset = c(0.01, 0.01),
    legend = c(
      "phastCons100way conservation",
      "gnomAD populational AF",
      "Homoplasmic populational AF",
      "Heteroplasmic populational AF",
      "Homoplasmic mean AF",
      "Heteroplasmic mean AF",
      "Somatic mean AF"
    ),
    col = c(
      color_circos_track["phastCons100way"],
      color_circos_track["gnomad"],
      color_circos_track["homo_paf"],
      color_circos_track["hete_paf"],
      color_circos_track["homo_af"],
      color_circos_track["hete_af"],
      color_circos_track["somatic_af"]
    ),
    pch = c(NA, 15, 15, 15, 2, 3, 3),
    lty = c(1, NA, NA, NA, NA, NA, NA),
    lwd = c(2, NA, NA, NA, NA, NA, NA),
    pt.cex = c(NA, 1.6, 1.6, 1.6, 1.2, 1.2, 1.2),
    x.intersp = 0.8,
    y.intersp = 1.4,
    cex = 1.0,
    bty = "n",
    title = "Track Legend",
    title.font = 2,
    title.cex = 1.1
  )
}


# body --------------------------------------------------------------------
{
  outdir <- Sys.getenv("OUTDIR")
  pdf(
    file = path(
      outdir,
      "CIRCOS-ALL-VARIANTS-HOMO1266-HETE703-COLOR-FINAL.pdf"
    ),
    width = 13,
    height = 10
  )
  # Original color scheme
  fn_plot_circos(
    color_circos_track = c(
      "phastCons100way" = "#FFD700",
      "gnomad" = "#000000",
      "homo_paf" = "#ae00ff",
      "hete_paf" = "#00b3ff",
      "somatic_paf" = "#FF0000",
      "homo_af" = "#ae00ff",
      "hete_af" = "#00b3ff",
      "somatic_af" = "#FF0000",
      "gene_name_bg" = "#EAF7FF",
      "coverage" = "#3FB5FF"
    )
  )
  dev.off()
}


\() {
  outdir <- Sys.getenv("OUTDIR")
  pdf(
    file = path(
      outdir,
      "CIRCOS-ALL-VARIANTS-HOMO1266-HETE703-COLOR-SELECT.pdf"
    ),
    width = 13,
    height = 10
  )

  # Original color scheme
  fn_plot_circos(
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
  )

  # Color Scheme 1: Tol Bright (color-blind safe)
  fn_plot_circos(
    color_circos_track = c(
      "phastCons100way" = "#DDAA33",
      "gnomad" = "#004488",
      "homo_paf" = "#44AA99",
      "hete_paf" = "#BB5566",
      "somatic_paf" = "#000000",
      "homo_af" = "#77AADD",
      "hete_af" = "#EE6677",
      "somatic_af" = "#333333",
      "gene_name_bg" = "#F0F0F0",
      "coverage" = "#228833"
    )
  )

  # Color Scheme 2: Nature-inspired earth tones
  fn_plot_circos(
    color_circos_track = c(
      "phastCons100way" = "#E69F00",
      "gnomad" = "#0072B2",
      "homo_paf" = "#009E73",
      "hete_paf" = "#D55E00",
      "somatic_paf" = "#000000",
      "homo_af" = "#56B4E9",
      "hete_af" = "#CC79A7",
      "somatic_af" = "#2D2D2D",
      "gene_name_bg" = "#EFEFEF",
      "coverage" = "#F0E442"
    )
  )

  # Color Scheme 3: Okabe-Ito palette (universally accessible)
  fn_plot_circos(
    color_circos_track = c(
      "phastCons100way" = "#F0E442",
      "gnomad" = "#0173B2",
      "homo_paf" = "#029E73",
      "hete_paf" = "#DE8F05",
      "somatic_paf" = "#000000",
      "homo_af" = "#56B4E9",
      "hete_af" = "#CC78BC",
      "somatic_af" = "#404040",
      "gene_name_bg" = "#ECE7F2",
      "coverage" = "#FBAFE4"
    )
  )

  # Color Scheme 4: Vibrant but accessible
  fn_plot_circos(
    color_circos_track = c(
      "phastCons100way" = "#FDB863",
      "gnomad" = "#5E3C99",
      "homo_paf" = "#1B7837",
      "hete_paf" = "#D6604D",
      "somatic_paf" = "#000000",
      "homo_af" = "#B2ABD2",
      "hete_af" = "#F4A582",
      "somatic_af" = "#1A1A1A",
      "gene_name_bg" = "#F7F7F7",
      "coverage" = "#66C2A5"
    )
  )

  # Color Scheme 5: Cool professional tones
  fn_plot_circos(
    color_circos_track = c(
      "phastCons100way" = "#FFD92F",
      "gnomad" = "#4575B4",
      "homo_paf" = "#74ADD1",
      "hete_paf" = "#D73027",
      "somatic_paf" = "#313695",
      "homo_af" = "#ABD9E9",
      "hete_af" = "#F46D43",
      "somatic_af" = "#252525",
      "gene_name_bg" = "#E0F3F8",
      "coverage" = "#A50026"
    )
  )

  # Color Scheme 6: Warm accent palette
  fn_plot_circos(
    color_circos_track = c(
      "phastCons100way" = "#FEE08B",
      "gnomad" = "#998EC3",
      "homo_paf" = "#01665E",
      "hete_paf" = "#BF812D",
      "somatic_paf" = "#000000",
      "homo_af" = "#C7EAE5",
      "hete_af" = "#F6E8C3",
      "somatic_af" = "#303030",
      "gene_name_bg" = "#F5F5F5",
      "coverage" = "#5AB4AC"
    )
  )

  # Color Scheme 7: High contrast safe palette
  fn_plot_circos(
    color_circos_track = c(
      "phastCons100way" = "#FFCC33",
      "gnomad" = "#003366",
      "homo_paf" = "#00AA88",
      "hete_paf" = "#AA3377",
      "somatic_paf" = "#000000",
      "homo_af" = "#88CCEE",
      "hete_af" = "#EE7733",
      "somatic_af" = "#222222",
      "gene_name_bg" = "#EAEAEA",
      "coverage" = "#33BBEE"
    )
  )

  # Color Scheme 8: Muted scientific palette
  fn_plot_circos(
    color_circos_track = c(
      "phastCons100way" = "#E1BE6A",
      "gnomad" = "#40B0A6",
      "homo_paf" = "#5AB4AC",
      "hete_paf" = "#D8B365",
      "somatic_paf" = "#000000",
      "homo_af" = "#C7EAE5",
      "hete_af" = "#F5DEB3",
      "somatic_af" = "#2B2B2B",
      "gene_name_bg" = "#F6F6F6",
      "coverage" = "#01665E"
    )
  )

  # Color Scheme 9: Nature journal inspired
  fn_plot_circos(
    color_circos_track = c(
      "phastCons100way" = "#FDC086",
      "gnomad" = "#386CB0",
      "homo_paf" = "#7FC97F",
      "hete_paf" = "#F0027F",
      "somatic_paf" = "#000000",
      "homo_af" = "#BEAED4",
      "hete_af" = "#FDC086",
      "somatic_af" = "#353535",
      "gene_name_bg" = "#EDF8E9",
      "coverage" = "#BF5B17"
    )
  )

  # Color Scheme 10: Balanced dichromatic
  fn_plot_circos(
    color_circos_track = c(
      "phastCons100way" = "#FC8D62",
      "gnomad" = "#8DA0CB",
      "homo_paf" = "#66C2A5",
      "hete_paf" = "#E78AC3",
      "somatic_paf" = "#000000",
      "homo_af" = "#A6D854",
      "hete_af" = "#FFD92F",
      "somatic_af" = "#292929",
      "gene_name_bg" = "#E5F5E0",
      "coverage" = "#E5C494"
    )
  )

  # Color Scheme 11: Strong contrast scientific
  fn_plot_circos(
    color_circos_track = c(
      "phastCons100way" = "#EDAD08",
      "gnomad" = "#1965B0",
      "homo_paf" = "#7BAFDE",
      "hete_paf" = "#DC050C",
      "somatic_paf" = "#000000",
      "homo_af" = "#B17BA6",
      "hete_af" = "#F7B6D2",
      "somatic_af" = "#3C3C3C",
      "gene_name_bg" = "#F2F2F2",
      "coverage" = "#4EB265"
    )
  )

  dev.off()
}
# footer ------------------------------------------------------------------

# save image --------------------------------------------------------------
