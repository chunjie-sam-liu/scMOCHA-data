#' Plot Mitochondrial DNA Circos Plot with Variants and Annotations
#'
#' Creates a comprehensive circular plot (circos plot) of mitochondrial DNA showing
#' various genomic features including gene annotations, conservation scores, population
#' frequencies, variant allele frequencies, and sequencing coverage.
#'
#' @param a Numeric. Starting degree for the circular plot (default: 90)
#' @param b Numeric vector. X-axis limits for the canvas (default: c(-1, 1))
#' @param cc Numeric vector. Y-axis limits for the canvas (default: c(-1, 1))
#' @param d Numeric. Gap degree between sectors (default: 1)
#' @param e Numeric. Length of mitochondrial genome in base pairs (default: 16569)
#' @param f Data frame. Gene annotation data imported from mtdna_genes_dloop.qs
#' @param g Data frame. Conservation scores from phastCons100way analysis
#' @param h Data frame. Population frequency data from gnomAD database
#' @param i Data frame. Sequencing coverage data across mitochondrial positions
#' @param j Data frame. All variant data including homoplasmic and heteroplasmic variants
#'
#' @details
#' The function creates multiple concentric tracks in the circos plot:
#' - Outer track: Gene names and annotations
#' - Gene structure track: Visual representation of gene locations
#' - Conservation track: PhastCons conservation scores
#' - Population frequency track: gnomAD allele frequencies
#' - Variant frequency tracks: Homoplasmic and heteroplasmic variant frequencies
#' - Allele frequency tracks: Individual variant allele frequencies (points)
#' - Coverage track: Sequencing depth across positions
#'
#' Special highlighting is applied to D-Loop and MT rRNA regions.
#' Only variants with population allele frequency > 1% from gnomAD are displayed.
#' Note: Function implementation uses compact coding style for performance.
#'
#' @return NULL (creates a circos plot as side effect)
#'
#' @examples
#' # Basic usage with default parameters
#' fn_plot_mtdna_circos()
#'
#' # Custom canvas size and starting position
#' fn_plot_mtdna_circos(
#'   a = 0,
#'   b = c(-1.2, 1.2),
#'   cc = c(-1.2, 1.2)
#' )
#'
#' @note
#' Requires the circlize package. Function uses optimized parameter naming scheme.
#'
fn_plot_mtdna_circos <- function(a = 90, b = c(-1, 1), cc = c(-1, 1), d = 1, e = 16569, f = import("/home/liuc9/github/scMOCHA-data/config/mtdna_genes_dloop.qs"), g = import("/home/liuc9/github/scMOCHA-data/config/chrM.phastCons100way.wigFix.qs"), h = import("/home/liuc9/github/scMOCHA-data/analysis/zzz/db/gnomad.qs"), i = import("/home/liuc9/github/scMOCHA-data/analysis/zzz/clean-data/gse_data_coverage.fst"), j = import("/home/liuc9/github/scMOCHA-data/analysis/zzz/clean-data/all_variant.qs")) {
  h <- h |>
    dplyr::filter(filters == "PASS") |>
    dplyr::select(position, af = af_hom, ac = ac_hom) |>
    dplyr::arrange(position) |>
    dplyr::mutate(seqnames = "MT", start = position, end = position) |>
    dplyr::select(seqnames, start, end, af, ac) |>
    dplyr::filter(af > 0.01)
  i <- i |>
    dplyr::mutate(seqnames = "MT", start = pos, end = pos) |>
    dplyr::select(seqnames, start, end, depth) |>
    dplyr::mutate(depth = log10(depth + 1))
  j <- j |>
    dplyr::mutate(paf = n / 577) |>
    dplyr::arrange(Position)
  k <- j |>
    dplyr::filter(issomatic == "homoplasmic") |>
    dplyr::mutate(seqnames = "MT", start = Position, end = Position) |>
    dplyr::select(seqnames, start, end, variant, paf, af, Disease)
  l <- j |>
    dplyr::filter(issomatic == "heteroplasmic") |>
    dplyr::mutate(seqnames = "MT", start = Position, end = Position) |>
    dplyr::select(seqnames, start, end, variant, paf, af, Disease)
  m <- l |>
    dplyr::arrange(-paf) |>
    head(5)
  source("/home/liuc9/github/scMOCHA-data/analysis/00-colors.R")
  library(circlize)
  circos.clear()
  circos.par(start.degree = a, canvas.xlim = b, canvas.ylim = cc, gap.degree = d)
  circos.genomicInitialize(data = i, plotType = "axis", axis.labels.cex = 0.8 * par("cex"))
  n <- f |> dplyr::filter(TYPE %in% c("D-Loop", "MT rRNA"))
  for (o in seq_len(nrow(n))) {
    p = circlize(c(n$start[o], n$end[o]), c(0, 1), sector.index = "MT")
    draw.sector(p[1, "theta"], p[2, "theta"], rou1 = 0.95, rou2 = 0.09, clock.wise = TRUE, col = prismatic::clr_alpha(n$COLOR[o], alpha = 0.3), border = NA)
  }
  circos.genomicTrack(g, track.height = 0.03, ylim = c(0, 1), track.margin = c(0, 0.005), cell.padding = c(0, 0, 0, 0), bg.border = NA, bg.col = NA, panel.fun = function(region, value, ...) circos.genomicLines(region = region, value = value, col = color_circos_track["phastCons100way"], lwd = 0.5))
  circos.genomicTrack(h, track.height = 0.1, ylim = c(0, 1), track.margin = c(0, 0.01), cell.padding = c(0, 0, 0, 0), bg.border = NA, bg.col = NA, panel.fun = function(region, value, ...) {
    pos = region$start
    val = value$af
    circos.barplot(value = val, pos = pos, col = color_circos_track["gnomad"], border = color_circos_track["gnomad"])
  })
  circos.genomicTrack(k, track.height = 0.1, ylim = c(0, 1), track.margin = c(0, 0.01), cell.padding = c(0, 0, 0, 0), bg.border = NA, bg.col = NA, panel.fun = function(region, value, ...) {
    pos = region$start
    val = value$paf
    circos.barplot(value = val, pos = pos, col = color_circos_track["homoplasmic_paf"], border = color_circos_track["homoplasmic_paf"])
  })
  circos.genomicTrack(l, track.height = 0.1, ylim = c(0, 1), track.margin = c(0, 0), cell.padding = c(0, 0, 0, 0), bg.border = NA, bg.col = NA, panel.fun = function(region, value, ...) {
    pos = region$start
    val = value$paf
    circos.barplot(value = val, pos = pos, col = color_circos_track["heteroplasmic_paf"], border = color_circos_track["heteroplasmic_paf"])
  })
  circos.genomicTrack(f, ylim = c(0, 1), track.height = 0.1, track.margin = c(0, 0.01), cell.padding = c(0, 0, 0, 0), bg.border = color_circos_track["gene_name_bg"], bg.col = color_circos_track["gene_name_bg"], panel.fun = function(region, value, ...) {
    q <- sort(unique(value$TYPE))
    for (r in q) {
      s <- region[value$TYPE == r, ]
      t <- value[value$TYPE == r, ]
      u <- gsub(pattern = "MT-|-", replacement = "", x = t$gene_name)
      circos.genomicText(region = s, value = t, y = 0, adj = c(0, 0.5), labels = u, facing = "clockwise", niceFacing = TRUE, cex = 0.9, col = "black")
    }
  })
  circos.genomicTrack(f, track.height = 0.05, ylim = c(0, 1), bg.border = NA, track.margin = c(0, 0), cell.padding = c(0, 0, 0, 0), panel.fun = function(region, value, ...) {
    q <- sort(unique(value$TYPE))
    for (r in q) {
      s <- region[value$TYPE == r, ]
      t <- value[value$TYPE == r, ]
      v <- value[value$TYPE == r, ]$COLOR
      circos.genomicRect(region = s, value = t, col = v, border = "white", lty = 1)
    }
  })
  circos.genomicTrack(k, track.height = 0.1, ylim = c(0, 1), track.margin = c(0, 0.01), cell.padding = c(0, 0, 0, 0), bg.border = NA, bg.col = NA, panel.fun = function(region, value, ...) {
    pos = region$start
    val = value$af
    circos.genomicPoints(region = region, value = value, pch = 2, cex = 0.5, col = color_circos_track["homoplasmic_af"], bg.border = color_circos_track["homoplasmic_af"])
  })
  circos.genomicTrack(l, track.height = 0.1, ylim = c(0, 1), track.margin = c(0, 0.01), cell.padding = c(0, 0, 0, 0), bg.border = NA, bg.col = NA, panel.fun = function(region, value, ...) {
    pos = region$start
    val = value$af
    circos.genomicPoints(region = region, value = value, pch = 3, cex = 0.5, col = color_circos_track["heteroplasmic_af"], bg.border = color_circos_track["heteroplasmic_af"])
  })
  circos.genomicTrack(i, track.height = 0.1, ylim = c(0, max(i$depth)), track.margin = c(0, 0.01), cell.padding = c(0, 0, 0, 0), bg.border = NA, bg.col = NA, panel.fun = function(region, value, ...) circos.genomicLines(region = region, value = value, col = color_circos_track["coverage"], border = color_circos_track["coverage"], area = TRUE))
  circos.clear()
}

fn_plot_mtdna_circos()
