fn_plot_mtdna_circos <- function(
    start.degree = 90,
    canvas.xlim = c(-1, 1),
    canvas.ylim = c(-1, 1),
    gap.degree = 1
    LENGTH <- 16569,
    gtf_gene_df = import("/home/liuc9/github/scMOCHA-data/config/mtdna_genes_dloop.qs"),
    phastCons100way = import("/home/liuc9/github/scMOCHA-data/config/chrM.phastCons100way.wigFix.qs"),
    gnomad = import("/home/liuc9/github/scMOCHA-data/analysis/zzz/db/gnomad.qs"),
    coverage = import("/home/liuc9/github/scMOCHA-data/analysis/zzz/clean-data/gse_data_coverage.fst"),
    all_variant = import("/home/liuc9/github/scMOCHA-data/analysis/zzz/clean-data/all_variant.qs")
    ) {


  gnomad <- gnomad |>
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

  coverage <- coverage |>
    dplyr::mutate(
      seqnames = "MT",
      start = pos,
      end = pos
    ) |>
    dplyr::select(
      seqnames,
      start,
      end,
      depth
    ) |>
    dplyr::mutate(
      depth = log10(depth + 1)
    )

  all_variant <- all_variant |>
    dplyr::mutate(
      paf = n / 577
    ) |>
    dplyr::arrange(Position)

  all_variant |>
    dplyr::filter(issomatic == "homoplasmic") |>
    dplyr::mutate(
      seqnames = "MT",
      start = Position,
      end = Position
    ) |>
    dplyr::select(
      seqnames,
      start,
      end,
      variant,
      paf,
      af,
      Disease
    ) ->
  homoplasmic_variant_af

  all_variant |>
    dplyr::filter(issomatic == "heteroplasmic") |>
    dplyr::mutate(
      seqnames = "MT",
      start = Position,
      end = Position
    ) |>
    dplyr::select(
      seqnames,
      start,
      end,
      variant,
      paf,
      af,
      Disease
    ) ->
  heteroplasmic_variant_af

  heteroplasmic_variant_af |>
    dplyr::arrange(-paf) |>
    head(5) ->
  top_variants

  source("/home/liuc9/github/scMOCHA-data/analysis/00-colors.R")

  library(circlize)

  circos.clear()
  circos.par(
    start.degree = start.degree,
    canvas.xlim = canvas.xlim,
    canvas.ylim = canvas.ylim,
    gap.degree = gap.degree
  )

  circos.genomicInitialize(
    data = coverage,
    plotType = "axis",
    axis.labels.cex = 0.8 * par("cex"),
  )

  gtf_gene_df |>
    dplyr::filter(TYPE %in% c("D-Loop", "MT rRNA")) ->
  highlight_df

  for (i in seq_len(nrow(highlight_df))) {
    pos = circlize(
      c(highlight_df$start[i], highlight_df$end[i]),
      c(0, 1),
      sector.index = "MT"
    )
    draw.sector(
      pos[1, "theta"],
      pos[2, "theta"],
      rou1 = 0.95,
      rou2 = 0.09,
      clock.wise = TRUE,
      col = prismatic::clr_alpha(
        highlight_df$COLOR[i],
        alpha = 0.3
      ),
      border = NA
    )
  }

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

  circos.genomicTrack(
    homoplasmic_variant_af,
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
        col = color_circos_track["homoplasmic_paf"],
        border = color_circos_track["homoplasmic_paf"]
      )
    }
  )

  circos.genomicTrack(
    heteroplasmic_variant_af,
    track.height = 0.1,
    ylim = c(0, 1),
    track.margin = c(0, 0),
    cell.padding = c(0, 0, 0, 0),
    bg.border = NA,
    bg.col = NA,
    panel.fun = function(region, value, ...) {
      pos = region$start
      val = value$paf
      circos.barplot(
        value = val,
        pos = pos,
        col = color_circos_track["heteroplasmic_paf"],
        border = color_circos_track["heteroplasmic_paf"]
      )
    }
  )

  circos.genomicTrack(
    gtf_gene_df,
    ylim = c(0, 1),
    track.height = 0.1,
    track.margin = c(0, 0.01),
    cell.padding = c(0, 0, 0, 0),
    bg.border = color_circos_track["gene_name_bg"],
    bg.col = color_circos_track["gene_name_bg"],
    panel.fun = function(region, value, ...) {
      genetypes <- sort(unique(value$TYPE))
      for (genetype in genetypes) {
        r <- region[value$TYPE == genetype, ]
        v <- value[value$TYPE == genetype, ]
        v_gene_name <- gsub(
          pattern = "MT-|-",
          replacement = "",
          x = v$gene_name
        )
        circos.genomicText(
          region = r,
          value = v,
          y = 0,
          adj = c(0, 0.5),
          labels = v_gene_name,
          facing = "clockwise",
          niceFacing = TRUE,
          cex = 0.9,
          col = "black",
        )
      }
    }
  )

  circos.genomicTrack(
    gtf_gene_df,
    track.height = 0.05,
    ylim = c(0, 1),
    bg.border = NA,
    track.margin = c(0, 0),
    cell.padding = c(0, 0, 0, 0),
    panel.fun = function(region, value, ...) {
      genetypes <- sort(unique(value$TYPE))
      for (genetype in genetypes) {
        r <- region[value$TYPE == genetype, ]
        v <- value[value$TYPE == genetype, ]
        vcol <- value[value$TYPE == genetype, ]$COLOR
        circos.genomicRect(
          region = r,
          value = v,
          col = vcol,
          border = "white",
          lty = 1
        )
      }
    },
  )

  circos.genomicTrack(
    homoplasmic_variant_af,
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
        col = color_circos_track["homoplasmic_af"],
        bg.border = color_circos_track["homoplasmic_af"]
      )
    }
  )

  circos.genomicTrack(
    heteroplasmic_variant_af,
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
        pch = 3,
        cex = 0.5,
        col = color_circos_track["heteroplasmic_af"],
        bg.border = color_circos_track["heteroplasmic_af"]
      )
    }
  )

  circos.genomicTrack(
    coverage,
    track.height = 0.1,
    ylim = c(0, max(coverage$depth)),
    track.margin = c(0, 0.01),
    cell.padding = c(0, 0, 0, 0),
    bg.border = NA,
    bg.col = NA,
    panel.fun = function(region, value, ...) {
      circos.genomicLines(
        region = region,
        value = value,
        col = color_circos_track["coverage"],
        border = color_circos_track["coverage"],
        area = TRUE
      )
    }
  )

  circos.clear()
}
