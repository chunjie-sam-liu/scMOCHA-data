fn_plot_mtdna_cov_circos <- function(
    start.degree = 80,
    canvas.xlim = c(-1, 1),
    canvas.ylim = c(-1, 1),
    gap.degree = 20,
    scaley = FALSE) {
  LENGTH <- 16569
  gtf_gene_df <- import("/home/liuc9/github/scMOCHA-data/config/mtdna_genes_dloop.qs")


  # coverage
  coverage <- import(
    file.path(
      "/home/liuc9/github/scMOCHA-data/analysis/zzz/clean-data/", "gse_data_coverage_chemistry.csv"
    )
  ) |>
    dplyr::mutate(
      seqnames = "MT",
      start = pos,
      end = pos
    ) |>
    dplyr::select(
      seqnames,
      start,
      end,
      depth,
      chemistry
    ) |>
    dplyr::mutate(
      depth = log2(depth + 1),
    )

  source("/home/liuc9/github/scMOCHA-data/analysis/00-colors.R")
  # conserve_rate

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

  circos.genomicInitialize(
    data = coverage,
    plotType = "axis",
    axis.labels.cex = 0.8
  )

  # ! highlights --------------------------------------------------------------------
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
      rou2 = 0.4,
      clock.wise = TRUE,
      col = prismatic::clr_alpha(
        highlight_df$COLOR[i],
        alpha = 0.3
      ),
      border = NA
    )
  }

  # ! coverage --------------------------------------------------------------------
  for (chem in names(color_chemistry)) {
    coverage |>
      dplyr::filter(chemistry == chem) ->
    coverage_chem

    maxy <- max(coverage$depth)

    if (scaley == TRUE) {
      maxy = max(coverage_chem$depth)
    }

    circos.genomicTrack(
      coverage_chem,
      track.height = 0.1,
      ylim = c(0, maxy),
      track.margin = c(0, 0.01),
      cell.padding = c(0, 0, 0, 0),
      bg.border = NA,
      bg.col = NA,
      panel.fun = function(region, value, ...) {
        circos.yaxis()
        circos.genomicLines(
          region = region,
          value = value,
          col = color_chemistry[chem],
          border = color_chemistry[chem],
          area = TRUE
        )
      }
    )
  }



  # ! gene name--------------------------------------------------------------------

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
      for (genotype in genetypes) {
        r <- region[value$TYPE == genotype, ]
        v <- value[value$TYPE == genotype, ]
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


  # ! gene region --------------------------------------------------------------------


  circos.genomicTrack(
    gtf_gene_df,
    track.height = 0.05,
    ylim = c(0, 1),
    bg.border = NA,
    track.margin = c(0, 0),
    cell.padding = c(0, 0, 0, 0),
    panel.fun = function(region, value, ...) {
      genetypes <- sort(unique(value$TYPE))
      for (genotype in genetypes) {
        r <- region[value$TYPE == genotype, ]
        v <- value[value$TYPE == genotype, ]
        vcol <- value[value$TYPE == genotype, ]$COLOR
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



  circos.clear()
}

{
  pdf(
    file = "/home/liuc9/github/scMOCHA-data/analysis/zzz/out-heteroplasmic/circos-coverage.pdf",
    width = 13,
    height = 10
  )
  fn_plot_mtdna_cov_circos()
  dev.off()
}
{
  pdf(
    file = "/home/liuc9/github/scMOCHA-data/analysis/zzz/out-heteroplasmic/circos-coverage-scaley.pdf",
    width = 13,
    height = 10
  )
  fn_plot_mtdna_cov_circos(scaley = TRUE)
  dev.off()
}
{
  pdf(
    file = "/home/liuc9/github/scMOCHA-data/analysis/zzz/out-heteroplasmic/circos-coverage-90.pdf",
    width = 13,
    height = 10
  )
  fn_plot_mtdna_cov_circos(
    start.degree = 90,
    gap.degree = 90
  )
  dev.off()
}
