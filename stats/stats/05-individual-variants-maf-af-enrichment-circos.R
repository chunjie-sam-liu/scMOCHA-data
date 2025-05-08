fn_plot_mtdna_circos <- function(
    start.degree = 90,
    canvas.xlim = c(-1, 1),
    canvas.ylim = c(-1, 1),
    gap.degree = 1) {
  LENGTH <- 16569
  gtf_gene_df <- readr::read_rds("/home/liuc9/github/scMOCHA-data/config/mtdna_genes_dloop.rds.gz")

  phastCons100way <- data.table::fread(
    "/home/liuc9/github/scMOCHA-data/config/chrM.phastCons100way.wigFix"
  ) |>
    tibble::rowid_to_column() |>
    tibble::add_column(
      seqnames = "MT",
      .before = 1
    ) |>
    dplyr::mutate(
      start = rowid,
      end = rowid
    ) |>
    dplyr::select(
      seqnames,
      start = rowid,
      end = rowid,
      phastCons100wayScore = "fixedStep chrom=chrM start=1 step=1"
    )

  # af_hom is gnomad AF
  gnomad <- data.table::fread(
    "/home/liuc9/github/scMOCHA-data/stats/stats/zzz/db/gnomad.csv"
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

  # coverage
  coverage <- data.table::fread(
    file.path(
      "/home/liuc9/github/scMOCHA-data/stats/stats/zzz/clean-data/", "gse_data_coverage.csv"
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
      depth
    )

  # all variants
  all_variant <- readr::read_rds(
    file.path(
      "/home/liuc9/github/scMOCHA-data/stats/stats/zzz/clean-data/", "all_variant.rds"
    )
  ) |>
    dplyr::mutate(
      paf = n / 577
    ) |>
    dplyr::arrange(Position)

  all_variant |>
    dplyr::filter(issomatic == "homoplasmic") |>
    # dplyr::select(Position, paf) |>
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
    # dplyr::select(Position, paf) |>
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

  circos_track_colors <- c(
    "phastCons100way" = "#FFD700",
    "gnomad" = "#0000FF",
    "homoplasmic_paf" = "#03DC62",
    "heteroplasmic_paf" = "#7A0202",
    "homoplasmic_af" = "#A6FCCB",
    "heteroplasmic_af" = "#FF0000",
    "gene_name_bg" = "#EAF7FF",
    "coverage" = "#3FB5FF"
  )

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
    axis.labels.cex = 0.8 * par("cex"),
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
      rou2 = 0.09,
      clock.wise = TRUE,
      col = prismatic::clr_alpha(
        highlight_df$COLOR[i],
        alpha = 0.3
      ),
      border = NA
    )
  }


  # ! phastCons100way --------------------------------------------------------------------

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
        col = circos_track_colors["phastCons100way"],
        lwd = 0.5
      )
    }
  )

  # ! gnomad --------------------------------------------------------------------


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
        col = circos_track_colors["gnomad"],
        border = circos_track_colors["gnomad"]
      )
    }
  )


  # ! homoplasmic paf--------------------------------------------------------------------


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
        col = circos_track_colors["homoplasmic_paf"],
        border = circos_track_colors["homoplasmic_paf"]
      )
    }
  )

  # ! heteroplasmic paf--------------------------------------------------------------------


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
        col = circos_track_colors["heteroplasmic_paf"],
        border = circos_track_colors["heteroplasmic_paf"]
      )
    }
  )

  # ! gene name--------------------------------------------------------------------

  circos.genomicTrack(
    gtf_gene_df,
    ylim = c(0, 1),
    track.height = 0.1,
    track.margin = c(0, 0.01),
    cell.padding = c(0, 0, 0, 0),
    bg.border = circos_track_colors["gene_name_bg"],
    bg.col = circos_track_colors["gene_name_bg"],
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



  # ! homoplasmic af--------------------------------------------------------------------


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
      # circos.barplot(
      #   value = val,
      #   pos = pos,
      #   col = "#a6fccb",
      #   border = "#a6fccb"
      # )
      circos.genomicPoints(
        region = region,
        value = value,
        pch = 2,
        cex = 0.5,
        col = circos_track_colors["homoplasmic_af"],
        bg.border = circos_track_colors["homoplasmic_af"]
      )
    }
  )


  # ! heteroplasmic af--------------------------------------------------------------------


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
        col = circos_track_colors["heteroplasmic_af"],
        bg.border = circos_track_colors["heteroplasmic_af"]
      )
    }
  )

  # ! coverage --------------------------------------------------------------------
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
        col = circos_track_colors["coverage"],
        border = circos_track_colors["coverage"],
        area = TRUE
      )
    }
  )



  # ! labels --------------------------------------------------------------------
  # circos.labels(
  #   rep("MT", 5),
  #   x = top_variants$start,
  #   labels = top_variants$variant,
  #   side = "inside",
  #   cex = 0.5,
  #   track.margin = c(0, 0.01),
  # )

  circos.clear()
}

{
  pdf(
    file = "/home/liuc9/github/scMOCHA-data/stats/stats/zzz/heteroplasmic/circos-homo-hetero.pdf",
    width = 13,
    height = 10
  )
  fn_plot_mtdna_circos()
  dev.off()
}


{
  pdf(
    file = "/home/liuc9/github/scMOCHA-data/stats/stats/zzz/heteroplasmic/circos-homo-hetero-90.pdf",
    width = 13,
    height = 10
  )
  fn_plot_mtdna_circos(
    start.degree = 90,
    canvas.xlim = c(-1, 1),
    canvas.ylim = c(-1, 1),
    gap.degree = 30
  )
  dev.off()
}
