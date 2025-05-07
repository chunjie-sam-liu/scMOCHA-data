fn_plot_mtdna_circle <- function() {
  LENGTH <- 16569
  # rCRS <- Biostrings::readDNAStringSet("/home/liuc9/github/scMOCHA-data/config/rCRS.MT.fasta")
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

  # conserve_rate


  library(circlize)


  # ! init --------------------------------------------------------------------


  circos.clear()
  circos.par(
    start.degree = 90
  )


  # ! axis --------------------------------------------------------------------

  circos.genomicInitialize(
    data = phastCons100way,
    plotType = "axis",
    axis.labels.cex = 0.8 * par("cex"),
  )

  # ! highlights --------------------------------------------------------------------
  gtf_gene_df |> dplyr::filter(TYPE %in% c("D-Loop", "MT rRNA")) -> highlight_df

  for (i in seq_len(nrow(highlight_df))) {
    pos = circlize(c(highlight_df$start[i], highlight_df$end[i]), c(0, 1), sector.index = "MT")
    draw.sector(
      pos[1, "theta"],
      pos[2, "theta"],
      rou1 = 0.95,
      rou2 = 0.2,
      clock.wise = TRUE,
      col = prismatic::clr_alpha(highlight_df$COLOR[i], alpha = 0.3),
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
        col = "gold",
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
        col = "blue",
        border = "blue"
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
        col = "#03DC62",
        border = "#03DC62"
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
        col = "#7a0202",
        border = "#7a0202"
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
    bg.border = "#EAF7FFFF",
    bg.col = "#EAF7FFFF",
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
        col = "#a6fccb",
        bg.border = "#a6fccb"
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
        col = "red",
        bg.border = "red"
      )
    }
  )


  # ! labels --------------------------------------------------------------------
  circos.labels(
    rep("MT", 5),
    x = top_variants$start,
    labels = top_variants$variant,
    side = "inside",
    cex = 0.5,
    track.margin = c(0, 0.01),
  )

  circos.clear()
}

{
  pdf(
    file = "/home/liuc9/github/scMOCHA-data/stats/stats/zzz/heteroplasmic/circos-homo-hetero.pdf",
    width = 13,
    height = 10
  )
  fn_plot_mtdna_circle()
  dev.off()
}
