fn_plot_mtdna <- function() {
  # mt_exons_df <- "/home/liuc9/github/scMOCHA/fasta/mt_exons.df.rds.gz"

  LENGTH <- 16569
  # rCRS <- Biostrings::readDNAStringSet("/home/liuc9/github/scMOCHA-data/config/rCRS.MT.fasta")
  gtf_gene_df <- readr::read_rds("/home/liuc9/github/scMOCHA-data/config/mtdna_genes_dloop.rds.gz")


  library(gggenes)
  ggplot(
    gtf_gene_df,
    aes(
      xmin = start,
      xmax = end,
      y = seqnames,
    )
  ) +
    # geom_gene_arrow() +
    geom_gene_arrow(
      aes(
        fill = COLOR
      ),
      arrowhead_height = unit(3, "mm"), arrowhead_width = unit(1, "mm"),
    ) +
    scale_fill_identity(
      name = "Gene type",
      guide = "legend",
      labels = c("MT rRNA", "Protein coding", "MT tRNA", "MT OLR", "D-Loop")
    ) +
    # scale_fill_brewer(
    #   palette = "Set1",
    #   name = "Gene type",
    #   labels = c("D-Loop", "MT rRNA", "MT tRNA", "Protein coding")
    # ) +
    ggrepel::geom_text_repel(
      aes(
        x = (start + end) / 2,
        label = gsub(
          pattern = "MT-",
          replacement = "",
          x = gene_name
        ),
      ),
      color = "black",
      # fill = "white",
      # nudge_x =1,
      # nudge_y =0.001,
      size = 3,
      show.legend = F,
      max.overlaps = Inf,
    ) +
    # scale_color_brewer(palette = "Set1") +
    scale_x_continuous(
      limits = c(0, LENGTH),
      breaks = c(seq(0, LENGTH, 1000), LENGTH),
      labels = c(seq(0, LENGTH, 1000), LENGTH),
      expand = expansion(mult = c(0, 0.01)),
    ) +
    scale_y_discrete(
      expand = expansion(mult = c(0, 0), add = c(0, 0))
    ) +
    # theme_genes() +
    theme(
      legend.position = "bottom",
      axis.title = element_blank(),
      axis.text.y = element_blank(),
      # axis.text.x = element_text(size = 14),
      # legend.text = element_text(size = 14),
      panel.background = element_blank(),
      panel.grid = element_blank(),
      axis.ticks.y = element_blank(),
      axis.ticks.x = element_line(color = "black"),
      axis.line.x = element_line(color = "black"),
      axis.text.x = element_text(
        vjust = -1,
      ),
    )
}
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


  library(circlize)

  circos.clear()
  circos.par(
    start.degree = 90
  )
  circos.genomicInitialize(
    data = phastCons100way,
    plotType = "axis",
    axis.labels.cex = 0.8 * par("cex"),
  )
  circos.genomicTrack(
    gnomad,
    ylim = c(0, 1),
    track.margin = c(0, 0.01),
    cell.padding = c(0, 0, 0, 0),
    track.height = 0.1,
    bg.border = NA,
    bg.col = NA,
    panel.fun = function(region, value, ...) {
      pos = region$start
      val = value$af
      circos.barplot(
        value = val,
        pos = pos,
        col = "red",
        border = "red"
      )
    }
  )
  circos.genomicTrack(
    phastCons100way,
    ylim = c(0, 1),
    track.margin = c(0, 0.01),
    cell.padding = c(0, 0, 0, 0),
    track.height = 0.1,
    bg.border = NA,
    bg.col = NA,
    panel.fun = function(region, value, ...) {
      circos.genomicLines(
        region = region,
        value = value,
        col = "#ADA32E",
        lwd = 1
      )
      # pos = region$start
      # val = value$phastCons100wayScore
      # circos.barplot(
      #   value = val,
      #   pos = pos,
      #   col = "#ADA32E",
      #   border = "#ADA32E"
      # )
    }
  )
  circos.genomicTrack(
    gtf_gene_df,
    ylim = c(0, 1),
    track.margin = c(0, 0.01),
    cell.padding = c(0, 0, 0, 0),
    track.height = 0.15,
    bg.border = "#E1F0DD",
    bg.col = "#E1F0DD",
    panel.fun = function(region, value, ...) {
      genetypes <- sort(unique(value$TYPE))
      for (genotype in genetypes) {
        r <- region[value$TYPE == genotype, ]
        v <- value[value$TYPE == genotype, ]
        v_gene_name <- gsub(
          pattern = "MT-",
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
    ylim = c(0, 1),
    bg.border = NA,
    track.height = 0.1,
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
