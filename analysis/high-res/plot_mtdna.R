fn_plot_mtdna <- function() {
  # mt_exons_df <- "/home/liuc9/github/scMOCHA/fasta/mt_exons.df.rds.gz"

  LENGTH <- 16569
  # rCRS <- Biostrings::readDNAStringSet("/home/liuc9/github/scMOCHA-data/config/rCRS.MT.fasta")
  gtf_gene_df <- import(
    "/home/liuc9/github/scMOCHA-data/config/mtdna_genes_dloop.qs"
  )

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
      arrowhead_height = unit(3, "mm"),
      arrowhead_width = unit(1, "mm"),
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
      plot.margin = margin(t = 0, b = 0, unit = "cm"),
      legend.position = "bottom",
      axis.title = element_blank(),
      axis.text.y = element_blank(),
      # axis.text.x = element_text(size = 14),
      # legend.text = element_text(size = 14),
      # panel.background = element_rect(
      #   color = "red"
      # ),
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
