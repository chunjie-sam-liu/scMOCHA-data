#!/usr/bin/env Rscript
# Metainfo ----------------------------------------------------------------
# @AUTHOR: Chun-Jie Liu
# @CONTACT: chunjie.sam.liu.at.gmail.com
# @DATE: 2025-12-18 14:01:53
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

outdir <- path("/home/liuc9/github/scMOCHA-data/analysis/zzz/MANUSCRIPTFIGURES")

METAFULL <- import(outdir / "SAMPLES-METADATA-FULL.xlsx")
ALLVARIANTS <- import(
  outdir / "SAMPLE-VARIANT-CLASSIFICATION-CLUSTER-BULK-AF.xlsx"
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

SOMATIC_VARIANTS <- ALLVARIANTS |>
  dplyr::filter(variant_type == "somatic")
HOMO_VARIANTS <- ALLVARIANTS |>
  dplyr::filter(variant_type == "homo")
HETE_VARIANTS <- ALLVARIANTS |>
  dplyr::filter(variant_type == "hete")

# load conn ---------------------------------------------------------------

# src ---------------------------------------------------------------------
source("/home/liuc9/github/scMOCHA-data/analysis/high-res/00-colors.R")
source("/home/liuc9/github/scMOCHA-data/analysis/high-res/plot_mtdna.R")
source(
  "/home/liuc9/github/scMOCHA-data/analysis/high-res/plot_variant_celltype_af_haplogroup.R"
)


# function ----------------------------------------------------------------
fn_plot_freq_mtdna <- function(df, label_variants = NULL) {
  fn_xy_breaks_limits(
    df$n_individuals,
  ) -> ybl_

  df |> dplyr::filter(variant %in% label_variants) -> forlabel
  df |>
    ggplot(aes(
      x = start,
      y = n_individuals
    )) +
    geom_segment(
      aes(x = start, xend = start, y = 0, yend = n_individuals)
    ) +
    geom_point(
      aes(size = n_individuals),
      color = "red",
      fill = "red",
      # alpha = 0.7,
      shape = 21,
      stroke = 1
    ) +
    ggrepel::geom_text_repel(
      data = forlabel,
      aes(label = variant),
      # size = 3,
      nudge_y = -0.1,
      nudge_x = 0.1,
      show.legend = FALSE,
      max.overlaps = Inf,
      segment.size = 0.2,
      segment.color = "black",
      box.padding = 0.5
    ) +
    scale_x_continuous(
      expand = expansion(mult = c(0, 0.01)),
      limits = c(1, 16569),
      breaks = c(seq(0, 17000, 1000), 16569),
      labels = c(seq(0, 17000, 1000), 16569),
    ) +
    scale_y_continuous(
      expand = expansion(add = c(.05, 0.05), mult = c(0.01, 0.01)),
      limits = ybl_$limits,
      breaks = ybl_$breaks,
      # labels = ybl_$limits,
    ) +
    theme(
      plot.margin = margin(t = 0, b = 0, unit = "cm"),
      # panel.background = element_rect(color = "red"),
      panel.background = element_blank(),
      panel.grid = element_blank(),
      axis.line.y.left = element_line(color = "black"),
      # axis.line.x.bottom = element_line(color = "black"),
      axis.ticks.x = element_blank(),
      axis.text.x = element_blank(),
      axis.line.x = element_blank(),
      axis.title.x = element_blank(),
      # legend.position = c(0.8, 0.5),
      legend.position = "none",
      legend.key = element_blank(),
      axis.title.y = element_text(size = 16, color = "black", face = "bold"),
      axis.text.y = element_text(color = "black", size = 12),
      legend.text = element_text(
        size = 14,
        color = "black"
      ),
      legend.title = element_text(
        size = 16,
        colour = "black"
      ),
      strip.background = element_blank(),
      strip.text = element_text(
        size = 8,
        color = "black",
        face = "bold"
      ),
      plot.title = element_text(
        size = 16,
        color = "black",
        hjust = 0.5,
        face = "bold"
      )
    ) +
    labs(
      title = "mtDNA Variant Hotspots",
      y = "# of Samples",
    )
}

fn_plot_somatic_dist <- function() {
  SOMATIC_VARIANTS |>
    dplyr::mutate(srridvariant = glue("{srrid}_{variant}")) |>
    dplyr::select(
      srridvariant,
      c(B, CD4_T, CD8_T, DC, Mono, NK, other, Bulk, other_T)
    ) |>
    tibble::column_to_rownames("srridvariant") |>
    as.matrix() -> mat_somatic_af

  # Identify dominant cell type for each variant (excluding Bulk)
  mat_specific <- mat_somatic_af[, c(
    "B",
    "CD4_T",
    "CD8_T",
    "other_T",
    "NK",
    "DC",
    "Mono",
    "other"
  )]

  # For each variant, find the cell type with maximum AF
  dominant_celltype <- apply(mat_specific, 1, function(x) {
    if (max(x, na.rm = TRUE) < 0.01) {
      return("Bulk")
    }
    colnames(mat_specific)[which.max(x)]
  })

  # Create ordering: group by dominant cell type, then by AF within each group
  celltype_order <- c(
    "Bulk",
    "B",
    "CD4_T",
    "CD8_T",
    "other_T",
    "NK",
    "DC",
    "Mono",
    "other"
  )

  variant_df <- data.frame(
    srridvariant = rownames(mat_somatic_af),
    dominant_celltype = dominant_celltype,
    max_af = apply(mat_specific, 1, max, na.rm = TRUE),
    stringsAsFactors = FALSE
  )

  variant_df$dominant_celltype <- factor(
    variant_df$dominant_celltype,
    levels = celltype_order
  )
  variant_ordered <- variant_df[
    order(variant_df$dominant_celltype, -variant_df$max_af),
  ]
  clustered_order <- variant_ordered$srridvariant

  SOMATIC_VARIANTS |>
    dplyr::mutate(srridvariant = glue("{srrid}_{variant}")) |>
    tidyr::pivot_longer(
      cols = c(B, CD4_T, CD8_T, DC, Mono, NK, other, Bulk, other_T),
      names_to = "celltype",
      values_to = "af"
    ) |>
    dplyr::mutate(
      celltype = factor(
        celltype,
        levels = c(
          "Bulk",
          "B",
          "CD4_T",
          "CD8_T",
          "other_T",
          "NK",
          "DC",
          "Mono",
          "other"
        )
      )
    ) |>
    dplyr::mutate(
      srridvariant = factor(
        srridvariant,
        levels = clustered_order
      )
    ) -> forplot

  forplot |>
    dplyr::select(Haplogroup, Verbose_haplogroup) |>
    dplyr::distinct() |>
    dplyr::mutate(
      Haplogroup_s = purrr::map_chr(
        .x = Haplogroup,
        .f = \(.x) {
          # if (stringr::str_starts(.x, "L")) {
          #   gsub("L", "L0", .x)
          # }
          gsub("\\d+.*", "", .x)
        }
      )
    ) |>
    dplyr::mutate(
      color_haplogroup = color(color_haplogroup[Haplogroup_s])
    ) |>
    dplyr::mutate(
      color_verbose_haplogroup = ifelse(
        Haplogroup == Verbose_haplogroup,
        color_haplogroup,
        prismatic::clr_lighten(color_haplogroup, 0.5)
      )
    ) |>
    dplyr::filter(!is.na(Haplogroup)) |>
    dplyr::select(-Haplogroup_s) -> haplo_colors

  c(
    haplo_colors$color_haplogroup,
    haplo_colors$color_verbose_haplogroup
  ) -> haplo_color_vector
  names(haplo_color_vector) <- c(
    haplo_colors$Haplogroup,
    haplo_colors$Verbose_haplogroup
  )

  forplot |>
    dplyr::select(srridvariant, srrid, Haplogroup, Verbose_haplogroup) |>
    dplyr::mutate(
      srridvariant = factor(
        srridvariant,
        levels = clustered_order
      )
    ) |>
    tidyr::pivot_longer(
      cols = c(Haplogroup, Verbose_haplogroup),
      names_to = "type",
      values_to = "Haplogroup"
    ) -> forplot_haplogroup

  forplot_haplogroup |>
    ggplot(aes(
      x = type,
      y = srridvariant,
      fill = Haplogroup
    )) +
    geom_tile(
      show.legend = FALSE
    ) +
    scale_fill_manual(
      values = haplo_color_vector
    ) +
    geom_text(
      data = forplot_haplogroup |> dplyr::distinct(),
      aes(
        label = Haplogroup,
      ),
      color = "black",
      fontface = "bold"
    ) +
    scale_x_discrete(
      expand = c(0, 0)
    ) +
    scale_y_discrete(
      labels = function(x) sub("_.*", "", x)
    ) +
    theme_classic() +
    theme(
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      axis.title.x = element_blank(),
      axis.line = element_blank(),
      axis.title.y = element_blank(),
      axis.text.y = element_text(
        face = "bold",
        size = 12
      ),
      axis.ticks.y = element_blank(),
    ) -> p_haplogroup

  gtf_gene_df <- import(
    "/home/liuc9/github/scMOCHA-data/config/mtdna_genes_dloop.qs"
  )

  forplot |>
    dplyr::select(srridvariant, variant) |>
    dplyr::distinct() |>
    dplyr::mutate(
      pos = as.integer(gsub("(\\d+)[A-Z]>[A-Z]", "\\1", variant))
    ) |>
    dplyr::mutate(
      gene_df = purrr::map(
        .x = pos,
        .f = \(.x) {
          gtf_gene_df |>
            dplyr::filter(start <= .x & end >= .x) |>
            dplyr::select(gene_name, TYPE, COLOR)
        }
      )
    ) |>
    tidyr::unnest(cols = gene_df) -> forplot_variant_gene

  # Create color vector: COLOR is for TYPE, lighten for different gene_name
  type_color_df <- forplot_variant_gene |>
    dplyr::select(TYPE, COLOR) |>
    dplyr::distinct()

  # For each TYPE, create lightened colors for each gene_name
  gene_colors_list <- forplot_variant_gene |>
    dplyr::select(TYPE, gene_name, COLOR) |>
    dplyr::distinct() |>
    dplyr::arrange(TYPE, gene_name) |>
    dplyr::group_by(TYPE) |>
    dplyr::mutate(
      gene_idx = dplyr::row_number(),
      n_genes = dplyr::n()
    ) |>
    dplyr::ungroup() |>
    dplyr::mutate(
      gene_color = purrr::pmap_chr(
        .l = list(
          .color = COLOR,
          .idx = gene_idx
        ),
        .f = \(.color, .idx) {
          if (.idx == 1) {
            .color
          } else {
            as.character(prismatic::clr_lighten(.color, shift = 0.05 * .idx))
          }
        }
      )
    )

  # Create named vectors
  type_color_vector <- type_color_df$COLOR
  names(type_color_vector) <- type_color_df$TYPE

  gene_color_vector <- gene_colors_list$gene_color
  names(gene_color_vector) <- gene_colors_list$gene_name

  all_color_vector <- c(type_color_vector, gene_color_vector)

  forplot_variant_gene |>
    dplyr::select(srridvariant, gene_name, TYPE, COLOR) |>
    dplyr::mutate(
      srridvariant = factor(srridvariant, levels = clustered_order)
    ) |>
    tidyr::pivot_longer(
      cols = c(TYPE, gene_name),
      names_to = "type",
      values_to = "label"
    ) |>
    dplyr::mutate(
      type = factor(type, levels = c("TYPE", "gene_name"))
    ) -> forplot_gene

  forplot_gene |>
    ggplot(aes(
      x = type,
      y = srridvariant,
      fill = label
    )) +
    geom_tile(
      show.legend = FALSE
    ) +
    scale_fill_manual(
      values = all_color_vector
    ) +
    geom_text(
      data = forplot_gene |> dplyr::distinct(),
      aes(label = label),
      color = "black",
      fontface = "bold",
      size = 3
    ) +
    scale_x_discrete(
      expand = c(0, 0)
    ) +
    scale_y_discrete(
      labels = function(x) sub(".*_", "", x)
    ) +
    theme_classic() +
    theme(
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      axis.title.x = element_blank(),
      axis.line = element_blank(),
      axis.title.y = element_blank(),
      axis.text.y = element_text(
        face = "bold",
        size = 10
      ),
      axis.ticks.y = element_blank(),
    ) -> p_gene

  forplot |>
    ggplot(aes(
      x = celltype,
      y = srridvariant,
      fill = af
    )) +
    geom_tile() +
    geom_text(
      aes(
        label = ifelse(af >= 0.01, sprintf("%.2f", af), "")
      ),
      color = "black",
      fontface = "bold"
    ) +
    scale_fill_gradient(
      name = "AF",
      low = "white",
      high = "red"
    ) +
    scale_x_discrete(
      expand = c(0, 0)
    ) +
    scale_y_discrete(
      labels = function(x) sub(".*_", "", x)
    ) +
    theme_classic() +
    labs(
      x = "Celltype",
    ) +
    theme(
      plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
      axis.text.x = element_text(
        face = "bold",
        size = 12,
        angle = 15,
        hjust = 1
      ),
      # axis.title.x = element_text(face = "bold", size = 12),
      axis.title = element_blank(),
      axis.line.y = element_blank(),
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
    ) -> p_af

  SOMATIC_VARIANTS$variant |> unique() |> length() -> n_variants
  SOMATIC_VARIANTS$srrid |> unique() |> length() -> n_samples

  wrap_plots(
    p_haplogroup,
    p_gene,
    p_af,
    ncol = 3,
    widths = c(0.2, 0.3, 1),
    guides = "collect"
  ) +
    plot_annotation(
      title = glue::glue(
        "Somatic variants distribution across cell types and haplogroups\n({n_variants} variants across {n_samples} samples)"
      ),
      theme = theme(
        plot.title = element_text(size = 14, face = "bold", hjust = 0.5)
      )
    ) -> p_collect

  pdf(
    file = path(
      "/home/liuc9/github/scMOCHA-data/analysis/high-res/MANUSCRIPTFIGURES-notuse"
    ) /
      "HOTSPOTS-SOMATIC-CLUSTER.pdf",
    width = 15,
    height = 35
  )
  print(p_collect)
  dev.off()
}

fn_plot_somatic_dist()
# body --------------------------------------------------------------------

#
#
# ? plot somatic --------------------------------------------------------------------
#
#
{
  SOMATIC_VARIANTS |>
    tidyr::nest(
      .by = c(variant, start)
    ) |>
    dplyr::mutate(
      n_individuals = purrr::map_int(
        .x = data,
        .f = nrow
      )
    ) |>
    dplyr::arrange(
      desc(n_individuals)
    ) |>
    dplyr::select(
      variant,
      start,
      n_individuals
    ) -> SOMATIC_VARIANT_N_INDIVIDUALS

  wrap_plots(
    fn_plot_freq_mtdna(
      SOMATIC_VARIANT_N_INDIVIDUALS,
      label_variants = c("14082C>G", "15169A>G", "3240C>G")
    ),
    plot_spacer(),
    fn_plot_mtdna(),
    ncol = 1,
    heights = c(15, -0.7, 1)
  ) -> p_somatic_hotspots

  pdf(
    file = outdir / "HOTSPOTS-SOMATIC.pdf",
    width = 11,
    height = 6
  )
  print(p_somatic_hotspots)
  dev.off()

  p_somatic_14082 <- fn_plot_variant_celltype_af_haplogroup(
    thevariant = "14082C>G",
    vtype = "somatic"
  )$p_collect_2

  pdf(
    file = outdir / "HOTSPOTS-SOMATIC-EXAMPLE-14082C>G.pdf",
    width = 10,
    height = 6
  )
  print(p_somatic_14082)
  dev.off()

  p_somatic_15169 <- fn_plot_variant_celltype_af_haplogroup(
    thevariant = "15169A>G",
    vtype = "somatic"
  )$p_collect_2

  pdf(
    file = outdir / "HOTSPOTS-SOMATIC-EXAMPLE-15169A>G.pdf",
    width = 10,
    height = 4
  )
  print(p_somatic_15169)
  dev.off()

  p_somatic_3240 <- fn_plot_variant_celltype_af_haplogroup(
    thevariant = "3240C>G",
    vtype = "somatic"
  )$p_collect_2

  pdf(
    file = outdir / "HOTSPOTS-SOMATIC-EXAMPLE-3240C>G.pdf",
    width = 10,
    height = 4
  )
  print(p_somatic_3240)
  dev.off()
}


#
#
# ? hete --------------------------------------------------------------------
#
#

{
  HETE_VARIANTS |>
    tidyr::nest(
      .by = c(variant, start)
    ) |>
    dplyr::mutate(
      n_individuals = purrr::map_int(
        .x = data,
        .f = nrow
      )
    ) |>
    dplyr::arrange(
      desc(n_individuals)
    ) |>
    dplyr::select(
      variant,
      start,
      n_individuals
    ) -> HETE_VARIANT_N_INDIVIDUALS

  wrap_plots(
    fn_plot_freq_mtdna(
      HETE_VARIANT_N_INDIVIDUALS,
      label_variants = c("3173G>A", "1670A>G", "15666T>C", "2442T>C", "7428G>A")
    ),
    plot_spacer(),
    fn_plot_mtdna(),
    ncol = 1,
    heights = c(15, -0.7, 1)
  ) -> p_hete_hotspots

  pdf(
    file = outdir / "HOTSPOTS-HETE.pdf",
    width = 11,
    height = 6
  )
  print(p_hete_hotspots)
  dev.off()

  p_hete_3173 <- fn_plot_variant_celltype_af_haplogroup(
    thevariant = "3173G>A",
    vtype = "hete"
  )$p_collect_2

  pdf(
    file = outdir / "HOTSPOTS-HETE-EXAMPLE-3173G>A.pdf",
    width = 10,
    height = 15
  )
  print(p_hete_3173)
  dev.off()

  p_hete_1670 <- fn_plot_variant_celltype_af_haplogroup(
    thevariant = "1670A>G",
    vtype = "hete"
  )$p_collect_2

  pdf(
    file = outdir / "HOTSPOTS-HETE-EXAMPLE-1670A>G.pdf",
    width = 10,
    height = 17
  )
  print(p_hete_1670)
  dev.off()

  p_hete_15666 <- fn_plot_variant_celltype_af_haplogroup(
    thevariant = "15666T>C",
    vtype = "hete"
  )$p_collect_2

  pdf(
    file = outdir / "HOTSPOTS-HETE-EXAMPLE-15666T>C.pdf",
    width = 10,
    height = 17
  )
  print(p_hete_15666)
  dev.off()

  p_hete_2442 <- fn_plot_variant_celltype_af_haplogroup(
    thevariant = "2442T>C",
    vtype = "hete"
  )$p_collect_2
  pdf(
    file = outdir / "HOTSPOTS-HETE-EXAMPLE-2442T>C.pdf",
    width = 10,
    height = 15
  )
  print(p_hete_2442)
  dev.off()

  p_hete_7428 <- fn_plot_variant_celltype_af_haplogroup(
    thevariant = "7428G>A",
    vtype = "hete"
  )$p_collect_2
  pdf(
    file = outdir / "HOTSPOTS-HETE-EXAMPLE-7428G>A.pdf",
    width = 10,
    height = 15
  )
  print(p_hete_7428)
  dev.off()
}


#
#
# ? homo --------------------------------------------------------------------
#
#

{
  HOMO_VARIANTS |>
    tidyr::nest(
      .by = c(variant, start)
    ) |>
    dplyr::mutate(
      n_individuals = purrr::map_int(
        .x = data,
        .f = nrow
      )
    ) |>
    dplyr::arrange(
      desc(n_individuals)
    ) |>
    dplyr::select(
      variant,
      start,
      n_individuals
    ) -> HOMO_VARIANT_N_INDIVIDUALS

  wrap_plots(
    fn_plot_freq_mtdna(
      HOMO_VARIANT_N_INDIVIDUALS,
      label_variants = c("750A>G", "15326A>G", "2706A>G")
    ),
    plot_spacer(),
    fn_plot_mtdna(),
    ncol = 1,
    heights = c(15, -0.7, 1)
  ) -> p_homo_hotspots

  pdf(
    file = outdir / "HOTSPOTS-HOMO.pdf",
    width = 11,
    height = 6
  )
  print(p_homo_hotspots)
  dev.off()

  p_homo_750 <- fn_plot_variant_celltype_af_haplogroup(
    thevariant = "750A>G",
    vtype = "homo"
  )$p_collect_2

  pdf(
    file = outdir / "HOTSPOTS-HOMO-EXAMPLE-750A>G.pdf",
    width = 10,
    height = 60
  )
  print(p_homo_750)
  dev.off()

  p_homo_15326 <- fn_plot_variant_celltype_af_haplogroup(
    thevariant = "15326A>G",
    vtype = "homo"
  )$p_collect_2

  pdf(
    file = outdir / "HOTSPOTS-HOMO-EXAMPLE-15326A>G.pdf",
    width = 10,
    height = 60
  )
  print(p_homo_15326)
  dev.off()

  p_homo_2706 <- fn_plot_variant_celltype_af_haplogroup(
    thevariant = "2706A>G",
    vtype = "homo"
  )$p_collect_2

  pdf(
    file = outdir / "HOTSPOTS-HOMO-EXAMPLE-2706A>G.pdf",
    width = 10,
    height = 60
  )
  print(p_homo_2706)
  dev.off()
}

# footer ------------------------------------------------------------------

# save image --------------------------------------------------------------
