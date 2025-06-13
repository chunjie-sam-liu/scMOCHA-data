#!/usr/bin/env Rscript
# Metainfo ----------------------------------------------------------------
# @AUTHOR: Chun-Jie Liu
# @CONTACT: chunjie.sam.liu.at.gmail.com
# @DATE: 2025-06-11 11:48:13
# @DESCRIPTION: filename
# @VERSION: v0.0.1



# Library -----------------------------------------------------------------

suppressPackageStartupMessages(library(magrittr))
library(ggplot2)
library(patchwork)
library(prismatic)
library(paletteer)
library(data.table)
# library(rlang)
library(glue)
library(parallel)
library(GetoptLong)
library(logger)

# args --------------------------------------------------------------------

# s: string, i: integer, f: float, !: boolean
# @: array
# %: hash
# default: default value specified here.
verbose <- FALSE
spec <- "
Usage: Rscript foorbar.R [options]
Options:

<verbose!> Print messages
"

GetoptLong.options(help_style = "two-column")
GetoptLong(spec, template_control = list(opt_width = 21))

# src ---------------------------------------------------------------------

# header ------------------------------------------------------------------


# future: :plan(future: :multisession, workers = 10)


# load data ---------------------------------------------------------------
conn <- DBI::dbConnect(
  duckdb::duckdb(),
  dbdir = "/home/liuc9/github/scMOCHA-data/analysis/zzz/clean-data/all_variant_cell.duckdb.1.2.1"
)

# function ----------------------------------------------------------------


# body --------------------------------------------------------------------
dplyr::tbl(conn, "all_variant_cell") |>
  dplyr::filter(
    variant_in_cell_cluster == "cell"
  ) |>
  dplyr::select(
    gseid, srrid, variant
  ) |>
  dplyr::distinct() |>
  as.data.table() ->
gseid_srrid_variant

all_variant_cell_table <- dplyr::tbl(conn, "all_variant_cell")

gseid_srrid_variant |>
  # head(100) |>
  dplyr::mutate(
    co = parallel::mcmapply(
      .x = srrid,
      .y = variant,
      FUN = \(.x, .y) {
        # .x <- "GSM4762179"
        # .y <- "11251A>G"

        log_trace(
          glue::glue(
            "Processing variant {.y} for srrid {.x}"
          )
        )
        all_variant_cell_table |>
          dplyr::filter(
            variant == .y,
            srrid == .x,
            variant_in_cell_cluster == "cell"
          ) |>
          dplyr::select(
            barcode, af, depth, variant_type, celltype
          ) |>
          as.data.table() ->
        .d
        .d |>
          dplyr::group_by(celltype) |>
          dplyr::summarise(
            sum_depth = sum(depth, na.rm = TRUE),
            mean_depth = mean(depth, na.rm = T)
          ) ->
        .dd
        log_trace("has data in database ", nrow(.d))
        .d |>
          dplyr::count(
            celltype, variant_type
          ) |>
          dplyr::left_join(
            .dd,
            by = "celltype"
          )
      },
      mc.cores = 20,
      SIMPLIFY = FALSE
    )
  ) ->
gseid_srrid_variant_co


gseid_srrid_variant_co |>
  dplyr::count(gseid, srrid) |>
  dplyr::arrange(-n)

gseid_srrid_variant_co |>
  tidyr::unnest(cols = co) |>
  tidyr::pivot_wider(
    names_from = variant_type,
    values_from = c(n),
  ) |>
  tidyr::nest(
    .by = c(gseid, srrid, variant),
    .key = "variant_celltype"
  ) ->
gseid_srrid_variant_celltype

# gseid_srrid_variant_celltype$variant_celltype[[1]] -> .x


gseid_srrid_variant_celltype |>
  dplyr::mutate(
    n_colorful = parallel::mcmapply(
      .x = variant_celltype,
      FUN = \(.x) {
        .x |>
          tidyr::pivot_longer(
            cols = -c(celltype, sum_depth, mean_depth),
            names_to = "group",
            values_to = "n",
          ) |>
          dplyr::mutate(
            n = ifelse(
              n >= 4,
              n,
              NA_real_
            )
          ) |>
          dplyr::filter(
            !is.na(n)
          ) |>
          dplyr::count(group) |>
          tidyr::pivot_wider(
            names_from = group,
            values_from = n,
            names_prefix = "n_"
          )
      },
      mc.cores = 20,
      SIMPLIFY = FALSE
    )
  ) |>
  tidyr::unnest(n_colorful) ->
gseid_srrid_variant_celltype_n


gseid_srrid_variant_celltype_n |>
  dplyr::filter(
    !is.na(n_black),
    n_black == 8,
    n_colorful < 2
  ) |>
  # dplyr::slice(6) |>
  dplyr::filter(
    srrid == "GSM7080031"
  ) |>
  tidyr::unnest(cols = variant_celltype)



# ? real somatic mutation --------------------------------------------------------------------


ALLVARIANTS <- import(file.path(
  "/home/liuc9/github/scMOCHA-data/analysis/zzz/clean-data/", "all_variant.qs"
)) |>
  dplyr::filter(
    issomatic == "heteroplasmic"
  )

META <- import("/home/liuc9/github/scMOCHA-data/analysis/zzz/clean-data/gse_dataset_metadata_full.qs") |>
  dplyr::select(gseid, srrid, Age_new, Age_group)


gseid_srrid_variant_celltype_n |>
  dplyr::filter(
    # srrid == "GSM7080031"
    variant %in% ALLVARIANTS$variant
  ) |>
  dplyr::filter(
    n_black >= 7,
    n_colorful < 6
  ) ->
somatic_variants

export(
  somatic_variants,
  file = file.path(
    "/home/liuc9/github/scMOCHA-data/analysis/zzz/clean-data/",
    "real_somatic_variant_celltype.qs"
  )
)

somatic_variants <- import(
  file.path(
    "/home/liuc9/github/scMOCHA-data/analysis/zzz/clean-data/",
    "real_somatic_variant_celltype.qs"
  )
)


# ? for plot --------------------------------------------------------------------

META |>
  dplyr::left_join(
    somatic_variants |>
      dplyr::count(gseid, srrid)
  ) |>
  dplyr::filter(
    Age_group != "Unknown"
  ) |>
  dplyr::mutate(
    n = ifelse(
      is.na(n),
      0,
      n
    )
  ) |>
  ggpubr::ggscatter(
    x = "Age_new", y = "n",
    color = "black", shape = 20, size = 3, # Points color, shape and size
    add = "loess", # Add regressin line
    add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
    conf.int = TRUE, # Add confidence interval
    cor.coef = TRUE, # Add correlation coefficient. see ?stat_cor
    cor.coeff.args = list(method = "pearson", label.x = 3, label.sep = "\n"),
    jitter = 0.2,
    xlab = "Age (years)", ylab = "Number of somatic variants",
  ) ->
p_real_somatic_variants_age

p_real_somatic_variants_age |>
  ggplot2::ggsave(
    filename = file.path(
      "/home/liuc9/github/scMOCHA-data/analysis/zzz/plot-real-somatic-variant",
      "real_somatic_variants_age.pdf"
    ),
    width = 8, height = 6
  )




# ? somatic mutation accumulation --------------------------------------------------------------------

source("analysis/00-colors.R")
gseid_srrid_variant_celltype_n |>
  dplyr::filter(
    # srrid == "GSM7080031"
    variant %in% ALLVARIANTS$variant
  ) |>
  dplyr::filter(
    n_black >= 6,
    n_colorful <= 2
  ) |>
  # dplyr::slice(1) |>
  tidyr::unnest(cols = variant_celltype) |>
  dplyr::mutate(
    issomatic = ifelse(
      colorful >= 4,
      "yes",
      "no"
    )
  ) |>
  dplyr::select(gseid, srrid, variant, celltype, issomatic) |>
  dplyr::filter(issomatic == "yes") |>
  dplyr::left_join(
    META,
    by = c("gseid", "srrid")
  ) |>
  dplyr::filter(
    Age_group != "Unknown"
  ) |>
  dplyr::mutate(
    Age_new_group = dplyr::case_when(
      Age_new < 30 ~ "30<",
      Age_new >= 30 & Age_new < 40 ~ "30~40",
      Age_new >= 40 & Age_new < 50 ~ "40~50",
      Age_new >= 50 & Age_new < 60 ~ "50~60",
      Age_new >= 60 & Age_new < 70 ~ "60~70",
      Age_new >= 70 & Age_new < 80 ~ "70~80",
      Age_new >= 80 ~ ">=80",
    )
  ) |>
  dplyr::mutate(
    Age_new_group = factor(
      Age_new_group,
      levels = c("30<", "30~40", "40~50", "50~60", "60~70", "70~80", ">=80")
    )
  ) ->
somatic_variants_age_group

somatic_variants_age_group |>
  dplyr::select(
    Age_new_group, celltype, variant
  ) |>
  dplyr::distinct() |>
  dplyr::count(
    Age_new_group,
    celltype
  ) |>
  dplyr::mutate(
    celltype = as.character(celltype),
  ) |>
  dplyr::mutate(
    celltype = gsub(
      pattern = "_",
      replacement = " ",
      x = celltype
    )
  ) |>
  dplyr::mutate(
    celltype = factor(
      celltype,
      levels = names(color_celltype)
    )
  ) ->
forplot_age_group

# forplot_age_group |>
# dplyr::filter(Age_group == "45~50") |>
forplot_age_group |>
  ggplot(aes(
    x = Age_new_group,
    y = n,
    fill = celltype
  )) +
  geom_col() +
  scale_fill_manual(
    name = "Cell Type",
    values = color_celltype,
  ) +
  scale_y_continuous(
    expand = expansion(mult = c(0, 0.01))
  ) +
  theme(
    panel.background = element_blank(),
    panel.grid = element_blank(),
    # axis.ticks = element_blank(),
    axis.text.x = element_text(size = 14, color = "black"),
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = 16, color = "black"),
    axis.text.y = element_text(size = 14, color = "black"),
    legend.position = "right",
    legend.key = element_blank(),
    legend.text = element_text(size = 14, color = "black"),
    legend.title = element_text(size = 16, colour = "black"),
    plot.margin = margin(t = 0, b = 0, unit = "cm"),
    axis.line = element_line(color = "black"),
  ) +
  labs(
    y = "Number of mutation",
  ) ->
p_somatic_variants_age_group
p_somatic_variants_age_group |>
  ggplot2::ggsave(
    filename = file.path(
      "/home/liuc9/github/scMOCHA-data/analysis/zzz/plot-real-somatic-variant",
      "real_somatic_variants_age_group.pdf"
    ),
    width = 8, height = 6
  )


# ? somatic variant example --------------------------------------------------------------------

somatic_variants |>
  dplyr::filter(
    n_colorful <= 2
  ) |>
  dplyr::slice(2) |>
  tidyr::unnest(cols = variant_celltype)

# thevariant <- "7757G>A"
# thesrrid <- "GSM7437874"

thevariant <- "6967G>A"
thesrrid <- "GSM7080026"
fn_plot_somatic_variant <- function(thevariant, thesrrid) {
  all_variant_cell_table |>
    dplyr::filter(
      srrid == thesrrid,
      variant == thevariant
    ) |>
    dplyr::collect() |>
    dplyr::mutate(
      variant_type = dplyr::case_match(
        variant_type,
        "colorful" ~ "red",
        "black" ~ "darkblue",
        "white" ~ "white",
        "grey" ~ "gray",
        NA ~ "white"
      )
    ) |>
    dplyr::mutate(
      variant_type = factor(
        variant_type,
        levels = c("red", "darkblue", "gray", "white")
      )
    ) |>
    dplyr::arrange(
      variant_type,
      -af
    ) ->
  forplot_



  forplot_ |>
    dplyr::mutate(
      barcode = factor(
        barcode,
        levels = forplot_$barcode
      )
    ) ->
  forplot
  source("analysis/00-colors.R")

  thetheme <- theme(
    panel.background = element_blank(),
    panel.grid = element_blank(),
    axis.ticks = element_blank(),
    axis.text = element_blank(),
    axis.title.x = element_blank(),
    plot.margin = margin(t = 0, b = 0, unit = "cm"),
  )


  forplot |>
    dplyr::mutate(
      celltype = gsub(
        "_",
        " ",
        celltype
      )
    ) |>
    dplyr::mutate(
      celltype = factor(
        celltype,
        names(color_celltype)
      )
    ) |>
    ggplot(aes(
      x = barcode,
      y = 1,
      fill = celltype
    )) +
    geom_col() +
    scale_fill_manual(
      name = "Cell Type",
      values = color_celltype,
    ) +
    scale_y_continuous(expand = expansion(mult = c(0, 0)), ) +
    thetheme +
    # theme(panel.background = element_rect(color = "red")) +
    labs(
      y = "Cell Type",
    ) ->
  p1_celltype

  forplot |>
    dplyr::mutate(
      af = ifelse(
        af < 0.01,
        NA_real_,
        af
      )
    ) |>
    ggplot(aes(
      x = barcode,
      y = af,
      fill = af
    )) +
    geom_col() +
    scale_fill_gradient(
      name = "Allele Frequency",
      high = "red",
      low = "white"
    ) +
    scale_y_continuous(
      expand = expansion(mult = c(0, 0)),
    ) +
    theme(
      panel.background = element_blank(),
      panel.grid = element_blank(),
      axis.line = element_line(colour = "black"),
      axis.ticks.x = element_blank(),
      axis.text.x = element_blank(),
      axis.title.x = element_blank(),
      plot.margin = margin(t = 0, b = 0, unit = "cm"),
    ) +
    labs(
      y = "Allele Frequency",
    ) ->
  p2_af



  forplot |>
    dplyr::mutate(
      variant_type = as.character(variant_type),
    ) |>
    ggplot(aes(
      x = barcode,
      y = 1,
      fill = variant_type
    )) +
    geom_col() +
    scale_y_continuous(expand = expansion(mult = c(0, 0)), ) +
    scale_fill_identity(
      guide = "legend",
      name = "Variant cell",
      breaks = c("red", "darkblue", "gray", "white"),
      labels = c("Heteroplasmy", "Suficcient reads", "No sufficient reads", "No reads")
    ) +
    thetheme +
    labs(
      y = "Variant cells",
    ) ->
  p3_variant_cells


  forplot |>
    dplyr::mutate(
      depth = log2(depth + 1) # log2 transform to reduce skewness
    ) |>
    ggplot(aes(
      x = barcode,
      y = depth,
      fill = depth
    )) +
    geom_col() +
    scale_fill_gradient(
      name = "log2(depth + 1)",
      high = "gold",
      low = "white"
    ) +
    scale_y_continuous(
      expand = expansion(mult = c(0, 0)),
    ) +
    theme(
      panel.background = element_blank(),
      panel.grid = element_blank(),
      axis.line = element_line(colour = "black"),
      axis.ticks.x = element_blank(),
      axis.text.x = element_blank(),
      axis.title.x = element_blank(),
      plot.margin = margin(t = 0, b = 0, unit = "cm"),
    ) +
    labs(
      y = "Log2(Depth + 1)",
    ) ->
  p4_depth

  wrap_plots(
    p2_af,
    plot_spacer(),
    p4_depth,
    plot_spacer(),
    p3_variant_cells,
    plot_spacer(),
    p1_celltype,
    ncol = 1,
    heights = c(15, -1.05, 15, -1.05, 10, -1.05, 10),
    guides = "collect"
  ) +
    plot_annotation(
      title =
        glue::glue(
          "Variant {thevariant} in {thesrrid}"
        ),
      theme = theme(
        plot.title = element_text(hjust = 0.5, size = 16, face = "bold")
      )
    ) ->
  p_all
  # p_all

  ggsave(
    p_all,
    filename = file.path(
      "/home/liuc9/github/scMOCHA-data/analysis/zzz/plot-real-somatic-variant",
      glue::glue("somatic_variant_{thevariant}_{thesrrid}.pdf")
    ),
    width = 13, height = 8
  )
}


fn_plot_somatic_variant("6967G>A", "GSM7080026")
fn_plot_somatic_variant("7757G>A", "GSM7437874")
fn_plot_somatic_variant("12501G>A", "GSM5227130")

srrid_variant_pairs <- data.frame(
  srrid = c("GSM5227130", "GSM4905211", "GSM4905214", "GSM5494119", "GSM7493839", "GSM4670210", "GSM4670211", "GSM7080026", "GSM7437874"),
  variant = c("12501G>A", "1314C>T", "1314C>T", "13271T>C", "13271T>C", "14530T>C", "14530T>C", "6967G>A", "7757G>A")
)
srrid_variant_pairs |>
  dplyr::mutate(
    a = parallel::mcmapply(
      thesrrid = srrid,
      thevariant = variant,
      FUN = fn_plot_somatic_variant,
      mc.cores = 5,
      SIMPLIFY = FALSE
    )
  )

# ? somatic variants hotspot ----------------------------------------------------

fn_plot_mtdna <- function() {
  # mt_exons_df <- "/home/liuc9/github/scMOCHA/fasta/mt_exons.df.rds.gz"

  LENGTH <- 16569
  # rCRS <- Biostrings::readDNAStringSet("/home/liuc9/github/scMOCHA-data/config/rCRS.MT.fasta")
  gtf_gene_df <- import("/home/liuc9/github/scMOCHA-data/config/mtdna_genes_dloop.qs")


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

somatic_variants |>
  dplyr::count(variant) |>
  dplyr::arrange(-n)

somatic_variants |>
  dplyr::filter(
    variant == "14530T>C"
  ) |>
  dplyr::slice(1) |>
  tidyr::unnest(cols = variant_celltype)

fn_plot_somatic_variant("14530T>C", "GSM4670211")

somatic_variants |>
  dplyr::count(variant) |>
  dplyr::arrange(-n) |>
  dplyr::mutate(
    pos = gsub(
      "[ATGC]+|>",
      "",
      variant
    ) |>
      as.integer()
  ) ->
forplot_somatic_variant_hotspot
forplot_somatic_variant_hotspot |>
  ggplot(aes(
    x = pos,
    y = n
  )) +
  geom_segment(
    aes(x = pos, xend = pos, y = 0, yend = n)
  ) +
  geom_point(
    aes(size = n),
    color = "red",
    fill = "red",
    # alpha = 0.7,
    shape = 21,
    stroke = 1
  ) +
  ggrepel::geom_text_repel(
    data = forplot_somatic_variant_hotspot |>
      dplyr::filter(
        n > 2
      ),
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
    expand = expansion(add = c(.05, 0.1)),
    limits = c(0, 6),
    breaks = seq(0, 6, 1),
    labels = seq(0, 6, 1)
  ) +
  # scale_fill_identity(
  #   name = "Sample"
  # ) +
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
    axis.title.y = element_text(size = 16, color = "black"),
    axis.text.y = element_text(color = "black"),
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
    )
  ) +
  labs(
    y = "Number of somatic variants",
  ) ->
p_somatic_variant_hotspot


wrap_plots(
  p_somatic_variant_hotspot,
  plot_spacer(),
  fn_plot_mtdna(),
  ncol = 1,
  heights = c(15, -0.7, 1)
)

ggsave(
  plot = wrap_plots(
    p_somatic_variant_hotspot,
    plot_spacer(),
    fn_plot_mtdna(),
    ncol = 1,
    heights = c(15, -0.7, 1)
  ),
  filename = file.path(
    "/home/liuc9/github/scMOCHA-data/analysis/zzz/plot-real-somatic-variant",
    "somatic_variant_hotspot.pdf"
  ),
  width = 12, height = 6
)
# footer ------------------------------------------------------------------

# future: :plan(future: :sequential)

# save image --------------------------------------------------------------
