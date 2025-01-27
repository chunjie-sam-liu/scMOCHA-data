#!/usr/bin/env Rscript
# Metainfo ----------------------------------------------------------------
# @AUTHOR: Chun-Jie Liu
# @CONTACT: chunjie.sam.liu.at.gmail.com
# @DATE: 2025-01-24 15:18:13
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
source("/home/liuc9/github/scMOCHA-data/src/check_position_variant.R")
# header ------------------------------------------------------------------
log_threshold(TRACE)
log_layout(layout_glue_colors)

# future: :plan(future: :multisession, workers = 10)

# function ----------------------------------------------------------------

fn_plot_all <- function(thepath, thevariants = thevariants, outdir = outdir) {
  log_info("Start to plot ", thepath)
  # ! parse --------------------------------------------------------------------
  if (!file.exists(outdir)) {
    dir.create(outdir, showWarnings = FALSE)
  }

  gsmid <- basename(thepath)
  gseid <- basename(dirname(dirname(thepath)))

  theposes <- thevariants |>
    purrr::map(~ gsub(pattern = "[>|AGCT]", "", x = .)) |>
    purrr::map_int(as.integer)

  # load data ---------------------------------------------------------------
  # load sc
  sc <- fn_load_by_path(thepath)
  # load count
  cluster_n_forplot <- fn_load_count(thepath, type = "cluster")


  # vaf cell umap -------------------------------------------------------------

  fn_plot_vaf_featureplot_multi(
    .thevariants = thevariants,
    sc = sc
  ) -> p_vaf_feature

  # p_vaf_feature

  ggsave(
    filename = "{gseid}-{gsmid}-selected_variants_vaf_featureplot.pdf" |> glue::glue(),
    path = outdir,
    plot = p_vaf_feature,
    width = 9,
    height = 4,
  )

  # read count -------------------------------------------------------------------
  fn_plot_count_multi(
    cluster_n_forplot,
    theposes = theposes
  ) -> p_count

  # p_count

  ggsave(
    filename = "{gseid}-{gsmid}-selected_variants_count.pdf" |> glue::glue(),
    path = outdir,
    plot = p_count,
    width = 15,
    height = 5,
  )

  # depth -------------------------------------------------------------------
  p_mtdna <- fn_plot_mtdna()
  p_depth <- fn_plot_coverage(thepath, theposes)

  ggsave(
    filename = "{gseid}-{gsmid}-depth-celltype.pdf" |> glue::glue(),
    path = outdir,
    plot = wrap_plots(
      p_depth$p_mt_depth_celltype,
      p_mtdna,
      ncol = 1,
      heights = c(0.7, 0.1)
    ),
    width = 17,
    height = 9,
  )


  # # hotspots ----------------------------------------------------------------
  # p_hotspots <- fn_plot_hotspots(thepath, thevariants)

  # ggsave(
  #   filename = "{gseid}-{gsmid}-hotspots_final_af_somatic.pdf" |> glue::glue(),
  #   path = outdir,
  #   plot = wrap_plots(
  #     p_hotspots,
  #     p_depth$p_mt_depth_allcell,
  #     p_mtdna,
  #     ncol = 1,
  #     heights = c(1.6, 0.4, 0.1),
  #     axes = "collect_x"
  #   ),
  #   device = "pdf",
  #   width = 24,
  #   height = 12
  # )
  log_success("Finish to plot ", thepath)
}

fn_get_count_all <- function(thepath, thevariants = thevariants, outdir = outdir) {
  log_info("Start to plot ", thepath)
  # ! parse --------------------------------------------------------------------


  gsmid <- basename(thepath)
  gseid <- basename(dirname(dirname(thepath)))

  theposes <- thevariants |>
    purrr::map(~ gsub(pattern = "[>|AGCT]", "", x = .)) |>
    purrr::map_int(as.integer)

  # load data ---------------------------------------------------------------
  # load sc
  sc <- fn_load_by_path(thepath)
  # load count
  cluster_n_forplot <- fn_load_count(thepath, type = "cluster")

  log_success("Finish to plot ", thepath)

  list(
    gseid = gseid,
    gsmid = gsmid,
    thevariants = list(thevariants),
    theposes = list(theposes),
    cluster_n_forplot = list(cluster_n_forplot)
  )
}


# load data ---------------------------------------------------------------

basedir <- "/home/liuc9/github/scMOCHA-data/data"
outdir <- "/home/liuc9/github/scMOCHA-data/data/out_variant_check"

thepaths <- c(
  "/mnt/isilon/u01_project/large-scale/liuc9/raw/GSE226602/final/GSM7080044",
  "/mnt/isilon/u01_project/large-scale/liuc9/raw/GSE163668/final/GSM4995425",
  "/mnt/isilon/u01_project/large-scale/liuc9/raw/GSE163668/final/GSM4995448",
  "/mnt/isilon/u01_project/large-scale/liuc9/raw/GSE279945/final/GSM8583898"
)

sc5p_pe_variant <- c(
  "1255T>C", "1314C>T", "1315G>A", "1380G>T", "1382A>T", "1397T>A", "1670A>G", "2191A>C", "2285T>C", "2289G>T", "2442T>C", "3173G>A", "3176A>T", "3178T>A", "3727T>C",
  "3728C>T", "3734A>G", "7428G>A", "11560A>G", "13752T>G", "13954C>A", "14082C>G", "15666T>C"
)
sc5p_r2_variant <- c(
  "153A>G", "195T>C", "827A>G", "1002C>T", "2352T>C", "3547A>G", "3766T>C", "4820G>A", "4977T>C", "6164C>T", "6473C>T", "8362T>G", "8598T>C", "8730A>G", "9196G>A",
  "9497T>C", "10604T>A", "10819A>G", "11177C>T", "14212T>C", "14905G>A", "15047G>A", "15535C>T", "15747T>C"
)

sample_chem <- tibble::tibble(
  gseid = c("GSE226602", "GSE163668", "GSE163668", "GSE279945"),
  gsmid = c("GSM7080044", "GSM4995425", "GSM4995448", "GSM8583898"),
  thepath = thepaths,
  thevariant = list(sc5p_pe_variant, sc5p_r2_variant, NULL, NULL),
  chemistry = c("SC5P-PE", "SC5P-R2", "SC5P-R2", "SC3Pv3"),
  color = c("red", "blue", "black", "black")
) |>
  dplyr::mutate(
    gsmid_label = glue::glue("{gsmid} ({chemistry})")
  ) |>
  dplyr::mutate(
    gsmid_label = factor(gsmid_label, levels = gsmid_label),
    gsmid = factor(gsmid, levels = gsmid),
    chemistry = factor(chemistry, levels = chemistry |> unique())
  )

pcc <- readr::read_tsv(file = "https://raw.githubusercontent.com/chunjie-sam-liu/chunjie-sam-liu.life/master/public/data/pcc.tsv") |>
  dplyr::arrange(cancer_types)

# body --------------------------------------------------------------------

thevariants <- c(
  sc5p_pe_variant, sc5p_r2_variant
)
theposes <- thevariants |>
  purrr::map(~ gsub(pattern = "[>|AGCT]", "", x = .)) |>
  purrr::map_int(as.integer)



# ! variant distribution --------------------------------------------------------------------
fn_plot_variant_depth_distribution <- function(thepath, thevariants, outdir) {
  if (!file.exists(outdir)) {
    dir.create(outdir, showWarnings = FALSE)
  }

  gsmid <- basename(thepath)
  gseid <- basename(dirname(dirname(thepath)))

  theposes <- thevariants |>
    purrr::map(~ gsub(pattern = "[>|AGCT]", "", x = .)) |>
    purrr::map_int(as.integer)

  sc <- fn_load_by_path(thepath)
  # load count
  cluster_n_forplot <- fn_load_count(thepath, type = "cluster")

  # depth -------------------------------------------------------------------
  p_mtdna <- fn_plot_mtdna()
  p_depth <- fn_plot_coverage(thepath, theposes)

  pp <- wrap_plots(
    p_depth$p_mt_depth_celltype,
    p_mtdna,
    ncol = 1,
    heights = c(0.7, 0.1)
  )

  ggsave(
    filename = "{gseid}-{gsmid}-depth-celltype.pdf" |> glue::glue(),
    path = outdir,
    plot = pp,
    width = 17,
    height = 9,
  )
  tibble::tibble(
    p_depth = list(p_depth),
    p_mtdna = list(p_mtdna)
  )
}

fn_plot_variant_depth_combined <- function(.d, .v) {
  .d |>
    dplyr::select(pos, depth, gsmid_label) |>
    ggplot(aes(x = pos, y = depth, fill = gsmid_label)) +
    geom_bar(stat = "identity", show.legend = FALSE) +
    scale_x_continuous(
      limits = c(0, 17000),
      breaks = seq(0, 17000, 1000),
      labels = seq(0, 17000, 1000),
      expand = expansion(mult = c(0, 0.01)),
    ) +
    # .plot_vline +
    geom_vline(
      xintercept = .v |>
        purrr::map(~ gsub(pattern = "[>|AGCT]", "", x = .)) |>
        purrr::map_int(as.integer),
      linetype = "dashed",
      color = "red",
      linetype = 21
    ) +
    scale_y_continuous(
      expand = c(0.01, 0),
      # limits = c(0, 520000),
      label = scales::label_number()
    ) +
    scale_fill_manual(
      name = "Cell type",
      values = pcc$color
    ) +
    ggh4x::facet_wrap2(
      ~gsmid_label,
      ncol = 1,
      strip.position = "right",
      strip = ggh4x::strip_themed(
        background_y = ggh4x::elem_list_rect(
          fill = pcc$color
        ),
        text_y = ggh4x::elem_list_text(
          colour = "white",
          face = c("bold")
        ),
        by_layer_y = FALSE,
      ),
      scales = "free_y",
    ) +
    theme(
      plot.margin = margin(t = 0, b = 0, unit = "cm"),
      panel.background = element_blank(),
      panel.grid = element_blank(),
      axis.line.y.left = element_line(color = "black"),
      # axis.line.x.bottom = element_line(color = "black"),
      axis.ticks.x = element_blank(),
      axis.text.x = element_blank(),
      axis.line.x = element_blank(),
      axis.title.x = element_blank(),
      legend.position = c(0.8, 0.5),
      legend.key = element_blank(),
      axis.title.y = element_text(color = "black"),
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
    coord_cartesian(xlim = c(0, 17000)) +
    labs(y = "Depth") ->
  .pd

  pp <- wrap_plots(
    .pd,
    .d$p_mtdna[[1]],
    ncol = 1,
    heights = c(0.7, 0.1)
  )

  pp
}

sample_chem |>
  dplyr::mutate(
    p = purrr::map(
      .x = thepath,
      fn_plot_variant_depth_distribution,
      thevariants = sc5p_pe_variant,
      outdir = file.path(
        outdir,
        "sc5p_pe_variant"
      )
    )
  ) |>
  tidyr::unnest(cols = p) |>
  dplyr::mutate(
    cov = purrr::map(
      .x = p_depth,
      .f = ~ {
        .x$p_mt_depth_allcell -> .xx
        ggplot_build(.xx)$plot$data
      }
    ),
  ) |>
  tidyr::unnest(cols = c(cov)) ->
sc5p_pe_variant_plot


ggsave(
  filename = "sc5p_pe_variant_depth_combined.pdf" |> glue::glue(),
  path = outdir,
  plot = fn_plot_variant_depth_combined(sc5p_pe_variant_plot, sc5p_pe_variant),
  width = 17,
  height = 9,
)


sample_chem |>
  dplyr::mutate(
    p = purrr::map(
      .x = thepath,
      fn_plot_variant_depth_distribution,
      thevariants = sc5p_r2_variant,
      outdir = file.path(
        outdir,
        "sc5p_r2_variant"
      )
    )
  ) |>
  tidyr::unnest(cols = p) |>
  dplyr::mutate(
    cov = purrr::map(
      .x = p_depth,
      .f = ~ {
        .x$p_mt_depth_allcell -> .xx
        ggplot_build(.xx)$plot$data
      }
    ),
  ) |>
  tidyr::unnest(cols = c(cov)) ->
sc5p_r2_variant_plot



ggsave(
  filename = "sc5p_r2_variant_depth_combined.pdf" |> glue::glue(),
  path = outdir,
  plot = fn_plot_variant_depth_combined(sc5p_r2_variant_plot, sc5p_r2_variant),
  width = 17,
  height = 9,
)


# ! read count --------------------------------------------------------------------

parallel::mclapply(thepaths, function(path) {
  fn_get_count_all(path, thevariants = thevariants, outdir = outdir)
}, mc.cores = length(thepaths)) ->
load_count

load_count |>
  purrr::map(
    ~ {
      .x |>
        tibble::enframe() |>
        tidyr::spread(key = name, value = value) |>
        tidyr::unnest(cols = c(thevariants, theposes, cluster_n_forplot, gseid, gsmid)) |>
        dplyr::mutate(
          cluster_n_forplot = purrr::map2(
            .x = cluster_n_forplot,
            .y = theposes,
            ~ {
              .x |>
                dplyr::filter(pos %in% .y)
            }
          )
        )
    }
  ) |>
  dplyr::bind_rows() ->
load_count_unnest


load_count_unnest |>
  dplyr::select(3, 2, 1) |>
  tidyr::unnest(cols = c(cluster_n_forplot)) ->
cluster_n_forplot_

gt <- factor(c("A", "G", "C", "T"), levels = c("A", "G", "C", "T"))
posref = cluster_n_forplot_$posref |> unique()
group = cluster_n_forplot_$group |> unique()

posref_df <- data.table::data.table(
  gt = rep(gt, each = length(posref)),
  posref = rep(posref, length(gt)),
  group = rep(group, each = length(gt) * length(posref))
)



tibble::tibble(
  pos = theposes,
  thevariant = thevariants,
) |>
  dplyr::mutate(
    posref = purrr::map_chr(
      .x = thevariants,
      ~ {
        gsub(pattern = ">[AGCT]", "", x = .x)
      }
    )
  ) |>
  dplyr::arrange(pos) |>
  dplyr::mutate(
    variant_group = purrr::map_chr(
      .x = posref,
      ~ {
        gsub(pattern = "[0-9]", "", x = .x)
      }
    )
  ) |>
  dplyr::mutate(
    color = ifelse(
      thevariant %in% sc5p_pe_variant,
      "red",
      ifelse(
        thevariant %in% sc5p_r2_variant,
        "blue",
        "black"
      )
    )
  ) ->
posref_df_rank

c("A", "G", "C", "T") |>
  purrr::map(
    ~ {
      cluster_n_forplot_ |>
        dplyr::left_join(sample_chem, by = c("gseid", "gsmid")) |>
        dplyr::mutate(
          gsmid_label = factor(gsmid_label, levels = sample_chem$gsmid_label)
        ) |>
        dplyr::filter(group == "B") |>
        dplyr::filter(pos %in% (posref_df_rank |>
          dplyr::filter(variant_group == .x) |>
          dplyr::pull(pos)
        )) |>
        dplyr::mutate(
          posref = factor(posref, levels = posref_df_rank |>
            dplyr::filter(variant_group == .x) |>
            dplyr::pull(posref))
        ) |>
        ggplot(aes(x = posref, y = gsmid_label)) +
        geom_tile(aes(fill = ratio)) +
        geom_text(aes(label = label), size = 5) +
        scale_fill_gradient(
          low = "white",
          high = "red",
          na.value = "white"
        ) +
        ggh4x::facet_wrap2(
          ~gt,
          ncol = 1,
          strip.position = "right",
          strip = ggh4x::strip_themed(
            background_y = ggh4x::elem_list_rect(
              fill = pcc$color
            ),
            text_y = ggh4x::elem_list_text(
              colour = "white",
              face = c("bold")
            ),
            by_layer_x = FALSE,
          )
        ) +
        scale_x_discrete(
          expand = expansion(mult = c(0.01, 0.01))
        ) +
        scale_y_discrete(
          expand = expansion(mult = c(0, 0))
        ) +
        theme(
          panel.background = element_rect(
            color = "black",
            fill = NA,
            linewidth = 0.2
          ),
          panel.grid = element_line(colour = "grey", linetype = "dashed"),
          panel.grid.major = element_line(
            colour = "grey",
            linetype = "dashed",
            size = 0.2
          ),
          # axis.ticks = element_blank(),
          axis.title = element_blank(),
          axis.text = element_text(
            color = "black",
            size = 18
          ),
          axis.text.y = element_text(
            color = sample_chem$color,
            size = 14,
          ),
          axis.text.x = element_text(
            color = posref_df_rank |>
              dplyr::filter(variant_group == .x) |>
              dplyr::pull(color),
            size = 14,
          ),
          legend.position = "none ",
          plot.title = element_text(
            size = 16,
            hjust = 0.5
          ),
          strip.background = element_rect(
            fill = NA,
            color = "black",
          ),
          strip.text = element_text(
            color = "black",
            size = 14,
            face = "bold"
          ),
          axis.line = element_line(
            color = "black"
          )
        ) ->
      pv

      ggsave(
        filename = "selected_variants_ratio_{.x}.pdf" |> glue::glue(),
        path = outdir,
        plot = pv,
        width = 26,
        height = 15,
      )
    }
  )


# footer ------------------------------------------------------------------

# future: :plan(future: :sequential)

# save image --------------------------------------------------------------
