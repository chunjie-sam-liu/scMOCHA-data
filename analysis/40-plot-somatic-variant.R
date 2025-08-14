#!/usr/bin/env Rscript
# Metainfo ----------------------------------------------------------------
# @AUTHOR: Chun-Jie Liu
# @CONTACT: chunjie.sam.liu.at.gmail.com
# @DATE: 2025-08-04 17:14:17
# @DESCRIPTION: filename
# @VERSION: v0.0.1

# Library -----------------------------------------------------------------

suppressPackageStartupMessages(library(magrittr))
library(ggplot2)
library(patchwork)
library(prismatic)
library(paletteer)
library(data.table)
#library(rlang)
library(glue)
library(parallel)
library(GetoptLong)
library(logger)
library(future.apply)

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

# header ------------------------------------------------------------------

# future: :plan(future: :multisession, workers = 10)

# load data ---------------------------------------------------------------

# load conn ---------------------------------------------------------------

conn_all_variant_cell <- DBI::dbConnect(
  duckdb::duckdb(),
  dbdir = "/home/liuc9/github/scMOCHA-data/analysis/zzz/clean-data/all_variant_cell.duckdb.1.2.1"
)
DBI::dbListTables(conn_all_variant_cell)

conn_all_hetero_af <- DBI::dbConnect(
  duckdb::duckdb(),
  "/home/liuc9/github/scMOCHA-data/analysis/zzz/clean-data/all_hetero_af.cell.duckdb.1.2.1"
)
DBI::dbListTables(conn_all_hetero_af)

# conn_all_hetero_af |>
#   dplyr::tbl("gseid_srrid_variant")

tbl_allvariants <- conn_all_hetero_af |>
  dplyr::tbl("allvariants")


tbl_gseid_srrid_variant <- conn_all_hetero_af |>
  dplyr::tbl("gseid_srrid_variant")

tbl_gseid_srrid_variant |>
  dplyr::collect() |>
  dplyr::mutate(
    a = purrr::map(
      .x = variant_alltype,
      ~ {
        jsonlite::fromJSON(.x) |>
          purrr::pluck("heteroplasmic_variant") -> .v
        if (length(.v) == 0) {
          return(NULL)
        } else {
          return(tibble::tibble(variant = .v))
        }
      }
    )
  ) |>
  dplyr::select(-variant_alltype) |>
  tidyr::unnest(cols = c(a)) -> gseid_srrid_variant_hetero

# src ---------------------------------------------------------------------

# function ----------------------------------------------------------------

fn_forplot <- function(thevariant, thesrrid) {
  source("analysis/00-colors.R")

  colorcode <- setNames(names(color_variantcell), color_variantcell)

  dplyr::tbl(
    conn_all_variant_cell,
    "all_variant_cell"
  ) |>
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
        levels = color_variantcell
      )
    ) |>
    dplyr::arrange(
      variant_type,
      -af
    ) -> forplot_

  forplot_ |>
    dplyr::mutate(
      barcode = factor(
        barcode,
        levels = forplot_$barcode
      )
    ) |>
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
    dplyr::mutate(
      af = ifelse(
        af < 0.01,
        NA_real_,
        af
      )
    ) |>
    # dplyr::mutate(
    #   variant_type = as.character(variant_type),
    # ) |>
    dplyr::mutate(
      depth = log2(depth + 1) # log2 transform to reduce skewness
    ) |>
    dplyr::mutate(
      cellvarianttype = colorcode[variant_type]
    ) |>
    dplyr::mutate(
      cellvarianttype = factor(
        cellvarianttype,
        levels = colorcode
      )
    ) -> forplot
  forplot
}
fn_plot_somatic_variant <- function(forplot_) {
  source("analysis/00-colors.R")

  colorcode <- setNames(names(color_variantcell), color_variantcell)
  thetheme <- theme(
    panel.background = element_blank(),
    panel.grid = element_blank(),
    axis.ticks = element_blank(),
    axis.text = element_blank(),
    axis.title.x = element_blank(),
    plot.margin = margin(t = 0, b = 0, unit = "cm"),
  )

  forplot_ |>
    dplyr::mutate(
      barcode = factor(
        barcode,
        levels = forplot_$barcode
      )
    ) |>
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
    dplyr::mutate(
      af = ifelse(
        af < 0.01,
        NA_real_,
        af
      )
    ) |>
    # dplyr::mutate(
    #   variant_type = as.character(variant_type),
    # ) |>
    dplyr::mutate(
      depth = log2(depth + 1) # log2 transform to reduce skewness
    ) |>
    dplyr::mutate(
      cellvarianttype = colorcode[variant_type]
    ) |>
    dplyr::mutate(
      cellvarianttype = factor(
        cellvarianttype,
        levels = colorcode
      )
    ) -> forplot

  forplot |>
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
    ) -> p1_celltype
  forplot |>
    ggplot(aes(
      x = barcode,
      y = af,
      fill = af
    )) +
    geom_col() +
    scale_fill_gradient2(
      name = "Allele Frequency",
      high = "#FDE725FF",
      mid = "#21908CFF",
      low = "#440154FF"
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
    ) -> p2_af

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
      labels = c(
        "Heteroplasmy",
        "Sufficient reads",
        "No sufficient reads",
        "No reads"
      )
    ) +
    thetheme +
    labs(
      y = "Variant cells",
    ) -> p3_variant_cells

  forplot |>
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
    ) -> p4_depth

  .gseid <- unique(forplot$gseid)
  .srrid <- unique(forplot$srrid)
  .variant <- unique(forplot$variant)

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
      title = glue::glue(
        "Variant {.variant} in {.gseid}-{.srrid}"
      ),
      theme = theme(
        plot.title = element_text(hjust = 0.5, size = 16, face = "bold")
      )
    ) -> p_all

  p_all
}
fn_cellvarianttype <- function(forplot) {
  source("analysis/00-colors.R")
  colorcode <- setNames(names(color_variantcell), color_variantcell)

  forplot |>
    dplyr::mutate(
      color = as.character(variant_type),
    ) |>
    dplyr::count(color) |>
    dplyr::mutate(
      varianttype = colorcode[color]
    ) |>
    dplyr::mutate(
      varianttype = factor(
        varianttype,
        levels = colorcode
      )
    ) |>
    dplyr::arrange(varianttype) -> cellvarianttype
}
fn_cellvarianttype_ratio <- function(.x) {
  .x |>
    dplyr::select(n, varianttype) |>
    tidyr::pivot_wider(
      names_from = varianttype,
      values_from = n,
      values_fill = 0,
      names_prefix = "count_"
    ) |>
    dplyr::mutate(
      count_total = sum(dplyr::across(starts_with("count_")))
    ) -> .xx
  .x |>
    dplyr::mutate(
      ratio = n / sum(n)
    ) |>
    dplyr::select(
      varianttype,
      ratio
    ) |>
    tidyr::pivot_wider(
      names_from = varianttype,
      values_from = ratio,
      values_fill = 0,
      names_prefix = "ratio_"
    ) -> .xxx
  .xx |> dplyr::bind_cols(.xxx)
}
fn_plot <- function(thevariant, thesrrid) {
  #
  forplot <- fn_forplot(thevariant, thesrrid)
  p_all <- fn_plot_somatic_variant(forplot)
  cellvarianttype <- fn_cellvarianttype(forplot)

  tibble::tibble(
    plot = list(p_all),
    cellvarianttype = list(cellvarianttype),
  )
}

fn_plot_variant_ratio <- function(.d) {
  source("analysis/00-colors.R")
  colorcode <- setNames(names(color_variantcell), color_variantcell)

  .n_gse <- unique(.d$gseid) |> length()
  .n_srr <- unique(.d$srrid) |> length()
  .n_cells <- sum(.d$count)
  .variant <- unique(.d$variant)
  # scales::label_comma()(.n_srr)

  .d |>
    dplyr::mutate(
      varianttype = factor(
        varianttype,
        levels = names(color_variantcell) |>
          rev()
      )
    ) |>
    dplyr::group_by(
      srrid,
    ) |>
    dplyr::mutate(
      ratio = count / sum(count)
    ) |>
    dplyr::ungroup() -> .d_forplot

  .d_forplot |>
    dplyr::filter(
      varianttype == "Heteroplasmy"
    ) |>
    dplyr::arrange(ratio) |>
    dplyr::pull(srrid) -> rank_srrid

  .d_forplot |>
    dplyr::group_by(srrid) |>
    dplyr::summarise(
      count = sum(count),
      .groups = "drop"
    ) |>
    dplyr::pull(count) |>
    quantile(probs = seq(0, 1, 0.1)) -> count_quantiles

  .d_forplot |>
    dplyr::mutate(
      srrid = factor(
        srrid,
        levels = rank_srrid
      )
    ) |>
    ggplot(aes(
      x = srrid,
      y = count,
    )) +
    geom_col(
      aes(
        fill = varianttype
      ),
      position = "stack"
    ) +
    scale_fill_manual(
      name = "Variant cell",
      values = color_variantcell
    ) +
    scale_x_discrete(
      limits = rank_srrid,
    ) +
    scale_y_continuous(
      limits = c(0, count_quantiles["100%"]),
      expand = expansion(mult = c(0.005, 0.03)),
      labels = scales::label_comma()
    ) +
    theme(
      # panel.background = element_blank(),
      panel.grid = element_blank(),
      # axis.ticks = element_blank(),
      axis.line = element_line(color = "black"),
      axis.text.x = element_text(angle = 45, hjust = 1),
      axis.title.x = element_blank(),
      plot.margin = margin(t = 0, b = 0, unit = "cm"),
    ) +
    labs(y = "# of cell") -> p_count

  .d_forplot |>
    dplyr::mutate(
      srrid = factor(
        srrid,
        levels = rank_srrid
      )
    ) |>
    ggplot(aes(
      x = srrid,
      y = ratio,
    )) +
    geom_col(
      aes(
        fill = varianttype
      ),
      position = "stack"
    ) +
    # scale_fill_manual(
    #   name = "Cell Type",
    #   values = RColorBrewer::brewer.pal(8, "Set2")
    # ) +
    scale_fill_manual(
      name = "Variant cell",
      values = color_variantcell
    ) +
    scale_y_continuous(
      expand = expansion(add = c(0.005, 0.01)),
      labels = scales::percent_format(accuracy = 1)
    ) +
    scale_x_discrete(
      limits = rank_srrid,
    ) +
    theme(
      panel.grid = element_blank(),
      # axis.ticks = element_blank(),
      # axis.text.x = element_text(angle = 45, hjust = 1),
      axis.text.x = element_blank(),
      axis.title.x = element_blank(),
      plot.margin = margin(t = 0, b = 0, unit = "cm"),
      axis.line = element_line(color = "black"),
    ) +
    labs(y = "Cell ratio") -> p_ratio

  wrap_plots(
    p_ratio,
    p_count,
    ncol = 1,
    heights = c(1, 1),
    guides = "collect"
  ) +
    plot_annotation(
      title = "{.variant} in {.n_gse} projects and {scales::label_comma()(.n_srr)} samples and {scales::label_comma()(.n_cells)} cells" |>
        glue::glue(),
      theme = theme(
        plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
        # legend.position = "bottom",
        # legend.direction = "horizontal",
      )
    )
}


# body --------------------------------------------------------------------
thevariant <- "3173G>A"
thesrrid <- "GSM7080053"
thevariants <- c(
  "3173G>A",
  "3176A>T",
  "3178T>A",
  "3727T>C",
  "3728C>T"
)

gseid_srrid_variant_hetero |>
  dplyr::filter(variant %in% thevariants) |>
  dplyr::mutate(
    p = parallel::mcmapply(
      FUN = fn_plot,
      thevariant = variant,
      thesrrid = srrid,
      SIMPLIFY = FALSE,
      mc.cores = 20
    )
  ) -> gseid_srrid_variant_hetero_plot

# save image
{
  gseid_srrid_variant_hetero_plot |>
    tidyr::unnest(cols = c(p)) |>
    dplyr::mutate(
      # save image to /home/liuc9/github/scMOCHA-data/analysis/zzz/plot-real-somatic-variant/main-variants
      a = parallel::mcmapply(
        FUN = \(.gseid, .srrid, .variant, .p) {
          .filename = "{.gseid}_{.srrid}_{.variant}.pdf" |> glue::glue()
          ggsave(
            plot = .p,
            filename = .filename,
            path = "/home/liuc9/github/scMOCHA-data/analysis/zzz/plot-real-somatic-variant/main-variants/individual_variant",
            width = 13,
            height = 8
          )
        },
        .gseid = gseid,
        .srrid = srrid,
        .variant = variant,
        .p = plot,
        SIMPLIFY = FALSE,
        mc.cores = 20
      )
    ) -> gseid_srrid_variant_hetero_plot_save
}


gseid_srrid_variant_hetero_plot |>
  tidyr::unnest(cols = c(p)) |>
  dplyr::mutate(
    ratio = purrr::map(
      .x = cellvarianttype,
      fn_cellvarianttype_ratio
    )
  ) |>
  tidyr::unnest(cols = c(ratio)) -> gseid_srrid_variant_hetero_plot_ratio

gseid_srrid_variant_hetero_plot_ratio |>
  dplyr::select(-c(plot, cellvarianttype)) |>
  dplyr::group_by(variant) |>
  dplyr::slice_max(order_by = ratio_Heteroplasmy, n = 1) |>
  dplyr::ungroup() |>
  dplyr::select(gseid, srrid, variant, ratio_Heteroplasmy, count_total)


gseid_srrid_variant_hetero_plot_ratio |>
  dplyr::select(-c(plot, cellvarianttype)) |>
  dplyr::filter(variant == "3727T>C", srrid == "GSM7493836") |>
  dplyr::glimpse()


gseid_srrid_variant_hetero_plot_ratio |>
  dplyr::select(-c(plot, cellvarianttype)) |>
  dplyr::select(gseid, srrid, variant, dplyr::contains("count_")) |>
  tidyr::pivot_longer(
    cols = dplyr::contains("count_"),
    names_to = "varianttype",
    values_to = "count",
    names_prefix = "count_"
  ) |>
  dplyr::filter(varianttype != "total") |>
  tidyr::nest(
    .by = variant,
    .key = "ratio"
  ) |>
  dplyr::mutate(
    p = purrr::map2(
      .x = variant,
      .y = ratio,
      ~ {
        .y |> dplyr::mutate(variant = .x) -> .d
        fn_plot_variant_ratio(.d) -> .p
        .p
      }
    )
  ) -> gseid_srrid_variant_hetero_plot_ratio_plot


# save plots
gseid_srrid_variant_hetero_plot_ratio_plot |>
  dplyr::mutate(
    a = purrr::map2(
      .x = variant,
      .y = p,
      ~ {
        .dir <- glue::glue(
          "/home/liuc9/github/scMOCHA-data/analysis/zzz/plot-real-somatic-variant/main-variants/{.x}"
        )
        if (dir.exists(.dir) == FALSE) {
          dir.create(.dir, recursive = TRUE)
        }
        ggsave(
          plot = .y,
          filename = glue::glue("{.x}-individual-cell-variant-proportion.pdf"),
          path = .dir,
          width = 10,
          height = 6
        )
      }
    )
  )

# footer ------------------------------------------------------------------

# future: :plan(future: :sequential)

# save image --------------------------------------------------------------
