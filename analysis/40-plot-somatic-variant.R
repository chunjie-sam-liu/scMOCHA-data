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

# conn_all_variant_cell <- DBI::dbConnect(
#   duckdb::duckdb(),
#   dbdir = "/home/liuc9/github/scMOCHA-data/analysis/zzz/clean-data/all_variant_cell.duckdb.1.2.1"
# )
# DBI::dbListTables(conn_all_variant_cell)

conn_all_hetero_af <- DBI::dbConnect(
  duckdb::duckdb(),
  "/home/liuc9/github/scMOCHA-data/analysis/zzz/clean-data/all_hetero_af.cell.duckdb.1.2.1"
)
DBI::dbListTables(conn_all_hetero_af)


tbl_allvariants <- conn_all_hetero_af |>
  dplyr::tbl("allvariants")


tbl_gseid_srrid_variant <- conn_all_hetero_af |>
  dplyr::tbl("gseid_srrid_variant")

conn_all_hetero_af |>
  dplyr::tbl("gseid_srrid_variant")


tbl_gseid_srrid_variant |>
  dplyr::collect() |>
  dplyr::mutate(
    a = purrr::map(
      .x = variant_alltype,
      ~ {
        # a$variant_alltype[[1]] -> .x
        jsonlite::fromJSON(.x) -> .j
        purrr::pluck(.j, "heteroplasmic_variant") |> unlist() -> .v_hete
        purrr::pluck(.j, "homoplasmic_variant") |> unlist() -> .v_homo

        tibble::tibble(
          variant = c(.v_hete, .v_homo),
          variant_type = c(
            rep("heteroplasmic", length(.v_hete)),
            rep("homoplasmic", length(.v_homo))
          )
        ) -> .d

        if (length(.v_hete) == 0 & length(.v_homo) == 0) {
          return(NULL)
        } else {
          return(.d)
        }
      }
    )
  ) |>
  dplyr::select(-variant_alltype) |>
  tidyr::unnest(cols = c(a)) -> gseid_srrid_variant_hetero

DBI::dbListTables(conn_all_hetero_af)

tbl_all_hetero_af_bulk <- dplyr::tbl(
  conn_all_hetero_af,
  # "all_hetero_af_bulk"
  "allvariants_af_bulk"
)
tbl_all_hetero_af_cluster <- dplyr::tbl(
  conn_all_hetero_af,
  # "all_hetero_af_cluster"
  "allvariants_af_cluster"
)
tbl_all_hetero_af_cell <- dplyr::tbl(
  conn_all_hetero_af,
  # "all_hetero_af_cell"
  "allvariants_af_cell"
)

# src ---------------------------------------------------------------------

# function ----------------------------------------------------------------

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

# fn_plot_cell_af_depth
fn_plot_cell_af_depth <- function(thevariant, thesrrid) {
  #
  forplot_ <- fn_plot_cell_af_depth_forplot(thevariant, thesrrid)
  p_all <- fn_plot_cell_af_somatic_variant(forplot_)
  cellvarianttype <- fn_plot_cell_af_cellvarianttype(forplot_)

  tibble::tibble(
    plot = list(p_all),
    cellvarianttype = list(cellvarianttype),
    forplot = list(forplot_)
  )
}

#' Very important function
fn_plot_cell_af_depth_forplot <- function(thevariant, thesrrid) {
  source("analysis/00-colors.R")

  colorcode <- setNames(names(color_variantcell), color_variantcell)

  dplyr::tbl(
    conn_all_hetero_af,
    "allvariants_cell"
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

fn_plot_cell_af_depth_forplot_fisher <- function(thevariant, thesrrid) {
  source("analysis/00-colors.R")

  colorcode <- setNames(names(color_variantcell), color_variantcell)

  dplyr::tbl(
    conn_all_hetero_af,
    # "allvariants_cell"
    "allvariants_cell_fishertest"
  ) |>
    dplyr::filter(
      srrid == thesrrid,
      variant == thevariant
    ) |>
    dplyr::collect() -> d

  ref <- gsub("\\d*|>.*", "", thevariant)
  alt <- gsub(".*>", "", thevariant)

  d |>
    dplyr::select(
      variant_type,
      reff = !!sym(paste0(ref, "FO")),
      refr = !!sym(paste0(ref, "RE")),
      altf = !!sym(paste0(alt, "FO")),
      altr = !!sym(paste0(alt, "RE")),
      fisher_test_pvalue,
      alt_strand_ratio
    ) |>
    dplyr::mutate(
      variant_type2 = purrr::pmap_chr(
        list(
          variant_type,
          reff,
          refr,
          altf,
          altr,
          fisher_test_pvalue,
          alt_strand_ratio
        ),
        \(
          variant_type,
          reff,
          refr,
          altf,
          altr,
          fisher_test_pvalue,
          alt_strand_ratio
        ) {
          if (variant_type != "colorful") {
            return(variant_type)
          } else {
            if (fisher_test_pvalue < 0.05) {
              return("black")
            } else {
              if (alt_strand_ratio < 0.1 | alt_strand_ratio > 0.9) {
                return("black")
              } else {
                if (altf >= 2 & altr >= 2) {
                  return("colorful")
                } else {
                  return("black")
                }
              }
            }
          }
        }
      )
    ) -> d_

  d |>
    dplyr::mutate(
      variant_type = d_$variant_type2
    ) |>
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

fn_plot_cell_af_somatic_variant <- function(forplot_) {
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
    dplyr::mutate(
      cellvarianttype = colorcode[variant_type]
    ) |>
    dplyr::mutate(
      cellvarianttype = factor(
        cellvarianttype,
        levels = colorcode
      )
    ) -> forplot

  fn_plot_cell_af_somatic_variant_af(forplot, thetheme) -> p_af
  fn_plot_cell_af_somatic_variant_depth(forplot, thetheme) -> p_depth
  fn_plot_cell_af_somatic_variant_cell(forplot, thetheme) -> p_variant_cells
  fn_plot_cell_af_somatic_variant_celltype(forplot, thetheme) -> p_celltype

  .gseid <- unique(forplot$gseid)
  .srrid <- unique(forplot$srrid)
  .variant <- unique(forplot$variant)

  wrap_plots(
    p_af,
    plot_spacer(),
    p_depth,
    plot_spacer(),
    p_variant_cells,
    plot_spacer(),
    p_celltype,
    ncol = 1,
    heights = c(15, -1.05, 15, -1.05, 10, -1.05, 10),
    guides = "collect"
  ) +
    plot_annotation(
      title = glue::glue(
        "Variant {.variant} in {.gseid}-{.srrid}"
      ),
      theme = theme(
        plot.title = element_text(
          hjust = 0.5,
          size = 16,
          face = "bold"
        )
      )
    ) -> p_all

  p_all
}

fn_plot_cell_af_somatic_variant_fisher <- function(forplot_) {
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
    dplyr::mutate(
      cellvarianttype = colorcode[variant_type]
    ) |>
    dplyr::mutate(
      cellvarianttype = factor(
        cellvarianttype,
        levels = colorcode
      )
    ) -> forplot

  fn_plot_cell_af_somatic_variant_af(forplot, thetheme) -> p_af
  fn_plot_cell_af_somatic_variant_depth(forplot, thetheme) -> p_depth
  fn_plot_cell_af_somatic_variant_cell(forplot, thetheme) -> p_variant_cells
  fn_plot_cell_af_somatic_variant_celltype(forplot, thetheme) -> p_celltype
  fn_plot_cell_af_somatic_variant_forwardreverse(forplot, thetheme) -> p_fr

  .gseid <- unique(forplot$gseid)
  .srrid <- unique(forplot$srrid)
  .variant <- unique(forplot$variant)

  wrap_plots(
    p_af,
    plot_spacer(),
    p_depth,
    plot_spacer(),
    p_fr,
    plot_spacer(),
    p_variant_cells,
    plot_spacer(),
    p_celltype,
    ncol = 1,
    heights = c(15, -1.05, 15, -1.05, 15, -1.05, 10, -1.05, 10),
    guides = "collect"
  ) +
    plot_annotation(
      title = glue::glue(
        "Variant {.variant} in {.gseid}-{.srrid}"
      ),
      theme = theme(
        plot.title = element_text(
          hjust = 0.5,
          size = 16,
          face = "bold"
        )
      )
    ) -> p_all

  p_all
}

fn_plot_cell_af_somatic_variant_celltype <- function(forplot, thetheme) {
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
    )
}
fn_plot_cell_af_somatic_variant_af <- function(forplot, thetheme) {
  fn_xy_breaks_limits(forplot$af, step = 0.2) -> .ybl

  forplot |>
    ggplot(aes(
      x = barcode,
      y = af,
      fill = af
    )) +
    geom_col() +
    geom_hline(
      aes(yintercept = 0.05),
      linetype = 20,
      color = "red"
    ) +
    geom_hline(
      aes(yintercept = 0.1),
      linetype = 21,
      color = "black"
    ) +
    scale_fill_gradient2(
      name = "Allele Frequency",
      high = "#FDE725FF",
      mid = "#21908CFF",
      low = "#440154FF"
    ) +
    scale_y_continuous(
      limits = .ybl$limits,
      breaks = c(.ybl$breaks, 0.05, 0.1) |> unique() |> sort(),
      labels = \(b) {
        dplyr::case_when(
          b == 0.1 ~ "gnomAD cutoff 10%",
          b == 0.05 ~ "our cutoff 5%",
          TRUE ~ scales::percent_format(accuracy = 1)(b)
        )
      },
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
    )
}
fn_plot_cell_af_somatic_variant_depth <- function(forplot, thetheme) {
  fn_xy_breaks_limits(forplot$depth, step = 1) -> .ybl_depth
  forplot |>
    ggplot(aes(
      x = barcode,
      y = depth,
      fill = depth
    )) +
    geom_col() +
    geom_hline(
      aes(yintercept = log2(10 + 1)),
      linetype = 21,
      color = "black"
    ) +
    scale_fill_gradient(
      name = "log2(depth + 1)",
      high = "gold",
      low = "white"
    ) +
    scale_y_continuous(
      limits = .ybl_depth$limits,
      breaks = c(.ybl_depth$breaks, log2(10 + 1)) |> unique() |> sort(),
      labels = \(b) {
        dplyr::case_when(
          b == log2(10 + 1) ~ "cutoff 10",
          TRUE ~ scales::label_number()(b)
        )
      },
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
    )
}

fn_plot_cell_af_somatic_variant_cell <- function(forplot, thetheme) {
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
    scale_y_continuous(
      expand = expansion(mult = c(0, 0)),
    ) +
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
    )
}

fn_plot_cell_af_cellvarianttype <- function(forplot) {
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

# pseudo-bulk
fn_plot_hetero_pseudo_bulk <- function(.d) {
  source("/home/liuc9/github/scMOCHA-data/analysis/00-colors.R")
  color_celltype_bulk <- c(
    "Pseudo-bulk" = "red",
    color_celltype
  )
  .variant <- unique(.d$variant)
  .n_srr <- unique(.d$srrid) |> length()
  .n_gse <- unique(.d$gseid) |> length()

  dplyr::bind_rows(
    tbl_all_hetero_af_bulk |>
      dplyr::filter(variant == .variant) |>
      dplyr::filter(srrid %in% .d$srrid) |>
      dplyr::collect(),
    tbl_all_hetero_af_cluster |>
      dplyr::filter(variant == .variant) |>
      dplyr::filter(srrid %in% .d$srrid) |>
      dplyr::collect()
  ) -> .all_hetero_af_thevariant

  .all_hetero_af_thevariant |>
    dplyr::mutate(
      barcode = gsub(barcode, pattern = "_", replacement = " "),
      barcode = ifelse(barcode == "bulk", "Pseudo-bulk", barcode),
    ) |>
    dplyr::mutate(
      barcode = factor(
        barcode,
        levels = names(color_celltype_bulk)
      ),
    ) -> .forplot

  .forplot |>
    dplyr::filter(barcode == "Pseudo-bulk") |>
    dplyr::arrange(-af) -> .rank_pseudo_bulk

  fn_xy_breaks_limits(.forplot$af, step = 0.1) -> .ybl

  .forplot |>
    dplyr::mutate(
      srrid = factor(srrid, levels = .rank_pseudo_bulk$srrid),
    ) |>
    ggplot(aes(x = srrid, y = af, fill = barcode)) +
    geom_col() +
    geom_hline(
      aes(yintercept = 0.05),
      linetype = 20,
      color = "red"
    ) +
    geom_hline(
      aes(yintercept = 0.1),
      linetype = 21,
      color = "black"
    ) +
    scale_fill_manual(
      name = "Cell type",
      values = color_celltype_bulk
    ) +
    scale_y_continuous(
      name = "Heteroplasmy frequency",
      limits = .ybl$limits,
      breaks = c(.ybl$breaks, 0.05, 0.1) |> unique() |> sort(),
      expand = expansion(mult = c(0.005, 0.03)),
      labels = \(b) {
        dplyr::case_when(
          b == 0.1 ~ "gnomAD cutoff 10%",
          b == 0.05 ~ "our cutoff 5%",
          TRUE ~ scales::percent_format(accuracy = 1)(b)
        )
      },
    ) +
    ggh4x::facet_grid2(
      ~barcode,
      strip = ggh4x::strip_themed(
        background_x = ggh4x::elem_list_rect(
          fill = color_celltype_bulk,
          color = NA
        ),
        text_x = ggh4x::elem_list_text(
          colour = "white",
          face = c("bold")
        )
      ),
      switch = "x",
    ) +
    theme(
      plot.margin = margin(t = 0.2, b = 0.1, l = 0.1, r = 0.2, unit = "cm"),
      panel.grid = element_blank(),
      panel.background = element_blank(),
      axis.text.x = element_blank(),
      axis.title = element_text(size = 12, color = "black", face = "bold"),
      axis.ticks.x = element_blank(),
      axis.line = element_line(color = "black"),
      # legend.position = c(0.2, 0.6),
      legend.position = "none",
      strip.placement = "outside",
      plot.title = element_text(hjust = 0.5, size = 16, face = "bold")
    ) +
    labs(
      x = "Individuals",
      title = glue::glue(
        "Heteroplasmy frequency of {.variant} in {.n_gse} projects and {.n_srr} samples"
      ),
    )
}

# variant ratio
fn_plot_variant_ratio <- function(.d) {
  source("analysis/00-colors.R")
  colorcode <- setNames(names(color_variantcell), color_variantcell)
  color_celltype_bulk <- c(
    "Pseudo-bulk" = "red",
    color_celltype
  )

  .n_gse <- unique(.d$gseid) |> length()
  .n_srr <- unique(.d$srrid) |> length()
  .n_cells <- sum(.d$count)
  .n_cells_hete <- sum(
    .d |>
      dplyr::filter(varianttype == "Heteroplasmy") |>
      dplyr::pull(count),
    na.rm = TRUE
  )
  .variant <- unique(.d$variant)

  .d |>
    dplyr::select(-forplot) |>
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
      ratio = count / sum(count, na.rm = TRUE)
    ) |>
    dplyr::ungroup() -> .dd

  .dd |>
    dplyr::filter(
      varianttype == "Heteroplasmy"
    ) |>
    dplyr::arrange(ratio) |>
    dplyr::pull(srrid) -> rank_srrid

  .dd |>
    dplyr::mutate(
      srrid = factor(
        srrid,
        levels = rank_srrid
      )
    ) -> .d_forplot

  #plot count
  fn_plot_variant_ratio_count(.d_forplot, rank_srrid) -> p_count
  # plot ratio
  fn_plot_variant_ratio_ratio(.d_forplot, rank_srrid) -> p_ratio
  # paf
  fn_plot_variant_ratio_paf(.d, rank_srrid) -> p_haf_list
  p_haf <- p_haf_list$p
  .mean_af <- p_haf_list$mean_af

  wrap_plots(
    p_haf,
    p_ratio,
    p_count,
    ncol = 1,
    heights = c(0.8, 1, 1),
    guides = "collect"
  ) +
    plot_annotation(
      title = "{.variant} in {.n_gse} projects and  {scales::label_comma()(.n_srr)} samples" |>
        glue::glue(),
      subtitle = "{scales::label_percent(accuracy = 0.01)(.n_cells_hete/.n_cells)} ({scales::label_comma()(.n_cells_hete)}/{scales::label_comma()(.n_cells)}) cells with average HAF {scales::label_number(accuracy=0.01)(.mean_af)}" |>
        glue::glue(),
      theme = theme(
        plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
        plot.subtitle = element_text(hjust = 0.5, size = 14, face = "bold"),
      )
    )
}
fn_plot_variant_ratio_count <- function(.d_forplot, rank_srrid) {
  .d_forplot |>
    dplyr::group_by(srrid) |>
    dplyr::summarise(
      count = sum(count),
      .groups = "drop"
    ) -> .m

  fn_xy_breaks_limits(.m$count, step = 2000) -> .count_ybl
  # count
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
      limits = .count_ybl$limits,
      breaks = .count_ybl$breaks,
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
    labs(y = "# of cell")
}
fn_plot_variant_ratio_ratio <- function(.d_forplot, rank_srrid) {
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
    labs(y = "Cell ratio")
}
fn_get_hetero_pseudo_bulk <- function(.d) {
  source("/home/liuc9/github/scMOCHA-data/analysis/00-colors.R")
  color_celltype_bulk <- c(
    "Pseudo-bulk" = "red",
    color_celltype
  )
  .variant <- unique(.d$variant)
  .n_srr <- unique(.d$srrid) |> length()
  .n_gse <- unique(.d$gseid) |> length()

  dplyr::bind_rows(
    tbl_all_hetero_af_bulk |>
      dplyr::filter(variant == .variant) |>
      dplyr::filter(srrid %in% .d$srrid) |>
      dplyr::collect(),
    tbl_all_hetero_af_cluster |>
      dplyr::filter(variant == .variant) |>
      dplyr::filter(srrid %in% .d$srrid) |>
      dplyr::collect()
  ) -> .all_hetero_af_thevariant

  .all_hetero_af_thevariant |>
    dplyr::mutate(
      barcode = gsub(barcode, pattern = "_", replacement = " "),
      barcode = ifelse(barcode == "bulk", "Pseudo-bulk", barcode),
    ) |>
    dplyr::mutate(
      barcode = factor(
        barcode,
        levels = names(color_celltype_bulk)
      ),
    ) -> .forplot
  .forplot
}
fn_plot_variant_ratio_paf <- function(.d, rank_srrid) {
  .d |>
    dplyr::select(gseid, srrid, variant) |>
    dplyr::distinct() |>
    fn_get_hetero_pseudo_bulk() |>
    dplyr::filter(barcode == "Pseudo-bulk") -> .bulk_forplot

  fn_xy_breaks_limits(.bulk_forplot$af, step = 0.1) -> .ybl

  mean(.bulk_forplot$af) -> .mean_af

  .bulk_forplot |>
    dplyr::mutate(
      srrid = factor(
        srrid,
        levels = rank_srrid
      )
    ) |>
    ggplot(aes(x = srrid, y = af, fill = barcode)) +
    geom_col() +
    geom_hline(
      aes(yintercept = 0.05),
      linetype = 20,
      color = "red"
    ) +
    geom_hline(
      aes(yintercept = 0.1),
      linetype = 21,
      color = "black"
    ) +
    scale_fill_manual(
      name = "Cell type",
      # values = color_celltype_bulk
      values = "gold"
    ) +
    scale_y_continuous(
      name = "Pseudo-bulk HAF",
      limits = .ybl$limits,
      breaks = c(.ybl$breaks, 0.05, 0.1) |> unique() |> sort(),
      labels = \(b) {
        dplyr::case_when(
          b == 0.1 ~ "gnomAD cutoff 10%",
          b == 0.05 ~ "our cutoff 5%",
          TRUE ~ scales::percent_format(accuracy = 1)(b)
        )
      },
      expand = expansion(mult = c(0.01, 0.01)),
    ) +
    theme(
      panel.grid = element_blank(),
      # axis.ticks = element_blank(),
      # axis.text.x = element_text(angle = 45, hjust = 1),
      axis.text.x = element_blank(),
      axis.title.x = element_blank(),
      plot.margin = margin(t = 0, b = 0, unit = "cm"),
      axis.line = element_line(color = "black"),
      panel.background = element_blank(),
      legend.position = "none",
    ) -> .p

  list(
    p = .p,
    mean_af = .mean_af
  )
}
fn_plot_variant_ratio_swarm <- function(.d, rank_srrid) {
  # Heteroplasmy
  .d |>
    dplyr::select(forplot) |>
    tidyr::unnest(cols = forplot) |>
    dplyr::filter(cellvarianttype == "Heteroplasmy") -> .dd

  fn_xy_breaks_limits(
    .dd$af,
    step = 0.1
  ) -> .ybl

  mean(.dd$af) -> .mean_af

  .dd |>
    dplyr::group_by(srrid) |>
    dplyr::summarise(
      meanaf = mean(af, na.rm = TRUE),
    ) -> .dd_meanaf

  .dd |>
    dplyr::mutate(
      srrid = factor(
        srrid,
        levels = rank_srrid
      )
    ) |>
    dplyr::left_join(
      .dd_meanaf,
      by = "srrid"
    ) |>
    ggplot(aes(x = srrid)) +
    geom_violin(
      aes(
        y = af,
        fill = meanaf,
      ),
      alpha = 0.5,
      size = 1,
      color = NA,
      show.legend = FALSE
    ) +
    scale_fill_gradient2(
      name = "AF",
      low = "white",
      mid = "red",
      high = "#3B0049",
      midpoint = 0.5,
    ) +
    ggbeeswarm::geom_quasirandom(
      aes(
        y = af,
        color = af
      ),
      size = 1,
      dodge.width = .75,
      alpha = .5,
    ) +
    scale_color_gradient2(
      name = "AF",
      low = "white",
      mid = "red",
      high = "#3B0049",
      midpoint = 0.5,
    ) +
    geom_hline(
      aes(yintercept = 0.05),
      linetype = 20,
      color = "red"
    ) +
    geom_hline(
      aes(yintercept = 0.1),
      linetype = 21,
      color = "black"
    ) +
    scale_y_continuous(
      name = "HAF",
      limits = .ybl$limits,
      breaks = c(.ybl$breaks, 0.05, 0.1) |> unique() |> sort(),
      labels = \(b) {
        dplyr::case_when(
          b == 0.1 ~ "gnomAD cutoff 10%",
          b == 0.05 ~ "our cutoff 5%",
          TRUE ~ scales::percent_format(accuracy = 1)(b)
        )
      },
      expand = expansion(mult = c(0.01, 0.01)),
    ) +
    theme(
      panel.grid = element_blank(),
      # axis.ticks = element_blank(),
      # axis.text.x = element_text(angle = 45, hjust = 1),
      axis.text.x = element_blank(),
      axis.title.x = element_blank(),
      plot.margin = margin(t = 0, b = 0, unit = "cm"),
      axis.line = element_line(color = "black"),
      panel.background = element_blank(),
      legend.position = "none",
    )
}

# body --------------------------------------------------------------------

thevariant <- "3173G>A"
thesrrid <- "GSM7493843"
thevariants <- c(
  "3173G>A",
  "3176A>T",
  "3178T>A",
  "3727T>C",
  "3728C>T",
  "13271T>C",
  "14063T>C",
  "14831G>A",
  "1643A>G",
  "3667T>G",
  "4175G>A",
  "5513G>A",
  "7065G>A",
  "9025G>A",
  "9237G>A",
  "10398A>G"
)

fn_plot_cell_af_depth(
  thevariant = thevariant,
  thesrrid = "GSM7080018"
) -> a
a$plot
# thevariants <- c("10398A>G")
# thevariant <- "10398A>G"
# thevariants <- c(thevariants, disease_variant$variant) |> unique()
gseid_srrid_variant_hetero |>
  # head(20) |>
  dplyr::filter(variant %in% thevariants) |>
  dplyr::mutate(
    p = parallel::mcmapply(
      FUN = fn_plot_cell_af_depth,
      thevariant = variant,
      thesrrid = srrid,
      SIMPLIFY = FALSE,
      mc.cores = 20
    )
  ) -> gseid_srrid_variant_hetero_plot

# save image
gseid_srrid_variant_hetero_plot |> dplyr::count(variant)
#
{
  gseid_srrid_variant_hetero_plot |>
    tidyr::unnest(cols = c(p)) |>
    dplyr::mutate(
      # save image to /home/liuc9/github/scMOCHA-data/analysis/zzz/plot-real-somatic-variant/main-variants
      a = parallel::mcmapply(
        FUN = \(.gseid, .srrid, .variant, .p) {
          .filename = "{.variant}_{.gseid}_{.srrid}.pdf" |> glue::glue()
          .dir <- glue::glue(
            "/home/liuc9/github/scMOCHA-data/analysis/zzz/plot-real-somatic-variant/main-variants/{.variant}/somatic_plot"
          )
          if (dir.exists(.dir) == FALSE) {
            dir.create(.dir, recursive = TRUE)
          }
          log_info(
            "Saving plot for {.variant} in {.gseid}-{.srrid} to {.dir}/{.filename}"
          )
          ggsave(
            plot = .p,
            filename = .filename,
            path = .dir,
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


#
#
# ? bulk pseudo-bulk-variant-proportion--------------------------------------------------------------------

gseid_srrid_variant_hetero |>
  dplyr::filter(variant %in% thevariants) |>
  tidyr::nest(
    .by = "variant",
    .key = "gse_srr"
  ) -> variant_individuals

# plot
variant_individuals |>
  dplyr::mutate(
    p = purrr::map2(
      variant,
      gse_srr,
      ~ {
        # .x <- variant_individuals$variant[[1]]
        # .y <- variant_individuals$gse_srr[[1]]

        .y |> dplyr::mutate(variant = .x) -> .d
        fn_plot_hetero_pseudo_bulk(.d) -> .p
        .p
      }
    )
  ) -> variant_individuals_plot

# save plot

variant_individuals_plot |>
  dplyr::mutate(
    a = parallel::mcmapply(
      FUN = \(.x, .y) {
        .dir <- glue::glue(
          "/home/liuc9/github/scMOCHA-data/analysis/zzz/plot-real-somatic-variant/main-variants/{.x}"
        )
        if (dir.exists(.dir) == FALSE) {
          dir.create(.dir, recursive = TRUE)
        }
        log_info(
          "Saving bulk variant proportion plot for {.x} to {.dir}/{.x}-bulk-variant-proportion.pdf"
        )
        ggsave(
          plot = .y,
          filename = glue::glue("{.x}-bulk-variant-proportion.pdf"),
          path = .dir,
          width = 15,
          height = 8
        )
      },
      variant,
      p,
      mc.cores = 20,
      SIMPLIFY = FALSE
    )
  )

#
#
# ? plot ratio variant individual-cell-variant-proportion--------------------------------------------------------------------

gseid_srrid_variant_hetero_plot |>
  tidyr::unnest(cols = c(p)) |>
  dplyr::mutate(
    ratio = parallel::mclapply(
      X = cellvarianttype,
      FUN = fn_cellvarianttype_ratio,
      mc.cores = 20
    )
  ) |>
  tidyr::unnest(cols = c(ratio)) -> gseid_srrid_variant_hetero_plot_ratio

export(
  gseid_srrid_variant_hetero_plot_ratio,
  "/home/liuc9/github/scMOCHA-data/analysis/zzz/plot-real-somatic-variant/main-variants/gseid_srrid_variant_hetero_plot_ratio.qs"
)


gseid_srrid_variant_hetero_plot_ratio |>
  dplyr::select(-c(plot, cellvarianttype)) |>
  dplyr::select(gseid, srrid, variant, forplot, dplyr::contains("count_")) |>
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
  ) -> gseid_srrid_variant_hetero_plot_ratio_

gseid_srrid_variant_hetero_plot_ratio_ |>
  dplyr::mutate(
    p = parallel::mcmapply(
      FUN = \(.x, .y) {
        .y |>
          dplyr::mutate(variant = .x) |>
          tidyr::replace_na(
            list(
              count = 0
            )
          ) -> .d
        log_info("Plotting variant ratio for variant {.x}")
        fn_plot_variant_ratio(.d) -> .p
        .p
      },
      .x = variant,
      .y = ratio,
      mc.cores = 20,
      SIMPLIFY = FALSE
    )
  ) -> gseid_srrid_variant_hetero_plot_ratio_plot


# save plots
gseid_srrid_variant_hetero_plot_ratio_plot |>
  dplyr::mutate(
    a = parallel::mcmapply(
      FUN = \(.x, .y) {
        .dir <- glue::glue(
          "/home/liuc9/github/scMOCHA-data/analysis/zzz/plot-real-somatic-variant/main-variants/{.x}"
        )
        if (dir.exists(.dir) == FALSE) {
          dir.create(.dir, recursive = TRUE)
        }
        log_info(
          "Saving individual cell variant proportion plot for {.x} to {.dir}/{.x}-individual-cell-variant-proportion.pdf"
        )
        ggsave(
          plot = .y,
          filename = glue::glue("{.x}-individual-cell-variant-proportion.pdf"),
          path = .dir,
          width = 15,
          height = 8
        )
      },
      .x = variant,
      .y = p,
      mc.cores = 20,
      SIMPLIFY = FALSE
    )
  )


#
# ! don't run below
# ? 3727 3728 --------------------------------------------------------------------
#
#

gseid_srrid_variant_hetero_plot_ratio_$ratio[[4]] |>
  dplyr::filter(varianttype == "Heteroplasmy") -> v3728

gseid_srrid_variant_hetero_plot_ratio_$ratio[[5]] |>
  dplyr::filter(varianttype == "Heteroplasmy") -> v3727

shared_individuals <- intersect(
  v3728$srrid,
  v3727$srrid
)

ggvenn::ggvenn(
  data = list(
    "m.3727T>C" = v3727$srrid,
    "m.3728C>T" = v3728$srrid
  ),
  fill_color = ggsci::pal_aaas()(2),
  stroke_color = "white"
) -> venn_3727_3728

ggsave(
  plot = venn_3727_3728,
  filename = "venn_3727_3728.pdf",
  path = "/home/liuc9/github/scMOCHA-data/analysis/zzz/plot-real-somatic-variant/main-variants/",
  width = 6,
  height = 5
)

v3728 |>
  dplyr::filter(srrid %in% shared_individuals) |>
  dplyr::select(forplot) |>
  tidyr::unnest(cols = forplot) |>
  dplyr::filter(
    cellvarianttype == "Heteroplasmy"
  ) |>
  dplyr::mutate(
    barcode_new = glue::glue("{gseid}_{srrid}_{barcode}")
  ) |>
  dplyr::select(
    gseid,
    srrid,
    barcode_new,
    v3728_af = af,
    celltype,
  ) -> v3728_shared

v3727 |>
  dplyr::filter(srrid %in% shared_individuals) |>
  dplyr::select(forplot) |>
  tidyr::unnest(cols = forplot) |>
  dplyr::filter(
    cellvarianttype == "Heteroplasmy"
  ) |>
  dplyr::mutate(
    barcode_new = glue::glue("{gseid}_{srrid}_{barcode}")
  ) |>
  dplyr::select(
    gseid,
    srrid,
    barcode_new,
    v3727_af = af,
    celltype,
  ) -> v3727_shared

v3728_shared |>
  dplyr::select(-c(gseid, srrid)) |>
  dplyr::full_join(
    v3727_shared |> dplyr::select(-c(gseid, srrid)),
    by = c("barcode_new", "celltype"),
    # suffix = c("_3728", "_3727")
  ) |>
  tidyr::replace_na(
    list(
      v3728_af = 0,
      v3727_af = 0
    )
  ) |>
  dplyr::mutate(
    celltype = factor(
      celltype,
      levels = names(color_celltype) |> rev()
    )
  ) -> v3728_3727_shared

v3728_3727_shared |>
  dplyr::filter(v3728_af > 0 & v3727_af > 0) |>
  nrow() -> n_cells_shared


v3728_3727_shared |>
  ggplot(aes(
    x = v3728_af,
    y = v3727_af,
    color = celltype
  )) +
  geom_point() +
  scale_color_manual(
    values = color_celltype,
    na.value = "grey50",
    name = "Celltype"
  ) +
  geom_vline(
    aes(xintercept = 0.05),
    linetype = 20,
    color = "red"
  ) +
  geom_hline(
    aes(yintercept = 0.05),
    linetype = 20,
    color = "red"
  ) +
  scale_x_continuous(
    name = "m.3728C>T",
    limits = c(0, 1),
    breaks = c(seq(0, 1, 0.1), 0.05, 0.1) |> unique() |> sort(),
    labels = \(b) {
      dplyr::case_when(
        b == 0.05 ~ "5%",
        TRUE ~ scales::percent_format(accuracy = 1)(b)
      )
    },
    expand = expansion(mult = c(0.01, 0.01)),
  ) +
  scale_y_continuous(
    name = "m.3727T>C",
    limits = c(0, 1),
    breaks = c(seq(0, 1, 0.1), 0.05, 0.1) |> unique() |> sort(),
    labels = \(b) {
      dplyr::case_when(
        b == 0.05 ~ "our cutoff 5%",
        TRUE ~ scales::percent_format(accuracy = 1)(b)
      )
    },
    expand = expansion(mult = c(0.01, 0.01)),
  ) +
  theme(
    panel.grid = element_blank(),
    # axis.ticks = element_blank(),
    # axis.text.x = element_text(angle = 45, hjust = 1),
    # axis.text.x = element_blank(),
    # axis.title.x = element_blank(),
    axis.title = element_text(
      size = 14,
      color = "black",
      face = "bold"
    ),
    # plot.margin = margin(t = 0, b = 0, unit = "cm"),
    axis.line = element_line(color = "black"),
    panel.background = element_blank(),
    legend.position = c(0.8, 0.7),
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5, size = 14),
  ) +
  labs(
    title = "Variant m.3727T>C and m.3728C>T mutate independently",
    subtitle = "Variant m.3727T>C in {scales::label_comma()(nrow(v3727_shared))} from {scales::label_comma()(length(unique(v3727_shared$gseid)))} projects and {scales::label_comma()(length(unique(v3727_shared$srrid)))} samples\nVariant m.3728C>T in {scales::label_comma()(nrow(v3728_shared))} from {scales::label_comma()(length(unique(v3728_shared$gseid)))} projects and {scales::label_comma()(length(unique(v3728_shared$srrid)))} samples\nOnly {scales::label_comma()(n_cells_shared)} cells both variants satisfied cutoff" |>
      glue::glue()
  ) -> p_3727_3728_corr

ggsave(
  plot = p_3727_3728_corr,
  filename = "p_3727_3728_corr.pdf",
  path = "/home/liuc9/github/scMOCHA-data/analysis/zzz/plot-real-somatic-variant/main-variants/",
  width = 13,
  height = 8
)

# footer ------------------------------------------------------------------

# future: :plan(future: :sequential)

# save image --------------------------------------------------------------
