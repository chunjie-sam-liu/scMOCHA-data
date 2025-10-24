#!/usr/bin/env Rscript
# Metainfo ----------------------------------------------------------------
# @AUTHOR: Chun-Jie Liu
# @CONTACT: chunjie.sam.liu.at.gmail.com
# @DATE: 2025-06-05 10:53:21
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

# header ------------------------------------------------------------------

# future: :plan(future: :multisession, workers = 10)

# load data ---------------------------------------------------------------

cleandatadir <- "/home/liuc9/github/scMOCHA-data/data/zzz/clean-data"

conn <- DBI::dbConnect(
  duckdb::duckdb(),
  "/home/liuc9/github/scMOCHA-data/analysis/zzz/clean-data/all_hetero_af.cell.duckdb.1.2.1",
  read_only = TRUE
)

DBI::dbListTables(conn)
tbl_allvariants <- dplyr::tbl(conn, "allvariants")
tbl_allvariants <- dplyr::tbl(
  conn,
  "allvariants_fisher"
)

# tbl_allvariants_fisher |>
#   data.table::as.data.table() |> dplyr::count(issomatic)

tbl_all_hetero_af_bulk <- dplyr::tbl(
  conn,
  "all_hetero_af_bulk"
)
DBI::dbListTables(conn)

tbl_gseid_srrid_variant <- dplyr::tbl(
  conn,
  "gseid_srrid_variant_fisher"
)

tbl_gseid_srrid_variant |>
  dplyr::collect() |>
  dplyr::mutate(
    a = purrr::map(
      .x = variant_alltype,
      ~ {
        # .x <- a$variant_alltype[[38]]
        jsonlite::fromJSON(.x) -> .v
        .v[["heteroplasmic_variant"]] |> as.character() -> .hete
        .v[["homoplasmic_variant"]] |> as.character() -> .homo

        if (length(.v) == 0) {
          return(NULL)
        } else {
          return(tibble::tibble(variant = c(.hete, .homo)))
        }
      }
    )
  ) |>
  dplyr::select(-variant_alltype) |>
  tidyr::unnest(cols = c(a)) -> gseid_srrid_variant_hetero

# all_variant_ <- import(
#   file.path(cleandatadir, "all_variant.qs")
# )
all_variant_ <- tbl_allvariants |> as.data.table()
dplyr::left_join(
  gseid_srrid_variant_hetero,
  all_variant_,
  by = "variant"
) -> all_variant

all_variant |>
  dplyr::filter(
    issomatic == "heteroplasmic"
  ) -> v_hete

all_variant |>
  dplyr::filter(
    issomatic == "homoplasmic"
  ) -> v_homo

all_variant |>
  dplyr::filter(
    issomatic %in% c("homoplasmic", "heteroplasmic")
  ) -> v_homo_hete

# function ----------------------------------------------------------------

fn_variant_L_H_strand <- function(.v) {
  # .v <- v_hete
  Ori <- c(
    "pos1" = 210,
    "pos2" = 16172
  )
  Other <- c(211, 16171)
  L_ref <- c("C", "T")
  L_variant <- c(
    "C>A",
    "C>G",
    "C>T",
    "T>A",
    "T>G",
    "T>C"
  )
  H_ref <- c("G", "A")
  H_variant <- c(
    "G>A",
    "G>C",
    "G>T",
    "A>C",
    "A>T",
    "A>G"
  )
  variant_reverse <- c(
    "G>A" = "C>T",
    "G>C" = "C>G",
    "G>T" = "C>A",
    "T>A" = "A>T",
    "T>C" = "A>G",
    "T>G" = "A>C",
    "C>A" = "G>T",
    "C>G" = "G>C",
    "C>T" = "G>A",
    "A>C" = "T>G",
    "A>G" = "T>C",
    "A>T" = "T>A"
  )
  spontaneous_deamination <- c(
    "A>G",
    "C>T",
    # complementary
    "T>C",
    "G>A"
  )
  ros_damage <- c(
    "G>T",
    "G>C",
    # complementary
    "C>A",
    "C>G"
  )
  other_damage <- c(
    "T>G",
    "T>A",
    # complementary
    "A>C",
    "A>T"
  )
  .v |>
    dplyr::mutate(
      variant_short = gsub(
        "[0-9]*",
        "",
        variant
      )
    ) |>
    tidyr::separate(
      variant_short,
      into = c("ref", "alt"),
      remove = FALSE,
    ) |>
    dplyr::mutate(
      variant_location = ifelse(
        dplyr::between(
          Position,
          Other[1],
          Other[2]
        ),
        "Other region",
        "Ori region"
      )
    ) |>
    dplyr::mutate(
      L_H_strand = ifelse(
        ref %in% L_ref,
        "L",
        "H"
      )
    ) |>
    dplyr::mutate(
      deamination_ros = dplyr::case_when(
        variant_short %in% spontaneous_deamination ~ "Spontaneous deamination",
        variant_short %in% ros_damage ~ "ROS damage",
        variant_short %in% other_damage ~ "Other damage",
        TRUE ~ "Unknown"
      )
    ) |>
    dplyr::mutate(
      variant_six = purrr::map2_chr(
        .x = variant_short,
        .y = L_H_strand,
        .f = function(x, y) {
          if (y == "L") {
            return(x)
          } else {
            return(variant_reverse[x])
          }
        }
      )
    )
}

fn_plot_pie <- function(.d, .colors = NULL) {
  .d |>
    dplyr::select(group = 1, n) |>
    dplyr::arrange(-n) |>
    dplyr::mutate(csum = rev(cumsum(rev(n)))) %>%
    dplyr::mutate(pos = n / 2 + dplyr::lead(csum, 1)) %>%
    dplyr::mutate(pos = dplyr::if_else(is.na(pos), n / 2, pos)) %>%
    dplyr::mutate(percentage = n / sum(n)) |>
    dplyr::mutate(group = factor(group, levels = group)) -> .dd

  .scalefill <- if (is.null(.colors)) {
    ggsci::scale_fill_aaas(
      name = NULL
    )
  } else {
    scale_fill_manual(
      name = NULL,
      values = .colors
    )
  }
  .scalecolor <- if (is.null(.colors)) {
    ggsci::scale_color_aaas(
      name = NULL
    )
  } else {
    scale_color_manual(
      name = NULL,
      values = .colors
    )
  }

  .dd |>
    ggplot(aes(
      x = "",
      y = n,
    )) +
    geom_bar(
      aes(fill = group),
      stat = "identity",
      width = 1,
      color = "white",
      show.legend = FALSE
    ) +
    .scalefill +
    ggrepel::geom_label_repel(
      aes(
        y = pos,
        label = glue::glue("{group}\n{n} ({scales::percent(percentage)})"),
        color = group,
      ),
      size = 6,
      nudge_x = 1,
      nudge_y = 0,
      show.legend = FALSE,
      max.overlaps = Inf,
    ) +
    .scalecolor +
    coord_polar(theta = "y", start = 0) +
    theme_void() +
    theme(
      plot.title = element_text(
        hjust = 0.5,
        size = 22,
      ),
      # legend.position = "none"
    )
}

fn_plot_bar <- function(.d) {
  .d |>
    dplyr::group_by(
      variant_six,
      L_H_strand
    ) |>
    dplyr::count() |>
    dplyr::ungroup() |>
    dplyr::mutate(percent = n / sum(n)) |>
    dplyr::rename(
      sigvar = variant_six,
      strand = L_H_strand,
    ) |>
    ggplot(aes(x = sigvar, y = percent, fill = sigvar, alpha = strand)) +
    geom_hline(yintercept = c(0, 0.2, 0.4), color = "gray90") +
    geom_bar(stat = "identity", position = "dodge") +
    # facet_wrap(~ factor(tissue2, levels = c("colon", "fibroblast", "blood"))) +
    theme_classic() +
    # geom_text(aes(label=round(percent,1)), position=position_dodge(width=1), vjust=-0.5) +
    ggtitle("mutational signature") +
    ggsci::scale_fill_cosmic(
      palette = "signature_substitutions",
      name = "Mutation"
    ) +
    scale_alpha_manual(values = c(1, 0.5)) +
    theme(
      axis.title.x = element_blank(),
      plot.title = element_text(hjust = 0.5),
      axis.text.x = element_text(
        angle = 45,
        hjust = 1
      ),
    ) +
    labs(
      y = "Proportion"
    )
}
# body --------------------------------------------------------------------

v_hete_L_H_strand <- fn_variant_L_H_strand(v_hete)
v_homo_hete_L_H_strand <- fn_variant_L_H_strand(v_homo_hete)

# fn_variant_L_H_strand(v_homo_hete) |> fn_plot_bar()
fn_variant_L_H_strand(v_homo_hete) |>
  dplyr::count(variant_six) |>
  fn_plot_pie()


fn_plot_bar(fn_variant_L_H_strand(
  v_hete
)) +
  labs(title = "495 Heteroplasmic variants") -> p_hete

fn_plot_bar(fn_variant_L_H_strand(
  v_homo
)) +
  labs(title = "1049 Homoplasmic variants") -> p_homo

fn_plot_bar(fn_variant_L_H_strand(v_homo_hete)) +
  labs(title = "1544 Homoplasmic and Heteroplasmic variants") -> p_homo_hete


wrap_plots(
  p_hete,
  plot_spacer(),
  p_homo,
  plot_spacer(),
  p_homo_hete,
  nrow = 1,
  widths = c(15, -1.05, 15, -1.05, 10),
  guides = "collect"
) -> p_heavy_light
p_heavy_light

outdir <- "/home/liuc9/github/scMOCHA-data/analysis/zzz/plot-heavy-light"
ggsave(
  path = outdir,
  filename = "heavy_light_variants-fisher.pdf",
  plot = p_heavy_light,
  width = 11,
  height = 5,
  dpi = 300
)

# footer ------------------------------------------------------------------

# future: :plan(future: :sequential)

# save image --------------------------------------------------------------
