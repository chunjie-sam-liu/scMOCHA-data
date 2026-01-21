#!/usr/bin/env Rscript
# Metainfo ----------------------------------------------------------------
# @AUTHOR: Chun-Jie Liu
# @CONTACT: chunjie.sam.liu.at.gmail.com
# @DATE: 2025-12-18 16:50:17
# @DESCRIPTION: this script is used for ...

# Library -----------------------------------------------------------------

load_pkg(jutils)

dotenv(".env")

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
dotenv(".env")
outdir <- path(Sys.getenv("OUTDIR"))

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
          Position = pos
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
HOMO_HETE_VARIANTS <- ALLVARIANTS |>
  dplyr::filter(variant_type %in% c("homo", "hete"))
# load conn ---------------------------------------------------------------

# src ---------------------------------------------------------------------

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
    dplyr::mutate(csum = rev(cumsum(rev(n)))) |>
    dplyr::mutate(pos = n / 2 + dplyr::lead(csum, 1)) |>
    dplyr::mutate(pos = dplyr::if_else(is.na(pos), n / 2, pos)) |>
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
    dplyr::select(variant, variant_six, L_H_strand) -> .m
  .m |> distinct() -> .mm
  .m |>
    # dplyr::distinct() |> # don't use distinct here, count all variants
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
    ) -> .dd

  .dd |>
    select(-n) |>
    pivot_wider(
      names_from = strand,
      values_from = percent,
      values_fill = 0
    ) |>
    pivot_longer(
      cols = c("L", "H"),
      names_to = "strand",
      values_to = "percent"
    ) |>
    ggplot(aes(
      x = sigvar,
      y = percent,
      fill = sigvar,
      alpha = strand
    )) +
    geom_hline(yintercept = c(0, 0.2, 0.4), color = "gray90") +
    geom_bar(
      stat = "identity",
      position = position_dodge(preserve = "single")
    ) +
    # facet_wrap(~ factor(tissue2, levels = c("colon", "fibroblast", "blood"))) +
    theme_classic() +
    # geom_text(aes(label=round(percent,1)), position=position_dodge(width=1), vjust=-0.5) +
    ggsci::scale_fill_cosmic(
      palette = "signature_substitutions",
      name = "Mutation"
    ) +
    scale_alpha_manual(values = c(1, 0.5)) +
    theme(
      axis.title.x = element_blank(),
      plot.title = element_text(hjust = 0.5),
      plot.subtitle = element_text(hjust = 0.5),
      axis.text.x = element_text(
        angle = 45,
        hjust = 1
      ),
    ) +
    labs(
      title = "Variant Light/Heavy Chain Distribution",
      subtitle = glue(
        "Total variants: {nrow(.m)}, Unique variants: {nrow(.mm)}"
      ),
      y = "Proportion"
    )
}

fn_plot_pie_nature_aging <- function(.d) {
  .d |>
    dplyr::select(
      variant,
      deamination_ros,
      variant_six
    ) |>
    # dplyr::distinct() |> # don't use distinct here, count all variants
    dplyr::count(
      deamination_ros,
      variant_six
    ) |>
    dplyr::mutate(
      variant = dplyr::recode(
        variant_six,
        "C>A" = "C:G>A:T",
        "C>G" = "C:G>G:C",
        "C>T" = "C:G>T:A",
        "T>A" = "T:A>A:T",
        "T>C" = "T:A>C:G",
        "T>G" = "T:A>G:C"
      )
    ) |>
    dplyr::mutate(
      variant = factor(
        variant,
        levels = c(
          "T:A>C:G",
          "C:G>T:A",
          "C:G>A:T",
          "C:G>G:C",
          "T:A>G:C",
          "T:A>A:T"
        )
      )
    ) |>
    dplyr::arrange(
      variant
    ) |>
    dplyr::mutate(
      pick_color = c(
        "#29833C",
        "#5BB446",
        "#F28107",
        "#FFB10A",
        "#01467B",
        "#25648C"
      )
    ) -> .dd

  .dd |>
    # dplyr::select(group = 1, n) |>
    # dplyr::arrange(-n) |>
    dplyr::mutate(group = variant) |>
    dplyr::mutate(csum = rev(cumsum(rev(n)))) |>
    dplyr::mutate(pos = n / 2 + dplyr::lead(csum, 1)) |>
    dplyr::mutate(pos = dplyr::if_else(is.na(pos), n / 2, pos)) |>
    dplyr::mutate(percentage = n / sum(n)) |>
    dplyr::mutate(
      group = factor(group, levels = group),
      pick_color = factor(pick_color, levels = pick_color)
    ) -> forplot_pie

  forplot_pie |>
    ggplot(aes(
      x = "",
      y = n,
    )) +
    geom_bar(
      aes(fill = group),
      stat = "identity",
      width = 1,
      color = "white",
      # show.legend = TRUE
    ) +
    # scale_fill_identity() +
    scale_fill_manual(
      name = "Base substitution",
      values = levels(forplot_pie$pick_color)
    ) +
    ggrepel::geom_label_repel(
      aes(
        y = pos,
        # label = glue::glue("{group}\n{n} ({scales::percent(percentage)})"),
        label = glue::glue("{scales::percent(percentage, accuracy = 0.01)}"),
        color = group,
      ),
      size = 6,
      # nudge_x = 1,
      # nudge_y = 0,
      show.legend = FALSE,
      max.overlaps = Inf,
    ) +
    # .scalecolor +
    # scale_color_identity() +
    scale_color_manual(
      name = NULL,
      values = levels(forplot_pie$pick_color)
    ) +
    coord_polar(theta = "y", start = 0) +
    theme_void() +
    theme(
      plot.title = element_text(
        hjust = 0.5,
        size = 22,
      ),
      legend.position = "right"
    )
}

# body --------------------------------------------------------------------
# fn_variant_L_H_strand(SOMATIC_VARIANTS) |>
#   dplyr::count(variant_six) |>
#   fn_plot_pie()

# fn_plot_bar(fn_variant_L_H_strand(SOMATIC_VARIANTS))
{
  # somatic variants
  pdf(
    file = outdir / "PBMC-VARIANT-LIGHT-HEAVY-CHAIN-SOMATIC.pdf",
    width = 8,
    height = 6
  )
  print(fn_plot_bar(fn_variant_L_H_strand(SOMATIC_VARIANTS)))
  print(fn_plot_pie_nature_aging(fn_variant_L_H_strand(SOMATIC_VARIANTS)))
  dev.off()
}
{
  # hete
  pdf(
    file = outdir / "PBMC-VARIANT-LIGHT-HEAVY-CHAIN-HETE.pdf",
    width = 8,
    height = 6
  )
  print(fn_plot_bar(fn_variant_L_H_strand(HETE_VARIANTS)))
  print(fn_plot_pie_nature_aging(fn_variant_L_H_strand(HETE_VARIANTS)))
  dev.off()
}
{
  # homo
  pdf(
    file = outdir / "PBMC-VARIANT-LIGHT-HEAVY-CHAIN-HOMO.pdf",
    width = 8,
    height = 6
  )
  print(fn_plot_bar(fn_variant_L_H_strand(HOMO_VARIANTS)))
  print(fn_plot_pie_nature_aging(fn_variant_L_H_strand(HOMO_VARIANTS)))
  dev.off()
}
{
  # homo_hete
  pdf(
    file = outdir / "PBMC-VARIANT-LIGHT-HEAVY-CHAIN-HOMO-HETE.pdf",
    width = 8,
    height = 6
  )
  print(fn_plot_bar(fn_variant_L_H_strand(HOMO_HETE_VARIANTS)))
  print(fn_plot_pie_nature_aging(fn_variant_L_H_strand(HOMO_HETE_VARIANTS)))
  dev.off()
}
# footer ------------------------------------------------------------------

# save image --------------------------------------------------------------
