#!/usr/bin/env Rscript
# Metainfo ----------------------------------------------------------------
# @AUTHOR: Chun-Jie Liu
# @CONTACT: chunjie.sam.liu.at.gmail.com
# @DATE: 2025-08-21 14:14:19
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

conn <- DBI::dbConnect(
  duckdb::duckdb(),
  "/home/liuc9/github/scMOCHA-data/analysis/zzz/clean-data/all_hetero_af.cell.duckdb.1.2.1"
)

DBI::dbListTables(conn)
tbl_allvariants <- dplyr::tbl(conn, "allvariants_fisher")

tbl_all_hetero_af_bulk <- dplyr::tbl(
  conn,
  "all_hetero_af_bulk_fisher"
)
DBI::dbListTables(conn)

tbl_gseid_srrid_variant <- dplyr::tbl(
  conn,
  "gseid_srrid_variant_fisher"
)

gse_data <- import(
  "/home/liuc9/github/scMOCHA-data/analysis/zzz/clean-data/gse_data_fisher.qs"
)

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

# body --------------------------------------------------------------------

#
#
# ? homo --------------------------------------------------------------------
#
#

tbl_allvariants |>
  dplyr::filter(issomatic == "homoplasmic") |>
  dplyr::group_by(Position) |>
  dplyr::filter(dplyr::n() == 1) |>
  as.data.table() |>
  dplyr::ungroup() |>
  dplyr::mutate(
    freq = n / nrow(gse_data),
    `Gnomad Frequency` = `Gnomad Frequency` / 100
  ) -> homo_variants

fn_xy_breaks_limits(
  homo_variants$freq,
  step = 0.2,
  max = FALSE
) -> xbl_homo
fn_xy_breaks_limits(
  homo_variants$`Gnomad Frequency`,
  step = 0.2,
  max = FALSE
) -> ybl_homo

homo_variants |>
  ggpubr::ggscatter(
    x = "freq",
    y = "Gnomad Frequency",
    # conf.int = TRUE,
    cor.coef = TRUE, # Add correlation coefficient. see ?stat_cor
    cor.coeff.args = list(
      method = "pearson",
      label.x = 0.1,
      label.sep = "\n",
      size = 5
    )
  ) +
  # geom_point(alpha = 0.7) +
  ggpointdensity::geom_pointdensity(
    adjust = 0.01,
    show.legend = FALSE,
  ) +
  viridis::scale_color_viridis() +
  scale_x_continuous(
    limits = xbl_homo$limits,
    breaks = xbl_homo$breaks,
    labels = scales::label_number(
      accuracy = 0.1
    ),
    expand = expansion(mult = c(0.01, 0.01))
  ) +
  scale_y_continuous(
    limits = ybl_homo$limits,
    breaks = ybl_homo$breaks,
    labels = scales::label_number(
      accuracy = 0.1
    ),
    expand = expansion(mult = c(0.01, 0.01))
  ) +
  labs(
    x = "Homoplasmic variant population frequency (n=577)",
    y = "gnomAD Frequency"
  ) +
  geom_abline(
    slope = 1,
    intercept = 0,
    linetype = "dashed",
    color = "darkgray"
  ) +
  theme(
    axis.title = element_text(size = 16),
  ) -> p_homo_corr_gnomad
p_homo_corr_gnomad

ggsave(
  filename = file.path(
    "/home/liuc9/github/scMOCHA-data/analysis/zzz/plot-heteroplasmic/gnomad",
    "homoplasmic_variant_correlates_with_gnomad-fisher.pdf"
  ),
  plot = p_homo_corr_gnomad,
  width = 6,
  height = 5,
  dpi = 300
)
#
#
# ? hete --------------------------------------------------------------------
#
#

tbl_allvariants |>
  dplyr::filter(issomatic == "heteroplasmic") |>
  dplyr::group_by(Position) |>
  dplyr::filter(dplyr::n() == 1) |>
  as.data.table() |>
  dplyr::ungroup() |>
  dplyr::mutate(
    freq = n / nrow(gse_data),
    `Gnomad Frequency` = `Gnomad Frequency` / 100
  ) -> hete_variants

fn_xy_breaks_limits(
  hete_variants$freq,
  step = 0.05,
  max = FALSE
) -> xbl_hete
fn_xy_breaks_limits(
  hete_variants$`Gnomad Frequency`,
  step = 0.05,
  max = FALSE
) -> ybl_hete

hete_variants |>
  ggpubr::ggscatter(
    x = "freq",
    y = "Gnomad Frequency",
    # conf.int = TRUE,
    cor.coef = TRUE, # Add correlation coefficient. see ?stat_cor
    cor.coeff.args = list(
      method = "pearson",
      label.x = 0.1,
      label.sep = "\n",
      size = 5
    )
  ) +
  # geom_point(alpha = 0.7) +
  ggpointdensity::geom_pointdensity(
    adjust = 0.01,
    show.legend = FALSE,
  ) +
  viridis::scale_color_viridis() +
  scale_x_continuous(
    limits = xbl_hete$limits,
    breaks = xbl_hete$breaks,
    labels = scales::label_number(
      accuracy = 0.1
    ),
    expand = expansion(mult = c(0.01, 0.01))
  ) +
  scale_y_continuous(
    limits = xbl_hete$limits,
    breaks = xbl_hete$breaks,
    labels = scales::label_number(
      accuracy = 0.1
    ),
    expand = expansion(mult = c(0.01, 0.01))
  ) +
  labs(
    x = "Heteroplasmic variant population frequency (n=577)",
    y = "gnomAD Frequency"
  ) +
  geom_abline(
    slope = 1,
    intercept = 0,
    linetype = "dashed",
    color = "darkgray"
  ) +
  theme(
    axis.title = element_text(size = 16),
  ) -> p_hete_corr_gnomad
p_hete_corr_gnomad
ggsave(
  filename = file.path(
    "/home/liuc9/github/scMOCHA-data/analysis/zzz/plot-heteroplasmic/gnomad",
    "heteroplasmic_variant_correlates_with_gnomad-fisher.pdf"
  ),
  plot = p_hete_corr_gnomad,
  width = 6,
  height = 5,
  dpi = 300
)


#
#
# ? new cutoff to test gnomad --------------------------------------------------------------------
#
#
mtdna <- import("/home/liuc9/github/scMOCHA-data/config/mtdna_genes_dloop.qs")
rnr1 <- mtdna[gene_name == "MT-RNR1", .(start, end)] |> as.numeric()
rnr2 <- mtdna[gene_name == "MT-RNR2", .(start, end)] |> as.numeric()


gseid_srrid_variant_hetero |>
  dplyr::left_join(
    tbl_all_hetero_af_bulk |>
      as.data.table(),
    by = c("gseid", "srrid", "variant")
  ) |>
  dplyr::filter(!is.na(barcode)) |>
  # dplyr::mutate(
  #   new_af = ifelse(af > 0.1, af, 0)
  # ) |>
  dplyr::mutate(
    gnomad_cutoff = ifelse(af > 0.1, TRUE, FALSE)
  ) |>
  tidyr::nest(
    .by = variant
  ) |>
  dplyr::mutate(
    a = purrr::map(
      .x = data,
      .f = ~ {
        .x |>
          dplyr::count(gnomad_cutoff) |>
          tibble::deframe() -> .n
        tibble::tibble(
          gnomad_gt_0.1 = .n["TRUE"],
          gnomad_lt_0.1 = .n["FALSE"]
        ) |>
          tidyr::replace_na(list(
            gnomad_gt_0.1 = 0,
            gnomad_lt_0.1 = 0
          ))
      }
    )
  ) |>
  dplyr::select(-data) |>
  tidyr::unnest(cols = c(a)) |>
  dplyr::left_join(
    hete_variants,
    by = "variant"
  ) |>
  dplyr::mutate(
    inrnr = ifelse(
      dplyr::between(Position, rnr1[1], rnr1[2]) |
        dplyr::between(Position, rnr2[1], rnr2[2]),
      TRUE,
      FALSE
    )
  ) -> gseid_srrid_variant_hetero_gnomad

sum(
  gseid_srrid_variant_hetero_gnomad$gnomad_gt_0.1 >= 0.1
) -> n_gnomad_gt_0.1
sum(
  gseid_srrid_variant_hetero_gnomad$gnomad_lt_0.1 < 0.1
) -> n_gnomad_lt_0.1

gseid_srrid_variant_hetero_gnomad |>
  dplyr::mutate(
    innoncoding = grepl(
      "rRNA|MitoTIP|non-coding",
      aachange,
    )
  ) |>
  dplyr::mutate(
    ingnomad = gnomad_gt_0.1 #/ gnomad_lt_0.1 > 1
  ) |>
  dplyr::mutate(
    variant_type = ifelse(
      ingnomad,
      "Psudo-bulk AF >= 0.1 (n={n_gnomad_gt_0.1})" |> glue::glue(),
      "Psudo-bulk AF < 0.1 (n={n_gnomad_lt_0.1})" |> glue::glue()
    )
  ) |>
  dplyr::mutate(
    variant_type = ifelse(
      innoncoding & ingnomad,
      "Psudo-bulk AF >= 0.1 & noncoding(n=14)",
      variant_type
    )
  ) |>
  dplyr::mutate(
    variant_type = ifelse(
      (!innoncoding) & ingnomad,
      "Psudo-bulk AF >= 0.1 & coding(n=40)",
      variant_type
    )
  ) -> forplot


fn_xy_breaks_limits(
  forplot$freq,
  step = 0.05,
  max = FALSE
) -> xbl_hete
fn_xy_breaks_limits(
  forplot$`Gnomad Frequency`,
  step = 0.05,
  max = FALSE
) -> ybl_hete


forplot |>
  dplyr::filter(ingnomad) |>
  dplyr::arrange(-freq) |>
  head(5) -> forlabel

forplot |>
  ggplot(aes(
    x = `freq`,
    y = `Gnomad Frequency`,
    color = variant_type,
  )) +
  geom_point(
    alpha = 0.7,
    size = 1.5
  ) +
  ggrepel::geom_label_repel(
    data = forlabel,
    aes(
      label = variant,
    ),
    size = 3.5,
    nudge_x = 0.01,
    nudge_y = 0.01,
    show.legend = FALSE
  ) +
  scale_x_continuous(
    limits = xbl_hete$limits,
    breaks = xbl_hete$breaks,
    labels = scales::label_number(
      accuracy = 0.01
    ),
    expand = expansion(mult = c(0.01, 0.01))
  ) +
  scale_y_continuous(
    limits = xbl_hete$limits,
    breaks = xbl_hete$breaks,
    labels = scales::label_number(
      accuracy = 0.01
    ),
    expand = expansion(mult = c(0.01, 0.01))
  ) +
  scale_color_brewer(
    palette = "Set1",
    name = "Variant Type"
  ) +
  labs(
    x = "Heteroplasmic variant population frequency (n=577)",
    y = "gnomAD Frequency"
  ) +
  geom_abline(
    slope = 1,
    intercept = 0,
    linetype = "dashed",
    color = "darkgray"
  ) +
  theme(
    axis.title = element_text(size = 16),
    panel.grid = element_blank(),
    panel.background = element_rect(fill = NA),
    axis.line = element_line(
      color = "black",
      linewidth = 0.5
    ),
    legend.position = c(0.5, 0.6)
  ) -> p_hete_corr_gnomad_newcutoff

ggsave(
  filename = file.path(
    "/home/liuc9/github/scMOCHA-data/analysis/zzz/plot-heteroplasmic/gnomad",
    "heteroplasmic_variant_correlates_with_gnomad_newcutoff-fisher.pdf"
  ),
  plot = p_hete_corr_gnomad_newcutoff,
  width = 6,
  height = 5,
  dpi = 300
)

# footer ------------------------------------------------------------------

# future: :plan(future: :sequential)

# save image --------------------------------------------------------------
