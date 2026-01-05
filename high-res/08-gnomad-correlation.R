#!/usr/bin/env Rscript
# Metainfo ----------------------------------------------------------------
# @AUTHOR: Chun-Jie Liu
# @CONTACT: chunjie.sam.liu.at.gmail.com
# @DATE: 2025-12-17 14:15:55
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
allvariants <- import(
  "/home/liuc9/github/scMOCHA-data/analysis/high-res/MANUSCRIPTFIGURES/SAMPLE-VARIANT-CLASSIFICATION-CLUSTER-BULK-AF.xlsx"
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

allvariants |>
  dplyr::filter(variant_type == "somatic") |>
  dplyr::pull(variant) |>
  sort() |>
  unique() -> somatic_variants_list

gnomad <- import(
  "/home/liuc9/github/scMOCHA-data/analysis/zzz/db/gnomad.qs"
) |>
  dplyr::filter(filters == "PASS") |>
  dplyr::mutate(
    variant = glue("{position}{refna}>{regna}")
  ) |>
  dplyr::select(position, variant, gnomad_paf = af_hom) |>
  dplyr::arrange(position)

# load conn ---------------------------------------------------------------

# src ---------------------------------------------------------------------

# function ----------------------------------------------------------------

# body --------------------------------------------------------------------
N_INDIVIDUALS <- length(unique(allvariants$srrid))


#
#
# ? homo --------------------------------------------------------------------
#
#

allvariants |>
  dplyr::filter(variant_type == "homo") |>
  tidyr::nest(
    .by = c(variant, start, ref, alt)
  ) |>
  dplyr::mutate(
    n = purrr::map_int(data, nrow),
    avg_af = purrr::map_dbl(data, \(.df) {
      mean(.df$Bulk, na.rm = TRUE)
    })
  ) |>
  dplyr::select(-data) |>
  dplyr::mutate(paf = n / N_INDIVIDUALS) |>
  dplyr::left_join(
    gnomad,
    by = c("start" = "position", "variant" = "variant")
  ) |>
  dplyr::mutate(
    gnomad_paf = ifelse(is.na(gnomad_paf), 0, gnomad_paf)
  ) -> homo_variants


fn_xy_breaks_limits(
  homo_variants$paf,
  step = 0.2,
  max = FALSE
) -> xbl_homo
fn_xy_breaks_limits(
  homo_variants$gnomad_paf,
  step = 0.2,
  max = FALSE
) -> ybl_homo


homo_variants |>
  ggpubr::ggscatter(
    x = "paf",
    y = "gnomad_paf",
    cor.coef = TRUE, # Add correlation coefficient. see ?stat_cor
    cor.coeff.args = list(
      method = "pearson",
      label.x = 0.1,
      label.sep = "\n",
      size = 5
    )
  ) +
  # geom_point(alpha = 0.7) +
  # ggpointdensity::geom_pointdensity(
  #   adjust = 0.01,
  #   show.legend = FALSE,
  # ) +
  # viridis::scale_color_viridis() +
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
    x = glue(
      "Homoplasmic variant population frequency\n(individual={N_INDIVIDUALS}, variant={nrow(homo_variants)})"
    ),
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
  ) +
  coord_fixed() -> p_homo_corr_gnomad


ggsave(
  filename = path(
    "/home/liuc9/github/scMOCHA-data/analysis/high-res/MANUSCRIPTFIGURES",
    "gnomAD-HOMOPLASMIC-VARIANT-CORRELATION.pdf"
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
mtdna <- import("/home/liuc9/github/scMOCHA-data/config/mtdna_genes_dloop.qs")


allvariants |>
  dplyr::filter(variant_type == "hete") |>
  tidyr::nest(
    .by = c(variant, start, ref, alt)
  ) |>
  dplyr::mutate(
    n = purrr::map_int(data, nrow),
    avg_af = purrr::map_dbl(data, \(.df) {
      mean(.df$Bulk, na.rm = TRUE)
    })
  ) |>
  dplyr::select(-data) |>
  dplyr::mutate(paf = n / N_INDIVIDUALS) |>
  dplyr::left_join(
    gnomad,
    by = c("start" = "position", "variant" = "variant")
  ) |>
  dplyr::mutate(
    gnomad_paf = ifelse(is.na(gnomad_paf), 0, gnomad_paf)
  ) |>
  dplyr::mutate(
    variant_type = ifelse(
      paf >= 0.1,
      "Psudo-bulk AF >= 10%",
      "Psudo-bulk AF < 10%"
    )
  ) -> hete_variants


fn_xy_breaks_limits(
  hete_variants$paf,
  step = 0.05,
  max = FALSE
) -> xbl_hete
fn_xy_breaks_limits(
  hete_variants$gnomad_paf,
  step = 0.05,
  max = FALSE
) -> ybl_hete


hete_variants |>
  ggplot(aes(
    x = paf,
    y = gnomad_paf,
    color = variant_type,
  )) +
  geom_point(
    alpha = 0.7,
    size = 1.5
  ) +
  # ggrepel::geom_label_repel(
  #   data = forlabel,
  #   aes(
  #     label = variant,
  #   ),
  #   size = 3.5,
  #   nudge_x = 0.01,
  #   nudge_y = 0.01,
  #   show.legend = FALSE
  # ) +
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
    x = glue(
      "Heteroplasmic variant population frequency\n(individual={N_INDIVIDUALS}, variant={nrow(hete_variants)})"
    ),
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
  ) +
  coord_fixed() -> p_hete_corr_gnomad


ggsave(
  filename = path(
    "/home/liuc9/github/scMOCHA-data/analysis/high-res/MANUSCRIPTFIGURES",
    "gnomAD-HETEROPLASMIC-VARIANT-CORRELATION.pdf"
  ),
  plot = p_hete_corr_gnomad,
  width = 6,
  height = 5,
  dpi = 300
)

# footer ------------------------------------------------------------------

# save image --------------------------------------------------------------
