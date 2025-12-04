#!/usr/bin/env Rscript
# Metainfo ----------------------------------------------------------------
# @AUTHOR: Chun-Jie Liu
# @CONTACT: chunjie.sam.liu.at.gmail.com
# @DATE: 2025-12-02 16:23:08
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
gse_data_variant_classification_clusteraf_bulkaf <- import(
  "/home/liuc9/github/scMOCHA-data/analysis/high-res/MANUSCRIPTFIGURES/SAMPLE-VARIANT-CLASSIFICATION-CLUSTER-BULK-AF.xlsx"
)
# load conn ---------------------------------------------------------------

# src ---------------------------------------------------------------------
source("/home/liuc9/github/scMOCHA-data/analysis/high-res/plot_somatic.R")

# function ----------------------------------------------------------------
fn_plot_af_variant_type <- function(
  thevariant,
  df = gse_data_variant_classification_clusteraf_bulkaf
) {
  df |>
    dplyr::filter(variant == thevariant) |>
    dplyr::arrange(variant_type) |>
    dplyr::mutate(
      srrid = forcats::fct_reorder(srrid, Bulk, .desc = TRUE)
    ) -> df_thevariant

  df_thevariant |>
    dplyr::count(variant_type) |>
    tibble::deframe() -> n_hete_homo

  df_thevariant |>
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
    ) -> forplot

  forplot |>
    dplyr::select(srrid, Haplogroup, Verbose_haplogroup) |>
    tidyr::pivot_longer(
      cols = c(Haplogroup, Verbose_haplogroup),
      names_to = "type",
      values_to = "Haplogroup"
    ) |>
    ggplot(aes(
      x = type,
      y = srrid,
      fill = Haplogroup
    )) +
    geom_tile(
      show.legend = FALSE
    ) +
    geom_text(
      aes(
        label = Haplogroup
      ),
      color = "black",
    ) +
    scale_x_discrete(
      expand = c(0, 0)
    ) +
    theme_classic() +
    theme(
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      axis.title.x = element_blank(),
      axis.line = element_blank(),
      axis.title.y = element_blank(),
    ) -> p_haplogroup

  forplot |>
    ggplot(aes(
      x = celltype,
      y = srrid,
      fill = af
    )) +
    geom_tile() +
    scale_fill_gradient(
      name = "AF",
      low = "white",
      high = "red"
    ) +
    scale_x_discrete(
      expand = c(0, 0)
    ) +
    theme_classic() +
    labs(
      x = "Celltype",
      y = "Sample"
    ) +
    theme(
      plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
      axis.text.x = element_text(face = "bold", size = 12),
      # axis.title.x = element_text(face = "bold", size = 12),
      axis.title = element_blank(),
      axis.line.y = element_blank(),
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank()
    ) -> p_af

  df_thevariant |>
    dplyr::select(srrid, variant_type) |>
    dplyr::mutate(v = 1) |>
    tidyr::pivot_wider(
      names_from = variant_type,
      values_from = v,
      values_fill = NA_integer_
    ) |>
    dplyr::mutate(
      srrid = factor(srrid, forplot$srrid |> levels())
    ) |>
    tidyr::pivot_longer(
      cols = -srrid,
      names_to = "variant_type",
      values_to = "value"
    ) |>
    ggplot(
      aes(
        x = variant_type,
        y = srrid,
        fill = value
      )
    ) +
    geom_tile() +
    scale_x_discrete(
      expand = c(0, 0),
      limits = c("homo", "haplo", "hete", "somatic"),
      labels = c(
        "homo" = "Homoplasmic",
        "haplo" = "Ethnicity",
        "hete" = "Heteroplasmic",
        "somatic" = "Somatic"
      )
    ) +
    scale_fill_gradient(
      name = "Presence",
      low = "white",
      high = "blue"
    ) +
    theme_classic() +
    labs(
      x = "Variant Type",
      y = "Sample"
    ) +
    theme(
      plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
      axis.text.x = element_text(face = "bold", size = 12),
      # axis.title.x = element_text(face = "bold", size = 12),
      axis.title = element_blank(),
      axis.line.y = element_blank(),
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank()
    ) -> p_variant_type

  wrap_plots(
    p_haplogroup,
    p_af,
    p_variant_type,
    ncol = 3,
    widths = c(0.2, 1, 0.8),
    guides = "collect"
  ) +
    plot_annotation(
      title = glue::glue(
        "Variant {thevariant} Allele Frequency\n{n_hete_homo['hete']} samples heteroplasmic, {n_hete_homo['homo']} samples homoplasmic and 2 samples somatic"
      ),
      theme = theme(
        plot.title = element_text(size = 14, face = "bold", hjust = 0.5)
      )
    ) -> p_collect
  p_collect
}


fn_plot_af_individual <- function(
  thesrrid,
  df = gse_data_variant_classification_clusteraf_bulkaf
) {
  df |>
    dplyr::filter(srrid == thesrrid) |>
    dplyr::arrange(variant_type) |>
    dplyr::mutate(
      srrid = forcats::fct_reorder(srrid, Bulk, .desc = TRUE),
      variant = forcats::fct_reorder(variant, Bulk, .desc = TRUE)
    ) -> df_thevariant

  df_thevariant |>
    dplyr::count(variant_type) |>
    tibble::deframe() -> n_hete_homo

  df_thevariant |>
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
    ) -> forplot

  # forplot |>
  #   ggplot(aes(
  #     x = 1,
  #     y = variant,
  #     fill = variant_type
  #   )) +
  #   geom_tile(
  #     show.legend = FALSE
  #   ) +
  #   geom_text(
  #     aes(
  #       label = variant_type
  #     ),
  #     color = "black",
  #   ) +
  #   scale_x_continuous(
  #     expand = c(0, 0)
  #   ) +
  #   theme_classic() +
  #   theme(
  #     axis.text.x = element_blank(),
  #     axis.ticks.x = element_blank(),
  #     axis.title.x = element_blank(),
  #     axis.line = element_blank(),
  #     axis.title.y = element_blank(),
  #   ) -> p_haplogroup

  forplot |>
    ggplot(aes(
      x = celltype,
      y = variant,
      fill = af
    )) +
    geom_tile() +
    scale_fill_gradient(
      name = "AF",
      low = "white",
      high = "red"
    ) +
    scale_x_discrete(
      expand = c(0, 0)
    ) +
    theme_classic() +
    labs(
      x = "Celltype",
      y = "Sample"
    ) +
    theme(
      plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
      axis.text.x = element_text(face = "bold", size = 12),
      # axis.title.x = element_text(face = "bold", size = 12),
      axis.title = element_blank(),
      # axis.line.y = element_blank(),
      # axis.text.y = element_blank(),
      # axis.ticks.y = element_blank()
    ) -> p_af

  df_thevariant |>
    dplyr::select(variant, variant_type) |>
    dplyr::mutate(v = 1) |>
    tidyr::pivot_wider(
      names_from = variant_type,
      values_from = v,
      values_fill = NA_integer_
    ) |>
    dplyr::mutate(
      variant = factor(variant, forplot$variant |> levels())
    ) |>
    tidyr::pivot_longer(
      cols = -variant,
      names_to = "variant_type",
      values_to = "value"
    ) |>
    ggplot(
      aes(
        x = variant_type,
        y = variant,
        fill = value
      )
    ) +
    geom_tile() +
    scale_x_discrete(
      expand = c(0, 0),
      limits = c("homo", "haplo", "hete", "somatic"),
      labels = c(
        "homo" = "Homoplasmic",
        "haplo" = "Ethnicity",
        "hete" = "Heteroplasmic",
        "somatic" = "Somatic"
      )
    ) +
    scale_fill_gradient(
      name = "Presence",
      low = "white",
      high = "blue"
    ) +
    theme_classic() +
    labs(
      x = "Variant Type",
      y = "Sample"
    ) +
    theme(
      plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
      axis.text.x = element_text(face = "bold", size = 12),
      # axis.title.x = element_text(face = "bold", size = 12),
      axis.title = element_blank(),
      axis.line.y = element_blank(),
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank()
    ) -> p_variant_type

  wrap_plots(
    p_af,
    p_variant_type,
    ncol = 2,
    widths = c(1, 0.8),
    guides = "collect"
  ) +
    plot_annotation(
      title = glue::glue(
        "Variant {thevariant} Allele Frequency\n{unique(df_thevariant$gseid)}-{unique(df_thevariant$srrid)}"
      ),
      theme = theme(
        plot.title = element_text(size = 14, face = "bold", hjust = 0.5)
      )
    ) -> p_collect
  p_collect
}


# body --------------------------------------------------------------------
gse_data_variant_classification_clusteraf_bulkaf |>
  dplyr::filter(srrid == "GSM4712887")
#
# Euler venn diagram--------------------------------------------------------------------
#
#

gse_data_variant_classification_clusteraf_bulkaf |>
  dplyr::select(variant_type, variant) |>
  dplyr::distinct() -> variants_classification_unique

variants_classification_unique |>
  dplyr::mutate(isvariant = TRUE) |>
  tidyr::pivot_wider(
    names_from = variant_type,
    values_from = isvariant,
    values_fill = FALSE
  ) |>
  dplyr::rename(
    Homoplasmic = homo,
    Heteroplasmic = hete,
    Ethnicity = haplo,
    Somatic = somatic
  ) -> dt_allvariants_euler

color_variant_type <- c(
  "haplo" = "#BEBADAFF",
  "homo" = "#8DD3C7FF",
  "hete" = "#FFFFB3FF",
  "somatic" = "red"
)

library(eulerr)
fit_euler <- euler(
  dt_allvariants_euler[, c(
    "Ethnicity",
    "Homoplasmic",
    "Heteroplasmic",
    "Somatic"
  )],
  # shape = "ellipse",
  control = list(extraopt = FALSE)
)
plot(fit_euler)
plot(
  fit_euler,
  quantities = list(type = c("counts"), font = 3),
  labels = list(fontfamily = "serif"),
  edges = list(lty = 3),
  # fills = c("#BEBADAFF", "#8DD3C7FF", "#FFFFB3FF", "red")
  fills = color_variant_type
) -> p_euler


{
  outdir = "/home/liuc9/github/scMOCHA-data/analysis/high-res/MANUSCRIPTFIGURES/notuse"
  pdf(
    path(
      outdir,
      "Variant-type-Euler.pdf"
    ),
    width = 8,
    height = 6
  )
  print(p_euler)
  dev.off()
}
#
#
# ? variant is both hete and homo --------------------------------------------------------------------
#
#

{
  outdir <- "/home/liuc9/github/scMOCHA-data/analysis/high-res/MANUSCRIPTFIGURES/notuse"
  thevariant <- "2442T>C"
  ggsave(
    filename = glue::glue("Example-variant-{thevariant}-Hete-Homo.pdf"),
    plot = fn_plot_af_variant_type(thevariant),
    path = outdir,
    width = 6,
    height = 11
  )
  thevariant <- "8362T>G"
  ggsave(
    filename = glue::glue("Example-variant-{thevariant}-Hete-Homo.pdf"),
    plot = fn_plot_af_variant_type(thevariant),
    path = outdir,
    width = 6,
    height = 11
  )
}


#
#
# ? plot intersection --------------------------------------------------------------------
#
#

fn_plot_af_variant_type("15244A>G")
fn_plot_af_individual("GSM6793473")

dt_allvariants_euler |>
  dplyr::mutate(
    ncat = rowSums(dplyr::across(c(
      Ethnicity,
      Homoplasmic,
      Heteroplasmic,
      Somatic
    )))
  ) |>
  dplyr::filter(
    ncat > 1
  ) |>
  dplyr::filter(
    !(Homoplasmic & Ethnicity & !Heteroplasmic & !Somatic)
  ) |>
  dplyr::filter(
    !(Heteroplasmic & Somatic & !Homoplasmic & !Ethnicity)
  ) |>
  dplyr::arrange(ncat) |>
  dplyr::select(-ncat) |>
  tidyr::pivot_longer(
    -variant,
    names_to = "variant_type",
    values_to = "value"
  ) |>
  dplyr::filter(value) |>
  tidyr::nest(
    .by = variant,
  ) |>
  dplyr::mutate(
    label = purrr::map_chr(
      .x = data,
      .f = \(.d) {
        .d$variant_type |>
          paste(collapse = "-")
      }
    )
  ) |>
  dplyr::mutate(
    p = parallel::mclapply(
      X = variant,
      FUN = \(.v) {
        fn_plot_af_variant_type(.v)
      },
      mc.cores = 10
    )
  ) |>
  dplyr::select(-data) -> df_variants_intersection_plot


{
  outdir <- path(
    "/home/liuc9/github/scMOCHA-data/analysis/high-res/MANUSCRIPTFIGURES/notuse/euler-intersection-variants"
  )
  dir_create(outdir)
  df_variants_intersection_plot |>
    dplyr::mutate(
      a = parallel::mcmapply(
        FUN = \(variant, label, p) {
          ggsave(
            filename = glue::glue("{label}-{variant}.pdf"),
            plot = p,
            path = outdir,
            width = 15,
            height = 8
          )
        },
        variant,
        label,
        p,
        mc.cores = 10,
        SIMPLIFY = FALSE
      )
    )
}

#
#
# ? somatic ethnicity --------------------------------------------------------------------
#
#

{
  outdir <- "/home/liuc9/github/scMOCHA-data/analysis/high-res/MANUSCRIPTFIGURES/notuse"
  thevariant <- "9449C>T"
  fn_plot_af_variant_type(thevariant)
  ggsave(
    filename = glue::glue("Example-variant-{thevariant}-somatic.pdf"),
    plot = fn_plot_af_variant_type(thevariant),
    path = outdir,
    width = 15,
    height = 8
  )
}


#
#
# ? variant multiple allele--------------------------------------------------------------------
#
#

gse_data_variant_classification_clusteraf_bulkaf |>
  dplyr::filter(variant_type %in% c("homo", "hete")) -> df_variants


df_variants |>
  tidyr::nest(.by = variant, .key = "sample") -> df_variants_sample


df_variants_sample |>
  dplyr::mutate(
    a = parallel::mclapply(
      sample,
      function(.x) {
        # .x <- df_variants_sample$sample[[1]]
        .x |>
          dplyr::count(variant_type) |>
          tidyr::pivot_wider(
            names_from = variant_type,
            values_from = n,
            values_fill = NA_integer_,
            names_prefix = "n_"
          )
      },
      mc.cores = 10
    )
  ) |>
  dplyr::select(-sample) |>
  tidyr::unnest(a) -> df_variants_sample_count
#

df_variants_sample_count |>
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
  ) -> df_variants_sample_count_coord


df_variants$gseid |>
  unique() |>
  length()
df_variants$srrid |>
  unique() |>
  length()
df_variants_sample_count_coord$variant |>
  unique() |>
  length()
df_variants_sample_count_coord$start |> unique() |> length()

df_variants_sample_count_coord |>
  dplyr::arrange(start) |>
  dplyr::group_by(start) |>
  dplyr::filter(dplyr::n() > 2)


# thevariant <- "2442T>C"
thevariants1 <- c("13761A>G", "13761A>C", "13761A>T")
thevariants2 <- c("7585A>T", "7585A>C", "7585A>G")
thevariants <- c(thevariants1, thevariants2)
df_variants |>
  dplyr::filter(variant %in% thevariants)

df_variants |>
  dplyr::filter(variant %in% thevariants) |>
  dplyr::arrange(variant_type) |>
  dplyr::mutate(
    srrid = forcats::fct_reorder(srrid, Bulk, .desc = TRUE)
  ) |>
  tidyr::pivot_longer(
    cols = c(B, CD4_T, CD8_T, DC, Mono, NK, other, Bulk, other_T),
    names_to = "celltype",
    values_to = "af"
  ) |>
  dplyr::filter(celltype == "Bulk") |>
  ggplot(aes(
    x = variant,
    y = srrid,
    fill = af
  )) +
  geom_tile() +
  scale_fill_gradient(
    name = "AF",
    low = "white",
    high = "red"
  ) +
  theme_classic() +
  labs(
    x = "Variant",
    y = "Sample"
  ) +
  labs(
    title = glue::glue(
      "Variant in Bulk AF"
    ),
  ) +
  theme(
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
  ) -> p_multiple_variants_bulk

{
  outdir <- "/home/liuc9/github/scMOCHA-data/analysis/high-res/MANUSCRIPTFIGURES/notuse"
  ggsave(
    filename = glue::glue("Example-variant-multi.pdf"),
    plot = p_multiple_variants_bulk,
    path = outdir,
    width = 6,
    height = 8
  )
}

#
#
# plot somatic --------------------------------------------------------------------
#
#

# footer ------------------------------------------------------------------

# save image --------------------------------------------------------------
