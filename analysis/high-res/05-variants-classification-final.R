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

# function ----------------------------------------------------------------

# body --------------------------------------------------------------------

gse_data_variant_classification_clusteraf_bulkaf |>
  dplyr::filter(variant_type %in% c("homo", "hete")) -> df_variants


#
#
# ? Variant --------------------------------------------------------------------
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


library(eulerr)
plot(
  euler(
    dt_allvariants_euler[, c(
      "Ethnicity",
      "Homoplasmic",
      "Heteroplasmic",
      "Somatic"
    )],
    # shape = "ellipse",
    control = list(extraopt = FALSE)
  ),
  quantities = list(type = c("counts"), font = 3),
  labels = list(fontfamily = "serif"),
  edges = list(lty = 3),
  fills = c("#BEBADAFF", "#8DD3C7FF", "#FFFFB3FF", "red")
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

df_variants_sample_count |>
  dplyr::filter(!is.na(n_hete), !is.na(n_homo)) |>
  dplyr::mutate(total = (n_hete + n_homo))

thevariant <- "2442T>C"
thevariant <- "8362T>G"
df_variants |>
  dplyr::filter(variant == thevariant) |>
  dplyr::arrange(variant_type) |>
  dplyr::count(variant_type) |>
  tibble::deframe() -> n_hete_homo

df_variants |>
  dplyr::filter(variant == thevariant) |>
  dplyr::arrange(variant_type) |>
  dplyr::mutate(
    srrid = forcats::fct_reorder(srrid, Bulk, .desc = TRUE)
  ) |>
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
        "DC",
        "Mono",
        "NK",
        "other",
        "other_T"
      )
    )
  ) |>
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
  theme_classic() +
  labs(
    x = "Celltype",
    y = "Sample"
  ) +
  labs(
    title = glue::glue(
      "Variant {thevariant} Allele Frequency\n{n_hete_homo['hete']} samples heteroplasmic and {n_hete_homo['homo']} samples homoplasmic"
    ),
  ) +
  theme(
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
  ) -> p_both_hete_homo

{
  outdir <- "/home/liuc9/github/scMOCHA-data/analysis/high-res/MANUSCRIPTFIGURES/notuse"
  ggsave(
    filename = glue::glue("Example-variant-{thevariant}-Hete-Homo.pdf"),
    plot = p_both_hete_homo,
    path = outdir,
    width = 6,
    height = 11
  )
}


#
#
# ? variant multiple --------------------------------------------------------------------
#
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

# footer ------------------------------------------------------------------

# save image --------------------------------------------------------------
