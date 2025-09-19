regions_missalignment_error <- c(
  66:71,
  300:316,
  513:525,
  3106:3107,
  12418:12425,
  16182:16194
)
regions_rare_heteroplasmic_variants <- c(499, 538, 545, 10953, 12684)
variants_tobe_excluded <- c(
  regions_missalignment_error,
  regions_rare_heteroplasmic_variants
)

gse_data <- import(
  "/home/liuc9/github/scMOCHA-data/analysis/zzz/clean-data/gse_data_fisher.qs"
)

gse_dataset_metadata_full <- import(
  "analysis/zzz/clean-data/gse_dataset_metadata_full.qs"
)
gse_data |>
  dplyr::select(
    gseid,
    srrid,
    chemistry,
    anno,
    hetero,
    haplo_variant,
    haplo_violin = haplo_violin2,
    somatic_variant,
    celltype_ratio,
    clusteraf,
    bulkaf
  ) |>
  dplyr::left_join(
    gse_dataset_metadata_full |> dplyr::select(-gseid),
    by = c("srrid" = "srrid")
  ) |>
  dplyr::arrange(disease, Chemistry) -> gse_data_haplo_variant

gse_data_haplo_variant |>
  dplyr::mutate(
    heteroplasmic = parallel::mcmapply(
      .x = somatic_variant,
      .y = haplo_violin,
      FUN = \(.x, .y) {
        # .x <- gse_data_haplo_variant$somatic_variant[[1]]
        # .y <- gse_data_haplo_variant$haplo_violin[[1]]

        tibble::tibble(
          variant = .x$somatic
        ) |>
          dplyr::mutate(
            pos = stringr::str_extract(variant, "\\d+") |> as.integer(),
          ) |>
          dplyr::filter(
            !pos %in% variants_tobe_excluded,
          ) -> .xx
        # .xx$variant -> heteroplasmic_variant

        .y |>
          dplyr::filter(variant %in% .xx$variant) |>
          dplyr::filter(variant_type_fisher_test == "colorful") |>
          dplyr::count(variant) |>
          dplyr::filter(n >= 3) |>
          dplyr::pull(variant) -> heteroplasmic_variant

        c(.x$high_af, .x$haplo) |> unique() -> homoplasmic_variant
        .x$heteroplasmic_variant <- heteroplasmic_variant
        .x$homoplasmic_variant <- homoplasmic_variant

        .x
      },
      mc.cores = 20,
      SIMPLIFY = FALSE
    )
  ) |>
  dplyr::mutate(
    n_heteroplasmic = purrr::map(
      heteroplasmic,
      .f = \(.x) {
        tibble::tibble(
          n_heteroplasmic = length(.x$heteroplasmic_variant),
          n_homoplasmic = length(.x$homoplasmic_variant),
        )
      }
    )
  ) |>
  tidyr::unnest(cols = n_heteroplasmic) -> gse_data_variant_heteroplasmic

export(
  gse_data_variant_heteroplasmic,
  file = file.path(
    "/home/liuc9/github/scMOCHA-data/analysis/zzz/clean-data/",
    "gse_data_variant_heteroplasmic_fisher.qs"
  ),
  format = "qs",
)


gse_data_variant_heteroplasmic$heteroplasmic |>
  purrr::map(
    .f = \(.x) {
      .x$heteroplasmic_variant
    }
  ) |>
  purrr::reduce(
    union
  ) -> heteroplasmic_variant

gse_data_variant_heteroplasmic$heteroplasmic |>
  purrr::map(
    .f = \(.x) {
      .x$homoplasmic_variant
    }
  ) |>
  purrr::reduce(
    union
  ) -> homoplasmic_variant


# gse_data_variant_heteroplasmic |>
#   dplyr::select(srrid, hetero) |>
#   tidyr::unnest(cols = hetero) |>
#   dplyr::group_by(srrid, variant) |>
#   dplyr::summarise(
#     af = mean(af, na.rm = TRUE),
#   ) |>
#   dplyr::ungroup() |>
#   dplyr::group_by(variant) |>
#   dplyr::summarise(
#     af = mean(af, na.rm = TRUE),
#   ) -> variant_mean_af

gse_data_variant_heteroplasmic |>
  dplyr::select(srrid, bulkaf) |>
  tidyr::unnest(cols = bulkaf) |>
  dplyr::group_by(variant) |>
  dplyr::summarise(
    af = mean(bulkaf, na.rm = TRUE),
  ) -> variant_mean_af

gse_data_haplo_variant |>
  dplyr::select(gseid, srrid, chemistry, haplo_variant) |>
  tidyr::unnest(cols = haplo_variant) -> all_variants

all_variants |>
  dplyr::select(Position, variant, aachange, Disease, `Gnomad Frequency`) |>
  dplyr::mutate(
    Disease = ifelse(is.na(Disease), "", Disease),
  ) |>
  dplyr::distinct() |>
  dplyr::arrange(Position) -> variant_type


all_variants |>
  dplyr::count(variant) |>
  dplyr::left_join(
    variant_type,
    by = "variant"
  ) |>
  dplyr::mutate(
    issomatic = ifelse(
      variant %in% heteroplasmic_variant,
      "heteroplasmic",
      "other"
    ),
  ) |>
  dplyr::mutate(
    issomatic = ifelse(
      variant %in% homoplasmic_variant,
      "homoplasmic",
      issomatic
    ),
  ) |>
  dplyr::arrange(
    desc(n)
  ) |>
  dplyr::group_by(Position) |>
  dplyr::mutate(
    issomatic = ifelse(
      dplyr::n() > 1,
      "multiple",
      issomatic
    )
  ) |>
  dplyr::ungroup() |>
  dplyr::left_join(
    variant_mean_af,
    by = "variant"
  ) -> variant_count


{
  export(
    variant_count,
    file = file.path(
      "/home/liuc9/github/scMOCHA-data/analysis/zzz/clean-data/",
      "all_variant_fisher.csv"
    ),
    format = "both",
    sep = ",",
    row.names = FALSE,
    col.names = TRUE,
  )
  export(
    variant_count,
    file = file.path(
      "/home/liuc9/github/scMOCHA-data/analysis/zzz/clean-data/",
      "all_variant_fisher.rds"
    )
  )
  export(
    variant_count,
    file = file.path(
      "/home/liuc9/github/scMOCHA-data/analysis/zzz/clean-data/",
      "all_variant_fisher.qs"
    )
  )
}
