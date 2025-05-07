regions_missalignment_error <- c(66:71, 300:316, 513:525, 3106:3107, 12418:12425, 16182:16194)
regions_rare_heteroplasmic_variants <- c(499, 538, 545, 10953, 12684)
variants_tobe_excluded <- c(
  regions_missalignment_error,
  regions_rare_heteroplasmic_variants
)

gse_data <- readr::read_rds(
  "/home/liuc9/github/scMOCHA-data/stats/stats/zzz/clean-data/gse_data.rds"
)


gse_data |>
  dplyr::select(gseid, srrid, chemistry, anno, hetero, haplo_variant, haplo_violin, somatic_variant, celltype_ratio) |>
  dplyr::left_join(
    gse_dataset_metadata_full |> dplyr::select(-gseid),
    by = c("srrid" = "srrid")
  ) |>
  dplyr::arrange(disease, Chemistry) ->
gse_data_haplo_variant

gse_data_haplo_variant |>
  dplyr::mutate(
    heteroplasmic = purrr::map(
      somatic_variant,
      .f = \(.x) {
        tibble::tibble(
          variant = .x$somatic
        ) |>
          dplyr::mutate(
            pos = stringr::str_extract(variant, "\\d+") |> as.integer(),
          ) |>
          dplyr::filter(
            !pos %in% variants_tobe_excluded,
          ) ->
        .xx
        .xx$variant -> heteroplasmic_variant
        c(.x$high_af, .x$haplo) |> unique() -> homoplasmic_variant
        .x$heteroplasmic_variant <- heteroplasmic_variant
        .x$homoplasmic_variant <- homoplasmic_variant

        .x
      }
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
  tidyr::unnest(cols = n_heteroplasmic) ->
gse_data_variant_heteroplasmic

gse_data_variant_heteroplasmic$heteroplasmic |>
  purrr::map(
    .f = \(.x) {
      .x$heteroplasmic_variant
    }
  ) |>
  purrr::reduce(
    union
  ) ->
heteroplasmic_variant

gse_data_variant_heteroplasmic$heteroplasmic |>
  purrr::map(
    .f = \(.x) {
      .x$homoplasmic_variant
    }
  ) |>
  purrr::reduce(
    union
  ) ->
homoplasmic_variant


gse_data_variant_heteroplasmic |>
  dplyr::select(srrid, hetero) |>
  tidyr::unnest(cols = hetero) |>
  dplyr::group_by(srrid, variant) |>
  dplyr::summarise(
    af = mean(af, na.rm = TRUE),
  ) |>
  dplyr::ungroup() |>
  dplyr::group_by(variant) |>
  dplyr::summarise(
    af = mean(af, na.rm = TRUE),
  ) ->
variant_mean_af

gse_data_haplo_variant |>
  dplyr::select(gseid, srrid, chemistry, haplo_variant) |>
  tidyr::unnest(cols = haplo_variant) ->
all_variants

all_variants |>
  dplyr::select(Position, variant, aachange, Disease, `Gnomad Frequency`) |>
  dplyr::mutate(
    Disease = ifelse(is.na(Disease), "", Disease),
  ) |>
  dplyr::distinct() |>
  dplyr::arrange(Position) ->
variant_type


all_variants |>
  dplyr::count(variant) |>
  dplyr::left_join(
    variant_type,
    by = "variant"
  ) |>
  dplyr::mutate(
    issomatic = ifelse(variant %in% heteroplasmic_variant, "heteroplasmic", "other"),
  ) |>
  dplyr::mutate(
    issomatic = ifelse(variant %in% homoplasmic_variant, "homoplasmic", issomatic),
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
  ) ->
variant_count


{
  data.table::fwrite(
    variant_count,
    file = file.path(
      "/home/liuc9/github/scMOCHA-data/stats/stats/zzz/clean-data/", "all_variant.csv"
    ),
    sep = ",",
    row.names = FALSE,
    col.names = TRUE,
  )
  readr::write_rds(
    variant_count,
    file = file.path(
      "/home/liuc9/github/scMOCHA-data/stats/stats/zzz/clean-data/", "all_variant.rds"
    )
  )
}
