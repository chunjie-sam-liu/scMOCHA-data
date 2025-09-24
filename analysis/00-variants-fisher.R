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

# 3. tRNA p9 and RNA editing position
.editing_pos <- c(
  585,
  1610,
  3238,
  4271,
  5520,
  7526,
  8303, # tRNA p9
  9999,
  10413,
  12146,
  12274,
  14734,
  15896, # tRNA p9
  295,
  2617,
  13710 # RNA editing
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
    haplo_variant = haplo_variant_fisher,
    haplo_violin = haplo_violin_fisher,
    somatic_variant = somatic_variant_fisher,
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

        .x |>
          tibble::enframe(name = "type", value = "variant") |>
          tidyr::unnest(cols = variant) |>
          dplyr::distinct() |>
          dplyr::mutate(
            pos = stringr::str_extract(variant, "\\d+") |> as.integer(),
          ) -> .d
        .d |>
          dplyr::filter(
            pos %in% variants_tobe_excluded
          ) -> .d_excluded
        excluding_pos <- .d_excluded$variant

        .d |>
          dplyr::filter(type == "somatic") |>
          dplyr::filter(
            !pos %in% variants_tobe_excluded,
          ) -> .xx
        .xx$variant -> heteroplasmic_variant

        .d |>
          dplyr::filter(type == "haplo") |>
          dplyr::filter(
            !pos %in% variants_tobe_excluded,
          ) -> .yy
        .yy$variant -> homoplasmic_variant

        .x$heteroplasmic_variant <- heteroplasmic_variant
        .x$homoplasmic_variant <- homoplasmic_variant
        .x$excluding_pos <- excluding_pos

        .x
      },
      mc.cores = 20,
      SIMPLIFY = FALSE
    )
  ) -> gse_data_variant_heteroplasmic_


gse_data_variant_heteroplasmic_$heteroplasmic |>
  purrr::map(
    .f = \(.x) {
      .x$homoplasmic_variant
    }
  ) |>
  purrr::reduce(
    union
  ) -> homoplasmic_variant

gse_data_variant_heteroplasmic_$heteroplasmic |>
  purrr::map(
    .f = \(.x) {
      .x$heteroplasmic_variant
    }
  ) |>
  purrr::reduce(
    union
  ) |>
  setdiff(homoplasmic_variant) -> heteroplasmic_variant


ggvenn::ggvenn(
  list(
    heteroplasmic = heteroplasmic_variant,
    homoplasmic = homoplasmic_variant
  ),
  fill_color = c("#0073C2FF", "#EFC000FF"),
  stroke_size = 0.5,
  set_name_size = 4
)

intersect(heteroplasmic_variant, homoplasmic_variant)


gse_data_variant_heteroplasmic_ |>
  dplyr::mutate(
    heteroplasmic = purrr::map(
      .x = heteroplasmic,
      .f = \(.x) {
        .x$heteroplasmic_variant <- intersect(
          .x$heteroplasmic_variant,
          heteroplasmic_variant
        )
        .x$homoplasmic_variant <- intersect(
          .x$homoplasmic_variant,
          homoplasmic_variant
        )
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
  tidyr::unnest(cols = n_heteroplasmic) -> gse_data_variant_heteroplasmic

# gse_data_variant_heteroplasmic |>
#   dplyr::select(1, 2, heteroplasmic) |>
#   dplyr::mutate(
#     a = purrr::map(heteroplasmic, .f = \(.x) {
#       .x$heteroplasmic_variant
#     })
#   ) |>
#   tidyr::unnest(cols = a) -> n
# n$a |> unique() |> length()

export(
  gse_data_variant_heteroplasmic,
  file = file.path(
    "/home/liuc9/github/scMOCHA-data/analysis/zzz/clean-data/",
    "gse_data_variant_heteroplasmic_fisher.qs"
  ),
  format = "qs",
)

# gse_data_variant_heteroplasmic <- import(
#   "/home/liuc9/github/scMOCHA-data/analysis/zzz/clean-data/gse_data_variant_heteroplasmic_fisher.qs"
# )

gse_data_variant_heteroplasmic |>
  dplyr::select(gseid, srrid, heteroplasmic) |>
  dplyr::mutate(
    a = purrr::map_chr(
      heteroplasmic,
      \(.x) {
        .x |> jsonlite::toJSON()
      }
    )
  ) |>
  dplyr::select(
    gseid,
    srrid,
    variant_alltype = a
  ) -> gseid_srrid_variant_fisher

conn <- DBI::dbConnect(
  duckdb::duckdb(),
  "/home/liuc9/github/scMOCHA-data/analysis/zzz/clean-data/all_hetero_af.cell.duckdb.1.2.1",
  read_only = FALSE
)

DBI::dbWriteTable(
  conn,
  "gseid_srrid_variant_fisher",
  gseid_srrid_variant_fisher,
  overwrite = TRUE,
  temporary = FALSE
)
# DBI::dbDisconnect(conn, shutdown = TRUE)
# DBI::dbListTables(conn)

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
  dplyr::filter(
    purrr::map_lgl(
      haplo_variant,
      .f = \(.x) {
        nrow(.x) > 0
      }
    )
  ) |>
  tidyr::unnest(cols = haplo_variant) -> all_variants

all_variants |>
  dplyr::select(Position, variant, aachange, Disease, `Gnomad Frequency`) |>
  dplyr::mutate(
    Disease = ifelse(is.na(Disease), "", Disease),
  ) |>
  dplyr::distinct() |>
  dplyr::arrange(Position) -> variant_type

# all_variants |> dplyr::filter(variant == "3010G>A")

all_variants |>
  dplyr::count(variant) |>
  dplyr::left_join(
    variant_type,
    by = "variant"
  ) |>
  dplyr::mutate(
    issomatic = ifelse(
      variant %in% homoplasmic_variant,
      "homoplasmic",
      "other"
    ),
  ) |>
  dplyr::mutate(
    issomatic = ifelse(
      variant %in% heteroplasmic_variant,
      "heteroplasmic",
      issomatic
    ),
  ) |>
  dplyr::arrange(
    desc(n)
  ) |>
  # dplyr::group_by(Position) |>
  # dplyr::filter(dplyr::n() > 1) |> dplyr::arrange(variant)
  # dplyr::mutate(
  #   issomatic = ifelse(
  #     dplyr::n() > 1,
  #     "multiple",
  #     issomatic
  #   )
  # ) |>
  # dplyr::ungroup() |>
  dplyr::left_join(
    variant_mean_af,
    by = "variant"
  ) -> variant_count
variant_count |> dplyr::count(issomatic)

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

conn <- DBI::dbConnect(
  duckdb::duckdb(),
  "/home/liuc9/github/scMOCHA-data/analysis/zzz/clean-data/all_hetero_af.cell.duckdb.1.2.1"
)
DBI::dbListTables(conn)
dplyr::tbl(conn, "allvariants")
variant_count
DBI::dbWriteTable(
  conn,
  "allvariants_fisher",
  variant_count,
  overwrite = TRUE
)
DBI::dbDisconnect(conn, shutdown = TRUE)
