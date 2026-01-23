#!/usr/bin/env Rscript
# Metainfo ----------------------------------------------------------------
# @AUTHOR: Chun-Jie Liu
# @CONTACT: chunjie.sam.liu.at.gmail.com
# @DATE: 2025-12-02 14:09:56
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
dotenv(".env")
basedir <- path(Sys.getenv("BASEDIR"))
outdir <- path(Sys.getenv("OUTDIR"))
cleandatadir <- path(Sys.getenv("CLEANDATADIR"))
gse_data <- import(
  cleandatadir / "gse_data.qs"
)

gse_dataset_metadata_full <- import(
  cleandatadir / "gse_dataset_metadata_full.qs"
)
# load conn ---------------------------------------------------------------

# src ---------------------------------------------------------------------

# function ----------------------------------------------------------------

# body --------------------------------------------------------------------
gse_data |>
  dplyr::select(
    gseid,
    srrid,
    chemistry,
    anno,
    hetero,
    haplo_variant,
    haplo_violin,
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
    heteroplasmic = purrr::map(
      somatic_variant,
      .f = \(.x) {
        # .x <- gse_data_haplo_variant$somatic_variant[[1]]
        # tibble::tibble(
        #   variant = .x$somatic
        # ) |>
        #   dplyr::mutate(
        #     pos = stringr::str_extract(variant, "\\d+") |> as.integer(),
        #   ) |>
        #   dplyr::filter(
        #     !pos %in% variants_tobe_excluded,
        #   ) -> .xx
        # .xx$variant -> heteroplasmic_variant
        # c(.x$high_af, .x$haplo) |> unique() -> homoplasmic_variant
        .x$heteroplasmic_variant <- .x$hete
        .x$homoplasmic_variant <- .x$homo

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


export(
  gse_data_variant_heteroplasmic,
  file = cleandatadir / "gse_data_variant_heteroplasmic.qs",
  format = "qs",
)


gse_data_variant_heteroplasmic |>
  select(anno) |>
  unnest(cols = anno) |>
  select(
    variant,
    Position,
    Ref,
    Alt,
    Locus,
    Disease,
    Status,
    ntchange,
    aachange,
    # `Mitomap Frequency`,
    `Gnomad Frequency`
  ) |>

  mutate(
    Disease = ifelse(Disease == "", NA_character_, Disease),
    Status = ifelse(Status == "", NA_character_, Status),
  ) |>
  distinct() |>
  dplyr::filter(
    !(variant == "1382A>C" & Locus == "12S,MOTS-C")
  ) |>
  dplyr::filter(
    !(variant == "9055G>A" & Status == "Reported [B*]")
  ) |>
  dplyr::filter(
    !(variant == "12811T>C" & Status == "Reported [B*]")
  ) |>
  arrange(Position) -> variant_annotation_table

{
  variant_annotation_table |>
    export(
      file = outdir / "VARIANT-ANNOTATION-TABLE.xlsx",
    )
}

#
#
# ? variant classification --------------------------------------------------------------------
#
#

gse_data_variant_heteroplasmic |>
  dplyr::select(gseid, srrid, anno, heteroplasmic) |>
  dplyr::mutate(
    hap = purrr::map(
      .x = anno,
      .f = \(.x) {
        # .x <- gse_data_variant_heteroplasmic$anno[[1]]

        .x |>
          dplyr::select(Haplogroup, Verbose_haplogroup) |>
          dplyr::distinct() |>
          dplyr::filter(Haplogroup != "") -> .hap_info
        if (nrow(.hap_info) == 0) {
          tibble::tibble(
            Haplogroup = NA_character_,
            Verbose_haplogroup = NA_character_
          )
        } else {
          .hap_info
        }
      }
    )
  ) |>
  tidyr::unnest(cols = hap) |>
  dplyr::mutate(
    hhs = purrr::map(
      .x = heteroplasmic,
      .f = \(.x) {
        # .x <- gse_data_variant_heteroplasmic$heteroplasmic[[1]]

        .x |>
          tibble::enframe(
            name = "variant_type",
            value = "variant"
          ) |>
          dplyr::filter(
            !variant_type %in% c("heteroplasmic_variant", "homoplasmic_variant")
          ) |>
          dplyr::mutate(
            nvariant = purrr::map_int(
              variant,
              .f = \(.v) {
                length(.v)
              }
            )
          )
      }
    )
  ) |>
  dplyr::select(-c(anno, heteroplasmic)) |>
  tidyr::unnest(cols = hhs) -> gse_data_variant_classification

gse_data_variant_classification |>
  dplyr::select(-variant) |>
  tidyr::pivot_wider(
    names_from = variant_type,
    values_from = nvariant
  ) -> gse_data_variant_classification_wide

{
  export(
    gse_data_variant_classification_wide,
    file = outdir / "SAMPLE-VARIANT-CLASSIFICATION-COUNT.xlsx",
  )
}


#
#
# ? clusteraf and bulkaf --------------------------------------------------------------------
#
#

gse_data_variant_heteroplasmic |>
  # dplyr::select(gseid, srrid, clusteraf, bulkaf) |>
  dplyr::mutate(
    variant_cluster_bulk_af = parallel::mcmapply(
      FUN = \(.clusteraf, .bulkaf) {
        # .clusteraf <- gse_data$clusteraf[[1]]
        # .bulkaf <- gse_data$bulkaf[[1]]

        .clusteraf |>
          dplyr::select(variant, celltype, clusteraf) |>
          tidyr::pivot_wider(
            names_from = celltype,
            values_from = clusteraf
          ) -> .clusteraf_
        .bulkaf |> dplyr::select(variant, Bulk = bulkaf) -> .bulkaf_
        dplyr::left_join(
          .clusteraf_,
          .bulkaf_,
          by = "variant"
        )
      },
      .clusteraf = clusteraf,
      .bulkaf = bulkaf,
      SIMPLIFY = FALSE,
      mc.cores = 10
    )
  ) |>
  dplyr::select(gseid, srrid, variant_cluster_bulk_af) |>
  tidyr::unnest(cols = variant_cluster_bulk_af) -> gse_data_clusteraf_bulkaf

gse_data_variant_classification |>
  dplyr::filter(variant_type %in% c("homo", "hete", "somatic", "haplo")) |>
  dplyr::select(-nvariant) |>
  tidyr::unnest(cols = variant) |>
  dplyr::left_join(
    gse_data_clusteraf_bulkaf,
    by = c("gseid", "srrid", "variant" = "variant")
  ) -> gse_data_variant_classification_clusteraf_bulkaf


{
  export(
    gse_data_variant_classification_clusteraf_bulkaf,
    file = outdir / "SAMPLE-VARIANT-CLASSIFICATION-CLUSTER-BULK-AF.xlsx",
  )
}

#
# footer ------------------------------------------------------------------

# save image --------------------------------------------------------------
