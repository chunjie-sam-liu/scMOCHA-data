load_pkg(jutils)
dotenv(".env")
outdir <- path(Sys.getenv("OUTDIR"))
ALLVARIANTS <- import(
  outdir / "SAMPLE-VARIANT-CLASSIFICATION-CLUSTER-BULK-AF.xlsx"
)


HOMOHETE <- ALLVARIANTS |>
  dplyr::filter(
    variant_type %in% c("homo", "hete")
  ) |>
  mutate(
    variant_type = ifelse(
      variant_type == "homo",
      "homoplasmic",
      "heteroplasmic"
    )
  )

HAPLO <- ALLVARIANTS |>
  dplyr::filter(
    variant_type %in% c("haplo")
  ) |>
  select(gseid, srrid, variant, ishaplo = variant_type) |>
  dplyr::mutate(ishaplo = TRUE)

SOMATIC <- ALLVARIANTS |>
  dplyr::filter(
    variant_type %in% c("somatic")
  ) |>
  select(gseid, srrid, variant, issomatic = variant_type) |>
  dplyr::mutate(issomatic = TRUE)


HOMOHETE |>
  left_join(
    SOMATIC,
    by = c("gseid", "srrid", "variant")
  ) |>
  dplyr::mutate(
    issomatic = ifelse(is.na(issomatic), FALSE, issomatic)
  ) |>
  left_join(
    HAPLO,
    by = c("gseid", "srrid", "variant")
  ) |>
  dplyr::mutate(
    ishaplo = ifelse(is.na(ishaplo), FALSE, ishaplo)
  ) |>
  select(
    `GSE ID` = gseid,
    `Sample ID` = srrid,
    Haplogroup,
    `Verbose Haplogroup` = Verbose_haplogroup,
    `Variant type` = variant_type,
    `Is Somatic Variant` = issomatic,
    `Is Haplogroup Variant` = ishaplo,
    `Variant` = variant,
    B,
    CD4_T,
    CD8_T,
    DC,
    Mono,
    NK,
    other,
    other_T,
    Bulk
  ) |>
  export(
    outdir / "SAMPLE-VARIANT-CLASSIFICATION-CLUSTER-BULK-AF-CLEAN.xlsx"
  )
