gwas <- import("/mnt/isilon/xing_lab/liuc9/refdata/gwas/gwas_catalog_v1.0.2-associations_e107_r2022-10-08-ftp-efo.fst")

thedisease <- " Alzheimer's disease"

# gwas |>
#   dplyr::count(MAPPED_TRAIT) |>
#   dplyr::filter(grepl("Alzheimer disease", MAPPED_TRAIT))

gwas |> dplyr::glimpse()
# gwas |> dplyr::count(CHR_ID) |> dplyr::pull(CHR_ID) |> unique() |> sort()
# gwas  |>
#   dplyr::filter(CHR_ID == )
gwas |>
  dplyr::filter(grepl(thedisease, `DISEASE/TRAIT`)) ->
gwas_ad


gwas_ad |>
  dplyr::filter(filesize > 0) |>
  dplyr::glimpse()
gwas_ad |>
  dplyr::count(CHR_ID)

gwas_ad |>
  dplyr::filter(
    !CHR_ID %in% c(1:22, "X", "Y")
  )


h <- import("/mnt/isilon/xing_lab/liuc9/refdata/gwas/summary_statistics/34737426-GCST90042678-EFO_0009268.h.tsv.gz")

h |> dplyr::count(chromosome)
