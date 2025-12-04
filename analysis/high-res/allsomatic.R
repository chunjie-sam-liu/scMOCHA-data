source("/home/liuc9/github/scMOCHA-data/analysis/high-res/plot_somatic.R")
gse_data_variant_classification_clusteraf_bulkaf <- import(
  "/home/liuc9/github/scMOCHA-data/analysis/high-res/MANUSCRIPTFIGURES/SAMPLE-VARIANT-CLASSIFICATION-CLUSTER-BULK-AF.xlsx"
)


gse_data_variant_classification_clusteraf_bulkaf |>
  dplyr::filter(
    variant_type == "somatic"
  ) |>
  dplyr::select(gseid, srrid, variant) |>
  dplyr::distinct() |>
  tidyr::nest(
    .by = c(gseid, srrid),
    .key = "variant"
  ) -> somatic_variants_per_sample

\() {
  outdir <- "/home/liuc9/github/scMOCHA-data/analysis/high-res/MANUSCRIPTFIGURES/notuse/somatic"
  cli_progress_bar(
    "Plotting somatic variants per sample",
    total = nrow(somatic_variants_per_sample)
  )
  somatic_variants_per_sample |>
    # head(4) |>
    dplyr::mutate(
      plot = parallel::mcmapply(
        FUN = \(gseid, srrid, variant) {
          # gseid <- somatic_variants_per_sample$gseid[[1]]
          # srrid <- somatic_variants_per_sample$srrid[[1]]
          # variant <- somatic_variants_per_sample$variant[[1]]
          cli_alert_info(glue::glue(
            "Processing sample {gseid}-{srrid} with {nrow(variant)} somatic variants"
          ))

          # cli_progress_update()

          pdf(
            path(
              outdir,
              glue::glue("Somatic-{nrow(variant)}-{gseid}-{srrid}.pdf")
            ),
            width = 16,
            height = 8
          )
          variant$variant |>
            purrr::walk(
              \(.v) {
                cli_alert_info(glue::glue(
                  "Plotting variant {.v} in {gseid}-{srrid}"
                ))
                p <- fn_plot_somatic(
                  thevariant = .v,
                  thesrrid = srrid
                )
                print(p)
              }
            )
          dev.off()
          cli_alert_success(glue::glue(
            "Finished plotting sample {gseid}-{srrid}"
          ))
          1
        },
        gseid,
        srrid,
        variant,
        SIMPLIFY = FALSE,
        mc.cores = 20
      )
    )
  cli_progress_done()
}
