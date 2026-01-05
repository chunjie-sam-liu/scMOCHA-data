#!/usr/bin/env Rscript
# Metainfo ----------------------------------------------------------------
# @AUTHOR: Chun-Jie Liu
# @CONTACT: chunjie.sam.liu.at.gmail.com
# @DATE: 2026-01-01 15:42:50
# @DESCRIPTION: this script is used for ...

# Library -----------------------------------------------------------------

load_pkg(jutils)

dotenv(".env")

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

outdir <- path(Sys.getenv("OUTDIR"))
outdirnotuse <- path(
  Sys.getenv("OUTDIRNOTUSE")
)
cleandatadir <- path(Sys.getenv("CLEANDATADIR"))

METAFULL <- import(outdir / "SAMPLES-METADATA-FULL.xlsx")

ALLVARIANTS_TEST <- import(
  outdir / "VARIANT-KRUSKAL-WALLIS-TEST.xlsx"
)

# load conn ---------------------------------------------------------------

# src ---------------------------------------------------------------------
source(
  path(Sys.getenv(\"HIGHRESDIR\"), \"plot_celltype_specific_variant.R\")
)
source(
  path(Sys.getenv(\"HIGHRESDIR\"), \"plot_individual_proportion.R\")
)

# function ----------------------------------------------------------------

# body --------------------------------------------------------------------
ALLVARIANTS_TEST |>
  # dplyr::filter(p.value < 0.05) |>
  dplyr::select(variant, gseid, srrid) |>
  tidyr::nest(.by = variant, .key = "gsesrrid") -> ALLVARIANTS_TEST_SIG


ALLVARIANTS_TEST_SIG |>
  dplyr::mutate(
    m = parallel::mcmapply(
      function(thevariant, gsesrrid) {
        # thevariant <- ALLVARIANTS_TEST_SIG$variant[1]
        # gsesrrid <- ALLVARIANTS_TEST_SIG$gsesrrid[[1]]
        # p_joy <- fn_plot_joy(
        #   thevariant = thevariant,
        #   thegseid = thegseid,
        #   thesrrid = thesrrid
        # )
        pdf(
          file = outdirnotuse /
            "celltype-specific" /
            glue("VARIANT-{thevariant}-CELLTYPE-SPECIFIC-JOY-PLOT.pdf"),
          width = 11,
          height = 6
        )
        gsesrrid |>
          dplyr::mutate(
            p_joy = purrr::pmap(
              .l = list(
                gsesrrid$gseid,
                gsesrrid$srrid,
                rep(thevariant, nrow(gsesrrid))
              ),
              .f = \(thegseid, thesrrid, thevariant) {
                print(fn_plot_joy(thevariant, thegseid, thesrrid))
              }
            )
          ) -> p_joy
        dev.off()

        # p_hist <- fn_plot_hist(
        #   thevariant = thevariant,
        #   thegseid = thegseid,
        #   thesrrid = thesrrid
        # )

        pdf(
          file = outdirnotuse /
            "celltype-specific" /
            glue("VARIANT-{thevariant}-CELLTYPE-SPECIFIC-HIST-PLOT.pdf"),
          width = 11,
          height = 6
        )
        gsesrrid |>
          dplyr::mutate(
            p_hist = purrr::pmap(
              .l = list(
                gsesrrid$gseid,
                gsesrrid$srrid,
                rep(thevariant, nrow(gsesrrid))
              ),
              .f = \(thegseid, thesrrid, thevariant) {
                print(fn_plot_hist(thevariant, thegseid, thesrrid))
              }
            )
          ) -> p_hist
        dev.off()

        # p_cumfrac <- fn_plot_cumulative_fraction(
        #   thevariant = thevariant,
        #   thegseid = thegseid,
        #   thesrrid = thesrrid
        # )

        pdf(
          file = outdirnotuse /
            "celltype-specific" /
            glue("VARIANT-{thevariant}-CELLTYPE-SPECIFIC-CUMFRAC-PLOT.pdf"),
          width = 11,
          height = 6
        )
        gsesrrid |>
          dplyr::mutate(
            p_cumfrac = purrr::pmap(
              .l = list(
                gsesrrid$gseid,
                gsesrrid$srrid,
                rep(thevariant, nrow(gsesrrid))
              ),
              .f = \(thegseid, thesrrid, thevariant) {
                print(fn_plot_cumulative_fraction(
                  thevariant,
                  thegseid,
                  thesrrid
                ))
              }
            )
          ) -> p_cumfrac
        dev.off()

        # p_celltype_detail <- fn_plot_joy_celltype_detail(
        #   thevariant = thevariant,
        #   thegseid = thegseid,
        #   thesrrid = thesrrid
        # )
        pdf(
          file = outdirnotuse /
            "celltype-specific" /
            glue("VARIANT-{thevariant}-CELLTYPE-SPECIFIC-DETAIL-PLOT.pdf"),
          width = 20,
          height = 12
        )
        gsesrrid |>
          dplyr::mutate(
            p_detail = purrr::pmap(
              .l = list(
                gsesrrid$gseid,
                gsesrrid$srrid,
                rep(thevariant, nrow(gsesrrid))
              ),
              .f = \(thegseid, thesrrid, thevariant) {
                print(fn_plot_joy_celltype_detail(
                  thevariant,
                  thegseid,
                  thesrrid
                ))
              }
            )
          ) -> p_detail
        dev.off()

        pdf(
          file = outdirnotuse /
            "celltype-specific" /
            glue("VARIANT-{thevariant}-CELLTYPE-SPECIFIC-PSEUDO-BULK-PLOT.pdf"),
          width = 15,
          height = 8
        )
        print(fn_plot_hetero_pseudo_bulk(thevariant))
        dev.off()

        pdf(
          file = outdirnotuse /
            "celltype-specific" /
            glue(
              "VARIANT-{thevariant}-CELLTYPE-SPECIFIC-PSEUDO-BULK-PROPORTION-PLOT.pdf"
            ),
          width = 15,
          height = 8
        )
        print(fn_plot_variant_ratio(thevariant))
        dev.off()
      },
      thevariant = variant,
      gsesrrid = gsesrrid,
      mc.cores = 20
    )
  )

# footer ------------------------------------------------------------------

# save image --------------------------------------------------------------
