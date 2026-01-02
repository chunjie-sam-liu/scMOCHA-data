#!/usr/bin/env Rscript
# Metainfo ----------------------------------------------------------------
# @AUTHOR: Chun-Jie Liu
# @CONTACT: chunjie.sam.liu.at.gmail.com
# @DATE: 2026-01-02 14:36:05
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
outdir <- path("/home/liuc9/github/scMOCHA-data/analysis/zzz/MANUSCRIPTFIGURES")
outdirnotuse <- path(
  "/home/liuc9/github/scMOCHA-data/analysis/high-res/MANUSCRIPTFIGURES-notuse"
)
scmergedir <- outdirnotuse / "scmerge"

cleandatadir <- path("/home/liuc9/github/scMOCHA-data/data/zzz/clean-data")

METAFULL <- import(outdir / "SAMPLES-METADATA-FULL.xlsx")
ALLVARIANTS <- import(
  outdir / "SAMPLE-VARIANT-CLASSIFICATION-CLUSTER-BULK-AF.xlsx"
)

HOMO_HETE_VARIANTS <- ALLVARIANTS |>
  dplyr::filter(variant_type %in% c("homo", "hete"))

# load conn ---------------------------------------------------------------

# src ---------------------------------------------------------------------

# function ----------------------------------------------------------------
load_pkg(Seurat)

fn_load_sc <- function(.filepath, thegseid, thesrrid) {
  .sc <- import(.filepath)
  sc_azimuth <- .sc$sc_azimuth
  rm(.sc)
  gc()

  sc_azimuth@meta.data |>
    as.data.table(keep.rownames = "barcode") |>
    dplyr::mutate(
      barcode_new = glue::glue("{thegseid}-{thesrrid}-{barcode}")
    ) |>
    as.data.frame() -> d_merge

  new_names <- setNames(
    d_merge$barcode_new,
    d_merge$barcode
  )

  sc_azimuth <- RenameCells(
    sc_azimuth,
    new.names = new_names
  )

  sc_azimuth@meta.data <- d_merge |>
    tibble::column_to_rownames("barcode_new")
  sc_azimuth
}

fn_norm <- function(thegseid, thesrrid) {
  # thegseid <- VARIANT_GSEID_SRRID$gseid[[1]]
  # thesrrid <- VARIANT_GSEID_SRRID$srrid[[1]]

  .dir <- path(
    "/home/liuc9/github/scMOCHA-data/data/",
    thegseid,
    "final",
    thesrrid
  )

  sc_azimuth <- fn_load_sc(
    .filepath = path(
      .dir,
      "sc_azimuth.rds.gz"
    ),
    thegseid = thegseid,
    thesrrid = thesrrid
  )

  .forintegrationdir <- .dir / "for_integration"
  dir_create(.forintegrationdir)

  export(
    sc_azimuth,
    path(
      .forintegrationdir,
      "sc_azimuth.qs"
    )
  )

  sc_azimuth |>
    Seurat::NormalizeData() |>
    Seurat::FindVariableFeatures() |>
    Seurat::ScaleData() |>
    Seurat::RunPCA() -> sc_azimuth

  export(
    sc_azimuth,
    path(
      .forintegrationdir,
      "sc_azimuth.norm.qs"
    )
  )
}

# body --------------------------------------------------------------------

HOMO_HETE_VARIANTS |>
  dplyr::select(gseid, srrid) |>
  dplyr::distinct() -> VARIANT_GSEID_SRRID


VARIANT_GSEID_SRRID |>
  dplyr::mutate(
    sct = parallel::mcmapply(
      FUN = fn_norm,
      thegseid = gseid,
      thesrrid = srrid,
      SIMPLIFY = FALSE,
      mc.cores = 20
    )
  ) -> gseid_srrid_srrdir_sct
# footer ------------------------------------------------------------------

# save image --------------------------------------------------------------
