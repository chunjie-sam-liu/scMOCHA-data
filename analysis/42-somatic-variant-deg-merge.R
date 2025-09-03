#!/usr/bin/env Rscript
# Metainfo ----------------------------------------------------------------
# @AUTHOR: Chun-Jie Liu
# @CONTACT: chunjie.sam.liu.at.gmail.com
# @DATE: 2025-09-03 10:43:09
# @DESCRIPTION: filename
# @VERSION: v0.0.1

# Library -----------------------------------------------------------------

suppressPackageStartupMessages(library(magrittr))
library(ggplot2)
library(patchwork)
library(prismatic)
library(paletteer)
library(data.table)
#library(rlang)
library(glue)
library(parallel)
library(GetoptLong)
library(logger)

# args --------------------------------------------------------------------

# s: string, i: integer, f: float, !: boolean
# @: array
# %: hash
# default: default value specified here.
verbose <- FALSE
spec <- "
Usage: Rscript foorbar.R [options]
Options:

<verbose!> Print messages
"

GetoptLong.options(help_style = "two-column")
GetoptLong(spec, template_control = list(opt_width = 21))

# header ------------------------------------------------------------------

# future: :plan(future: :multisession, workers = 10)

# load data ---------------------------------------------------------------
gseid_srrid_variant <- import(
  "/home/liuc9/github/scMOCHA-data/analysis/zzz/plot-real-somatic-variant/main-variants/gseid_srrid_variant.fst"
)
# load conn ---------------------------------------------------------------

# src ---------------------------------------------------------------------

# function ----------------------------------------------------------------

# body --------------------------------------------------------------------

gseid_srrid_variant |>
  dplyr::filter(
    variant %in% c("3727T>C", "3728C>T")
  ) |>
  dplyr::mutate(
    sc_file = file.path(
      "/home/liuc9/github/scMOCHA-data/data/",
      gseid,
      "final",
      srrid,
      "de",
      "sc_azimuth.sct.qs"
    )
  ) |>
  dplyr::mutate(
    file_exists = file.exists(sc_file)
  ) -> gseid_srrid_variant_sc


#
#
# ? 3727T>C --------------------------------------------------------------------
#
#
library(Seurat)
gseid_srrid_variant_sc |>
  dplyr::filter(variant == "3727T>C") |>
  dplyr::mutate(
    load = parallel::mclapply(
      sc_file,
      function(f) {
        .sc <- import(f)
        .sc[["SCT"]]@scale.data <- matrix()
        .sc
      },
      mc.cores = 10
    )
  ) -> gseid_srrid_variant_sc_filtered


sc_list <- gseid_srrid_variant_sc_filtered |> dplyr::pull(load)

lapply(
  sc_list,
  Seurat::VariableFeatures
) |>
  unlist() |>
  unique() -> var_features

sc_merge <- merge(
  x = sc_list[[1]],
  y = sc_list[2:length(sc_list)],
  merge.data = FALSE # not merge the scale.data, for memory sake
)

Seurat::VariableFeatures(sc_merge) <- var_features

DefaultAssay(sc_merge)
Assays(sc_merge)
Layers(sc_merge[["SCT"]])

sc_merge <- Seurat::PrepSCTFindMarkers(
  sc_merge,
  # features = Seurat::VariableFeatures(sc_merge)
)

export(
  sc_merge,
  "/home/liuc9/github/scMOCHA-data/analysis/zzz/plot-real-somatic-variant/main-variants/sc_merge.sct.3727T>C.qs"
)

markers <- Seurat::FindMarkers(
  object = sc_merge,
  ident.1 = "Heteroplasmy",
  ident.2 = "Sufficient reads",
  assay = "SCT",
  slot = "data",
  test.use = "wilcox",
  group.by = "cellvarianttype",
  latent.vars = "srrid",
  features = Seurat::VariableFeatures(sc_merge)
)
export(
  markers,
  "/home/liuc9/github/scMOCHA-data/analysis/zzz/plot-real-somatic-variant/main-variants/markers.hetero_vs_sufficient.3727T>C.qs"
)

#
#
# ? 3728C>T --------------------------------------------------------------------
#
#
gseid_srrid_variant_sc |>
  dplyr::filter(variant == "3728C>T") |>
  dplyr::mutate(
    load = parallel::mclapply(
      sc_file,
      function(f) {
        .sc <- import(f)
        .sc[["SCT"]]@scale.data <- matrix()
        .sc
      },
      mc.cores = 10,
    )
  ) -> gseid_srrid_variant_sc_filtered


sc_list <- gseid_srrid_variant_sc_filtered |> dplyr::pull(load)

lapply(
  sc_list,
  Seurat::VariableFeatures
) |>
  unlist() |>
  unique() -> var_features

sc_merge <- merge(
  x = sc_list[[1]],
  y = sc_list[2:length(sc_list)],
  merge.data = FALSE # not merge the scale.data, for memory sake
)

Seurat::VariableFeatures(sc_merge) <- var_features

DefaultAssay(sc_merge)
Assays(sc_merge)
Layers(sc_merge[["SCT"]])

sc_merge <- Seurat::PrepSCTFindMarkers(
  sc_merge,
  # features = Seurat::VariableFeatures(sc_merge)
)

export(
  sc_merge,
  "/home/liuc9/github/scMOCHA-data/analysis/zzz/plot-real-somatic-variant/main-variants/sc_merge.sct.3728C>T.qs"
)

markers <- Seurat::FindMarkers(
  object = sc_merge,
  ident.1 = "Heteroplasmy",
  ident.2 = "Sufficient reads",
  assay = "SCT",
  slot = "data",
  test.use = "wilcox",
  group.by = "cellvarianttype",
  latent.vars = "srrid",
  features = Seurat::VariableFeatures(sc_merge)
)
export(
  markers,
  "/home/liuc9/github/scMOCHA-data/analysis/zzz/plot-real-somatic-variant/main-variants/markers.hetero_vs_sufficient.3728C>T.qs"
)

# footer ------------------------------------------------------------------

# future: :plan(future: :sequential)

# save image --------------------------------------------------------------
