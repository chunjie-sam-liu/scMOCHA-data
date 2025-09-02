#!/usr/bin/env Rscript
# Metainfo ----------------------------------------------------------------
# @AUTHOR: Chun-Jie Liu
# @CONTACT: chunjie.sam.liu.at.gmail.com
# @DATE: 2025-09-02 12:01:16
# @DESCRIPTION: filename
# @VERSION: v0.0.1

# Library -----------------------------------------------------------------

suppressPackageStartupMessages(library(magrittr))
library(ggplot2)
library(patchwork)
library(prismatic)
library(paletteer)
library(data.table)
# library(rlang)
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
gseid_srrid_variant_hetero_plot_ratio <- import(
  "/home/liuc9/github/scMOCHA-data/analysis/zzz/plot-real-somatic-variant/main-variants/gseid_srrid_variant_hetero_plot_ratio.qs"
)
# load conn ---------------------------------------------------------------

# src ---------------------------------------------------------------------

# function ----------------------------------------------------------------
fn_de <- function(thegseid, thesrrid, thevariant, forplot_) {
  library(Seurat)
  sc <- import(
    file.path(
      "/home/liuc9/github/scMOCHA-data/data/",
      thegseid,
      "final",
      thesrrid,
      "sc_azimuth.rds.gz"
    )
  )
  sc_azimuth <- sc$sc_azimuth
  rm(sc)
  gc()

  sc_azimuth@meta.data |>
    tibble::rownames_to_column("barcode") |>
    as.data.table() |>
    dplyr::left_join(
      forplot_ |>
        as.data.table() |>
        dplyr::mutate(barcode = as.character(barcode)),
    ) |>
    dplyr::mutate(
      barcode_new = glue::glue("{gseid}-{srrid}-{barcode}")
    ) |>
    as.data.frame() -> d_merge
  rownames(d_merge) <- d_merge$barcode_new
  sc_azimuth@meta.data <- d_merge
}

# body --------------------------------------------------------------------

thevariant <- "3727T>C"
thesrrid <- ""

gseid_srrid_variant_hetero_plot_ratio |>
  # dplyr::filter(variant == thevariant) |>
  dplyr::filter(variant %in% c("3727T>C", "3728C>T")) |>
  # dplyr::select(tidyselect::contains("ratio"))
  dplyr::select(
    gseid,
    srrid,
    variant,
    forplot,
    # tidyselect::contains("ratio")
  ) -> filtered_data


thevariant <- "3727T>C"
thesrrid <- "GSM7493836"
thegseid <- "GSE235050"

filtered_data |>
  dplyr::filter(
    variant == thevariant,
    srrid == thesrrid
  ) |>
  dplyr::select(forplot) |>
  tidyr::unnest(cols = c(forplot)) -> variant_cell_barcode

#
#
# ? load expression --------------------------------------------------------------------
#
#
library(Seurat)
sc <- import(
  file.path(
    "/home/liuc9/github/scMOCHA-data/data/",
    thegseid,
    "final",
    thesrrid,
    "sc_azimuth.rds.gz"
  )
)

sc_azimuth <- sc$sc_azimuth


# variant_cell_barcode |>
#   dplyr::slice(
#     match(variant_cell_barcode$barcode, colnames(sc_azimuth))
#   ) -> variant_cell_barcode_matched

sc_azimuth@meta.data |>
  tibble::rownames_to_column("barcode") |>
  dplyr::left_join(
    variant_cell_barcode |>
      dplyr::mutate(barcode = as.character(barcode)),
  ) |>
  as.data.frame() -> d_merge

rownames(d_merge) <- d_merge$barcode
sc_azimuth@meta.data <- d_merge


# DefaultAssay(sc_azimuth)
# sc_azimuth[["RNA"]]
# dim(sc_azimuth@assays$RNA@counts)
sc_azimuth[["RNA"]] <- as(object = sc_azimuth[["RNA"]], Class = "Assay")
markers_deseq2 <- FindMarkers(
  object = sc_azimuth,
  ident.1 = "Heteroplasmy",
  ident.2 = "Sufficient reads",
  test.use = "DESeq2",
  group.by = "cellvarianttype"
)
# min.pct = 0.1
# logfc.threshold = 0.25

dim(markers_deseq2)

markers_deseq2 |>
  tibble::rownames_to_column("gene") |>
  dplyr::arrange(p_val) -> markers_deseq2_s

markers_deseq2_s |>
  dplyr::mutate(
    fdr = -log10(p_val_adj)
  ) |>
  dplyr::mutate(
    avg_log2FC = ifelse(
      abs(avg_log2FC) > 100,
      sign(avg_log2FC) * 100,
      avg_log2FC
    )
  ) |>
  dplyr::mutate(
    color = dplyr::case_when(
      p_val_adj < 0.05 & avg_log2FC > 1 ~ "red",
      p_val_adj < 0.05 & avg_log2FC < 1 ~ "blue",
      TRUE ~ "grey"
    )
  ) -> forplot

forplot |>
  ggplot(aes(
    x = avg_log2FC,
    y = fdr,
    color = color
  )) +
  geom_point(aes()) +
  ggrepel::geom_text_repel(
    data = forplot |>
      dplyr::filter(
        (p_val_adj < 0.05 & avg_log2FC > 1) |
          (p_val_adj < 0.05 & avg_log2FC < -1)
      ) |>
      dplyr::slice(1:10),
    aes(label = gene),
    size = 3,
    max.overlaps = 20
  ) +
  scale_color_identity() +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  theme_classic() +
  labs(
    x = "Fold change(Heteroplasmy/Sufficient reads)",
    y = "FDR"
  )

# footer------------------------------------------------------------------

# future: :plan(future: :sequential)

# save image --------------------------------------------------------------
