#!/usr/bin/env Rscript
# Metainfo ----------------------------------------------------------------
# @AUTHOR: Chun-Jie Liu
# @CONTACT: chunjie.sam.liu.at.gmail.com
# @DATE: 2025-05-22 12:39:40
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

# src ---------------------------------------------------------------------

# header ------------------------------------------------------------------

# future: :plan(future: :multisession, workers = 10)

# function ----------------------------------------------------------------

# load data ---------------------------------------------------------------

# load data ---------------------------------------------------------------
cleandatadir <- "/home/liuc9/github/scMOCHA-data/data/zzz/clean-data"

gse_data <- import(
  "/home/liuc9/github/scMOCHA-data/analysis/zzz/clean-data/gse_data.qs"
)

sex_pred <- import(
  "/home/liuc9/github/scMOCHA-data/analysis/zzz/clean-data/gse_srrid_srrdir_sex.qs"
) |>
  dplyr::select(
    srrid,
    Sex = sex
  )


METADATA <- import(
  file.path(cleandatadir, "gse_dataset_metadata_full.qs")
) |>
  dplyr::filter(
    disease %in% c("Alzheimer's Disease", "Healthy"),
    Chemistry == "SC5P-PE"
  ) |>
  dplyr::mutate(
    Haplogroup_s = purrr::map_chr(
      .x = Haplogroup,
      .f = \(.x) {
        # if (stringr::str_starts(.x, "L")) {
        #   gsub("L", "L0", .x)
        # }
        gsub("\\d+.*", "", .x)
      }
    )
  ) |>
  dplyr::left_join(
    sex_pred,
    by = "srrid"
  )

METADATA |>
  dplyr::pull(srrid) -> srrids

all_variant <- import(
  file.path(cleandatadir, "all_variant.fst")
)

all_heteroplasmic_af <- import(
  file.path(cleandatadir, "all_hetero_af.cluster.fst")
)

all_variant |>
  dplyr::filter(issomatic == "heteroplasmic") |>
  dplyr::arrange(Position) -> heteroplasmic

gse_data |>
  dplyr::filter(srrid %in% srrids) |>
  dplyr::select(srrid, haplo_variant) |>
  tidyr::unnest(cols = haplo_variant) |>
  dplyr::filter(variant %in% heteroplasmic$variant) -> gse_variant_het

VARIANTS <- gse_variant_het$variant |>
  sort() |>
  unique()
length(VARIANTS)

# ! body --------------------------------------------------------------------

source("/home/liuc9/github/scMOCHA-data/analysis/00-colors.R")

all_heteroplasmic_af |>
  # dtplyr::lazy_dt() |>
  dplyr::filter(
    num_variants > 0,
    srrid %in% srrids
  ) -> all_heteroplasmic_af_1
dim(all_heteroplasmic_af_1)

all_heteroplasmic_af_1 |>
  dplyr::select(dplyr::any_of(VARIANTS)) |>
  as.matrix() -> all_heteroplasmic_af_1_mat


all_heteroplasmic_af_1 |>
  dplyr::select(gseid, srrid, barcode) |>
  dplyr::mutate(
    colname = paste0(gseid, "_", srrid, "_", barcode)
  ) |>
  dplyr::mutate(
    barcode = gsub("_", " ", barcode)
  ) |>
  dplyr::mutate(
    Cluster = factor(barcode, levels = names(color_celltype))
  ) |>
  dplyr::left_join(
    METADATA |>
      dplyr::select(
        gseid,
        srrid,
        Chemistry,
        Sex,
        Age_new,
        disease,
        Haplogroup_s
      ),
    by = c("gseid", "srrid")
  ) |>
  dplyr::select(-c(gseid, srrid, barcode)) |>
  dplyr::mutate(
    Sex = factor(Sex, levels = names(color_gender))
  ) |>
  dplyr::rename(
    Age = Age_new,
    Disease = disease,
    Haplogroup = Haplogroup_s
  ) -> .af_cluster_before


.af_cluster_before |>
  tibble::column_to_rownames(var = "colname") -> .af_cluster

suppressPackageStartupMessages(library(ComplexHeatmap))
library(circlize)
pcc <- import(
  file = "https://raw.githubusercontent.com/chunjie-sam-liu/chunjie-sam-liu.life/master/public/data/pcc.tsv"
) |>
  dplyr::arrange(cancer_types)


colSums(all_heteroplasmic_af_1_mat) |>
  as.data.frame() |>
  tibble::rownames_to_column(var = "variant") |>
  dplyr::select(variant, freq = `colSums(all_heteroplasmic_af_1_mat)`) |>
  dplyr::arrange(desc(freq)) -> sort_variants

dplyr::bind_rows(
  dplyr::slice_head(sort_variants, n = 10)
) |>
  dplyr::distinct() |>
  dplyr::filter(variant != "8545G>A") |>
  dplyr::left_join(
    heteroplasmic,
    by = "variant"
  ) |>
  dplyr::mutate(
    label = glue::glue("{variant};{aachange}\n{Disease}")
  ) -> top_variants


.af_mtx <- all_heteroplasmic_af_1_mat |> t()
colnames(.af_mtx) <- .af_cluster$colname
dim(.af_mtx)


# chm top color setting
Cluster_ <- levels(.af_cluster$Cluster)
CLUSTER_ <- color_celltype
names(CLUSTER_) <- Cluster_
Chemistry_ <- levels(.af_cluster$Chemistry)[c(1)]
CHEMISTRY_ <- color_chemistry[c(1)]
names(CHEMISTRY_) <- Chemistry_
Sex_ <- levels(.af_cluster$Sex)[c(1, 2)]
SEX_ <- color_gender[c(1, 2)]
names(SEX_) <- Sex_
Disease_ <- levels(.af_cluster$Disease)[c(1, 3)]
DISEASE_ <- color_disease[c(1, 3)]
names(DISEASE_) <- Disease_
Haplogroup_ <- sort(unique(.af_cluster$Haplogroup))
HAPLOGROUP_ <- pcc$color[1:length(Haplogroup_)]
names(HAPLOGROUP_) <- Haplogroup_

chm_top <- ComplexHeatmap::HeatmapAnnotation(
  df = .af_cluster,
  # gap = unit(c(2, 2), "mm"),
  col = list(
    Cluster = CLUSTER_,
    Chemistry = CHEMISTRY_,
    Sex = SEX_,
    Disease = DISEASE_,
    Haplogroup = HAPLOGROUP_,
    `Age` = circlize::colorRamp2(
      breaks = c(2, 92),
      colors = c("white", "red"),
      space = "RGB"
    )
  ),
  which = "column"
)

hma_index_left <- match(top_variants$variant, rownames(.af_mtx))
hma_left <- ComplexHeatmap::rowAnnotation(
  link = anno_mark(
    at = hma_index_left,
    labels = top_variants$label,
    which = "row",
    side = "left",
    lines_gp = gpar(
      lwd = 0.5,
      col = "black"
    ),
    labels_gp = gpar(
      fontsize = 10,
      col = "black"
    ),
    padding = unit(0.5, "mm"),
    link_width = unit(5, "mm"),
  )
)

sort(unique(as.numeric(.af_mtx))) |> quantile(seq(0, 1, by = 0.01)) -> .seq_af
col_option <- "turbo"
col_pick = c(
  "white",
  viridis::viridis_pal(
    alpha = 1,
    begin = 0,
    end = 1,
    direction = 1,
    option = col_option
  )(
    # n_break
    length(.seq_af) - 1
  )
)
col_fun = circlize::colorRamp2(
  .seq_af,
  col_pick
)

# ! cluster --------------------------------------------------------------------

ComplexHeatmap::Heatmap(
  matrix = .af_mtx,
  col = col_fun,
  name = "Allele Freq",
  na_col = "white",
  color_space = "LAB",
  rect_gp = gpar(col = NA),
  border = NA,
  border_gp = gpar(col = NA),
  cell_fun = NULL,
  layer_fun = NULL,
  jitter = FALSE,
  # row
  cluster_rows = TRUE,
  cluster_row_slices = TRUE,
  clustering_distance_rows = "pearson",
  clustering_method_rows = "ward.D",
  show_row_dend = FALSE,
  show_row_names = FALSE,
  row_dend_reorder = TRUE,
  row_dend_gp = gpar(),
  row_split = 2,
  row_title = NULL,
  # column
  cluster_columns = TRUE,
  cluster_column_slices = TRUE,
  clustering_distance_columns = "pearson",
  clustering_method_columns = "ward.D",
  # column_title_gp = gpar(fontsize = 40),
  show_column_dend = FALSE,
  show_column_names = FALSE,
  row_names_side = "left",
  column_split = 3,
  column_title = NULL,

  # annotation
  top_annotation = chm_top,
  left_annotation = hma_left,
  heatmap_legend_param = list(
    title = "Allele Freq",
    at = c(0, 0.5, 1),
    labels = c("0", "0.5", "1"),
    legend_direction = "vertical",
    title_gp = gpar(fontsize = 10)
  )
) -> ch_af

{
  pdf(
    file = "/home/liuc9/github/scMOCHA-data/analysis/zzz/plot-heteroplasmic/heatmap_cluster_af_SC5P-PE_AD.pdf",
    width = 20,
    height = 9
  )
  ComplexHeatmap::draw(object = ch_af)
  dev.off()
}

# footer ------------------------------------------------------------------

# future: :plan(future: :sequential)

# save image --------------------------------------------------------------
