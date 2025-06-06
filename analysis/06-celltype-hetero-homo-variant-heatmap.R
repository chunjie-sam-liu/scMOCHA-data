#!/usr/bin/env Rscript
# Metainfo ----------------------------------------------------------------
# @AUTHOR: Chun-Jie Liu
# @CONTACT: chunjie.sam.liu.at.gmail.com
# @DATE: 2025-04-24 22:03:34
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
log_threshold(TRACE)
log_layout(layout_glue_colors)

# future: :plan(future: :multisession, workers = 10)

# function ----------------------------------------------------------------


# load data ---------------------------------------------------------------
cleandatadir <- "/home/liuc9/github/scMOCHA-data/data/zzz/clean-data"

pcc <- import(file = "https://raw.githubusercontent.com/chunjie-sam-liu/chunjie-sam-liu.life/master/public/data/pcc.tsv") |>
  dplyr::arrange(cancer_types)

all_homo_af <- import(file.path(cleandatadir, "all_hetero_af.cluster.fst"))
all_heteroplasmic_af <- import(
  file.path(cleandatadir, "all_homo_af.cluster.fst"),
) |>
  dplyr::left_join(all_homo_af, by = c("gseid", "srrid", "barcode")) |>
  dplyr::mutate(
    num_variants = num_variants.x + num_variants.y
  ) |>
  dplyr::select(-c(num_variants.x, num_variants.y))

# all_homo_af <- data.table::fread(
#   file.path(cleandatadir, "all_homo_af.cluster.csv"),
#   header = TRUE,
#   sep = ","
# )


all_variant <- import(
  file.path(cleandatadir, "all_variant.qs")
)

all_variant |> dplyr::count(issomatic)

all_variant |>
  dplyr::filter(issomatic %in% c("homoplasmic", "heteroplasmic")) |>
  dplyr::arrange(Position) ->
heteroplasmic

heteroplasmic |>
  dplyr::filter(Disease != "") |>
  dplyr::arrange(-n) ->
heteroplasmic_disease


# #
# all_heteroplasmic_af |>
#   dplyr::select(1, 2, 3, `11467A>G`, `15301G>A`) |>
#   dplyr::mutate(
#     colname = glue::glue("{gseid}_{srrid}_{barcode}"),
#   ) |>
#   dplyr::select(-c(gseid, srrid)) |>
#   tidyr::pivot_longer(
#     cols = -c(colname, barcode),
#     names_to = "variant",
#     values_to = "af"
#   ) ->
# forplot
# forplot |>
#   dplyr::filter(variant == "11467A>G") |>
#   ggplot(aes(
#     x = barcode,
#     y = af,
#   )) +
#   geom_point() +
#   theme(
#     axis.text = element_blank(),
#     axis.ticks = element_blank(),
#   )
# forplot |>
#   ggstatsplot::ggbetweenstats(
#     x = barcode,
#     y = af,
#     xlab = "Cluster",
#     ylab = "Allele Frequency",
#     # title = glue::glue("{the_srrid} - {one_variant}"),
#     pairwise.display = "significant",
#     p.adjust.method = "BH"
#   ) +
#   scale_y_continuous(
#     limits = c(0, 1)
#   )



VARIANTS <- heteroplasmic$variant
length(VARIANTS)

CELLBARCODE <- import(
  file.path(cleandatadir, "barcode_celltype.fst")
)

sex_pred <- import("/home/liuc9/github/scMOCHA-data/analysis/zzz/clean-data/gse_srrid_srrdir_sex.qs") |>
  dplyr::select(
    srrid,
    Sex = sex
  )

METADATA <- import(
  file.path(cleandatadir, "gse_dataset_metadata_full.qs")
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
# METADATA |>
#   # dplyr::count(Chemistry)
#   dplyr::filter(Chemistry == "SC5P-PE") |>
#   dplyr::pull(srrid) ->
# srrids

# body --------------------------------------------------------------------

source("/home/liuc9/github/scMOCHA-data/analysis/00-colors.R")
all_heteroplasmic_af |>
  # dtplyr::lazy_dt() |>
  dplyr::filter(
    num_variants > 0,
    # srrid %in% srrids
  ) ->
all_heteroplasmic_af_1
nrow(all_heteroplasmic_af_1)

all_heteroplasmic_af_1 |>
  dplyr::select(dplyr::any_of(VARIANTS)) |>
  as.matrix() ->
all_heteroplasmic_af_1_mat

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
    METADATA |> dplyr::select(gseid, srrid, Chemistry, Sex, Age_new, disease, Haplogroup_s),
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
  ) ->
.af_cluster_before

.af_cluster_before |>
  tibble::column_to_rownames(var = "colname") ->
.af_cluster



suppressPackageStartupMessages(library(ComplexHeatmap))
library(circlize)

colSums(all_heteroplasmic_af_1_mat) |>
  as.data.frame() |>
  tibble::rownames_to_column(var = "variant") |>
  dplyr::select(variant, freq = `colSums(all_heteroplasmic_af_1_mat)`) |>
  dplyr::arrange(desc(freq)) ->
sort_variants

dplyr::bind_rows(
  dplyr::slice_head(sort_variants, n = 10),
  dplyr::slice_tail(sort_variants, n = 10)
  # sort_variants |>
  #   dplyr::filter(variant %in% heteroplasmic_disease$variant) |>
  #   head(10)
) |>
  dplyr::distinct() |>
  dplyr::filter(variant != "8545G>A") |>
  dplyr::left_join(
    heteroplasmic,
    by = "variant"
  ) |>
  dplyr::mutate(
    label = glue::glue("{variant};{aachange}\n{Disease}")
  ) ->
top_variants



.af_mtx <- all_heteroplasmic_af_1_mat |> t()
colnames(.af_mtx) <- .af_cluster$colname
dim(.af_mtx)

# chm top color setting
Cluster_ <- levels(.af_cluster$Cluster)
CLUSTER_ <- color_celltype
names(CLUSTER_) <- Cluster_
Chemistry_ <- levels(.af_cluster$Chemistry)
CHEMISTRY_ <- color_chemistry
names(CHEMISTRY_) <- Chemistry_
Sex_ <- levels(.af_cluster$Sex)[c(1, 2)]
SEX_ <- color_gender[c(1, 2)]
names(SEX_) <- Sex_
Disease_ <- levels(.af_cluster$Disease)
DISEASE_ <- color_disease
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

hma_index_right <- match(top_variants$variant, rownames(.af_mtx))
hma_left <- ComplexHeatmap::rowAnnotation(
  link = anno_mark(
    at = hma_index_right,
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
# seq(0, 1, 0.01) -> .seq_af
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
  # seq(col_start,
  #   col_end,
  #   length.out = n_break
  # ),
  .seq_af,
  col_pick
)


ComplexHeatmap::Heatmap(
  matrix = .af_mtx,
  # col = circlize::colorRamp2(
  #   breaks = c(-0.1, 0, 1),
  #   colors = c("lightgrey", "gold", "blue"),
  #   space = "RGB"
  # ),
  col = col_fun,
  name = "Allele Freq",
  na_col = "white",
  color_space = "LAB",
  rect_gp = gpar(col = NA),
  border = NA,
  cell_fun = NULL,
  layer_fun = NULL,
  jitter = FALSE,
  # row
  cluster_rows = TRUE,
  cluster_row_slices = T,
  show_row_names = FALSE,
  show_row_dend = FALSE,
  clustering_distance_rows = "pearson",
  clustering_method_rows = "ward.D",
  # row_names_gp = gpar(
  #   # fontsize = 20,
  #   col = .gcol$cell_variants
  # ),
  # column
  # column_title = paste0("palette = '", col_option, "'"),
  # column_title = .column_title,
  column_title_gp = gpar(fontsize = 40),
  cluster_columns = TRUE,
  # cluster_columns = cluster_within_group(
  #   mat = .af_mtx,
  #   factor = .af_cluster_before |>
  #     dplyr::select(colname, Cluster) |>
  #     tibble::deframe()
  # ),
  cluster_column_slices = T,
  show_column_dend = FALSE,
  clustering_distance_columns = "pearson",
  clustering_method_columns = "ward.D",
  show_column_names = FALSE,
  row_names_side = "left",
  top_annotation = chm_top,
  # left_annotation = hma_left,
  # right_annotation = hma_right,
  heatmap_legend_param = list(
    title = "Allele Freq",
    at = c(0, 0.5, 1),
    labels = c("0", "0.5", "1"),
    legend_direction = "vertical",
    title_gp = gpar(fontsize = 10)
  )
) ->
ch_af
# ch_af


{
  pdf(
    file = "/home/liuc9/github/scMOCHA-data/analysis/zzz/plot-heteroplasmic/heatmap_cluster_af.homo+hetero.pdf",
    width = 22,
    height = 9
  )
  ComplexHeatmap::draw(object = ch_af)
  dev.off()
}

# footer ------------------------------------------------------------------

# future: :plan(future: :sequential)

# save image --------------------------------------------------------------
