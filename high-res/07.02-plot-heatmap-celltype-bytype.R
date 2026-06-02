#!/usr/bin/env Rscript
# Metainfo ----------------------------------------------------------------
# @AUTHOR: Chun-Jie Liu
# @CONTACT: chunjie.sam.liu.at.gmail.com
# @DATE: 2026-06-02
# @DESCRIPTION: Heatmaps for all variants + per-type (Somatic/Heteroplasmic/Homoplasmic)
#               AF < 0.05 shown as white; columns clustered within Haplogroup groups.
# @VERSION: v0.0.2

# Library -----------------------------------------------------------------

load_pkg(jutils)
dotenv(".env")

GetoptLong.options(help_style = "two-column")
verbose = TRUE
GetoptLong("verbose!", "print messages")

logger::log_threshold(logger::TRACE)
logger::log_layout(logger::layout_glue_colors)

load_pkg(ComplexHeatmap, circlize)

# load data ---------------------------------------------------------------

pcc <- import(
  file = "https://raw.githubusercontent.com/chunjie-sam-liu/chunjie-sam-liu.life/master/public/data/pcc.tsv"
) |>
  dplyr::arrange(cancer_types)

allvariants <- import(
  path(Sys.getenv("OUTDIR")) /
    "SAMPLE-VARIANT-CLASSIFICATION-CLUSTER-BULK-AF.xlsx"
)

# VARIANT_TYPE wide -------------------------------------------------------

allvariants |>
  dplyr::select(variant, variant_type) |>
  dplyr::distinct() |>
  dplyr::filter(variant_type %in% c("homo", "hete", "somatic")) |>
  dplyr::mutate(
    variant_type = plyr::revalue(
      variant_type,
      replace = c(
        "homo" = "Homoplasmic",
        "hete" = "Heteroplasmic",
        "somatic" = "Somatic"
      )
    )
  ) |>
  dplyr::mutate(presence = 1) |>
  tidyr::pivot_wider(
    names_from = variant_type,
    values_from = presence,
    values_fill = 0
  ) -> VARIANT_TYPE

# VARIANT_TYPE long -------------------------------------------------------

allvariants |>
  dplyr::select(variant, variant_type) |>
  dplyr::distinct() |>
  dplyr::filter(variant_type %in% c("homo", "hete", "somatic")) |>
  dplyr::mutate(
    variant_type = plyr::revalue(
      variant_type,
      replace = c(
        "homo" = "Homoplasmic",
        "hete" = "Heteroplasmic",
        "somatic" = "Somatic"
      )
    )
  ) -> VARIANT_TYPE_LONG

# METADATA and colors -----------------------------------------------------

METADATA <- import(
  path(Sys.getenv("CLEANDATADIR"), "gse_dataset_metadata_full.qs"),
  lazy = FALSE
) |>
  dplyr::mutate(
    Haplogroup_s = purrr::map_chr(.x = Haplogroup, .f = \(.x) {
      gsub("\\d+.*", "", .x)
    })
  ) |>
  dplyr::mutate(Sex = SEXPRED)

source(path(Sys.getenv("HIGHRESDIR"), "00-colors.R"))

# load cluster AF data ----------------------------------------------------

variant_clusteraf_unlist <- import(
  path(Sys.getenv("OUTDIR")) / "ALLVARIANT-ALLSAMPLES-CLUSTERAF.qs"
) |>
  filter(celltype != "Bulk")

# build PLOTCOLMETADATA (no Chemistry, no Disease) ------------------------

variant_clusteraf_unlist |>
  dplyr::select(c(gseid, srrid, barcode = celltype)) |>
  dplyr::mutate(colname = paste0(gseid, "_", srrid, "_", barcode)) |>
  dplyr::mutate(Cluster = factor(barcode, levels = names(color_celltype))) |>
  dplyr::left_join(
    METADATA |> dplyr::select(gseid, srrid, Sex, Age_new, Haplogroup_s),
    by = c("gseid", "srrid")
  ) |>
  dplyr::select(-c(gseid, srrid, barcode)) |>
  dplyr::mutate(Sex = factor(Sex, levels = names(color_gender))) |>
  dplyr::rename(Age = Age_new, Haplogroup = Haplogroup_s) |>
  as.data.frame() |>
  tibble::column_to_rownames(var = "colname") -> PLOTCOLMETADATA

# build matrix and PLOTROWDATA --------------------------------------------

variant_clusteraf_unlist |>
  dplyr::select(-c(gseid, srrid, celltype)) |>
  as.matrix() -> PLOTDATAMATRIX
rownames(PLOTDATAMATRIX) <- rownames(PLOTCOLMETADATA)

AFMTX_ALL <- PLOTDATAMATRIX |> t()

colMeans(PLOTDATAMATRIX) |>
  as.data.table(keep.rownames = TRUE) |>
  dplyr::select(variant = V1, AVGAF = V2) |>
  dplyr::left_join(VARIANT_TYPE, by = "variant") |>
  as.data.frame() |>
  tibble::column_to_rownames(var = "variant") -> PLOTROWDATA

PLOTROWDATA |>
  dplyr::arrange(desc(Homoplasmic), desc(Heteroplasmic), desc(Somatic)) |>
  rownames() -> .rownames_sorted
AFMTX_ALL <- AFMTX_ALL[.rownames_sorted, ]
PLOTROWDATA <- PLOTROWDATA[.rownames_sorted, ]

# top annotation (no Chemistry, no Disease) --------------------------------

Haplogroup_ <- sort(unique(PLOTCOLMETADATA$Haplogroup))
HAPLOGROUP_ <- pcc$color[seq_along(Haplogroup_)]
names(HAPLOGROUP_) <- Haplogroup_

CHM_TOP <- ComplexHeatmap::HeatmapAnnotation(
  df = PLOTCOLMETADATA,
  col = list(
    Cluster = color_celltype[levels(PLOTCOLMETADATA$Cluster)],
    Sex = color_gender[levels(PLOTCOLMETADATA$Sex)[c(1, 2)]],
    Haplogroup = HAPLOGROUP_,
    Age = circlize::colorRamp2(
      breaks = c(2, 92),
      colors = c("white", "red"),
      space = "RGB"
    )
  ),
  which = "column"
)

# left annotation for all-variants plot ------------------------------------

CHM_LEFT_ALL <- ComplexHeatmap::rowAnnotation(
  df = PLOTROWDATA |> dplyr::select(-AVGAF),
  col = list(
    Homoplasmic = c("0" = "white", "1" = "#ae00ff"),
    Heteroplasmic = c("0" = "white", "1" = "#00b3ff"),
    Somatic = c("0" = "white", "1" = "#FF0000")
  )
)

# haplogroup factor for cluster_within_group --------------------------------

haplogroup_factor <- setNames(
  PLOTCOLMETADATA$Haplogroup,
  rownames(PLOTCOLMETADATA)
)

# color scale: AF < 0.05 = white -------------------------------------------

col_option <- "turbo"
n_turbo <- 99
col_breaks_05white <- c(0, 0.0499, seq(0.05, 1, length.out = n_turbo))
col_colors_05white <- c(
  "white",
  "white",
  viridis::viridis_pal(
    alpha = 1,
    begin = 0,
    end = 1,
    direction = 1,
    option = col_option
  )(n_turbo)
)
col_fun_05white <- circlize::colorRamp2(col_breaks_05white, col_colors_05white)
target_values_05 <- c(0, 0.05, 0.1, 0.5, 1.0)
legend_positions_05 <- sapply(target_values_05, \(x) {
  which.min(abs(col_breaks_05white - x))
})

# output directory ---------------------------------------------------------

OUTDIR_HEATMAP <- path(Sys.getenv("OUTDIR"), "heatmap_0602")
fs::dir_create(OUTDIR_HEATMAP)

# helper: draw and save one heatmap ----------------------------------------

.draw_heatmap <- function(
  mat,
  row_label,
  filename,
  left_annotation = NULL,
  width = 22,
  height = 9
) {
  set.seed(123)
  # Build within-Haplogroup column order manually (avoids cluster_within_group bugs)
  hap_fac <- haplogroup_factor[colnames(mat)]
  hap_fac[is.na(hap_fac)] <- "Unknown"
  col_ord <- integer(0)
  for (.hg in sort(unique(hap_fac))) {
    .idx <- which(hap_fac == .hg)
    if (length(.idx) == 1) {
      col_ord <- c(col_ord, .idx)
    } else {
      .hc <- hclust(dist(t(mat[, .idx, drop = FALSE])), method = "ward.D2")
      col_ord <- c(col_ord, .idx[.hc$order])
    }
  }

  # Build within-type row order: Homoplasmic -> Heteroplasmic -> Somatic -> Other
  .rowdata <- PLOTROWDATA[rownames(mat), , drop = FALSE]
  .rowdata$priority_group <- dplyr::case_when(
    !is.na(.rowdata$Homoplasmic) & .rowdata$Homoplasmic == 1 ~ "1_Homoplasmic",
    !is.na(.rowdata$Heteroplasmic) & .rowdata$Heteroplasmic == 1 ~ "2_Heteroplasmic",
    !is.na(.rowdata$Somatic) & .rowdata$Somatic == 1 ~ "3_Somatic",
    TRUE ~ "4_Other"
  )
  row_ord <- integer(0)
  for (.grp in sort(unique(.rowdata$priority_group))) {
    .ridx <- which(.rowdata$priority_group == .grp)
    if (length(.ridx) == 1) {
      row_ord <- c(row_ord, .ridx)
    } else {
      .hc <- hclust(dist(mat[.ridx, , drop = FALSE]), method = "ward.D2")
      row_ord <- c(row_ord, .ridx[.hc$order])
    }
  }

  ComplexHeatmap::Heatmap(
    matrix = mat,
    col = col_fun_05white,
    name = "Allele Freq",
    na_col = "white",
    rect_gp = gpar(col = NA),
    border = NA,
    # rows: pre-computed within-type order (Homoplasmic -> Heteroplasmic -> Somatic)
    cluster_rows = FALSE,
    row_order = row_ord,
    show_row_names = nrow(mat) <= 80,
    row_names_gp = gpar(fontsize = 7),
    show_row_dend = FALSE,
    row_title = row_label,
    row_title_gp = gpar(fontsize = 12, fontface = "bold"),
    # columns: pre-computed within-Haplogroup order
    cluster_columns = FALSE,
    column_order = col_ord,
    show_column_dend = FALSE,
    show_column_names = FALSE,
    row_names_side = "left",
    top_annotation = CHM_TOP,
    left_annotation = left_annotation,
    heatmap_legend_param = list(
      title = "Allele Freq",
      at = col_breaks_05white[legend_positions_05],
      labels = c("0", "0.05", "0.1", "0.5", "1"),
      legend_direction = "vertical",
      title_gp = gpar(fontsize = 10)
    )
  ) -> chm

  pdf(file = OUTDIR_HEATMAP / filename, width = width, height = height)
  ComplexHeatmap::draw(object = chm)
  dev.off()
  logger::log_info("Saved: {filename}")
}

# body: all-variants heatmap -----------------------------------------------

.draw_heatmap(
  mat = AFMTX_ALL,
  row_label = paste0("All variants (n=", nrow(AFMTX_ALL), ")"),
  filename = "HEATMAP-COL-CELLTYPE-ROW-ALLVARIANT-AF05WHITE.pdf",
  left_annotation = CHM_LEFT_ALL
)

# body: per-type heatmaps --------------------------------------------------

for (.vtype in c("Somatic", "Heteroplasmic", "Homoplasmic")) {
  logger::log_info("Processing variant type: {.vtype}")

  .variants <- VARIANT_TYPE_LONG |>
    dplyr::filter(variant_type == .vtype) |>
    dplyr::pull(variant)
  .variants <- intersect(.variants, rownames(AFMTX_ALL))

  if (length(.variants) == 0) {
    logger::log_warn("No variants found for {.vtype}, skipping.")
    next
  }

  .mat <- AFMTX_ALL[.variants, , drop = FALSE]

  .draw_heatmap(
    mat = .mat,
    row_label = paste0(.vtype, " variants (n=", length(.variants), ")"),
    filename = paste0(
      "HEATMAP-COL-CELLTYPE-ROW-",
      toupper(.vtype),
      "-AF05WHITE.pdf"
    ),
    left_annotation = NULL
  )
}

logger::log_info("All heatmaps done.")
