#!/usr/bin/env Rscript
# Metainfo ----------------------------------------------------------------
# @AUTHOR: Chun-Jie Liu
# @CONTACT: chunjie.sam.liu.at.gmail.com
# @DATE: 2025-12-16 16:13:47
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

pcc <- import(
  file = "https://raw.githubusercontent.com/chunjie-sam-liu/chunjie-sam-liu.life/master/public/data/pcc.tsv"
) |>
  dplyr::arrange(cancer_types)

allvariants <- import(
  path(Sys.getenv("HIGHRESDIR"), "MANUSCRIPTFIGURES/SAMPLE-VARIANT-CLASSIFICATION-CLUSTER-BULK-AF.xlsx")
) |>
  dplyr::mutate(
    coord = parallel::mclapply(
      X = variant,
      FUN = \(.v) {
        # .v <- gse_data_variant_classification_clusteraf_bulkaf$variant[[1]]
        pos <- gsub(
          pattern = "(\\d+)([A-Z])>([A-Z])",
          replacement = "\\1",
          x = .v
        ) |>
          as.integer()
        ref <- gsub(
          pattern = "(\\d+)([A-Z])>([A-Z])",
          replacement = "\\2",
          x = .v
        )
        alt <- gsub(
          pattern = "(\\d+)([A-Z])>([A-Z])",
          replacement = "\\3",
          x = .v
        )
        data.table(
          seqnames = "MT",
          start = pos,
          end = pos,
          ref = ref,
          alt = alt
        )
      },
      mc.cores = 10
    )
  ) |>
  tidyr::unnest(
    cols = coord
  )

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

clusterlevel_dir <- path(
  Sys.getenv("ZZZDIR"), "new-variant-cell/homo-hete/clusterlevel"
)


METADATA <- import(
  path(
    Sys.getenv("CLEANDATADIR"),
    "gse_dataset_metadata_full.qs"
  )
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
  dplyr::mutate(Sex = SEXPRED)

source(path(
  Sys.getenv("HIGHRESDIR"),
  "00-colors.R"
))

# load conn ---------------------------------------------------------------

# src ---------------------------------------------------------------------

# function ----------------------------------------------------------------

# body --------------------------------------------------------------------
\() {
  # only run once.
  variant_clusteraf <- dir_ls(clusterlevel_dir) |>
    # head(20) |>
    purrr::map(
      .f = \(.x) {
        import(.x)
      }
    ) |>
    purrr::reduce(
      dplyr::left_join,
      by = c("gseid", "srrid", "celltype")
    )

  lobstr::obj_size(variant_clusteraf)

  # Convert list columns to character strings to avoid writexl warnings
  variant_clusteraf |>
    # dplyr::select(1, 2, 3, 4, 5) |>
    dplyr::mutate_all(.funs = \(x) {
      if (is.list(x)) as.double(x) else x
    }) -> variant_clusteraf_unlist

  export(
    variant_clusteraf_unlist,
    path(Sys.getenv("HIGHRESDIR"), "MANUSCRIPTFIGURES/ALLVARIANT-ALLSAMPLES-CLUSTERAF.xlsx")
  )
}
variant_clusteraf_unlist <- import(
  path(Sys.getenv("HIGHRESDIR"), "MANUSCRIPTFIGURES/ALLVARIANT-ALLSAMPLES-CLUSTERAF.xlsx")
)
#
#
# ? plot --------------------------------------------------------------------
#
#

variant_clusteraf_unlist[1:5, 1:6]
variant_clusteraf_unlist |>
  dplyr::select(c(gseid, srrid, barcode = celltype)) |>
  dplyr::mutate(
    colname = paste0(gseid, "_", srrid, "_", barcode)
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
  ) |>
  as.data.frame() |>
  tibble::column_to_rownames(var = "colname") -> PLOTCOLMETADATA

variant_clusteraf_unlist |>
  dplyr::select(-c(gseid, srrid, celltype)) |>
  as.matrix() -> PLOTDATAMATRIX
rownames(PLOTDATAMATRIX) <- rownames(PLOTCOLMETADATA)


colMeans(PLOTDATAMATRIX) |>
  as.data.table(keep.rownames = TRUE) |>
  dplyr::select(
    variant = V1,
    AVGAF = V2
  ) |>
  dplyr::left_join(VARIANT_TYPE, by = "variant") |>
  as.data.frame() |>
  tibble::column_to_rownames(var = "variant") -> PLOTROWDATA

load_pkg(ComplexHeatmap, circlize)

AFMTX <- PLOTDATAMATRIX |> t()


# chm top color setting
Cluster_ <- levels(PLOTCOLMETADATA$Cluster)
CLUSTER_ <- color_celltype
names(CLUSTER_) <- Cluster_
Chemistry_ <- levels(PLOTCOLMETADATA$Chemistry)
CHEMISTRY_ <- color_chemistry
names(CHEMISTRY_) <- Chemistry_
Sex_ <- levels(PLOTCOLMETADATA$Sex)[c(1, 2)]
SEX_ <- color_gender[c(1, 2)]
names(SEX_) <- Sex_
Disease_ <- levels(PLOTCOLMETADATA$Disease)
DISEASE_ <- color_disease
names(DISEASE_) <- Disease_
Haplogroup_ <- sort(unique(PLOTCOLMETADATA$Haplogroup))
HAPLOGROUP_ <- pcc$color[1:length(Haplogroup_)]
names(HAPLOGROUP_) <- Haplogroup_

CHM_TOP <- ComplexHeatmap::HeatmapAnnotation(
  df = PLOTCOLMETADATA,
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

PLOTROWDATA |>
  dplyr::arrange(Homoplasmic, Heteroplasmic, Somatic) |>
  rownames() -> .rownames_sorted

AFMTX <- AFMTX[.rownames_sorted, ]
PLOTROWDATA <- PLOTROWDATA[.rownames_sorted, ]

CHM_LEFT <- ComplexHeatmap::rowAnnotation(
  df = PLOTROWDATA |> dplyr::select(-AVGAF),
  col = list(
    Homoplasmic = c("0" = "white", "1" = "#ae00ff"),
    Heteroplasmic = c("0" = "white", "1" = "#00b3ff"),
    Somatic = c("0" = "white", "1" = "#FF0000")
  )
)


sort(unique(as.numeric(AFMTX))) |> quantile(seq(0, 1, by = 0.01)) -> SEQAF
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
    length(SEQAF) - 1
  )
)
col_fun = circlize::colorRamp2(
  # seq(col_start,
  #   col_end,
  #   length.out = n_break
  # ),
  SEQAF,
  col_pick
)


{
  ComplexHeatmap::Heatmap(
    matrix = AFMTX,
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
    # cluster_rows = FALSE, # disable row clustering
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
    top_annotation = CHM_TOP,
    left_annotation = CHM_LEFT,
    # right_annotation = hma_right,
    heatmap_legend_param = list(
      title = "Allele Freq",
      at = SEQAF[closest_positions],
      labels = c("0", "0.01", "0.05", "0.1", "0.5", "1"),
      legend_direction = "vertical",
      title_gp = gpar(fontsize = 10)
    )
  ) -> COMPLEXHEATMAP_AF

  OUTDIR <- path(\n    Sys.getenv(\"HIGHRESDIR\"), \"MANUSCRIPTFIGURES\"\n  )
  pdf(
    file = OUTDIR / "HEATMAP-COL-CELLTYPE-ROW-VARIANT-CLUSTERROW.pdf",
    width = 22,
    height = 9
  )
  ComplexHeatmap::draw(object = COMPLEXHEATMAP_AF)
  dev.off()
}


{
  library(ggplot2)

  legend_df <- data.frame(
    value = SEQAF,
    col = col_fun(SEQAF)
  ) |>
    tibble::rowid_to_column(var = "id")

  # Find closest SEQAF positions for desired labels
  target_values <- c(0, 0.01, 0.05, 0.1, 0.5, 1.0)
  closest_positions <- sapply(target_values, function(x) {
    which.min(abs(SEQAF - x))
  })

  p_legend <- ggplot(
    legend_df,
    aes(y = id, x = 1, fill = col)
  ) +
    geom_tile(width = 0.8, height = 1) +
    scale_fill_identity() +
    scale_y_continuous(
      breaks = closest_positions,
      labels = c("0", "0.01", "0.05", "0.1", "0.5", "1"),
      expand = c(0, 0),
      position = "right"
    ) +
    scale_x_continuous(expand = c(0, 0)) +
    coord_cartesian(xlim = c(0.6, 1.4)) +
    theme_void() +
    theme(
      axis.text.y.right = element_text(
        size = 10,
        hjust = 0,
        margin = margin(l = 5),
        color = "black"
      ),
      axis.title.y.right = element_text(
        size = 12,
        angle = 270,
        vjust = 0.5,
        margin = margin(l = 10),
        color = "black"
      ),
      plot.margin = margin(5, 25, 5, 5),
      panel.border = element_rect(
        color = "black",
        fill = NA,
        size = 0.5
      )
    ) +
    labs(y = "Allele Frequency")

  # Save legend as PDF
  OUTDIR <- path(
    Sys.getenv("HIGHRESDIR"), "MANUSCRIPTFIGURES"
  )
  ggsave(
    filename = OUTDIR /
      "HEATMAP-COL-CELLTYPE-ROW-VARIANT-CLUSTERROW-LEGEND.pdf",
    plot = p_legend,
    width = 2,
    height = 6,
    units = "in",
    dpi = 300
  )

  print(p_legend)
  cli_alert_success("Legend saved successfully.")
}
# footer ------------------------------------------------------------------

# save image --------------------------------------------------------------
