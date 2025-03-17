#!/usr/bin/env Rscript
# Metainfo ----------------------------------------------------------------
# @AUTHOR: Chun-Jie Liu
# @CONTACT: chunjie.sam.liu.at.gmail.com
# @DATE: 2025-03-07 16:21:52
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

fn_load_hetero <- function(.filename) {
  # .filename <- file.path(
  #   sc_dir,
  #   "mgatk_out/final",
  #   "sc.cell_heteroplasmic_df.tsv.gz"
  # )

  data.table::fread(input = .filename) -> .d

  data.table::setnames(.d, "V1", "barcode")
  .d <- data.table::melt(.d, id.vars = "barcode", variable.name = "variant", value.name = "af")
  .d[, pos := gsub(pattern = ">|[AGCT]", replacement = "", x = variant)]
  .d[, pos := as.integer(pos)]
  .d
}

fn_load_cluster <- function(.filename) {
  data.table::fread(
    input = .filename,
    sep = "\t",
    col.names = c("barcode", "tag", "celltype")
  ) |>
    dplyr::arrange(celltype) |>
    dplyr::mutate(celltype = factor(celltype)) |>
    dplyr::select(-tag)
}

fn_load_coverage <- function(.filename) {
  data.table::fread(
    input = .filename,
    sep = ",",
    col.names = c("pos", "barcode", "depth")
  )
}

fn_load_meta <- function(.filename) {
  data.table::fread(
    input = .filename,
    sep = "\t"
  ) |>
    dplyr::rename(
      barcode = cellbarcode
    ) |>
    dplyr::select(-orig.ident)
}

fn_forplot <- function(.af, .coverage, .meta) {
  # print(.meta)
  .af <- as.data.table(.af)

  # Ensure .af is a data.table
  setDT(.af)

  # Get columns containing ">"
  variant_cols <- grep(">", names(.af), value = TRUE)

  # Melt data.table to long format and summarize in one chain
  .rank <- melt(.af,
    id.vars = c("barcode", "cluster"),
    measure.vars = variant_cols,
    variable.name = "variant",
    value.name = "af"
  )[, .(s_af = sum(af, na.rm = TRUE)),
    by = .(barcode, cluster)
  ]

  # Sort by cluster and -s_af (modifies in place)
  setorder(.rank, cluster, -s_af)

  # Select barcode and variant columns, then melt to long format
  variant_cols <- grep(">", names(.af), value = TRUE)
  .forplot <- melt(.af[, c("barcode", variant_cols), with = FALSE],
    id.vars = "barcode",
    measure.vars = variant_cols,
    variable.name = "variant",
    value.name = "af"
  )

  # Extract position from variant
  .forplot[, pos := as.numeric(gsub(pattern = "([[:digit:]]*).*", "\\1", variant))]

  # Convert coverage to data.table if not already
  setDT(.coverage)

  # Perform left join with coverage
  .forplot <- merge(.forplot, .coverage, by = c("barcode", "pos"), all.x = TRUE)

  # Handle NAs and apply conditions on af
  .forplot[is.na(af), af := 0]
  .forplot[is.na(depth), af := NA]
  .forplot[depth < 10, af := -0.1]

  # Sort by position
  setorder(.forplot, pos)

  .coverage |>
    dplyr::group_by(barcode) |>
    dplyr::summarise(sum_depth = sum(depth, na.rm = TRUE)) ->
  .coverage_cell

  list(
    rank = .rank,
    forplot = .forplot,
    meta = .meta,
    coverage_cell = .coverage_cell
  )
}

# load data ---------------------------------------------------------------

basedir <- "/home/liuc9/github/scMOCHA-data/data"
foundation_out <- file.path(basedir, "scfoundation/out")

gse_dataset_metadata_full <- readr::read_rds(
  file.path(foundation_out, "gse_dataset_metadata_full.rds")
)
gseids <- c(
  "GSE155673",
  "GSE157344",
  "GSE149689",
  "GSE171555",
  "GSE155223",
  "GSE163668",
  "GSE175524",
  "GSE206283",
  "GSE226598",
  "GSE261140",
  "GSE279945",
  "GSE214865",
  "GSE220189",
  "GSE233844",
  "GSE175499",
  "GSE149313",
  "GSE154386",
  "GSE159117",
  "GSE188632",
  "GSE166992",
  "GSE162117",
  "GSE226602",
  "GSE161354",
  "GSE235050",
  "GSE181279",
  # scfoundation2
  "GSE143353",
  "GSE148215",
  "GSE163314",
  "GSE163633",
  "GSE164690",
  "GSE167825",
  "GSE174125",
  "GSE184703",
  "GSE153421",
  "GSE147794",
  "GSE168453"
)


# body --------------------------------------------------------------------
tibble::tibble(
  gseid = gseids
) |>
  dplyr::mutate(
    anno = purrr::map(
      .x = gseid,
      .f = \(.gseid) {
        .anno <- readr::read_rds(
          file.path(basedir, .gseid, "out", glue::glue("{.gseid}.scmocha.out.rds.gz"))
        )
      }
    )
  ) ->
gse_data_loaded


gse_data_loaded |>
  tidyr::unnest(cols = anno) ->
gse_data

gse_data |>
  dplyr::select(somatic_variant) |>
  dplyr::mutate(
    somatic_variant_new = purrr::map(
      .x = somatic_variant,
      .f = \(.x) {
        .x$somatic
      }
    )
  ) |>
  dplyr::pull(somatic_variant_new) |>
  purrr::reduce(union) ->
somatic_variants

gse_data |>
  dplyr::select(gseid, srrid, srrdir) |>
  dplyr::mutate(
    raw = parallel::mclapply(
      X = srrdir,
      FUN = function(.srrdir) {
        log_info(glue::glue("Processing {basename(.srrdir)}"))

        barcode_cluster_file <- "barcode_cluster.tsv"
        cluster_umap <- fn_load_cluster(
          .filename = file.path(.srrdir, barcode_cluster_file)
        )

        cell_hetero_raw_file <- "cell.cell_heteroplasmic_df_raw.tsv.gz"
        cell_hetero_raw <- fn_load_hetero(
          .filename = file.path(.srrdir, cell_hetero_raw_file)
        ) |>
          dplyr::filter(variant %in% somatic_variants)

        # cell_coverage_file <- "cell.coverage.txt.gz"
        # cell_coverage <- fn_load_coverage(
        #   .filename = file.path(.srrdir, cell_coverage_file)
        # )

        cell_meta_data_file <- "cell_meta_data.tsv"
        metadata <- fn_load_meta(
          .filename = file.path(.srrdir, cell_meta_data_file)
        )


        # cell_raw_cluster_af <- cluster_umap |>
        #   dplyr::left_join(cell_hetero_raw, by = "barcode") |>
        #   dplyr::rename(cluster = celltype) |>
        #   tidyr::pivot_wider(
        #     names_from = variant,
        #     values_from = af
        #   )

        # cell_raw_cluster_forplot <- fn_forplot(
        #   .af = cell_raw_cluster_af,
        #   .coverage = cell_coverage,
        #   .meta = metadata
        # )
        # log_success("Done")
        log_success(glue::glue("Processing {basename(.srrdir)}"))

        tibble::tibble(
          cluster_umap = list(cluster_umap),
          cell_hetero_raw = list(cell_hetero_raw),
          # cell_coverage = list(cell_coverage),
          metadata = list(metadata),
          # cell_raw_cluster_af = list(cell_raw_cluster_af),
          # cell_raw_cluster_forplot = list(cell_raw_cluster_forplot)
        )
      },
      mc.cores = 20
    )
  ) ->
gse_data_af

# gse_data_af |>
#   dplyr::filter(purrr::map_lgl(
#     .x = raw,
#     .f = \(.x) {
#       is.null(nrow(.x))
#     }
#   ))

gse_data_af |>
  # head(11) |>
  dplyr::mutate(
    raw = purrr::map2(
      .x = raw,
      .y = srrdir,
      .f = function(.raw, .srrdir) {
        if (!is.null(nrow(.raw))) {
          return(.raw)
        }
        log_info(glue::glue("Processing {basename(.srrdir)}"))

        tryCatch(
          {
            barcode_cluster_file <- "barcode_cluster.tsv"
            cluster_umap <- fn_load_cluster(
              .filename = file.path(.srrdir, barcode_cluster_file)
            )

            cell_hetero_raw_file <- "cell.cell_heteroplasmic_df_raw.tsv.gz"
            cell_hetero_raw <- fn_load_hetero(
              .filename = file.path(.srrdir, cell_hetero_raw_file)
            ) |>
              dplyr::filter(variant %in% somatic_variants)

            cell_meta_data_file <- "cell_meta_data.tsv"
            metadata <- fn_load_meta(
              .filename = file.path(.srrdir, cell_meta_data_file)
            )

            log_success(glue::glue("Processing {basename(.srrdir)}"))

            tibble::tibble(
              cluster_umap = list(cluster_umap),
              cell_hetero_raw = list(cell_hetero_raw),
              metadata = list(metadata)
            )
          },
          error = function(e) {
            log_error(glue::glue("Error processing {basename(.srrdir)}: {e$message}"))
            NULL
          }
        )
      }
    )
  ) ->
gse_data_af_new


# readr::write_rds(
#   gse_data_af_new,
#   file.path(
#     foundation_out,
#     "GSE226602-and-other-datasets-scfoundation-plottable-metadata-variant.R.gse_data_af.rds.gz"
#   )
# )

gse_data_af_new <- readr::read_rds(
  file.path(
    foundation_out,
    "GSE226602-and-other-datasets-scfoundation-plottable-metadata-variant.R.gse_data_af.rds.gz"
  )
)

# merge --------------------------------------------------------------------

gse_data
gse_data_af_new
gse_dataset_metadata_full

celltypes <- c("B", "CD4_T", "CD8_T", "DC", "Mono", "NK", "other")
gse_data_af_new |>
  dplyr::select(gseid, srrid, raw) |>
  tidyr::unnest(cols = raw) |>
  dplyr::select(-metadata) |>
  # head(20) |>
  dplyr::mutate(
    cell_hetero_raw_anno = parallel::mclapply(
      X = seq_along(cluster_umap),
      FUN = function(i) {
        .cluster_umap <- cluster_umap[[i]]
        .cell_hetero_raw <- cell_hetero_raw[[i]]

        .cell_hetero_raw |>
          dplyr::filter(variant %in% somatic_variants) |>
          dplyr::left_join(.cluster_umap, by = "barcode") |>
          dplyr::filter(celltype %in% celltypes)
      },
      mc.cores = 20
    )
  ) |>
  dplyr::select(gseid, srrid, cell_hetero_raw_anno) ->
gse_data_af_new_anno

gse_data_af_new_anno |>
  tidyr::unnest(cols = cell_hetero_raw_anno) |>
  as.data.table() ->
gse_data_af_new_anno_dt

gse_data_af_new |>
  dplyr::select(gseid, srrid, raw) |>
  tidyr::unnest(cols = raw) |>
  dplyr::select(-metadata) |>
  # head(20) |>
  dplyr::mutate(
    cell_hetero_raw_anno = parallel::mclapply(
      X = seq_along(cluster_umap),
      FUN = function(i) {
        .cluster_umap <- cluster_umap[[i]]
        .cell_hetero_raw <- cell_hetero_raw[[i]]

        .cell_hetero_raw |>
          dplyr::filter(variant %in% somatic_variants) |>
          dplyr::left_join(.cluster_umap, by = "barcode") |>
          dplyr::filter(celltype %in% celltypes) |>
          dplyr::select(celltype, variant, af) |>
          dplyr::group_by(celltype, variant) |>
          dplyr::summarize(af = mean(af, na.rm = TRUE)) |>
          dplyr::ungroup()
      },
      mc.cores = 20
    )
  ) |>
  dplyr::select(gseid, srrid, cell_hetero_raw_anno) ->
gse_data_af_new_celltype


gse_data_af_new_celltype |>
  tidyr::unnest(cols = cell_hetero_raw_anno) |>
  dplyr::mutate(
    celltype_new = glue::glue("{srrid}_{celltype}")
  ) ->
gse_data_af_new_celltype_new

gse_data_af_new_celltype_new |>
  dplyr::select(celltype_new, variant, af) |>
  tidyr::spread(key = celltype_new, value = af) ->
gse_data_af_new_celltype_new_wide



library(ComplexHeatmap)
library(circlize)


# No need to add variant column as it already exists
gse_data_af_new_celltype_new_wide |>
  as.data.frame() |>
  tibble::column_to_rownames(var = "variant") |>
  as.matrix() ->
af_mtx



gse_data_af_new_celltype_new |>
  dplyr::select(celltype_new, celltype, gseid, srrid) |>
  dplyr::mutate(celltype = factor(celltype, levels = celltypes)) |>
  dplyr::distinct() ->
af_cluster
pcc <- readr::read_tsv(file = "https://raw.githubusercontent.com/chunjie-sam-liu/chunjie-sam-liu.life/master/public/data/pcc.tsv") |>
  dplyr::arrange(cancer_types)
col_clusters <- levels(af_cluster$celltype)
col_colors <- pcc$color[1:length(levels(af_cluster$celltype))]
names(col_colors) <- col_clusters
chm_top <- ComplexHeatmap::HeatmapAnnotation(
  df = af_cluster,
  # gap = unit(c(2, 2), "mm"),
  col = list(
    celltype = col_colors
  ),
  which = "column"
)

sort(unique(as.numeric(af_mtx))) -> .seq_af

.seq_af


# .seq_af |> hist()
# .seq_af_sub <- c(
#   -0.1, 0, 0.05, 0.1, 0.15, 0.2,
#   .seq_af[(.seq_af > 0.2 & .seq_af < 0.8)],
#   0.8, 0.85, 0.9, 0.95, 1
# )
# .seq_af_sub |> hist()
col_option <- "turbo"
col_pick = c(
  "grey",
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
# (a) change the “0” value from dark blue to dark green;
# col_pick |> color()
#  (b) reverse the color, use red as “0” and violet as “1”
col_fun = circlize::colorRamp2(
  # seq(col_start,
  #   col_end,
  #   length.out = n_break
  # ),
  .seq_af,
  col_pick
)

ComplexHeatmap::Heatmap(
  matrix = af_mtx,
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
  show_row_names = F,
  cluster_rows = F,
  cluster_row_slices = T,
  clustering_distance_rows = "pearson",
  clustering_method_rows = "ward.D",
  row_names_gp = gpar(
    # fontsize = 20,
    # col = .gcol$cell_variants
  ),
  # column
  # column_title = paste0("palette = '", col_option, "'"),
  column_title = NULL,
  column_title_gp = gpar(fontsize = 40),
  cluster_columns = FALSE,
  cluster_column_slices = T,
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
)

# readr::write_rds(
#   gse_data_af_new_anno_dt,
#   file.path(
#     foundation_out,
#     "GSE226602-and-other-datasets-scfoundation-plottable-metadata-variant.R.gse_data_af_new_anno_dt.rds.gz"
#   )
# )

# footer ------------------------------------------------------------------

# future: :plan(future: :sequential)

# save image --------------------------------------------------------------
save.image(
  file = file.path(
    foundation_out,
    "/home/liuc9/github/scMOCHA-data/stats/scfoundation/GSE226602-and-other-datasets-scfoundation-plottable-metadata-variant.R.RData"
  )
)
