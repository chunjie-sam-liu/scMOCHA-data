#!/usr/bin/env Rscript
# Metainfo ----------------------------------------------------------------
# @AUTHOR: Chun-Jie Liu
# @CONTACT: chunjie.sam.liu.at.gmail.com
# @DATE: 2026-01-02 13:29:35
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
dims = 1:30

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
scintegrateddir <- outdirnotuse / "scintegrated"

# Ensure output directories exist
fs::dir_create(scmergedir)
fs::dir_create(scintegrateddir)

cleandatadir <- path("/home/liuc9/github/scMOCHA-data/data/zzz/clean-data")

METAFULL <- import(outdir / "SAMPLES-METADATA-FULL.xlsx")
ALLVARIANTS <- import(
  outdir / "SAMPLE-VARIANT-CLASSIFICATION-CLUSTER-BULK-AF.xlsx"
)

HOMO_HETE_VARIANTS <- ALLVARIANTS |>
  dplyr::filter(variant_type %in% c("homo", "hete"))

# load conn ---------------------------------------------------------------
dotenv(".env")

conn <- conn_db(
  Sys.getenv("DUCKDB_PATH"),
  read_only = TRUE
)
dplyr::tbl(
  conn,
  "allvariants_cell"
) -> tbl_allvariants_cell

# src ---------------------------------------------------------------------

# function ----------------------------------------------------------------

#' Very important function
#' @example fn_plot_cell_af_depth_forplot("10398A>G", "GSM5494107")
fn_plot_cell_af_depth_forplot <- function(thevariant, thesrrid) {
  source(path(
    Sys.getenv("HIGHRESDIR"),
    "00-colors.R"
  ))

  colorcode <- setNames(names(color_variantcell), color_variantcell)

  tbl_allvariants_cell |>
    dplyr::filter(
      srrid == thesrrid,
      variant == thevariant
    ) |>
    as.data.table() |>
    dplyr::mutate(
      variant_type = dplyr::case_match(
        variant_type,
        "colorful" ~ "red",
        "black" ~ "darkblue",
        "white" ~ "white",
        "grey" ~ "gray",
        NA ~ "white"
      )
    ) |>
    dplyr::mutate(
      variant_type = factor(
        variant_type,
        levels = color_variantcell
      )
    ) |>
    dplyr::arrange(
      variant_type,
      -af
    ) -> forplot_

  forplot_ |>
    dplyr::mutate(
      barcode = factor(
        barcode,
        levels = forplot_$barcode
      )
    ) |>
    dplyr::mutate(
      celltype = gsub(
        "_",
        " ",
        celltype
      )
    ) |>
    dplyr::mutate(
      celltype = factor(
        celltype,
        names(color_celltype)
      )
    ) |>
    dplyr::mutate(
      af = ifelse(
        af < 0.01,
        NA_real_,
        af
      )
    ) |>
    # dplyr::mutate(
    #   variant_type = as.character(variant_type),
    # ) |>
    dplyr::mutate(
      depth = log2(depth + 1) # log2 transform to reduce skewness
    ) |>
    dplyr::mutate(
      cellvarianttype = colorcode[variant_type]
    ) |>
    dplyr::mutate(
      cellvarianttype = factor(
        cellvarianttype,
        levels = colorcode
      )
    ) -> forplot
  forplot
}

library(Seurat)

fn_get_forplot_list <- function(df, thevariant) {
  lapply(
    X = df$srrid,
    FUN = function(thesrrid) {
      fn_plot_cell_af_depth_forplot(
        thevariant = thevariant,
        thesrrid = thesrrid
      )
    }
  ) |>
    dplyr::bind_rows() |>
    dplyr::mutate(
      barcode_new = glue::glue("{gseid}-{srrid}-{barcode}")
    ) |>
    as.data.table()
}
fn_load_sc_list <- function(df) {
  lapply(
    df$sc_file,
    function(f) {
      tryCatch(
        {
          .sc <- qs::qread(f)
          if (is.null(.sc)) {
            log_error("Failed to load file: {f}")
            return(NULL)
          }
          # .sc[["SCT"]]@scale.data <- matrix()
          .sc
        },
        error = function(e) {
          log_error("Error loading file {f}: {e$message}")
          return(NULL)
        }
      )
    }
  )
}
fn_get_var_features <- function(sc_list_loaded) {
  lapply(
    sc_list_loaded,
    Seurat::VariableFeatures
  ) |>
    unlist() |>
    unique()
}
fn_merge_with_progress <- function(sc_list_loaded, ...) {
  n <- length(sc_list_loaded)
  out <- sc_list_loaded[[1]]
  cli::cli_progress_bar(
    "Merging Seurat objects",
    total = n - 1,
    format = "{cli::pb_bar} {cli::pb_percent} ({cli::pb_current}/{cli::pb_total})"
  )
  for (i in 2:n) {
    cli::cli_progress_update()
    out <- merge(
      out,
      sc_list_loaded[[i]],
      ...
    )
  }
  cli::cli_progress_done()
  return(out)
}
fn_merge <- function(sc_list_loaded, forplot_list, thevariant, var_features) {
  if (length(sc_list_loaded) < 2) {
    log_info("Less than 2 sc objects for variant {thevariant}, cannot merge.")

    sc_merge <- sc_list_loaded[[1]]
    sc_merge@meta.data |>
      tibble::rownames_to_column("barcode_new") |>
      as.data.table() |>
      dplyr::left_join(
        forplot_list |>
          dplyr::select(
            -c(barcode, celltype)
          ),
        by = "barcode_new"
      ) |>
      as.data.frame() |>
      tibble::column_to_rownames("barcode_new") -> .d_merge

    sc_merge@meta.data <- .d_merge
  } else {
    sc_merge <- fn_merge_with_progress(
      sc_list_loaded
    ) # not merge the scale.data, for memory sake

    sc_merge@meta.data |>
      tibble::rownames_to_column("barcode_new") |>
      as.data.table() |>
      dplyr::left_join(
        forplot_list |>
          dplyr::select(
            -c(barcode, celltype)
          ),
        by = "barcode_new"
      ) |>
      as.data.frame() |>
      tibble::column_to_rownames("barcode_new") -> .d_merge

    sc_merge@meta.data <- .d_merge
  }
  Seurat::VariableFeatures(sc_merge) <- var_features
  sc_merge
}

fn_norm <- function(sc_merge) {
  sc_merge |>
    Seurat::NormalizeData() |>
    Seurat::FindVariableFeatures() |>
    Seurat::ScaleData() |>
    Seurat::RunPCA() |>
    Seurat::FindNeighbors(
      dims = dims,
      reduction = "pca"
    ) -> sc_merge
  sc_merge
}

fn_integrated <- function(sc_merge) {
  obj <- fn_norm(sc_merge)

  obj <- Seurat::IntegrateLayers(
    object = obj,
    method = CCAIntegration,
    orig.reduction = "pca",
    new.reduction = "integrated.cca",
    verbose = FALSE
  )

  cli_alert_success("Finished integration layers")
  obj[["RNA"]] <- SeuratObject::JoinLayers(obj[["RNA"]])
  cli_alert_success("Finished joining layers")

  # obj <- obj |>
  #   Seurat::FindNeighbors(
  #     reduction = "integrated.cca",
  #     dims = dims
  #   ) |>
  #   Seurat::FindClusters(
  #     resolution = resolution,
  #     cluster.name = "clusters_integrated"
  #   ) |>
  #   Seurat::RunTSNE(
  #     reduction = "integrated.cca",
  #     dims = dims,
  #     reduction.name = "tsne"
  #   ) |>
  #   Seurat::RunUMAP(
  #     reduction = "integrated.cca",
  #     dims = dims,
  #     reduction.name = "umap"
  #   )

  obj
}
fn_merge_sc_list_variant <- function(df, thevariant) {
  # VARIANT_GSEID_SRRID_SCFILE |>
  #   tibble::rowid_to_column("idx") |>
  #   dplyr::filter(variant == "3240C>G")
  # df <- VARIANT_GSEID_SRRID_SCFILE$gseid_srrid[[1438]]
  # thevariant <- VARIANT_GSEID_SRRID_SCFILE$variant[[1438]]

  log_info("Merging sc objects for variant {thevariant}")

  # step 1, if get forplot_list
  fn_get_forplot_list(df, thevariant) -> forplot_list
  log_success("Step 1: forplot_list for variant {thevariant} loaded.")

  # step 2, load sc_list
  fn_load_sc_list(df) -> sc_list_loaded

  # Filter out NULL entries from sc_list_loaded
  sc_list_loaded <- Filter(Negate(is.null), sc_list_loaded)

  if (length(sc_list_loaded) == 0) {
    log_warn(
      "No valid Seurat objects loaded for variant {thevariant}, skipping."
    )
    return(NULL)
  }

  n_samples <- length(sc_list_loaded)
  log_success(
    "Step 2: sc_list for variant {thevariant} loaded. {n_samples} valid objects."
  )

  # # step 3, get var features
  fn_get_var_features(sc_list_loaded) -> var_features
  log_success("Step 3: var_features for variant {thevariant} obtained.")

  # step 4, merge
  fn_merge(
    sc_list_loaded,
    forplot_list,
    thevariant,
    var_features
  ) -> sc_merge
  log_success("Step 4: sc_merge for variant {thevariant} obtained.")

  # clean up
  rm(forplot_list)
  rm(sc_list_loaded)
  gc()
  log_success("Clean up done for variant sc_list_loaded.")

  # step 5, export sc_merge
  export(
    sc_merge,
    scmergedir / glue::glue("sc.{thevariant}.merge.qs")
  )
  log_success("Step 5: sc_merge for variant {thevariant} exported.")

  # step 6, integrated analysis
  # sc_merge <- fn_integrated(sc_merge)
  # to reduce mem, dont use function.
  sc_merge |>
    Seurat::NormalizeData() |>
    Seurat::FindVariableFeatures() |>
    Seurat::ScaleData() |>
    Seurat::RunPCA() |>
    Seurat::FindNeighbors(
      dims = 1:30,
      reduction = "pca"
    ) -> sc_merge

  # Only integrate if there are multiple samples (layers)
  n_layers <- length(SeuratObject::Layers(
    sc_merge,
    assay = "RNA",
    search = "counts"
  ))
  if (n_layers >= 2) {
    sc_merge <- Seurat::IntegrateLayers(
      object = sc_merge,
      method = CCAIntegration,
      orig.reduction = "pca",
      new.reduction = "integrated.cca",
      verbose = FALSE
    )
    cli_alert_success("Finished integration layers for variant {thevariant}")
    sc_merge[["RNA"]] <- SeuratObject::JoinLayers(sc_merge[["RNA"]])
    cli_alert_success("Finished joining layers for variant {thevariant}")
  } else {
    log_info(
      "Only {n_layers} layer(s) for variant {thevariant}, skipping integration."
    )
  }

  export(
    sc_merge,
    scintegrateddir / glue::glue("sc.{thevariant}.integrated.qs")
  )

  # sc_merge
  rm(sc_merge)
  gc()
  1
}

# body --------------------------------------------------------------------
HOMO_HETE_VARIANTS |>
  dplyr::filter(variant_type == "hete") |>
  dplyr::count(variant) |>
  dplyr::arrange(n) |>
  dplyr::filter(n > 30) |>
  dplyr::pull(variant) -> variant_tobe_run

HOMO_HETE_VARIANTS |>
  dplyr::select(gseid, srrid, variant, variant_type) |>
  dplyr::distinct() |>
  dplyr::mutate(
    sc_file = path(
      "/home/liuc9/github/scMOCHA-data/data/",
      gseid,
      "final",
      srrid,
      "for_integration",
      "sc_azimuth.qs"
    )
  ) |>
  dplyr::mutate(
    file_exists = file.exists(sc_file)
  ) |>
  dplyr::filter(file_exists) |>
  tidyr::nest(.by = variant, .key = "gseid_srrid") |>
  dplyr::mutate(
    n = purrr::map_int(gseid_srrid, nrow)
  ) |>
  dplyr::arrange(n) -> VARIANT_GSEID_SRRID_SCFILE


VARIANT_GSEID_SRRID_SCFILE |>
  tibble::rowid_to_column("idx") |>
  dplyr::filter(variant %in% variant_tobe_run) |>
  dplyr::filter(
    idx > which(VARIANT_GSEID_SRRID_SCFILE$variant == "1380G>T")
  )

# dplyr::filter(n == 3)

# fn_merge_sc_list_variant(
#   VARIANT_GSEID_SRRID_SCFILE$gseid_srrid[[1]],
#   VARIANT_GSEID_SRRID_SCFILE$variant[[1]]
# )

# fn_merge_sc_list_variant(
#   VARIANT_GSEID_SRRID_SCFILE$gseid_srrid[[1250]],
#   VARIANT_GSEID_SRRID_SCFILE$variant[[1250]]
# )

VARIANT_GSEID_SRRID_SCFILE |>
  dplyr::filter(variant %in% variant_tobe_run) |>
  dplyr::mutate(
    a = parallel::mcmapply(
      FUN = function(df, thevariant) {
        fn_merge_sc_list_variant(
          df,
          thevariant
        )
      },
      df = gseid_srrid,
      thevariant = variant,
      mc.cores = 10
    )
  )

# footer ------------------------------------------------------------------

# save image --------------------------------------------------------------
