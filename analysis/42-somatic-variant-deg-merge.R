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
# gseid_srrid_variant <- import(
#   "/home/liuc9/github/scMOCHA-data/analysis/zzz/plot-real-somatic-variant/main-variants/gseid_srrid_variant.fst"
# )
gseid_srrid_variant_hetero_plot_ratio <- import(
  "/home/liuc9/github/scMOCHA-data/analysis/zzz/plot-real-somatic-variant/main-variants/gseid_srrid_variant_hetero_plot_ratio.qs"
)


gseid_srrid_variant_hetero_plot_ratio |>
  dplyr::select(gseid, srrid, variant) -> gseid_srrid_variant

# load conn ---------------------------------------------------------------

# src ---------------------------------------------------------------------

# function ----------------------------------------------------------------

fn_de_plot <- function(
  markers,
  .cutoff_pval = 0.05,
  .cutoff_log2fc = 0.25,
  .pct = 0.05
) {
  markers |>
    tibble::rownames_to_column("gene") |>
    dplyr::mutate(
      fdr = -log10(p_val_adj)
    ) |>
    dplyr::mutate(
      fdr = ifelse(
        fdr > -log10(1e-100),
        -log10(1e-100),
        fdr
      )
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
        p_val_adj < .cutoff_pval &
          (pct.1 >= .pct | pct.2 >= .pct) &
          avg_log2FC > .cutoff_log2fc ~
          "red",
        p_val_adj < .cutoff_pval &
          (pct.1 >= .pct | pct.2 >= .pct) &
          avg_log2FC < -.cutoff_log2fc ~
          "blue",
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
        dplyr::filter(color != "grey") |>
        dplyr::group_by(color) |>
        dplyr::slice_head(n = 20) |>
        dplyr::ungroup(),
      aes(label = gene),
      size = 3,
      max.overlaps = 20
    ) +
    scale_color_identity() +
    geom_vline(
      xintercept = c(
        -.cutoff_log2fc,
        .cutoff_log2fc
      ),
      linetype = "dashed"
    ) +
    geom_hline(
      yintercept = -log10(.cutoff_pval),
      linetype = "dashed"
    ) +
    theme_classic() +
    labs(
      x = "Fold change(Heteroplasmy/Sufficient reads)",
      y = "FDR",
      # title = "{thegseid}-{thesrrid}-m.{thevariant}" |> glue::glue()
    ) +
    theme(
      plot.title = element_text(
        hjust = 0.5,
        face = "bold",
        color = "black",
        size = 16
      )
    ) -> p
  list(
    p = p,
    markers = forplot
  )
}
fn_go_enrich <- function(cancer_sgene, ont = c("BP", "CC", "MF")) {
  .ont <- match.arg(ont)

  gobp <- clusterProfiler::enrichGO(
    gene = cancer_sgene |> unique(),
    OrgDb = "org.Hs.eg.db",
    keyType = "SYMBOL",
    # keyType = "ENSEMBL",
    ont = .ont
  )

  gobp
}

fn_plot_go <- function(.go, .topn = Inf, .ont = c("BP", "CC", "MF")) {
  .ont <- match.arg(.ont)

  base_fill <- c("BP" = "#AE1700", "CC" = "#DF8F44FF", "MF" = "#00A1D5FF")
  ont_fullname <- c(
    "BP" = "Biological Process",
    "CC" = "Cellular Component",
    "MF" = "Molecular Function"
  )

  .ont_fill <- base_fill[.ont]
  x_label <- ont_fullname[.ont]

  .go %>%
    tibble::as_tibble() %>%
    dplyr::mutate(
      Description = stringr::str_wrap(
        stringr::str_to_sentence(string = Description),
        width = 60
      )
    ) %>%
    dplyr::mutate(adjp = -log10(p.adjust)) %>%
    dplyr::select(ID, Description, adjp, Count, geneID) %>%
    dplyr::arrange(adjp, Count) %>%
    dplyr::mutate(
      Description = factor(Description, levels = Description)
    ) -> .go_bp_for_plot

  if (!is.infinite(.topn)) {
    .go_bp_for_plot |>
      tail(.topn) -> .go_bp_for_plot
  }

  .go_bp_for_plot |>
    ggplot(aes(x = Description, y = adjp)) +
    geom_col(fill = .ont_fill, color = NA, width = 0.7) +
    geom_text(aes(label = Count), hjust = 1, color = "white", size = 5) +
    labs(y = "-log10(Adj. P value)", x = x_label) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
    coord_flip() +
    theme(
      panel.background = element_rect(fill = NA),
      panel.grid = element_blank(),
      axis.line.x = element_line(color = "black"),
      axis.line.y = element_line(color = "black"),
      # axis.title.y = element_blank(),
      axis.text.y = element_text(color = "black", size = 13, hjust = 1),
      axis.ticks.length.y = unit(3, units = "mm"),
      axis.text.x = element_text(color = "black", size = 12),
      axis.title = element_text(colour = "black", size = 16, face = "bold")
    )
}

fn_variant_go <- function(markers, .variant) {
  # .variant <- variants[[1]]
  .corr_pos <- markers |>
    dplyr::filter(color == "red") |>
    dplyr::mutate(
      genename = gene
    )
  .corr_neg <- markers |>
    dplyr::filter(color == "blue") |>
    dplyr::mutate(
      genename = gene
    )

  .pos_bp <- fn_go_enrich(cancer_sgene = unique(.corr_pos$genename), "BP")
  .pos_cc <- fn_go_enrich(cancer_sgene = unique(.corr_pos$genename), "CC")
  .pos_mf <- fn_go_enrich(cancer_sgene = unique(.corr_pos$genename), "MF")

  .neg_bp <- fn_go_enrich(cancer_sgene = unique(.corr_neg$genename), "BP")
  .neg_cc <- fn_go_enrich(cancer_sgene = unique(.corr_neg$genename), "CC")
  .neg_mf <- fn_go_enrich(cancer_sgene = unique(.corr_neg$genename), "MF")

  tibble::tibble(
    # variant = .variant,
    pos_bp = list(.pos_bp),
    pos_cc = list(.pos_cc),
    pos_mf = list(.pos_mf),
    neg_bp = list(.neg_bp),
    neg_cc = list(.neg_cc),
    neg_mf = list(.neg_mf),
    pos_bp_plot = list(
      fn_plot_go(.pos_bp, 20, "BP") +
        labs(title = .variant) +
        theme(
          plot.title = element_text(size = 20)
        )
    ),
    pos_cc_plot = list(
      fn_plot_go(.pos_cc, 20, "CC") +
        labs(title = .variant) +
        theme(
          plot.title = element_text(size = 20)
        )
    ),
    pos_mf_plot = list(
      fn_plot_go(.pos_mf, 20, "MF") +
        labs(title = .variant) +
        theme(
          plot.title = element_text(size = 20)
        )
    ),
    neg_bp_plot = list(
      fn_plot_go(.neg_bp, 20, "BP") +
        labs(title = .variant) +
        theme(
          plot.title = element_text(size = 20)
        )
    ),
    neg_cc_plot = list(
      fn_plot_go(.neg_cc, 20, "CC") +
        labs(title = .variant) +
        theme(
          plot.title = element_text(size = 20)
        )
    ),
    neg_mf_plot = list(
      fn_plot_go(.neg_mf, 20, "MF") +
        labs(title = .variant) +
        theme(
          plot.title = element_text(size = 20)
        )
    )
  )
}


conn <- DBI::dbConnect(
  duckdb::duckdb(),
  "/home/liuc9/github/scMOCHA-data/analysis/zzz/clean-data/all_hetero_af.cell.duckdb.1.2.1",
  read_only = TRUE
)
dplyr::tbl(
  conn,
  "allvariants_cell"
) -> tbl_allvariants_cell

#' Very important function
#' @example fn_plot_cell_af_depth_forplot("10398A>G", "GSM5494107")
fn_plot_cell_af_depth_forplot <- function(thevariant, thesrrid) {
  source("analysis/00-colors.R")

  colorcode <- setNames(names(color_variantcell), color_variantcell)

  tbl_allvariants_cell |>
    dplyr::filter(
      srrid == thesrrid,
      variant == thevariant
    ) |>
    dplyr::collect() |>
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

fn_merge_sc_list_variant <- function(df, thevariant) {
  # df <- a$data[[5]]
  # thevariant <- a$variant[[5]]
  library(Seurat)
  log_info("Merging sc objects for variant {thevariant}")
  sc_list <- df |> dplyr::pull(sc_file)

  .outdir <- path(
    "/home/liuc9/github/scMOCHA-data/analysis/zzz/plot-real-somatic-variant/main-variants/",
    thevariant,
    "/deg_merge/"
  )
  dir_create(.outdir)

  sc_merge_path <- path(
    .outdir,
    glue::glue("sc_merge.sct.{thevariant}.qs")
  )

  # if (file_exists(sc_merge_path)) {
  #   glue::glue(
  #     "sc merge file for variant {thevariant} already exists, skip merging."
  #   ) |>
  #     log_debug()
  #   return(1)
  # }
  log_info("loading forplot for variant {thevariant}")
  parallel::mclapply(
    X = df$srrid,
    FUN = function(thesrrid) {
      fn_plot_cell_af_depth_forplot(
        thevariant = thevariant,
        thesrrid = thesrrid
      )
    },
    mc.cores = 10
  ) |>
    dplyr::bind_rows() |>
    dplyr::mutate(
      barcode_new = glue::glue("{gseid}-{srrid}-{barcode}")
    ) |>
    as.data.table() -> forplot_list
  # forplot_ <- fn_plot_cell_af_depth_forplot(
  #   thevariant = thevariant,
  #   thesrrid = thesrrid
  # )
  log_info("Loading sc objects for variant {thevariant}")
  parallel::mclapply(
    sc_list,
    function(f) {
      .sc <- import(f)
      .sc[["SCT"]]@scale.data <- matrix()
      .sc
    },
    mc.cores = 10
  ) -> sc_list_loaded

  log_fatal(
    "Loaded {length(sc_list_loaded)} sc objects for variant {thevariant}"
  )

  parallel::mclapply(
    sc_list_loaded,
    Seurat::VariableFeatures,
    mc.cores = 10
  ) |>
    unlist() |>
    unique() -> var_features

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
        merge.data = FALSE,
        ...
      )
    }
    cli::cli_progress_done()
    return(out)
  }

  log_info("Merging sc objects for variant {thevariant}")

  if (length(sc_list_loaded) < 2) {
    glue::glue(
      "Less than 2 sc objects for variant {thevariant}, cannot merge."
    )
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
    # sc_merge <- merge(
    #   x = sc_list_loaded[[1]],
    #   y = sc_list_loaded[2:length(sc_list_loaded)],
    #   merge.data = FALSE # not merge the scale.data, for memory sake
    # )
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
  log_success(
    "Merged sc object for variant {thevariant} with {ncol(sc_merge)} cells."
  )
  rm(sc_list_loaded)
  gc()

  Seurat::VariableFeatures(sc_merge) <- var_features

  log_info("Exporting merged sc object for variant {thevariant}")
  export(
    sc_merge,
    sc_merge_path
  )

  # sc_merge
  rm(sc_merge)
  gc()
  1
}


# body --------------------------------------------------------------------
thevariants <- c(
  "3173G>A",
  "3176A>T",
  "3178T>A",
  "3727T>C",
  "3728C>T",
  "13271T>C",
  "14063T>C",
  "14831G>A",
  "1643A>G",
  "3667T>G",
  "4175G>A",
  "5513G>A",
  "7065G>A",
  "9025G>A",
  "9237G>A",
  "10398A>G"
)

gseid_srrid_variant |>
  dplyr::filter(
    variant %in% thevariants
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

gseid_srrid_variant_sc |> dplyr::count(variant)

gseid_srrid_variant_sc |>
  tidyr::nest(.by = "variant") |>
  # dplyr::filter(
  #   !variant %in% c("3176A>T", "3173G>A", "3178T>A", "3728C>T", "3727T>C")
  # ) |>
  dplyr::mutate(
    data_merge = purrr::map2(
      .x = data,
      .y = variant,
      .f = fn_merge_sc_list_variant
    )
  ) -> gseid_srrid_variant_sc_merged


#
#
# ? check merge file --------------------------------------------------------------------
#
#

gseid_srrid_variant_sc |>
  dplyr::select(variant) |>
  dplyr::distinct() |>
  dplyr::mutate(
    variant_dir = path(
      "/home/liuc9/github/scMOCHA-data/analysis/zzz/plot-real-somatic-variant/main-variants/",
      variant,
      "/deg_merge/"
    )
  ) |>
  dplyr::mutate(
    merge_sct = path(
      variant_dir,
      glue::glue("sc_merge.sct.{variant}.qs")
    )
  ) |>
  dplyr::mutate(
    file_exists = file_exists(merge_sct)
  ) -> gseid_srrid_variant_sc_merged_check

gseid_srrid_variant_sc_merged_check |>
  dplyr::filter(!variant %in% c("3173G>A", "3176A>T")) |>
  # head(2) |>
  dplyr::mutate(
    variant_load = parallel::mcmapply(
      FUN = \(merge_sct) {
        # merge_sct <- gseid_srrid_variant_sc_merged_check$merge_sct[[1]]
        log_info("Start processing {merge_sct}")
        if (!file_exists(merge_sct)) {
          return(NULL)
        } else {
          log_info("Importing merged sc for {merge_sct}")
          .sc <- import(merge_sct)
          .variant <- .sc$variant |>
            unique() |>
            na.omit() |>
            as.character()
          rm(.sc)
          gc()
          return(.variant)
        }
      },
      merge_sct,
      mc.cores = 10,
      SIMPLIFY = FALSE
    )
  ) -> gseid_srrid_variant_sc_merged_check_load


gseid_srrid_variant_sc_merged_check_load |>
  dplyr::mutate(
    thesamevariant = purrr::map2_lgl(
      .x = variant,
      .y = variant_load,
      .f = \(v1, v2) {
        # v1 <- gseid_srrid_variant_sc_merged_check_load$variant[[1]]
        # v2 <- gseid_srrid_variant_sc_merged_check_load$variant_load[[1]]
        if (is.null(v2)) {
          return(FALSE)
        }
        if (all(v1 == v2)) {
          return(TRUE)
        } else {
          return(FALSE)
        }
      }
    )
  ) |>
  dplyr::filter(
    !thesamevariant
  ) -> gseid_srrid_variant_sc_merged_check_load_issue

gseid_srrid_variant_sc |>
  tidyr::nest(.by = "variant") |>
  dplyr::filter(
    variant %in% gseid_srrid_variant_sc_merged_check_load_issue$variant
  ) |>
  dplyr::mutate(
    data_merge = purrr::map2(
      .x = data,
      .y = variant,
      .f = fn_merge_sc_list_variant
    )
  ) -> gseid_srrid_variant_sc_merged_fixed


# export(
#   gseid_srrid_variant_sc_merged,
#   "/home/liuc9/github/scMOCHA-data/analysis/zzz/plot-real-somatic-variant/main-variants/gseid_srrid_variant_sc_merged.qs"
# )
#
#
# ! don't run below --------------------------------------------------------------------
# ! don't run below --------------------------------------------------------------------
# ! don't run below --------------------------------------------------------------------
# ! don't run below --------------------------------------------------------------------
#
#

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


export(
  sc_merge,
  "/home/liuc9/github/scMOCHA-data/analysis/zzz/plot-real-somatic-variant/main-variants/3727T>C/deg_merge/sc_merge.sct.3727T>C.qs"
)

#
# ? 3728C>T --------------------------------------------------------------------
#
#
gseid_srrid_variant_sc |>
  dplyr::filter(variant == "3728C>T") |>
  dplyr::mutate(
    # load = parallel::mclapply(
    load = lapply(
      sc_file,
      function(f) {
        .sc <- import(f)
        .sc[["SCT"]]@scale.data <- matrix()
        .sc
      }
      # mc.cores = 12,
    )
  ) -> gseid_srrid_variant_sc_filtered

# GSM7080030

# gseid_srrid_variant_sc_filtered |> dplyr::filter(srrid == "GSM7080030")

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


export(
  sc_merge,
  "/home/liuc9/github/scMOCHA-data/analysis/zzz/plot-real-somatic-variant/main-variants/3728C>T/deg_merge/sc_merge.sct.3728C>T.qs"
)

# footer ------------------------------------------------------------------

# future: :plan(future: :sequential)

# save image --------------------------------------------------------------
