#!/usr/bin/env Rscript
# Metainfo ----------------------------------------------------------------
# @AUTHOR: Chun-Jie Liu
# @CONTACT: chunjie.sam.liu.at.gmail.com
# @DATE: 2025-09-03 20:11:27
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
library(scales)

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

# load conn ---------------------------------------------------------------
dir_main_variant <- "/home/liuc9/github/scMOCHA-data/analysis/zzz/plot-real-somatic-variant/main-variants"

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
        fdr > -log10(1e-300),
        -log10(1e-300),
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

  # forplot |>
  #   dplyr::count(color) |>
  #   print()

  forplot |>
    dplyr::count(color) |>
    tibble::deframe() -> n_color

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
      subtitle = glue::glue(
        "Up:{
          scales::label_comma()(n_color[[3]])
        }, down: {
          scales::label_comma()(n_color[[1]])
        };(Cutoff: FDR<{.cutoff_pval}, log2FC>{.cutoff_log2fc}, Pct>{.pct})"
      )
    ) +
    theme(
      plot.title = element_text(
        hjust = 0.5,
        face = "bold",
        color = "black",
        size = 16
      ),
      plot.subtitle = element_text(
        hjust = 0.5,
        # face = "bold",
        color = "black",
        size = 12
      ),
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

fn_de_ <- function(
  thevariant,
  sc,
  .ident.1,
  .ident.2,
  .group.by,
  .prefix,
  .labs,
  .celltype = NULL
) {
  # .ident.1 <- glue::glue("{hetero_label} high")
  # .ident.2 <- glue::glue("{hetero_label} low")
  # .group.by <- "cellvarianttype2"
  # .prefix <- "hetero_high_vs_low"
  # .labs <- labs(
  #   x = "Fold change {hetero_label} High vs Low" |> glue::glue(),
  #   y = "FDR",
  #   # title = "m.{thevariant}" |> glue::glue()
  #   title = "Markers: {hetero_label} High vs Low (m.{thevariant})" |>
  #     glue::glue()
  # )

  markers_hetero_high_vs_low <- Seurat::FindMarkers(
    object = sc,
    ident.1 = .ident.1,
    ident.2 = .ident.2,
    assay = "SCT",
    slot = "data",
    test.use = "wilcox",
    group.by = .group.by,
    latent.vars = "srrid",
    features = Seurat::VariableFeatures(sc)
  )

  .outdir <- file.path(
    dir_main_variant,
    thevariant,
    "deg_merge"
  )

  if (!is.null(.celltype)) {
    .outdir <- file.path(
      dir_main_variant,
      thevariant,
      "deg_merge",
      .celltype
    )
  }

  dir.create(
    .outdir,
    showWarnings = FALSE,
    recursive = TRUE
  )

  export(
    markers_hetero_high_vs_low,
    file.path(
      .outdir,
      "markers.{.prefix}.{thevariant}.qs" |>
        glue::glue()
    )
  )

  fn_de_plot(
    markers_hetero_high_vs_low,
    .cutoff_pval = 0.05,
    .cutoff_log2fc = 0.25,
    .pct = 0.05
  ) -> p_hetero_high_vs_low
  p_hetero_high_vs_low$p + .labs -> p_hetero_high_vs_low$p

  ggsave(
    filename = "markers.{.prefix}.{thevariant}.pdf" |> glue::glue(),
    plot = p_hetero_high_vs_low$p,
    path = .outdir,
    device = "pdf",
    width = 10,
    height = 6
  )
  p_hetero_high_vs_low
}
fn_go_ <- function(
  thevariant,
  p_hetero_high_vs_low,
  .prefix,
  .celltype = NULL
) {
  # .prefix <- "hetero_high_vs_low"

  fn_variant_go(
    p_hetero_high_vs_low$markers,
    thevariant
  ) -> p_go_hetero_high_vs_low
  .outdir <- file.path(
    dir_main_variant,
    thevariant,
    "go_merge"
  )
  if (!is.null(.celltype)) {
    .outdir <- file.path(
      dir_main_variant,
      thevariant,
      "go_merge",
      .celltype
    )
  }
  dir.create(
    .outdir,
    showWarnings = FALSE,
    recursive = TRUE
  )

  export(
    p_go_hetero_high_vs_low,
    file.path(
      .outdir,
      "markers.{.prefix}.{thevariant}.go.qs" |>
        glue::glue()
    )
  )

  tibble::tibble(
    pn = c("pos", "neg") |> rep(each = 3),
    t = c("bp", "cc", "mf") |> rep(each = 2)
  ) |>
    dplyr::mutate(
      saveimage = purrr::map2(
        .x = pn,
        .y = t,
        .f = \(.x, .y) {
          .p <- p_go_hetero_high_vs_low[[glue::glue("{.x}_{.y}_plot")]][[1]]
          .filename <- "markers.{.prefix}.{thevariant}.go.{.x}_{.y}_plot.pdf" |>
            glue::glue()
          ggsave(
            filename = .filename,
            plot = .p +
              labs(
                title = glue::glue(
                  "m.{thevariant} {ifelse(is.null(.celltype), '', .celltype)}"
                )
              ),
            path = .outdir,
            device = "pdf",
            width = 10,
            height = 6
          )
        }
      )
    )
}

fn_load_sc <- function(thevariant) {
  library(Seurat)
  sc <- import(
    file.path(
      dir_main_variant,
      thevariant,
      "deg_merge",
      glue::glue("sc_merge.sct.{thevariant}.qs")
    )
  )
  sc
}

fn_variant_ <- function(
  thevariant,
  sc,
  .celltype = NULL
) {
  # thevariant <- "3727T>C"

  sc@meta.data |>
    as.data.table() |>
    dplyr::filter(cellvarianttype == "Heteroplasmy") |>
    dplyr::pull(af) |>
    median(na.rm = FALSE) -> median_af

  # scales::label_number(accuracy = 0.01)(median_af)
  hetero_label <- glue::glue(
    "Heteroplasmy (median {scales::label_number(accuracy = 0.01)(median_af)})"
  )

  sc@meta.data |>
    dplyr::mutate(
      cellvarianttype2 = dplyr::case_when(
        cellvarianttype == "Heteroplasmy" &
          af >= median_af ~
          glue::glue("{hetero_label} high"),
        cellvarianttype == "Heteroplasmy" &
          af < median_af ~
          glue::glue("{hetero_label} low"),
        TRUE ~ cellvarianttype
      )
    ) -> sc@meta.data

  DefaultAssay(sc) <- "SCT"

  sc@meta.data |>
    dplyr::count(cellvarianttype2) |>
    dplyr::arrange(cellvarianttype2) |>
    tibble::deframe() -> n_cellvarianttype2

  p_hetero_high_vs_low <- fn_de_(
    thevariant = thevariant,
    sc = sc,
    .ident.1 = glue::glue("{hetero_label} high"),
    .ident.2 = glue::glue("{hetero_label} low"),
    .group.by = "cellvarianttype2",
    .prefix = "hetero_high_vs_low",
    .labs <- labs(
      x = glue::glue(
        "Fold change {hetero_label} High (n={
          scales::label_comma()(n_cellvarianttype2[[1]])
        }) vs Low (n={
          scales::label_comma()(n_cellvarianttype2[[2]])
        })"
      ),
      y = "FDR",
      # title = "m.{thevariant}" |> glue::glue()
      title = "Markers: {hetero_label} High vs Low (m.{thevariant}) {ifelse(is.null(.celltype), '', .celltype)}" |>
        glue::glue()
    ),
    .celltype = .celltype
  )

  fn_go_(
    thevariant = thevariant,
    p_hetero_high_vs_low = p_hetero_high_vs_low,
    .prefix = "hetero_high_vs_low",
    .celltype = .celltype
  )

  p_hetero_vs_sufficient <- fn_de_(
    thevariant = thevariant,
    sc = sc,
    .ident.1 = "Heteroplasmy",
    .ident.2 = "Sufficient reads",
    .group.by = "cellvarianttype",
    .prefix = "hetero_vs_sufficient",
    .labs <- labs(
      x = "Fold change Heteroplasmy (n={
          scales::label_comma()(n_cellvarianttype2[[1]] + n_cellvarianttype2[[2]])
        }) vs Sufficient reads (n={
          scales::label_comma()(n_cellvarianttype2[[5]])
        })",
      y = "FDR",
      title = "Markers: Heteroplasmy vs Sufficient Reads (m.{thevariant}) {ifelse(is.null(.celltype), '', .celltype)}" |>
        glue::glue()
    ),
    .celltype = .celltype
  )

  fn_go_(
    thevariant = thevariant,
    p_hetero_high_vs_low = p_hetero_vs_sufficient,
    .prefix = "hetero_vs_sufficient",
    .celltype = .celltype
  )
}


fn_variant_cell_ <- function(thevariant) {
  library(Seurat)
  sc <- fn_load_sc(thevariant = thevariant)

  sc$predicted.celltype.l1 |> unique() -> celltypes

  parallel::mclapply(
    celltypes,
    function(.celltype) {
      sc_sub <- Seurat::subset(
        sc,
        subset = predicted.celltype.l1 == .celltype
      )
      fn_variant_(
        thevariant = thevariant,
        sc = sc_sub,
        .celltype = .celltype
      )
    },
    mc.cores = 10
  )
}
# body --------------------------------------------------------------------

#
#
# ? 3727T>C --------------------------------------------------------------------
#
#

fn_variant_(
  thevariant = "3727T>C",
  sc = fn_load_sc(
    thevariant = "3727T>C"
  )
)

fn_variant_cell_(
  thevariant = "3727T>C"
)
#
#
# ? 3728C>T --------------------------------------------------------------------
#
#
fn_variant_(
  thevariant = "3728C>T",
  sc = fn_load_sc(
    thevariant = "3728C>T"
  )
)
fn_variant_cell_(
  thevariant = "3728C>T"
)

# footer ------------------------------------------------------------------

# future: :plan(future: :sequential)

# save image --------------------------------------------------------------
