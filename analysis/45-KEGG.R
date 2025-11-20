#!/usr/bin/env Rscript
# Metainfo ----------------------------------------------------------------
# @AUTHOR: Chun-Jie Liu
# @CONTACT: chunjie.sam.liu.at.gmail.com
# @DATE: 2025-11-19 15:10:42
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
library(fs)

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
dir_main_variant <- "/home/liuc9/github/scMOCHA-data/analysis/zzz/plot-real-somatic-variant/main-variants"

msig_df <- msigdbr::msigdbr(species = "Homo sapiens") |> as.data.table()
msigdbr::msigdbr_species()
msigdbr::msigdbr_collections()
msig_df |>
  dplyr::mutate(gs_name = glue::glue("{gs_collection}#{gs_name}")) |>
  dplyr::select(
    gs_name,
    gs_collection,
    gs_subcollection,
    gs_collection_name,
    gene_symbol
  ) |>
  dplyr::select(gs_name, gene_symbol) -> msig_df_s


msig_df |>
  dplyr::filter(
    grepl("OXPHOS", gs_collection_name, ignore.case = TRUE)
  )
# load conn ---------------------------------------------------------------

# src ---------------------------------------------------------------------

# function ----------------------------------------------------------------

fn_gsea <- function(markers) {
  markers |>
    tibble::rownames_to_column(var = "GeneName") |>
    dplyr::rename(log2FC = avg_log2FC) -> .de

  .de %>%
    dplyr::arrange(-log2FC) %>%
    dplyr::select(GeneName, log2FC) %>%
    tidyr::drop_na() %>%
    tibble::deframe() -> geneList

  .GSEA <- clusterProfiler::GSEA(
    geneList = geneList,
    pAdjustMethod = "BH",
    TERM2GENE = msig_df_s,
    verbose = FALSE
  )

  .GSEA |>
    tibble::as_tibble() |>
    tidyr::separate(
      Description,
      into = c("Cat", "Description"),
      sep = "#"
    ) -> .GSEA_sep

  .GSEA_sep |>
    dplyr::filter(
      grepl("Atp", Description, ignore.case = TRUE)
    )
  # .GSEA_sep$Cat %>% table()

  .GSEA_sep |>
    dplyr::filter(
      grepl("OX", Description, ignore.case = TRUE)
    )

  .GSEA_sep %>%
    dplyr::filter(Cat == "C5") -> .GSEA_sep_C5

  .GSEA_sep_C5 %>%
    dplyr::mutate(adjp = -log10(p.adjust)) %>%
    dplyr::select(ID, Description, adjp, NES) %>%
    dplyr::filter(grepl(pattern = "gobp", Description, ignore.case = TRUE)) %>%
    dplyr::mutate(
      Description = gsub(pattern = "GOBP_", replacement = "", x = Description)
    ) %>%
    dplyr::mutate(
      Description = gsub(pattern = "_", replacement = " ", x = Description)
    ) %>%
    dplyr::mutate(
      Description = stringr::str_wrap(
        stringr::str_to_sentence(string = Description),
        width = 60
      )
    ) %>%
    dplyr::arrange(NES) %>%
    dplyr::mutate(color = ifelse(NES > 0, "P", "N")) %>%
    dplyr::slice(1:10, (dplyr::n() - 9):dplyr::n()) %>%
    dplyr::distinct() %>%
    dplyr::arrange(NES, -adjp) %>%
    dplyr::mutate(
      Description = factor(Description, levels = Description)
    ) -> .GSEA_sep_C5_for_plot

  .GSEA_sep_C5_for_plot %>%
    ggplot(aes(x = Description, y = NES)) +
    geom_col(aes(fill = color)) +
    scale_fill_manual(values = c("#54AE59", "#9F82B5")) +
    scale_x_discrete(
      limit = .GSEA_sep_C5_for_plot$Description,
      labels = stringr::str_wrap(
        stringr::str_replace_all(
          string = .GSEA_sep_C5_for_plot$Description,
          pattern = "_",
          replacement = " "
        ),
        width = 50
      )
    ) +
    labs(y = "Normalized Enrichment Score") +
    coord_flip() +
    theme(
      panel.background = element_rect(fill = NA),
      panel.grid = element_blank(),

      axis.line.x = element_line(color = "black"),
      axis.text.x = element_text(color = "black"),

      axis.title.y = element_blank(),
      axis.line.y = element_line(color = "black"),
      legend.position = "none"
    ) -> gseaplotc5

  .GSEA_sep %>%
    dplyr::filter(Cat == "H") -> .GSEA_sep_H

  .GSEA_sep_H %>%
    dplyr::mutate(adjp = -log10(p.adjust)) %>%
    dplyr::select(ID, Description, adjp, NES) %>%
    # dplyr::filter(grepl(pattern = "gobp", Description, ignore.case = TRUE)) %>%
    # dplyr::mutate(Description = gsub(pattern = "GOBP_", replacement = "", x = Description)) %>%
    dplyr::mutate(
      Description = gsub(pattern = "_", replacement = " ", x = Description)
    ) %>%
    dplyr::mutate(
      Description = stringr::str_wrap(
        Description,
        width = 60
      )
    ) %>%
    dplyr::arrange(NES) %>%
    dplyr::mutate(color = ifelse(NES > 0, "P", "N")) -> .GSEA_sep_H_d

  .GSEA_sep_H_for_plot <- if (nrow(.GSEA_sep_H) > 20) {
    .GSEA_sep_H_d %>%
      dplyr::slice(1:10, (dplyr::n() - 9):dplyr::n()) %>%
      dplyr::distinct() %>%
      dplyr::arrange(NES, -adjp) %>%
      dplyr::mutate(Description = factor(Description, levels = Description))
  } else {
    .GSEA_sep_H_d %>%
      dplyr::arrange(NES, -adjp) %>%
      dplyr::mutate(Description = factor(Description, levels = Description))
  }

  .GSEA_sep_H_for_plot %>%
    ggplot(aes(x = Description, y = NES)) +
    geom_col(aes(fill = color)) +
    scale_fill_manual(values = c("#54AE59", "#9F82B5")) +
    scale_x_discrete(
      limit = .GSEA_sep_H_for_plot$Description,
      labels = stringr::str_wrap(
        stringr::str_replace_all(
          string = .GSEA_sep_H_for_plot$Description,
          pattern = "_",
          replacement = " "
        ),
        width = 50
      )
    ) +
    labs(y = "Normalized Enrichment Score") +
    coord_flip() +
    theme(
      panel.background = element_rect(fill = NA),
      panel.grid = element_blank(),

      axis.line.x = element_line(color = "black"),
      axis.text.x = element_text(color = "black"),

      axis.title.y = element_blank(),
      axis.line.y = element_line(color = "black"),
      legend.position = "none"
    ) -> gseaploth

  tibble::tibble(
    gsea = list(.GSEA_sep),
    gseaplotc5 = list(gseaplotc5),
    gseaploth = list(gseaploth)
  )
}

fn_kegg <- function(
  markers,
  .cutoff_pval = 0.05,
  .cutoff_log2fc = 0.25,
  .pct = 0.05
) {
  markers |>
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
    ) -> markers_de

  markers_de |>
    tibble::rownames_to_column(var = "GeneName") |>
    dplyr::filter(color != "grey") -> .gs
  # .gs

  .gs_id <- clusterProfiler::bitr(
    geneID = .gs$GeneName,
    fromType = "SYMBOL",
    toType = "ENTREZID",
    OrgDb = org.Hs.eg.db::org.Hs.eg.db
  )

  .kegg <- clusterProfiler::enrichKEGG(
    gene = .gs_id$ENTREZID,
    organism = "hsa",
    pvalueCutoff = 0.5
  )

  .kegg |>
    as.data.table()
}

fn_gseGO <- function(.markers) {
  # fn_markers_update(markers = markers) -> .markers
  .markers |>
    dplyr::select(ENTREZID, avg_log2FC) |>
    dplyr::arrange(-avg_log2FC) |>
    tibble::deframe() -> geneList

  .onts <- c("BP", "CC", "MF")
  parallel::mclapply(
    X = .onts,
    FUN = function(.ont) {
      clusterProfiler::gseGO(
        geneList = geneList,
        OrgDb = org.Hs.eg.db::org.Hs.eg.db,
        keyType = "ENTREZID",
        ont = .ont,
        verbose = FALSE,
      )
    },
    mc.cores = 3
  ) -> .gsego_list
  names(.gsego_list) <- .onts
  list(
    geneList = geneList,
    gsego = .gsego_list
  )
}

fn_gseKEGG <- function(.markers) {
  .markers |>
    dplyr::select(ENTREZID, avg_log2FC) |>
    dplyr::arrange(-avg_log2FC) |>
    tibble::deframe() -> geneList

  clusterProfiler::gseKEGG(
    geneList = geneList,
    organism = "hsa",
    verbose = FALSE
  ) -> .gsekegg
  list(
    geneList = geneList,
    gsekegg = .gsekegg
  )
}

fn_gseWP <- function(.markers) {
  .markers |>
    dplyr::select(ENTREZID, avg_log2FC) |>
    dplyr::arrange(-avg_log2FC) |>
    tibble::deframe() -> geneList

  clusterProfiler::gseWP(
    geneList = geneList,
    organism = "Homo sapiens",
  )
  list(
    geneList = geneList,
    gsewp = .gsewp
  )
}

fn_gsePathway <- function(.markers) {
  .markers |>
    dplyr::select(ENTREZID, avg_log2FC) |>
    dplyr::arrange(-avg_log2FC) |>
    tibble::deframe() -> geneList
  ReactomePA::gsePathway(
    geneList = geneList,
    pvalueCutoff = 0.2,
    pAdjustMethod = "BH",
    verbose = FALSE
  ) -> .gsepathway
  list(
    geneList = geneList,
    gsepathway = .gsepathway
  )
}


fn_enrichKEGG <- function(.markers) {
  .markers |>
    dplyr::filter(color != "grey") |>
    dplyr::select(ENTREZID, avg_log2FC) |>
    dplyr::arrange(-avg_log2FC) |>
    tibble::deframe() -> gene

  clusterProfiler::enrichKEGG(
    gene = names(gene),
    organism = "hsa"
  ) -> .enrichkegg
  list(
    gene = gene,
    enrichkegg = .enrichkegg
  )
}

fn_enrichWP <- function(.markers) {
  .markers |>
    dplyr::filter(color != "grey") |>
    dplyr::select(ENTREZID, avg_log2FC) |>
    dplyr::arrange(-avg_log2FC) |>
    tibble::deframe() -> gene

  clusterProfiler::enrichWP(
    gene = names(gene),
    organism = "Homo sapiens"
  ) -> .enrichwp

  list(
    gene = gene,
    enrichwp = .enrichwp
  )
}

fn_enrichPathway <- function(.markers) {
  .markers |>
    dplyr::filter(color != "grey") |>
    dplyr::select(ENTREZID, avg_log2FC) |>
    dplyr::arrange(-avg_log2FC) |>
    tibble::deframe() -> gene

  ReactomePA::enrichPathway(
    gene = names(gene),
    readable = TRUE
  ) -> .enrichpathway
  list(
    gene = gene,
    enrichpathway = .enrichpathway
  )
}

ek <- fn_enrichKEGG(.markers = fn_markers_update(markers = markers))
clusterProfiler::browseKEGG(ek$enrichkegg, "hsa04010")
# library(pathview)
# hsa04210 <- pathview::pathview(
#   gene.data = geneList,
#   pathway.id = "hsa04210",
#   species = "hsa",
#   limit = list(gene = max(abs(geneList)), cpd = 1)
# )

fn_markers_update <- function(
  markers,
  .cutoff_pval = 0.05,
  .cutoff_log2fc = 0.25,
  .pct = 0.05
) {
  markers |>
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
    ) |>
    tibble::rownames_to_column(
      var = "SYMBOL"
    ) |>
    as.data.table() |>
    dplyr::arrange(-avg_log2FC) -> .d

  .gs_id <- clusterProfiler::bitr(
    geneID = .d$SYMBOL,
    fromType = "SYMBOL",
    toType = c("ENTREZID", "ENSEMBL"),
    OrgDb = org.Hs.eg.db::org.Hs.eg.db
  ) |>
    as.data.table()

  .d |>
    dplyr::left_join(
      .gs_id,
      by = "SYMBOL"
    ) |>
    dplyr::mutate(
      ENSEMBL = ifelse(
        stringr::str_detect(SYMBOL, "ENSG"),
        SYMBOL,
        ENSEMBL
      )
    )
}

# body --------------------------------------------------------------------
thevariant <- "4175G>A"
markers <- import(
  fs::path(
    dir_main_variant,
    thevariant,
    "deg_merge_new",
    "markers.hetero_vs_sufficient.{thevariant}.qs" |> glue::glue()
  )
)


dplyr::filter(
  ID == "hsa00190"
)

# footer ------------------------------------------------------------------

# future: :plan(future: :sequential)

# save image --------------------------------------------------------------
