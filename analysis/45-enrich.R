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

# msig_df <- msigdbr::msigdbr(species = "Homo sapiens") |> as.data.table()

# msigdbr::msigdbr_species()
# msigdbr::msigdbr_collections()
# msig_df |>
#   dplyr::mutate(gs_name = glue::glue("{gs_collection}#{gs_name}")) |>
#   dplyr::select(
#     gs_name,
#     gs_collection,
#     gs_subcollection,
#     gs_collection_name,
#     gene_symbol
#   ) |>
#   dplyr::select(gs_name, gene_symbol) -> msig_df_s

# msig_df |>
#   dplyr::filter(
#     grepl("OXPHOS", gs_collection_name, ignore.case = TRUE)
#   )
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

fn_gseGO <- function(geneList) {
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
    mc.cores = length(.onts)
  ) -> .gsego_list
  names(.gsego_list) <- .onts
  .gsego_list
}

fn_gseKEGG <- function(geneList) {
  clusterProfiler::gseKEGG(
    geneList = geneList,
    organism = "hsa",
    verbose = FALSE
  ) -> .gsekegg
}

fn_gseWP <- function(geneList) {
  clusterProfiler::gseWP(
    geneList = geneList,
    organism = "Homo sapiens",
  )
}

fn_gsePathway <- function(geneList) {
  ReactomePA::gsePathway(
    geneList = geneList,
    pvalueCutoff = 0.2,
    pAdjustMethod = "BH",
    verbose = FALSE
  ) -> .gsepathway
}

fn_enrichGO <- function(gene, geneList) {
  .onts <- c("BP", "CC", "MF")
  parallel::mclapply(
    X = .onts,
    FUN = function(.ont) {
      clusterProfiler::enrichGO(
        gene = names(gene),
        universe = names(geneList),
        OrgDb = org.Hs.eg.db::org.Hs.eg.db,
        keyType = "ENTREZID",
        ont = .ont,
        pvalueCutoff = 0.05
      )
    },
    mc.cores = length(.onts)
  ) -> .gsego_list
  names(.gsego_list) <- .onts
  .gsego_list
}

fn_enrichKEGG <- function(gene) {
  clusterProfiler::enrichKEGG(
    gene = names(gene),
    organism = "hsa"
  ) -> .enrichkegg
}

fn_enrichWP <- function(gene) {
  clusterProfiler::enrichWP(
    gene = names(gene),
    organism = "Homo sapiens"
  ) -> .enrichwp
}

fn_enrichPathway <- function(gene) {
  ReactomePA::enrichPathway(
    gene = names(gene),
    readable = TRUE
  ) -> .enrichpathway
}


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

fn_enrich <- function(markers) {
  ls_gses <- list(
    fn_gseGO,
    fn_gseKEGG,
    fn_gseWP,
    fn_gsePathway
  )
  ls_enrich <- list(
    fn_enrichKEGG,
    fn_enrichWP,
    fn_enrichPathway
  )

  markers |>
    dplyr::select(ENTREZID, avg_log2FC) |>
    dplyr::arrange(-avg_log2FC) |>
    tibble::deframe() -> geneList

  markers |>
    dplyr::filter(color != "grey") |>
    dplyr::select(ENTREZID, avg_log2FC) |>
    dplyr::arrange(-avg_log2FC) |>
    tibble::deframe() -> gene

  tibble::tibble(
    fn = c(ls_gses, ls_enrich),
    type = c(
      rep("gse", length(ls_gses)),
      rep("enrich", length(ls_enrich))
    ),
    name = c(
      "gseGO",
      "gseKEGG",
      "gseWP",
      "gsePathway",
      "enrichKEGG",
      "enrichWP",
      "enrichPathway"
    ),
    gene = c(
      rep(list(geneList), length(ls_gses)),
      rep(list(gene), length(ls_enrich))
    )
  ) -> .df

  .df |>
    dplyr::mutate(
      res = parallel::mcmapply(
        FUN = function(.fn, .gene) {
          tryCatch(
            expr = {
              .fn(.gene)
            },
            error = function(e) {
              log_info(glue::glue("Error in {.fn}: {e$message}"))
              return(NULL)
            }
          )
        },
        .fn = fn,
        .gene = gene,
        mc.cores = nrow(.df),
        SIMPLIFY = FALSE
      )
    ) -> .df_res

  list(
    df_res = .df_res,
    geneList = geneList,
    gene = gene
  )
}

fn_enrich_kegg_only <- function(markers) {
  markers |>
    dplyr::select(ENTREZID, avg_log2FC) |>
    dplyr::arrange(-avg_log2FC) |>
    tibble::deframe() -> geneList

  fn_gseKEGG(geneList) -> gse_kegg

  gse_kegg
}

sanitize_filename <- function(x) {
  x |>
    gsub(">", "GT", x = _, fixed = TRUE) |>
    gsub("<", "LT", x = _, fixed = TRUE) |>
    gsub("%", "pct", x = _, fixed = TRUE) |>
    gsub("=", "-", x = _, fixed = TRUE) |>
    gsub("[()]", "", x = _) |>
    gsub("[[:space:]]+", "_", x = _) |>
    trimws()
}

fn_variant_kegg <- function(thevariant) {
  variant_dir <- fs::path(
    dir_main_variant,
    thevariant
  )

  base_dir <- fs::path(
    variant_dir,
    "deg_merge_new"
  )
  markers_list <- dir_ls(
    path = base_dir,
    recurse = TRUE,
    regexp = "markers.*{thevariant}.qs" |> glue::glue()
  )
  tibble::tibble(
    marker_path = markers_list
  ) |>
    dplyr::mutate(
      celltype = fs::path_dir(marker_path) |> fs::path_file()
    ) |>
    dplyr::mutate(
      celltype = ifelse(
        celltype == "deg_merge_new",
        "all_cells",
        celltype
      )
    ) |>
    dplyr::mutate(
      filename = fs::path_file(marker_path)
    ) |>
    dplyr::mutate(
      filename = gsub(
        pattern = "markers.|.qs|.{thevariant}" |> glue::glue(),
        replacement = "",
        x = filename
      )
    ) |>
    # dplyr::select(filename) |>
    dplyr::mutate(
      markers = parallel::mclapply(
        X = marker_path,
        FUN = function(.path) {
          import(.path) |>
            fn_markers_update()
        },
        mc.cores = length(marker_path)
      )
    ) -> .df_variant

  .df_variant |>
    dplyr::mutate(
      kegg = parallel::mcmapply(
        FUN = function(.markers) {
          fn_enrich_kegg_only(markers = .markers)
        },
        .markers = markers,
        mc.cores = nrow(.df_variant),
        SIMPLIFY = FALSE
      )
    ) -> .df_variant_res
  outdir <- fs::path(
    variant_dir,
    "kegg"
  )
  dir_create(outdir)

  export(
    .df_variant_res,
    fs::path(
      outdir,
      "kegg_enrich.{thevariant}.qs" |> glue::glue()
    )
  )

  fs::file_delete(
    dir_ls(
      outdir,
      regexp = ".pdf$"
    )
  )

  .df_variant_res |>
    # dplyr::select(celltype, filename, kegg) |>
    dplyr::mutate(
      plot = parallel::mcmapply(
        FUN = function(.kegg) {
          as.data.table(.kegg) |>
            dplyr::filter(p.adjust < 0.05) |>
            dplyr::mutate(FDR = -log10(qvalue)) |>
            dplyr::mutate(
              y = glue::glue("{ID}_{Description}")
            ) -> .kegg_filtered

          .kegg_filtered |>
            ggplot(aes(
              x = NES,
              y = reorder(y, NES),
              size = FDR,
              color = NES
            )) +
            geom_point() +
            scale_color_gradient2(
              low = "blue",
              mid = "white",
              high = "red"
            ) +
            theme(
              panel.background = element_rect(fill = NA),
              panel.grid = element_blank(),
              axis.line.x = element_line(color = "black"),
              axis.text.x = element_text(color = "black"),
              axis.title.y = element_blank(),
              axis.line.y = element_line(color = "black"),
              legend.position = "right"
            )
        },
        .kegg = kegg,
        mc.cores = nrow(.df_variant_res),
        SIMPLIFY = FALSE
      )
    ) -> .df_variant_plots

  export(
    .df_variant_plots,
    fs::path(
      outdir,
      "kegg_enrich_plots.{thevariant}.qs" |> glue::glue()
    )
  )
  sanitize_filename <- function(x) {
    x |>
      gsub(">", "GT", x = _, fixed = TRUE) |>
      gsub("<", "LT", x = _, fixed = TRUE) |>
      gsub("%", "pct", x = _, fixed = TRUE) |>
      gsub("=", "-", x = _, fixed = TRUE) |>
      gsub("[()]", "", x = _) |>
      gsub("[[:space:]]+", "_", x = _) |>
      trimws()
  }

  .df_variant_plots |>
    dplyr::mutate(
      ggsave_path = fs::path(
        outdir,
        "kegg_enrich_plot.{celltype}.{filename}.{thevariant}.pdf" |>
          glue::glue() |>
          fs::path_sanitize() |>
          sanitize_filename()
      )
    ) |>
    dplyr::mutate(
      a = parallel::mcmapply(
        FUN = function(.plot, .ggsave_path) {
          ggsave(
            filename = .ggsave_path,
            plot = .plot,
            width = 8,
            height = 6
          )
        },
        .plot = plot,
        .ggsave_path = ggsave_path,
        mc.cores = nrow(.df_variant_plots),
        SIMPLIFY = FALSE
      )
    ) -> m
}

# body --------------------------------------------------------------------
thevariant <- "4175G>A"


fn_variant_kegg(thevariant = "4175G>A")
fn_variant_kegg(thevariant = "9025G>A")
fn_variant_kegg(thevariant = "13271T>C")

m <- import(
  "/home/liuc9/github/scMOCHA-data/analysis/zzz/plot-real-somatic-variant/main-variants/4175G>A/kegg/kegg_enrich.4175G>A.qs"
)

m |> dplyr::filter(celltype == "all_cells") -> mm
mm$plot[[1]]

library(clusterProfiler)
mm

gseaplot2(
  mm$kegg[[1]],
  geneSetID = "hsa00190",
  title = "Oxidative phosphorylation"
) -> p_Oxidative_phosphorylation

ggsave(
  filename = "Oxidative_phosphorylation.all_cells.{mm$filename[[1]]}.pdf" |>
    glue::glue() |>
    fs::path_sanitize() |>
    sanitize_filename(),
  plot = p_Oxidative_phosphorylation,
  path = "/home/liuc9/github/scMOCHA-data/analysis/zzz/plot-real-somatic-variant/main-variants/4175G>A/kegg/",
  width = 6,
  height = 4
)


mm$markers[[1]] |>
  dplyr::select(ENTREZID, avg_log2FC) |>
  dplyr::arrange(-avg_log2FC) |>
  tibble::deframe() -> geneList
library(pathview)
pathview::pathview(
  gene.data = geneList,
  pathway.id = "hsa00190",
  species = "hsa",
  limit = list(gene = max(abs(geneList)), cpd = 1),
  out.suffix = "4175G>A_Oxidative_phosphorylation",
  kegg.native = TRUE
)

mm$kegg[[1]] |> as.data.table() |> dplyr::filter(ID == "hsa00190")

mm$markers[[1]] |>
  dplyr::filter(color != "grey") |>
  dplyr::select(ENTREZID, avg_log2FC) |>
  dplyr::arrange(-avg_log2FC) |>
  tibble::deframe() -> gene

fn_enrichKEGG(gene) -> kk

browseKEGG(
  kk,
  pathID = "hsa00190"
)

# footer ------------------------------------------------------------------

# future: :plan(future: :sequential)

# save image --------------------------------------------------------------
