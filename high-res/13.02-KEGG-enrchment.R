#!/usr/bin/env Rscript
# Author: Chunjie Liu
# Contact: chunjie.sam.liu.at.gmail.com
# Date: 2026-02-02
# Description: KEGG enrichment analysis for DEGs
# Version: 0.3

# Library -----------------------------------------------------------------

library(conflicted)
conflict_prefer("filter", "dplyr")
conflict_prefer("path", "fs")

# Use load_pkg for the rest
if (!requireNamespace("jutils", quietly = TRUE)) {
  # If jutils is not installed, we might need to source it or install it
}
library(magrittr)
library(clusterProfiler)
library(org.Hs.eg.db)
library(ReactomePA)
library(writexl)
library(jutils)

# args --------------------------------------------------------------------

GetoptLong.options(help_style = "two-column")
VERSION = "v0.0.1"

verbose = TRUE
GetoptLong("verbose!", "print messages")

logger::log_threshold(logger::TRACE)
logger::log_layout(logger::layout_glue_colors)

# header ------------------------------------------------------------------

dotenv(".env")
outdir <- fs::path(Sys.getenv("OUTDIR"))
outdirnotuse <- fs::path(Sys.getenv("OUTDIRNOTUSE"))
cleandatadir <- fs::path(Sys.getenv("CLEANDATADIR"))
highresdir <- fs::path(Sys.getenv("HIGHRESDIR"))

# function ----------------------------------------------------------------

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
    ) -> .d

  if (!"rn" %in% colnames(.d)) {
    .d <- .d |> tibble::rownames_to_column(var = "SYMBOL")
  } else {
    .d <- .d |> dplyr::rename(SYMBOL = rn)
  }

  .d |>
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

fn_gseKEGG <- function(geneList) {
  tryCatch(
    clusterProfiler::gseKEGG(
      geneList = geneList,
      organism = "hsa",
      verbose = FALSE
    ),
    error = function(e) {
      log_warn(glue::glue("Error in gseKEGG: {e$message}"))
      return(NULL)
    }
  )
}

fn_enrich_kegg_only <- function(markers) {
  markers |>
    dplyr::filter(!is.na(ENTREZID)) |>
    dplyr::select(ENTREZID, avg_log2FC) |>
    # Handle duplicates by taking the absolute maximum log2FC
    dplyr::group_by(ENTREZID) |>
    dplyr::summarize(
      avg_log2FC = avg_log2FC[which.max(abs(avg_log2FC))],
      .groups = "drop"
    ) |>
    dplyr::arrange(-avg_log2FC) |>
    tibble::deframe() -> geneList

  if (length(geneList) == 0) {
    log_warn("No genes mapped to ENTREZID.")
    return(NULL)
  }

  log_info(glue::glue("Running gseKEGG with {length(geneList)} genes."))
  res <- fn_gseKEGG(geneList)
  if (is.null(res)) {
    log_warn("gseKEGG returned NULL.")
  } else {
    log_info(glue::glue("gseKEGG found {nrow(as.data.table(res))} pathways."))
  }
  return(res)
}

fn_variant_kegg <- function(thevariant, tmpdir = "deg_merge_vaf") {
  variant_dir <- fs::path(outdirnotuse, "deg", thevariant)
  base_dir <- fs::path(variant_dir, tmpdir)

  if (!fs::dir_exists(base_dir)) {
    log_warn(glue::glue(
      "Directory {base_dir} does not exist. Skipping variant {thevariant}."
    ))
    return(NULL)
  }

  markers_list <- dir_ls(
    path = base_dir,
    recurse = TRUE,
    regexp = glue::glue("markers.*{thevariant}.qs")
  )

  if (length(markers_list) == 0) {
    log_warn(glue::glue("No markers files found for {thevariant} in {tmpdir}."))
    return(NULL)
  }

  tibble::tibble(marker_path = markers_list) |>
    dplyr::mutate(
      celltype = fs::path_dir(marker_path) |> fs::path_file()
    ) |>
    dplyr::mutate(
      celltype = ifelse(celltype == tmpdir, "all_cells", celltype)
    ) |>
    dplyr::mutate(
      filename = fs::path_file(marker_path)
    ) |>
    dplyr::mutate(
      filename = gsub(
        pattern = glue::glue("markers.|.qs|.{thevariant}"),
        replacement = "",
        x = filename
      )
    ) |>
    dplyr::mutate(
      markers = lapply(
        X = marker_path,
        FUN = function(.path) {
          tryCatch(
            {
              import(.path) |> fn_markers_update()
            },
            error = function(e) {
              log_error(glue::glue(
                "Error reading or updating markers from {.path}: {e$message}"
              ))
              return(NULL)
            }
          )
        }
      )
    ) -> .df_variant

  .df_variant |>
    dplyr::mutate(
      kegg = mapply(
        FUN = function(.markers) {
          if (is.null(.markers)) {
            return(NULL)
          }
          fn_enrich_kegg_only(markers = .markers)
        },
        .markers = markers,
        SIMPLIFY = FALSE
      )
    ) -> .df_variant_res

  kegg_tmpdir <- gsub("deg", "kegg", tmpdir)
  outdir_kegg <- fs::path(variant_dir, kegg_tmpdir)
  dir_create(outdir_kegg)

  .df_variant_res |>
    dplyr::mutate(
      plot = mapply(
        FUN = function(.kegg) {
          if (is.null(.kegg)) {
            return(NULL)
          }
          .kegg_dt <- as.data.table(.kegg)
          if (nrow(.kegg_dt) == 0) {
            return(NULL)
          }

          .kegg_dt |>
            dplyr::filter(p.adjust < 0.05) |>
            dplyr::mutate(FDR = -log10(qvalue)) |>
            dplyr::mutate(
              y = glue::glue("{ID}_{Description}")
            ) -> .kegg_filtered

          if (nrow(.kegg_filtered) == 0) {
            return(NULL)
          }

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
        SIMPLIFY = FALSE
      )
    ) -> .df_variant_plots

  # Export plots object
  export(
    .df_variant_plots,
    fs::path(outdir_kegg, glue::glue("kegg_enrich_plots.{thevariant}.qs"))
  )

  # Save PDFs
  .df_variant_plots |>
    dplyr::mutate(
      ggsave_path = fs::path(
        outdir_kegg,
        glue::glue("kegg_enrich_plot.{celltype}.{filename}.{thevariant}.pdf") |>
          fs::path_sanitize() |>
          sanitize_filename()
      )
    ) |>
    dplyr::mutate(
      a = mapply(
        FUN = function(.plot, .ggsave_path) {
          if (is.null(.plot)) {
            return(NULL)
          }
          ggsave(
            filename = .ggsave_path,
            plot = .plot,
            width = 8,
            height = 6
          )
        },
        .plot = plot,
        .ggsave_path = ggsave_path,
        SIMPLIFY = FALSE
      )
    )

  return(.df_variant_plots)
}

fn_main_kegg <- function(thevariant, tmpdir = "deg_merge_vaf") {
  log_info(glue::glue(
    "Processing KEGG enrichment for variant {thevariant} in {tmpdir}"
  ))

  kegg_tmpdir <- gsub("deg", "kegg", tmpdir)
  excel_path <- fs::path(
    outdirnotuse,
    "deg",
    thevariant,
    kegg_tmpdir,
    glue::glue("kegg_enrich_details.{thevariant}.vaf.xlsx")
  )

  if (fs::file_exists(excel_path)) {
    log_info(glue::glue(
      "KEGG enrichment already processed for {thevariant}. Skipping."
    ))
    return(NULL)
  }

  vaf_kegg <- fn_variant_kegg(thevariant = thevariant, tmpdir = tmpdir)
  if (is.null(vaf_kegg)) {
    return(NULL)
  }

  kegg_url <- "https://www.kegg.jp/kegg-bin/show_pathway"

  colors_file <- fs::path(highresdir, "00-colors.R")
  if (fs::file_exists(colors_file)) {
    source(colors_file)
  } else {
    color_celltype <- NULL
  }

  vaf_kegg |>
    dplyr::select(celltype, filename, kegg) |>
    dplyr::mutate(
      a = lapply(
        X = kegg,
        FUN = function(.kegg) {
          if (is.null(.kegg) || inherits(.kegg, "try-error")) {
            return(NULL)
          }

          .kegg_dt <- as.data.table(.kegg)
          if (nrow(.kegg_dt) == 0) {
            return(NULL)
          }

          .kegg_dt |>
            dplyr::filter(p.adjust < 0.05) -> .kegg_filtered

          if (nrow(.kegg_filtered) == 0) {
            return(NULL)
          }

          .kegg_filtered |>
            dplyr::mutate(
              KEGG_URL = glue::glue("{kegg_url}?{ID}/{core_enrichment}")
            ) |>
            dplyr::mutate(
              Symbol = purrr::map_chr(
                .x = core_enrichment,
                .f = function(.genes) {
                  .genes_split <- strsplit(.genes, "/")[[1]]
                  tryCatch(
                    {
                      clusterProfiler::bitr(
                        geneID = .genes_split,
                        fromType = "ENTREZID",
                        toType = "SYMBOL",
                        OrgDb = org.Hs.eg.db::org.Hs.eg.db
                      ) |>
                        dplyr::pull(SYMBOL) |>
                        paste(collapse = "/")
                    },
                    error = function(e) {
                      return("")
                    }
                  )
                }
              )
            ) -> .d
          return(.d)
        }
      )
    ) -> vaf_kegg_enrich_details

  vaf_kegg_enrich_details_unnested <- vaf_kegg_enrich_details |>
    dplyr::select(-kegg) |>
    tidyr::unnest(cols = a, keep_empty = FALSE)

  if (nrow(vaf_kegg_enrich_details_unnested) == 0) {
    log_info(glue::glue("No significant KEGG pathways found for {thevariant}."))
    # Create empty excel with headers
    empty_df <- tibble::tibble(
      celltype = character(),
      compare = character(),
      kegg_id = character(),
      description = character(),
      nes = numeric(),
      pvalue = numeric(),
      p.adjust = numeric(),
      KEGG_URL = character(),
      Symbol = character(),
      core_enrichment = character(),
      setSize = integer(),
      rank = integer(),
      leading_edge = character()
    )
    writexl::write_xlsx(empty_df, path = excel_path)
    return(NULL)
  }

  if (!is.null(color_celltype)) {
    vaf_kegg_enrich_details_unnested |>
      dplyr::mutate(
        celltype = factor(
          celltype,
          levels = c("all_cells", names(color_celltype))
        )
      ) -> vaf_kegg_enrich_details_unnested
  }

  vaf_kegg_enrich_details_final <- vaf_kegg_enrich_details_unnested |>
    dplyr::select(
      celltype,
      compare = filename,
      kegg_id = ID,
      description = Description,
      nes = NES,
      pvalue = pvalue,
      p.adjust,
      KEGG_URL,
      Symbol,
      core_enrichment,
      setSize,
      rank,
      leading_edge
    ) |>
    dplyr::arrange(celltype, compare, p.adjust) |>
    dplyr::mutate(
      kegg_id = writexl::xl_hyperlink(
        url = KEGG_URL,
        name = kegg_id
      )
    )

  writexl::write_xlsx(
    vaf_kegg_enrich_details_final,
    path = excel_path
  )

  log_success(glue::glue(
    "Finished KEGG enrichment for {thevariant}. Results in {excel_path}"
  ))
}

# body --------------------------------------------------------------------

# List of variants to process
thevariants <- c(
  "14082C>G",
  "15169A>G",
  "3240C>G",
  "7757G>A",
  "3173G>A",
  "3176A>T",
  "3178T>A",
  "9025G>A",
  "9237G>A",
  "10398A>G",
  "4175G>A",
  "6227T>C",
  "9006A>G"
)


fn_main_kegg(
  thevariant = "6967G>A",
  tmpdir = "deg_merge_vaf"
)
fn_main_kegg(
  thevariant = "14831G>A",
  tmpdir = "deg_merge_vaf"
)

fn_main_kegg(
  thevariant = "7833T>C",
  tmpdir = "deg_merge_vaf"
)

fn_main_kegg(
  thevariant = "10500G>A",
  tmpdir = "deg_merge_vaf"
)

fn_main_kegg(
  thevariant = "4175G>A",
  tmpdir = "deg_merge_vaf"
)

# Process each variant
# for (thevariant in thevariants) {
#   tryCatch({
#     fn_main_kegg(thevariant = thevariant, tmpdir = "deg_merge_vaf")
#   }, error = function(e) {
#     log_error(glue::glue("Failed processing {thevariant}: {e$message}"))
#   })
# }

# footer ------------------------------------------------------------------
