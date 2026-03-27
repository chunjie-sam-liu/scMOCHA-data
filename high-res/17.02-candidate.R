#!/usr/bin/env Rscript
# Metainfo ----------------------------------------------------------------
# @AUTHOR: Chun-Jie Liu
# @CONTACT: chunjie.sam.liu.at.gmail.com
# @DATE: 2026-03-25 13:27:43
# @DESCRIPTION: this script is used for ...

# Reproducibility ----------------------------------------------------------
set.seed(1)
# Library -----------------------------------------------------------------

suppressMessages({
  load_pkg(jutils)
})

# Args --------------------------------------------------------------------

# s: string, i: integer, f: float, !: boolean, @: array, %: list
VERSION = "v0.0.1"

GetoptLong.options(help_style = "two-column")

# default: default value specified here.

nthread = 8
GetoptLong(
  "nthread=i",
  "Number of threads to use",
  "verbose",
  "Enable verbose logging"
)


# Logger ------------------------------------------------------------------

log_layout(layout_glue_colors)

if (isTRUE(verbose)) {
  log_threshold(TRACE)
  log_info("Verbose mode enabled")
} else {
  log_threshold(INFO)
}


# Load data ---------------------------------------------------------------
load_pkg(jutils)
dotenv()
suppressMessages({
  conflicted::conflict_prefer("filter", "dplyr")
})

outdir <- path(Sys.getenv("OUTDIR"))
outdirnotuse <- path(Sys.getenv("OUTDIRNOTUSE"))

ALLVARIANTS <- import(
  outdir / "SAMPLE-VARIANT-CLASSIFICATION-CLUSTER-BULK-AF.xlsx"
) |>
  dplyr::filter(variant_type %in% c("hete", "homo"))
METAFULL <- import(outdir / "SAMPLES-METADATA-FULL.xlsx")
variant_annotation <- import(outdir / "VARIANT-ANNOTATION-TABLE-APOGEE2.xlsx")

variant_list <-
  import(
    outdirnotuse /
      "allvariants-prioritize" /
      "allvariants-prioritize-variant-annotation-with-samples-n-clean-group.qs"
  )

# Source ---------------------------------------------------------------------
source(
  path(Sys.getenv("HIGHRESDIR"), "00-colors.R")
)
# Conn ---------------------------------------------------------------

# Function ----------------------------------------------------------------
fn_plot_cor_vaf_gene_sc <- function(
  thevariant,
  thegseid,
  thesrrid,
  thegene
) {
  load_pkg(Seurat)

  conn <- db_conn(
    Sys.getenv("DUCKDB_PATH"),
    readonly = TRUE
  )

  .dt <- tbl(conn, "allvariants_cell") |>
    dplyr::filter(
      srrid == thesrrid,
      variant == thevariant,
      variant_type %in% c("colorful", "black")
    ) |>
    as.data.table() |>
    dplyr::mutate(
      barcode = glue("{thegseid}-{thesrrid}-{barcode}")
    )

  sc <- import(
    "/mnt/isilon/u01_project/large-scale/liuc9/raw/{thegseid}/final/{thesrrid}/for_integration/sc_azimuth.qs" |>
      glue()
  )
  sc <- NormalizeData(
    sc,
    assay = "RNA",
    normalization.method = "LogNormalize",
    scale.factor = 10000
  )

  meta <- sc@meta.data |>
    as.data.table(keep.rownames = "barcode")
  celltype_col <- if ("celltype" %in% colnames(meta)) {
    "celltype"
  } else if ("predicted.celltype.l1" %in% colnames(meta)) {
    "predicted.celltype.l1"
  } else {
    stop("No celltype column found in Seurat metadata")
  }
  meta <- meta[, .(barcode, celltype = as.character(get(celltype_col)))]

  expr <- LayerData(sc, assay = "RNA", layer = "data")

  if (!(thegene %in% rownames(expr))) {
    stop(glue("Gene not found in expression matrix: {thegene}"))
  }

  col_idx <- match(.dt$barcode, colnames(expr))
  keep_cells <- !is.na(col_idx)
  if (!any(keep_cells)) {
    stop(
      "No overlapping barcodes between variant table and Seurat expression matrix"
    )
  }
  .dt_filt <- .dt[keep_cells]
  col_idx <- col_idx[keep_cells]

  plot_dt <- data.table(
    barcode = .dt_filt$barcode,
    af = .dt_filt$af,
    expr = as.numeric(expr[thegene, col_idx])
  ) |>
    dplyr::left_join(meta, by = "barcode") |>
    dplyr::filter(!is.na(af), !is.na(expr), expr > 0) |>
    dplyr::mutate(celltype = ifelse(is.na(celltype), "Unknown", celltype)) |>
    as.data.table()

  if (nrow(plot_dt) < 3) {
    stop(glue("Insufficient cells ({nrow(plot_dt)}) for plotting"))
  }

  # Combine ALL cells + per-celltype rows for faceting
  ct_levels <- c("ALL", sort(unique(plot_dt$celltype)))
  plot_combined <- rbind(
    plot_dt[, .(barcode, af, expr, celltype = "ALL")],
    plot_dt[, .(barcode, af, expr, celltype)]
  )
  plot_combined[, celltype := factor(celltype, levels = ct_levels)]

  # Compute Spearman correlation stats per facet
  cor_stats <- plot_combined[,
    {
      if (.N >= 3) {
        ct_res <- cor.test(af, expr, method = "spearman", exact = FALSE)
        .(
          label = glue(
            "rho={round(unname(ct_res$estimate), 3)}\np={signif(ct_res$p.value, 3)}\nn={.N}"
          )
        )
      } else {
        .(label = glue("n={.N}"))
      }
    },
    by = celltype
  ]

  p <- ggplot(plot_combined, aes(x = af, y = expr)) +
    geom_point(alpha = 0.4, size = 0.9, color = "#2B6CB0") +
    geom_smooth(
      method = "lm",
      formula = y ~ x,
      se = TRUE,
      color = "#D94841",
      linewidth = 0.8
    ) +
    geom_text(
      data = cor_stats,
      aes(label = label),
      x = Inf,
      y = Inf,
      hjust = 1.05,
      vjust = 1.3,
      size = 2.8,
      inherit.aes = FALSE
    ) +
    facet_wrap(~celltype, nrow = 1) +
    labs(
      title = glue("{thevariant} vs {thegene}"),
      subtitle = glue("{thegseid} | {thesrrid}"),
      x = glue("{thevariant} AF"),
      y = glue("{thegene} expression (log-normalized)")
    ) +
    theme_bw(base_size = 11) +
    theme(
      strip.text = element_text(size = 9),
      panel.spacing = unit(0.3, "lines")
    )

  return(p)
}


fn_cor_vaf_gene_sc <- function(
  thevariant,
  thegseid,
  thesrrid,
  min_cells = 10
) {
  # thevariant <- "961T>C"
  # thegseid <- "GSE235050"
  # thesrrid <- "GSM7493836"
  # min_cells <- 10
  load_pkg(Seurat)

  conn <- db_conn(
    Sys.getenv("DUCKDB_PATH"),
    readonly = TRUE
  )

  .dt <- tbl(conn, "allvariants_cell") |>
    dplyr::filter(
      srrid == thesrrid,
      variant == thevariant,
      variant_type %in% c("colorful", "black")
    ) |>
    as.data.table() |>
    dplyr::mutate(
      barcode = glue("{thegseid}-{thesrrid}-{barcode}")
    )

  sc <- import(
    "/mnt/isilon/u01_project/large-scale/liuc9/raw/{thegseid}/final/{thesrrid}/for_integration/sc_azimuth.qs" |>
      glue()
  )
  sc <- NormalizeData(
    sc,
    assay = "RNA",
    normalization.method = "LogNormalize",
    scale.factor = 10000
  )

  # --- celltype annotation ---
  meta <- sc@meta.data |>
    as.data.table(keep.rownames = "barcode")
  celltype_col <- if ("celltype" %in% colnames(meta)) {
    "celltype"
  } else if ("predicted.celltype.l1" %in% colnames(meta)) {
    "predicted.celltype.l1"
  } else {
    stop("No celltype column found in Seurat metadata")
  }
  meta <- meta[, .(barcode, celltype = as.character(get(celltype_col)))]

  # --- expression matrix + gene filter (applied on full matrix) ---
  expr <- LayerData(sc, assay = "RNA", layer = "data")
  keep_genes <- rowSums(expr > 0) >= ncol(expr) * 0.05 & rowMeans(expr) > 0.05
  expr <- expr[keep_genes, , drop = FALSE]

  # --- match barcodes from variant table to expression matrix ---
  col_idx <- match(.dt$barcode, colnames(expr))
  keep_match <- !is.na(col_idx)
  if (!any(keep_match)) {
    log_warn("No matched barcodes between variant table and Seurat object")
    return(data.table(
      celltype = character(),
      gene = character(),
      n_cells = integer(),
      cor.rho = numeric(),
      pval = numeric(),
      padj = numeric()
    ))
  }
  .dt_filt <- .dt[keep_match]
  col_idx <- col_idx[keep_match]

  # --- cell-level table with af + celltype ---
  dt_cell <- data.table(
    barcode = .dt_filt$barcode,
    af = .dt_filt$af
  ) |>
    dplyr::left_join(meta, by = "barcode") |>
    dplyr::mutate(celltype = ifelse(is.na(celltype), "Unknown", celltype)) |>
    as.data.table()

  # expression submatrix: genes x matched cells (same order as dt_cell)
  expr_sub <- as.matrix(expr[, col_idx, drop = FALSE])

  # drop cells with missing AF
  valid <- !is.na(dt_cell$af)
  dt_cell <- dt_cell[valid]
  expr_sub <- expr_sub[, valid, drop = FALSE]

  if (ncol(expr_sub) < min_cells) {
    log_warn("Insufficient valid cells ({ncol(expr_sub)}) for correlation")
    return(data.table(
      celltype = character(),
      gene = character(),
      n_cells = integer(),
      cor.rho = numeric(),
      pval = numeric(),
      padj = numeric()
    ))
  }

  # --- helper: Spearman cor for one group across all genes (genes in parallel) ---
  run_cor_group <- function(group_name, cell_idx) {
    n <- length(cell_idx)
    if (n < min_cells) {
      return(NULL)
    }

    vaf <- dt_cell$af[cell_idx]
    mat <- expr_sub[, cell_idx, drop = FALSE]

    if (stats::var(vaf) == 0) {
      log_warn("Zero variance in AF for group {group_name}, skipping")
      return(NULL)
    }

    # drop genes with zero variance in this subset
    gene_var <- apply(mat, 1, stats::var)
    mat <- mat[gene_var > 0, , drop = FALSE]
    if (nrow(mat) == 0) {
      return(NULL)
    }

    gene_names <- rownames(mat)

    # correlate each gene in parallel
    cor_list <- pbmclapply(
      gene_names,
      function(g) {
        gene_expr <- mat[g, ]
        ct_res <- cor.test(vaf, gene_expr, method = "spearman", exact = FALSE)
        data.table(
          gene = g,
          cor.rho = unname(ct_res$estimate),
          pval = ct_res$p.value
        )
      },
      mc.cores = nthread
    )

    res_mat <- rbindlist(cor_list)
    res_mat[, `:=`(
      celltype = group_name,
      n_cells = n,
      padj = p.adjust(pval, method = "BH")
    )]
    res_mat[, .(celltype, gene, n_cells, cor.rho, pval, padj)]
  }

  # --- run ALL cells first, then per-celltype sequentially ---
  groups <- c(
    list(ALL = seq_len(nrow(dt_cell))),
    setNames(
      lapply(sort(unique(dt_cell$celltype)), function(ct) {
        which(dt_cell$celltype == ct)
      }),
      sort(unique(dt_cell$celltype))
    )
  )

  res_list <- lapply(names(groups), function(grp) {
    run_cor_group(grp, groups[[grp]])
  })

  res <- rbindlist(res_list, use.names = TRUE, fill = TRUE)
  res <- res[!is.na(gene)]
  setcolorder(res, c("celltype", "gene", "n_cells", "cor.rho", "pval", "padj"))
  res <- res |>
    filter(abs(cor.rho) >= 0.3 & padj < 0.05) |>
    arrange(desc(abs(cor.rho)))
  return(res)
}

# GO enrichment for a character vector of gene SYMBOLs (over-representation)
fn_enrichGO_symbols <- function(gene_symbols, universe_symbols = NULL) {
  if (length(gene_symbols) < 3) {
    return(NULL)
  }
  load_pkg(clusterProfiler, org.Hs.eg.db)

  gene_ids <- suppressMessages(
    clusterProfiler::bitr(
      geneID = gene_symbols,
      fromType = "SYMBOL",
      toType = "ENTREZID",
      OrgDb = org.Hs.eg.db::org.Hs.eg.db
    )
  )
  if (nrow(gene_ids) == 0) {
    return(NULL)
  }

  universe_ids <- if (
    !is.null(universe_symbols) && length(universe_symbols) >= 10
  ) {
    suppressMessages(
      clusterProfiler::bitr(
        geneID = universe_symbols,
        fromType = "SYMBOL",
        toType = "ENTREZID",
        OrgDb = org.Hs.eg.db::org.Hs.eg.db
      )$ENTREZID
    )
  } else {
    NULL
  }

  onts <- c("BP", "CC", "MF")
  setNames(
    lapply(onts, function(.ont) {
      tryCatch(
        clusterProfiler::enrichGO(
          gene = gene_ids$ENTREZID,
          universe = universe_ids,
          OrgDb = org.Hs.eg.db::org.Hs.eg.db,
          keyType = "ENTREZID",
          ont = .ont,
          pvalueCutoff = 0.05,
          pAdjustMethod = "BH",
          readable = TRUE
        ),
        error = function(e) NULL
      )
    }),
    onts
  )
}

# Plot a single enrichGO result as a horizontal bar chart
fn_plot_go <- function(
  .go,
  .topn = 20,
  .ont = c("BP", "CC", "MF"),
  .title = NULL
) {
  .ont <- match.arg(.ont)

  if (is.null(.go) || nrow(as.data.frame(.go)) == 0) {
    return(NULL)
  }

  base_fill <- c("BP" = "#AE1700", "CC" = "#DF8F44FF", "MF" = "#00A1D5FF")
  ont_fullname <- c(
    "BP" = "Biological Process",
    "CC" = "Cellular Component",
    "MF" = "Molecular Function"
  )
  .ont_fill <- base_fill[.ont]
  x_label <- ont_fullname[.ont]

  tryCatch(
    {
      .go |>
        tibble::as_tibble() |>
        dplyr::mutate(
          Description = stringr::str_wrap(
            stringr::str_to_sentence(string = Description),
            width = 60
          )
        ) |>
        dplyr::mutate(adjp = -log10(p.adjust)) |>
        dplyr::select(ID, Description, adjp, Count, geneID) |>
        dplyr::arrange(adjp, Count) |>
        dplyr::mutate(
          Description = factor(Description, levels = Description)
        ) -> .go_for_plot

      if (!is.infinite(.topn)) {
        .go_for_plot <- tail(.go_for_plot, .topn)
      }
      if (nrow(.go_for_plot) == 0) {
        return(NULL)
      }

      p <- .go_for_plot |>
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
          axis.text.y = element_text(color = "black", size = 13, hjust = 1),
          axis.ticks.length.y = unit(3, units = "mm"),
          axis.text.x = element_text(color = "black", size = 12),
          axis.title = element_text(colour = "black", size = 16, face = "bold")
        )

      if (!is.null(.title)) {
        p <- p +
          labs(title = .title) +
          theme(plot.title = element_text(size = 12, face = "bold"))
      }
      p
    },
    error = function(e) NULL
  )
}

# Main --------------------------------------------------------------------
variant_annotation <- "MITOMAP reported disease"

variant_df <- variant_list[[variant_annotation]]

variant_disease_reported <- c(
  "961T>C",
  "13271T>C",
  "4175G>A",
  "14831G>A",
  "9804G>A",
  "1382A>C",
  "9025G>A",
  "8362T>G",
  "3667T>G",
  "7065G>A"
)


variant_df |>
  left_join(
    ALLVARIANTS |>
      left_join(
        METAFULL |>
          select(srrid, disease),
        by = c("srrid")
      ) |>
      mutate(
        variant_type = factor(variant_type, c("hete", "homo")),
        disease = factor(
          disease,
          levels = names(color_disease)
        )
      ) |>
      arrange(variant_type, disease) |>
      nest(.by = "variant", .key = "meta"),
    by = "variant"
  ) -> variant_df_meta


# variant_df |>
#   filter(variant == "13271T>C") |>
#   left_join(
#     ALLVARIANTS |>
#       left_join(
#         METAFULL |>
#           select(srrid, disease),
#         by = c("srrid")
#       ) |>
#       mutate(
#         variant_type = factor(variant_type, c("hete", "homo")),
#         disease = factor(
#           disease,
#           levels = names(color_disease)
#         )
#       ) |>
#       arrange(variant_type, disease) |>
#       nest(.by = "variant", .key = "meta"),
#     by = "variant"
#   ) -> a

# a$meta[[1]]

# fn_cor_vaf_gene_sc(
#   thevariant = "13271T>C",
#   thegseid = "GSE226602",
#   thesrrid = "GSM7080038"
# ) -> res_cor

# res_cor |> filter(celltype != "other")

# fn_plot_cor_vaf_gene_sc(
#   thevariant = "13271T>C",
#   thegseid = "GSE226602",
#   thesrrid = "GSM7080038",
#   thegene = "C21orf91"
# )

variant_pathogenic <- c(
  "3572G>G",
  "3727T>C",
  "3728C>T",
  "5343T>C",
  "6967G>A",
  "7583T>G",
  "1111T>G",
  "14520C>A",
  "15666T>C"
)

variant_vus <- c(
  "1670A>G",
  "2636G>A",
  "3734A>G",
  "4864C>T",
  "6190G>A",
  "6667C>A",
  "6668C>G",
  "6997T>G",
  "7847G>A",
  "8361G>T",
  "8849T>C",
  "11372G>A",
  "11708A>T"
)

thevariants <- c(
  variant_disease_reported,
  variant_pathogenic,
  variant_vus
)

# labeled data.table of all variants with group
thevariants_labeled <- rbindlist(list(
  data.table(variant = variant_disease_reported, group = "disease_reported"),
  data.table(variant = variant_pathogenic, group = "pathogenic"),
  data.table(variant = variant_vus, group = "vus")
))

# Build variant × sample map from ALLVARIANTS (has gseid + srrid columns)
variant_sample_map <- ALLVARIANTS |>
  dplyr::filter(variant %in% thevariants) |>
  dplyr::select(variant, gseid, srrid) |>
  dplyr::distinct() |>
  dplyr::left_join(thevariants_labeled, by = "variant") |>
  as.data.table()

log_info("Total variant-sample combinations: {nrow(variant_sample_map)}")

# Output directory
cor_outdir <- outdirnotuse /
  "allvariants-prioritize" /
  "17.02-variant-gene-correlation"
fs::dir_create(cor_outdir)

# Run correlation for every (variant, gseid, srrid) -----------------------
all_cor_results <- lapply(
  seq_len(nrow(variant_sample_map)),
  function(i) {
    .row <- variant_sample_map[i]
    log_info(
      "[{i}/{nrow(variant_sample_map)}] {.row$variant} | {.row$gseid} | {.row$srrid}"
    )
    tryCatch(
      {
        res <- fn_cor_vaf_gene_sc(
          thevariant = .row$variant,
          thegseid = .row$gseid,
          thesrrid = .row$srrid
        )
        if (nrow(res) > 0) {
          res[, `:=`(
            variant = .row$variant,
            gseid = .row$gseid,
            srrid = .row$srrid,
            group = .row$group
          )]
        }
        res
      },
      error = function(e) {
        log_error(
          "Correlation failed for {.row$variant} | {.row$srrid}: {e$message}"
        )
        NULL
      }
    )
  }
)

all_cor_dt <- rbindlist(all_cor_results, use.names = TRUE, fill = TRUE)

export(
  all_cor_dt,
  cor_outdir / "all-variant-gene-correlation.qs"
)

all_cor_dt |>
  nest(.by = c("variant", "group", "gseid", "srrid"), .key = "cor_results") |>
  arrange(group, variant) |>
  left_join(
    METAFULL |>
      select(srrid, disease, status, SEXPRED, Chemistry, Haplogroup) |>
      distinct() |>
      mutate(
        disease = factor(
          disease,
          levels = names(color_disease)
        )
      ),
  ) |>
  mutate(
    n = pbmclapply(
      cor_results,
      \(.x) {
        .xpos <- .x |> filter(cor.rho > 0.3 & padj < 0.05)
        .xneg <- .x |> filter(cor.rho < -0.3 & padj < 0.05)
        tibble(
          n_total = nrow(.x),
          n_pos = nrow(.xpos),
          n_neg = nrow(.xneg),
          gene_pos = list(.xpos$gene),
          gene_neg = list(.xneg$gene)
        )
      },
      mc.cores = 8
    )
  ) |>
  unnest(n) -> all_cor_dt_nested
# all_cor_dt_nested <- import(
#   cor_outdir / "all-variant-gene-correlation-nested.qs"
# )
export(
  all_cor_dt_nested,
  cor_outdir / "all-variant-gene-correlation-nested.qs"
)

# GO enrichment for pos and neg correlated genes per (variant, gseid, srrid) -----
log_info("Running GO enrichment for pos/neg correlated genes")

# universe = all genes tested in that sample (from cor_results, before significance filter)
all_cor_dt_universe <- all_cor_dt_nested |>
  dplyr::mutate(
    universe_genes = lapply(cor_results, function(.x) unique(.x$gene))
  )

all_cor_dt_enrich <- all_cor_dt_universe |>
  dplyr::mutate(
    enrich_pos = pbmclapply(
      seq_len(nrow(all_cor_dt_universe)),
      function(i) {
        fn_enrichGO_symbols(
          gene_symbols = all_cor_dt_universe$gene_pos[[i]],
          universe_symbols = all_cor_dt_universe$universe_genes[[i]]
        )
      },
      mc.cores = nthread
    ),
    enrich_neg = pbmclapply(
      seq_len(nrow(all_cor_dt_universe)),
      function(i) {
        fn_enrichGO_symbols(
          gene_symbols = all_cor_dt_universe$gene_neg[[i]],
          universe_symbols = all_cor_dt_universe$universe_genes[[i]]
        )
      },
      mc.cores = nthread
    )
  )

export(
  all_cor_dt_enrich,
  cor_outdir / "all-variant-gene-correlation-enrich.qs"
)
log_info("GO enrichment done, saved to {cor_outdir}")

# Build GO plot list column: pos/neg × BP, CC, MF per (variant, srrid) --------
log_info("Building GO enrichment plots")

all_cor_dt_enrich_plots <- all_cor_dt_enrich |>
  dplyr::mutate(
    go_plots = pbmclapply(
      seq_len(nrow(all_cor_dt_enrich)),
      function(i) {
        .epos <- all_cor_dt_enrich$enrich_pos[[i]]
        .eneg <- all_cor_dt_enrich$enrich_neg[[i]]
        .variant <- all_cor_dt_enrich$variant[[i]]
        .srrid <- all_cor_dt_enrich$srrid[[i]]
        .label <- glue("{.variant} | {.srrid}")

        plots <- list(
          pos_BP = fn_plot_go(
            .epos[["BP"]],
            20,
            "BP",
            glue("Pos corr \u2013 Biological Process: {.label}")
          ),
          pos_CC = fn_plot_go(
            .epos[["CC"]],
            20,
            "CC",
            glue("Pos corr \u2013 Cellular Component: {.label}")
          ),
          pos_MF = fn_plot_go(
            .epos[["MF"]],
            20,
            "MF",
            glue("Pos corr \u2013 Molecular Function: {.label}")
          ),
          neg_BP = fn_plot_go(
            .eneg[["BP"]],
            20,
            "BP",
            glue("Neg corr \u2013 Biological Process: {.label}")
          ),
          neg_CC = fn_plot_go(
            .eneg[["CC"]],
            20,
            "CC",
            glue("Neg corr \u2013 Cellular Component: {.label}")
          ),
          neg_MF = fn_plot_go(
            .eneg[["MF"]],
            20,
            "MF",
            glue("Neg corr \u2013 Molecular Function: {.label}")
          )
        )
        # keep only non-NULL panels
        Filter(Negate(is.null), plots)
      },
      mc.cores = nthread
    )
  )

export(
  all_cor_dt_enrich_plots,
  cor_outdir / "all-variant-gene-correlation-enrich-plots.qs"
)
log_info("GO enrichment plot object saved to {cor_outdir}")

# Save GO plot PDFs -------------------------------------------------------
# Directory layout: go-plots/{group}/{safe_variant}/{safe_variant}-{srrid}-GO.pdf
go_plot_dir <- cor_outdir / "go-plots"
log_info("Saving GO plot PDFs to {go_plot_dir}")

pbmclapply(
  seq_len(nrow(all_cor_dt_enrich_plots)),
  function(i) {
    .plots <- all_cor_dt_enrich_plots$go_plots[[i]]
    if (length(.plots) == 0) {
      return(invisible(NULL))
    }

    .variant <- all_cor_dt_enrich_plots$variant[[i]]
    .srrid <- all_cor_dt_enrich_plots$srrid[[i]]
    .group <- all_cor_dt_enrich_plots$group[[i]]
    safe_variant <- gsub("[^A-Za-z0-9_-]", "_", .variant)

    tryCatch(
      saveplot(
        plot = .plots,
        filename = go_plot_dir /
          .group /
          safe_variant /
          glue("{safe_variant}-{.srrid}-GO.pdf"),
        device = "pdf",
        width = 10,
        height = 7,
        create.dir = TRUE
      ),
      error = function(e) {
        log_error("GO plot save failed for {.variant} | {.srrid}: {e$message}")
      }
    )
  },
  mc.cores = nthread
)
log_info("GO plot PDFs saved to {go_plot_dir}")


#
#
# ? don't run below --------------------------------------------------------------------
#
#

all_cor_dt_nested |>
  arrange(-n_total) |>
  filter(variant == "8362T>G") |>
  arrange(-n_total) -> m

m$gene_neg[[1]]
m$gene_neg[[2]]
ggvenn::ggvenn(
  list(
    "8362T>G neg 1" = m$gene_pos[[1]],
    "8362T>G neg 2" = m$gene_pos[[2]],
    "8362T>G neg 3" = m$gene_pos[[3]],
    "8362T>G neg 4" = m$gene_pos[[4]]
  ),
  fill_color = c("#D94841", "#2B6CB0", "#D94841", "#2B6CB0"),
  stroke_size = 0.5,
  set_name_size = 4
) +
  theme(legend.position = "none")

#
#
# ? below is  --------------------------------------------------------------------
#
#

log_info(
  "Saved {nrow(all_cor_dt_nested )} significant correlations to {cor_outdir}"
)

# Scatter plots: one plot per (variant, gseid, srrid) using the top gene ---
plot_jobs <- all_cor_dt |>
  dplyr::filter(celltype == "ALL", !is.na(gene)) |>
  dplyr::group_by(variant, gseid, srrid, group) |>
  dplyr::slice_max(order_by = abs(cor.rho), n = 1, with_ties = FALSE) |>
  dplyr::ungroup() |>
  as.data.table()

log_info("Generating {nrow(plot_jobs)} scatter plots")

lapply(seq_len(nrow(plot_jobs)), function(i) {
  .row <- plot_jobs[i]
  log_info(
    "[{i}/{nrow(plot_jobs)}] variant={.row$variant} srrid={.row$srrid} gene={.row$gene}"
  )
  tryCatch(
    {
      p <- fn_plot_cor_vaf_gene_sc(
        thevariant = .row$variant,
        thegseid = .row$gseid,
        thesrrid = .row$srrid,
        thegene = .row$gene
      )
      safe_variant <- gsub("[^A-Za-z0-9_-]", "_", .row$variant)
      saveplot(
        plot = p,
        filename = cor_outdir /
          .row$group /
          glue("{safe_variant}-{.row$srrid}-{.row$gene}-scatter.pdf"),
        device = "pdf",
        width = 14,
        height = 4,
        create.dir = TRUE
      )
    },
    error = function(e) {
      log_error(
        "Plot failed for {.row$variant} | {.row$srrid}: {e$message}"
      )
    }
  )
})

# Save  --------------------------------------------------------------

# Session info -------------------------------------------------------------
if (isTRUE(verbose)) {
  sessionInfo()
}
