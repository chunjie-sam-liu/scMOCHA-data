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
  min_cells = 10,
  nthread = 8
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


fn_load_vaf <- function(
  thevariant,
  thegseid,
  thesrrid,
  min_cells = 10,
  nthread = 8
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
    ) |>
    select(
      barcode,
      af_cell = af,
      depth,
      variant_type_cell = variant_type,
      celltype
    )
  return(.dt)
}


# volcano plot for Seurat FindMarkers output (rownames = gene)
fn_de_plot <- function(
  markers,
  .cutoff_pval = 0.05,
  .cutoff_log2fc = 0.25,
  .pct = 0.05
) {
  markers |>
    tibble::rownames_to_column("gene") |>
    dplyr::mutate(fdr = -log10(p_val_adj)) |>
    dplyr::mutate(
      fdr = ifelse(fdr > -log10(1e-300), -log10(1e-300), fdr)
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
          avg_log2FC > .cutoff_log2fc ~ "red",
        p_val_adj < .cutoff_pval &
          (pct.1 >= .pct | pct.2 >= .pct) &
          avg_log2FC < -.cutoff_log2fc ~ "blue",
        TRUE ~ "grey"
      )
    ) -> forplot

  forplot |>
    dplyr::count(color) |>
    tibble::deframe() -> n_color

  p <- forplot |>
    ggplot(aes(x = avg_log2FC, y = fdr, color = color)) +
    geom_point(alpha = 0.7, size = 0.9) +
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
      xintercept = c(-.cutoff_log2fc, .cutoff_log2fc),
      linetype = "dashed",
      color = "grey50"
    ) +
    geom_hline(
      yintercept = -log10(.cutoff_pval),
      linetype = "dashed",
      color = "grey50"
    ) +
    theme_classic(base_size = 12) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 13),
      plot.subtitle = element_text(hjust = 0.5, color = "black", size = 11)
    )

  list(p = p, markers = forplot)
}


# GO enrichment for a character vector of gene SYMBOLs (over-representation)
fn_enrichGO_symbols <- function(gene_symbols, universe_symbols = NULL) {
  if (length(gene_symbols) < 3) {
    return(NULL)
  }
  load_pkg(clusterProfiler, org.Hs.eg.db)
  # resolve path conflict introduced by BiocGenerics (loaded via org.Hs.eg.db)
  conflicted::conflicts_prefer(fs::path, dplyr::filter)

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


# horizontal bar chart for a single enrichGO result
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


# DEG between two disease groups in variant-carrying cells
# a: sample-level data frame with af_cell list column (output of mutate/pbmcmapply)
# disease1: reference (e.g. "Healthy"), disease2: comparison (e.g. "Alzheimer's Disease")
# positive log2FC = higher in disease2
fn_deg_disease_compare <- function(
  a,
  disease1 = "Healthy",
  disease2 = "Alzheimer's Disease",
  the_variant = "8362T>G",
  min_cells = 20,
  nthread = 8
) {
  load_pkg(Seurat, Matrix)
  conflicted::conflict_prefer("filter", "dplyr")

  a_sub <- a |>
    dplyr::filter(disease %in% c(disease1, disease2))

  log_info(
    "Loading Seurat for DEG: {disease1} vs {disease2}, {nrow(a_sub)} samples"
  )

  sc_list <- list()

  for (i in seq_len(nrow(a_sub))) {
    gseid_i <- a_sub$gseid[i]
    srrid_i <- a_sub$srrid[i]
    disease_i <- as.character(a_sub$disease[i])
    cell_dt <- a_sub$af_cell[[i]]

    if (is.null(cell_dt) || nrow(cell_dt) == 0) {
      next
    }

    sc_path <- glue(
      "/mnt/isilon/u01_project/large-scale/liuc9/raw/{gseid_i}/final/{srrid_i}/for_integration/sc_azimuth.qs"
    )
    if (!file.exists(sc_path)) {
      log_warn("Not found: {sc_path}")
      next
    }

    sc <- tryCatch(import(sc_path), error = function(e) {
      log_warn("Failed to load {srrid_i}: {e$message}")
      NULL
    })
    if (is.null(sc)) {
      next
    }

    sc <- NormalizeData(
      sc,
      assay = "RNA",
      normalization.method = "LogNormalize",
      scale.factor = 10000
    )

    # barcodes in sc_azimuth already carry the gseid-srrid- prefix
    valid_bc <- intersect(cell_dt$barcode, colnames(sc))

    if (length(valid_bc) < 3) {
      log_warn("Too few variant cells in {srrid_i}: {length(valid_bc)}")
      next
    }

    sc_sub <- sc[, valid_bc]
    sc_sub$disease_group <- disease_i
    sc_sub$srrid_orig <- srrid_i
    sc_list[[length(sc_list) + 1]] <- sc_sub
    log_info("{srrid_i}: {length(valid_bc)} variant cells ({disease_i})")
  }

  if (length(sc_list) == 0) {
    log_warn("No Seurat objects loaded")
    return(NULL)
  }

  # check cell counts per disease
  n_d1 <- sum(
    sapply(sc_list, function(sc) sum(sc$disease_group == disease1))
  )
  n_d2 <- sum(
    sapply(sc_list, function(sc) sum(sc$disease_group == disease2))
  )
  log_info(
    "Variant-carrying cells: {disease1}={n_d1}, {disease2}={n_d2}"
  )

  if (n_d1 < min_cells || n_d2 < min_cells) {
    log_warn(
      "Insufficient cells for DEG (min={min_cells}): {disease1}={n_d1}, {disease2}={n_d2}"
    )
    return(NULL)
  }

  # merge all Seurat subsets
  merged_sc <- if (length(sc_list) == 1) {
    sc_list[[1]]
  } else {
    merge(sc_list[[1]], y = sc_list[-1], merge.data = TRUE) |>
      JoinLayers()
  }

  merged_sc <- FindVariableFeatures(merged_sc, nfeatures = 3000)
  Idents(merged_sc) <- merged_sc$disease_group

  markers <- tryCatch(
    Seurat::FindMarkers(
      object = merged_sc,
      ident.1 = disease2,
      ident.2 = disease1,
      assay = "RNA",
      test.use = "wilcox",
      min.pct = 0.1,
      logfc.threshold = 0.1
    ),
    error = function(e) {
      log_error("FindMarkers failed: {e$message}")
      NULL
    }
  )

  if (is.null(markers) || nrow(markers) == 0) {
    log_warn("No markers found for {disease2} vs {disease1}")
    return(NULL)
  }

  n_sig <- sum(
    markers$p_val_adj < 0.05 & abs(markers$avg_log2FC) >= 0.25
  )
  log_info(
    "FindMarkers: {nrow(markers)} genes tested, {n_sig} significant (padj<0.05, |log2FC|>=0.25)"
  )

  list(
    markers = markers,
    n_cells = c(setNames(n_d1, disease1), setNames(n_d2, disease2)),
    disease1 = disease1,
    disease2 = disease2,
    sc = merged_sc
  )
}


plot_af_violin_cell <- function(cell_data, the_variant, label = "cell") {
  suppressMessages({
    load_pkg(ggh4x, ggbeeswarm, ggsignif)
  })

  disease_order <- c("Healthy", "COVID-19", "Alzheimer's Disease")

  forplot <- cell_data |>
    dplyr::filter(variant == the_variant) |>
    dplyr::mutate(
      celltype = gsub(celltype, pattern = "_", replacement = " "),
      celltype = factor(celltype, levels = names(color_celltype))
    ) |>
    dplyr::filter(
      !is.na(celltype),
      !is.na(af_cell),
      disease %in% disease_order
    ) |>
    dplyr::mutate(
      disease = factor(disease, levels = disease_order)
    )

  # strip fill only for celltypes present in the data
  ct_present <- levels(droplevels(forplot$celltype))
  ct_colors <- color_celltype[ct_present]

  ymax <- max(forplot$af_cell, na.rm = TRUE)

  ggplot(forplot, aes(x = disease)) +
    ggh4x::facet_grid2(
      ~celltype,
      strip = ggh4x::strip_themed(
        background_x = ggh4x::elem_list_rect(
          fill = ct_colors,
          color = NA
        ),
        text_x = ggh4x::elem_list_text(colour = "white", face = "bold")
      )
    ) +
    geom_violin(
      aes(y = af_cell, fill = disease),
      alpha = 0.7,
      color = NA
    ) +
    scale_fill_manual(values = color_disease, name = "Disease") +
    ggbeeswarm::geom_quasirandom(
      aes(y = af_cell, color = disease),
      size = 0.8,
      dodge.width = 0.75,
      alpha = 0.8,
      show.legend = FALSE
    ) +
    scale_color_manual(values = color_disease, name = "Disease") +
    ggsignif::geom_signif(
      aes(y = af_cell),
      comparisons = list(
        c("Healthy", "COVID-19"),
        c("Healthy", "Alzheimer's Disease")
      ),
      y_position = c(ymax * 1.05, ymax * 1.15),
      tip_length = 0.01,
      textsize = 3
    ) +
    theme(
      plot.margin = margin(t = 0.2, b = 0.1, l = 0.1, r = 0.2, unit = "cm"),
      panel.background = element_rect(
        fill = NA,
        color = "black",
        linewidth = 0.5
      ),
      panel.grid = element_blank(),
      axis.line.y.left = element_line(color = "black"),
      legend.position = "top",
      legend.key = element_blank(),
      axis.title.y = element_text(color = "black", size = 14),
      axis.text.y = element_text(color = "black"),
      legend.text = element_text(size = 12, color = "black"),
      legend.title = element_text(size = 14, colour = "black"),
      strip.background = element_blank(),
      strip.text = element_text(size = 8, color = "black", face = "bold"),
      axis.ticks.x = element_blank(),
      axis.text.x = element_blank(),
      axis.title.x = element_blank()
    ) +
    labs(
      title = paste(the_variant, "AF distribution —", label, "level"),
      y = "Allele Frequency (cell)",
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


variant_df |>
  filter(variant == "8362T>G") |>
  left_join(
    ALLVARIANTS |>
      left_join(
        METAFULL |>
          select(srrid, disease, Chemistry),
        by = c("srrid")
      ) |>
      mutate(
        variant_type = factor(variant_type, c("hete", "homo")),
        disease = factor(
          disease,
          levels = names(color_disease)
        ),
        Chemistry = factor(
          Chemistry,
          levels = names(color_chemistry)
        )
      ) |>
      arrange(variant_type, disease),
    by = "variant"
  ) |>
  mutate(
    af_cell = pbmcmapply(
      FUN = fn_load_vaf,
      thevariant = variant,
      thegseid = gseid,
      thesrrid = srrid,
      mc.cores = 8,
      SIMPLIFY = FALSE
    )
  ) -> a

a |> count(variant_type, disease, Haplogroup) |> print(n = Inf)


a |>
  unnest(af_cell) -> a_cell


plot_af_violin_cell(a_cell, "8362T>G", "cell") +
  labs(
    title = "8362T>G AF distribution in single cells",
    subtitle = "Faceted by cell type, colored by disease"
  ) -> p_m.8362T_G_cell

saveplot(
  p_m.8362T_G_cell,
  filename = "/home/liuc9/github/scMOCHA-data/high-res-MANUSCRIPTFIGURES-notuse/allvariants-prioritize/17.02-variant-gene-correlation/AD-8362T-G-af-violin-cell.pdf",
  width = 12,
  height = 6
)


# helper: compute or load cached DEG + GO for one disease comparison ----------
fn_run_deg_and_go <- function(
  a,
  disease1,
  disease2,
  cache_path,
  the_variant = "8362T>G",
  min_cells = 20,
  nthread = 8
) {
  if (file.exists(cache_path)) {
    log_info("Loading cached DEG+GO from {cache_path}")
    cached <- import(cache_path)
    return(cached)
  }

  deg <- fn_deg_disease_compare(
    a,
    disease1 = disease1,
    disease2 = disease2,
    the_variant = the_variant,
    min_cells = min_cells,
    nthread = nthread
  )

  if (is.null(deg)) {
    log_warn("DEG returned NULL for {disease2} vs {disease1}, skipping")
    return(NULL)
  }

  markers <- deg$markers
  n_cells <- deg$n_cells
  universe <- rownames(markers)

  sig <- markers |>
    tibble::rownames_to_column("gene") |>
    dplyr::filter(
      p_val_adj < 0.05,
      abs(avg_log2FC) >= 0.25,
      pct.1 >= 0.05 | pct.2 >= 0.05
    )

  pos_genes <- sig |> dplyr::filter(avg_log2FC > 0) |> dplyr::pull(gene)
  neg_genes <- sig |> dplyr::filter(avg_log2FC < 0) |> dplyr::pull(gene)

  log_info(
    "DEGs {disease2} vs {disease1}: {length(pos_genes)} up, {length(neg_genes)} down"
  )

  go_pos <- fn_enrichGO_symbols(pos_genes, universe)
  go_neg <- fn_enrichGO_symbols(neg_genes, universe)

  result <- list(
    markers = markers,
    n_cells = n_cells,
    go_pos = go_pos,
    go_neg = go_neg,
    disease1 = disease1,
    disease2 = disease2
  )
  export(result, cache_path)
  log_info("DEG+GO cached to {cache_path}")
  result
}


# helper: save volcano + GO bar charts for one DEG result --------------------
fn_plot_deg_result <- function(
  result,
  label, # short label, e.g. "AD" or "COVID"
  outdirnotuse,
  variant_label = "8362T>G"
) {
  if (is.null(result)) {
    return(invisible(NULL))
  }

  d1 <- result$disease1
  d2 <- result$disease2
  n1 <- result$n_cells[d1]
  n2 <- result$n_cells[d2]
  n1_str <- if (!is.na(n1)) as.character(n1) else "N/A"
  n2_str <- if (!is.na(n2)) as.character(n2) else "N/A"

  outdir_corr <- fs::path(
    outdirnotuse,
    "allvariants-prioritize/17.02-variant-gene-correlation"
  )
  fs::dir_create(outdir_corr)

  # volcano
  p_volcano <- fn_de_plot(result$markers)
  p_volcano$p <- p_volcano$p +
    labs(
      title = glue("{variant_label}: {d2} vs {d1} (variant-carrying cells)"),
      subtitle = glue("{d2} cells: {n2_str},  {d1} cells: {n1_str}"),
      x = glue("avg log\u2082FC  ({d2} vs {d1})"),
      y = "-log\u2081\u2080(FDR)"
    )

  saveplot(
    p_volcano$p,
    filename = fs::path(
      outdir_corr,
      glue("AD-8362T-G-deg-{label}-vs-Healthy-volcano.pdf")
    ),
    width = 10,
    height = 6,
    device = "pdf"
  )

  # GO plots: up and down × 3 ontologies
  purrr::walk(c("pos", "neg"), function(.dir) {
    go_list <- result[[glue("go_{.dir}")]]
    dir_label <- if (.dir == "pos") {
      glue("UP in {d2}")
    } else {
      glue("DOWN in {d2} (up in {d1})")
    }
    dir_file <- if (.dir == "pos") "up" else "down"

    purrr::walk(c("BP", "CC", "MF"), function(.ont) {
      p_go <- fn_plot_go(
        go_list[[.ont]],
        .topn = 20,
        .ont = .ont,
        .title = glue("{variant_label} {dir_label} \u2013 {.ont}")
      )
      if (!is.null(p_go)) {
        saveplot(
          p_go,
          filename = fs::path(
            outdir_corr,
            glue("AD-8362T-G-go-{label}-{dir_file}-{tolower(.ont)}.pdf")
          ),
          width = 10,
          height = 8,
          device = "pdf"
        )
      }
    })
  })
  invisible(NULL)
}


# DEG analysis: Healthy vs Alzheimer's Disease in 8362T>G variant cells --------
.cache_deg_AD <- path(
  outdirnotuse,
  "allvariants-prioritize/17.02-variant-gene-correlation/AD-8362T-G-deg-AD-vs-Healthy.qs"
)
result_deg_AD <- fn_run_deg_and_go(
  a,
  disease1 = "Healthy",
  disease2 = "Alzheimer's Disease",
  cache_path = .cache_deg_AD,
  the_variant = "8362T>G",
  min_cells = 20,
  nthread = nthread
)
fn_plot_deg_result(result_deg_AD, label = "AD", outdirnotuse = outdirnotuse)


# DEG analysis: Healthy vs COVID-19 in 8362T>G variant cells ------------------
.cache_deg_COVID <- path(
  outdirnotuse,
  "allvariants-prioritize/17.02-variant-gene-correlation/AD-8362T-G-deg-COVID-vs-Healthy.qs"
)
result_deg_COVID <- fn_run_deg_and_go(
  a,
  disease1 = "Healthy",
  disease2 = "COVID-19",
  cache_path = .cache_deg_COVID,
  the_variant = "8362T>G",
  min_cells = 20,
  nthread = nthread
)
fn_plot_deg_result(
  result_deg_COVID,
  label = "COVID",
  outdirnotuse = outdirnotuse
)


# Save  --------------------------------------------------------------

# Session info -------------------------------------------------------------
if (isTRUE(verbose)) {
  sessionInfo()
}
