#!/usr/bin/env Rscript
# Metainfo ----------------------------------------------------------------
# @AUTHOR: Chun-Jie Liu
# @CONTACT: chunjie.sam.liu.at.gmail.com
# @DATE: 2026-04-01 15:31:27
# @DESCRIPTION: Compare high-AF versus low-AF variant-carrying cells by celltype across multiple AF cutoffs.
# @VERSION: v0.0.1

# Reproducibility ----------------------------------------------------------
set.seed(1)

# Library ------------------------------------------------------------------
suppressMessages({
  load_pkg(jutils)
})


# Helpers ------------------------------------------------------------------
fn_assign_af_group <- function(af_values, cutoff) {
  dplyr::case_when(
    af_values > cutoff ~ "high_af",
    af_values < cutoff ~ "low_af",
    TRUE ~ NA_character_
  )
}


fn_cutoff_dir <- function(cutoff) {
  as.character(glue("cutoff-{format(cutoff, nsmall = 1, trim = TRUE)}"))
}


fn_safe_name <- function(x) {
  gsub("[^A-Za-z0-9]+", "_", x)
}


fn_load_vaf <- function(thevariant, thegseid, thesrrid) {
  conn <- db_conn(
    Sys.getenv("DUCKDB_PATH"),
    readonly = TRUE
  )
  on.exit(DBI::dbDisconnect(conn), add = TRUE)

  tbl(conn, "allvariants_cell") |>
    dplyr::filter(
      srrid == thesrrid,
      variant == thevariant,
      variant_type %in% c("colorful", "black")
    ) |>
    as.data.table() |>
    dplyr::mutate(
      barcode = glue("{thegseid}-{thesrrid}-{barcode}")
    ) |>
    dplyr::select(
      barcode,
      af_cell = af,
      depth,
      variant_type_cell = variant_type,
      celltype
    ) |>
    as.data.table()
}


fn_enrichGO_symbols <- function(gene_symbols, universe_symbols = NULL) {
  if (length(gene_symbols) < 3) {
    return(NULL)
  }

  load_pkg(clusterProfiler, org.Hs.eg.db)
  conflicted::conflicts_prefer(fs::path, dplyr::filter)

  gene_symbols_for_mapping <- gsub("^MT-", "", gene_symbols)
  universe_symbols_for_mapping <- gsub("^MT-", "", universe_symbols)

  gene_ids <- suppressMessages(tryCatch(
    clusterProfiler::bitr(
      geneID = gene_symbols_for_mapping,
      fromType = "SYMBOL",
      toType = "ENTREZID",
      OrgDb = org.Hs.eg.db::org.Hs.eg.db
    ),
    error = function(e) data.frame(SYMBOL = character(), ENTREZID = character())
  ))
  if (is.null(gene_ids) || nrow(gene_ids) == 0) {
    return(NULL)
  }

  universe_ids <- if (
    !is.null(universe_symbols) && length(universe_symbols) >= 10
  ) {
    suppressMessages(tryCatch(
      clusterProfiler::bitr(
        geneID = universe_symbols_for_mapping,
        fromType = "SYMBOL",
        toType = "ENTREZID",
        OrgDb = org.Hs.eg.db::org.Hs.eg.db
      )$ENTREZID,
      error = function(e) NULL
    ))
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

  tryCatch(
    {
      .go |>
        tibble::as_tibble() |>
        dplyr::mutate(
          Description = stringr::str_wrap(
            stringr::str_to_sentence(Description),
            width = 60
          ),
          adjp = -log10(p.adjust)
        ) |>
        dplyr::select(ID, Description, adjp, Count, geneID) |>
        dplyr::arrange(adjp, Count) |>
        dplyr::mutate(Description = factor(Description, levels = Description)) -> .go_for_plot

      if (!is.infinite(.topn)) {
        .go_for_plot <- tail(.go_for_plot, .topn)
      }
      if (nrow(.go_for_plot) == 0) {
        return(NULL)
      }

      p <- .go_for_plot |>
        ggplot(aes(x = Description, y = adjp)) +
        geom_col(fill = base_fill[.ont], color = NA, width = 0.7) +
        geom_text(aes(label = Count), hjust = 1, color = "white", size = 5) +
        labs(y = "-log10(adj p)", x = ont_fullname[.ont]) +
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
      fdr = ifelse(fdr > -log10(1e-300), -log10(1e-300), fdr),
      avg_log2FC = ifelse(
        abs(avg_log2FC) > 100,
        sign(avg_log2FC) * 100,
        avg_log2FC
      ),
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

  n_up <- as.integer(dplyr::coalesce(n_color["red"], 0L))
  n_down <- as.integer(dplyr::coalesce(n_color["blue"], 0L))

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
    annotate(
      "text",
      x = Inf,
      y = Inf,
      label = paste0("Up: ", n_up),
      hjust = 1.1,
      vjust = 1.5,
      color = "red",
      size = 4,
      fontface = "bold"
    ) +
    annotate(
      "text",
      x = -Inf,
      y = Inf,
      label = paste0("Down: ", n_down),
      hjust = -0.1,
      vjust = 1.5,
      color = "blue",
      size = 4,
      fontface = "bold"
    ) +
    theme_classic(base_size = 12) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 13),
      plot.subtitle = element_text(hjust = 0.5, color = "black", size = 11)
    )

  list(p = p, markers = forplot)
}


fn_empty_go_table <- function() {
  data.table(
    ID = character(),
    Description = character(),
    GeneRatio = character(),
    BgRatio = character(),
    pvalue = numeric(),
    p.adjust = numeric(),
    qvalue = numeric(),
    geneID = character(),
    Count = integer()
  )
}


fn_go_to_table <- function(.go) {
  if (is.null(.go)) {
    return(fn_empty_go_table())
  }

  .go_df <- tryCatch(
    as.data.frame(.go),
    error = function(e) NULL
  )

  if (is.null(.go_df) || nrow(.go_df) == 0) {
    return(fn_empty_go_table())
  }

  as.data.table(.go_df)
}


fn_export_deg_excel <- function(deg_dt, outpath) {
  load_pkg(writexl)

  deg_export <- deg_dt |>
    dplyr::rename(neg_log10_fdr = fdr, volcano_color = color) |>
    dplyr::mutate(
      deg_direction = dplyr::case_when(
        volcano_color == "red" ~ "up_high_af",
        volcano_color == "blue" ~ "up_low_af",
        TRUE ~ "not_significant"
      )
    ) |>
    as.data.frame()

  writexl::write_xlsx(
    x = list(DEG = deg_export),
    path = outpath
  )
}


fn_export_go_excel <- function(result, outpath) {
  load_pkg(writexl)

  go_sheets <- list(
    up_high_af_BP = fn_go_to_table(result$go_pos$BP),
    up_high_af_CC = fn_go_to_table(result$go_pos$CC),
    up_high_af_MF = fn_go_to_table(result$go_pos$MF),
    up_low_af_BP = fn_go_to_table(result$go_neg$BP),
    up_low_af_CC = fn_go_to_table(result$go_neg$CC),
    up_low_af_MF = fn_go_to_table(result$go_neg$MF)
  )

  writexl::write_xlsx(
    x = lapply(go_sheets, as.data.frame),
    path = outpath
  )
}


fn_build_variant_sc <- function(
  candidate_dt,
  cache_path,
  min_variant_cells = 3
) {
  load_pkg(Seurat, Matrix)
  conflicted::conflicts_prefer(fs::path)
  conflicted::conflicts_prefer(base::intersect)

  if (file.exists(cache_path)) {
    log_info("Loading cached merged Seurat from {cache_path}")
    return(import(cache_path))
  }

  sc_list <- list()

  for (i in seq_len(nrow(candidate_dt))) {
    gseid_i <- candidate_dt$gseid[i]
    srrid_i <- candidate_dt$srrid[i]
    disease_i <- candidate_dt$disease[i]
    cell_dt <- candidate_dt$af_cell[[i]]

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

    sc <- tryCatch(
      import(sc_path),
      error = function(e) {
        log_warn("Failed {srrid_i}: {e$message}")
        NULL
      }
    )
    if (is.null(sc)) {
      next
    }

    sc <- NormalizeData(
      sc,
      assay = "RNA",
      normalization.method = "LogNormalize",
      scale.factor = 10000
    )

    valid_bc <- intersect(cell_dt$barcode, colnames(sc))
    if (length(valid_bc) < min_variant_cells) {
      log_warn("Too few variant cells in {srrid_i}: {length(valid_bc)}")
      next
    }

    match_idx <- match(valid_bc, cell_dt$barcode)
    sc_sub <- sc[, valid_bc]
    sc_sub$af_cell <- cell_dt$af_cell[match_idx]
    sc_sub$depth_cell <- cell_dt$depth[match_idx]
    sc_sub$srrid_orig <- srrid_i
    sc_sub$disease <- disease_i

    meta_cols <- colnames(sc_sub@meta.data)
    if (!"celltype_l1" %in% meta_cols) {
      if ("celltype" %in% meta_cols) {
        sc_sub$celltype_l1 <- sc_sub@meta.data$celltype
      } else if ("predicted.celltype.l1" %in% meta_cols) {
        sc_sub$celltype_l1 <- sc_sub@meta.data$predicted.celltype.l1
      } else {
        sc_sub$celltype_l1 <- cell_dt$celltype[match_idx]
      }
    }
    if (!"celltype_l2" %in% meta_cols) {
      if ("predicted.celltype.l2" %in% meta_cols) {
        sc_sub$celltype_l2 <- sc_sub@meta.data$predicted.celltype.l2
      } else {
        sc_sub$celltype_l2 <- NA_character_
      }
    }

    sc_list[[length(sc_list) + 1]] <- sc_sub
    log_info("{srrid_i}: {length(valid_bc)} variant cells loaded")
  }

  if (length(sc_list) == 0) {
    log_warn("No variant-carrying cells loaded for merged Seurat")
    return(NULL)
  }

  merged_sc <- if (length(sc_list) == 1) {
    sc_list[[1]]
  } else {
    merge(sc_list[[1]], y = sc_list[-1], merge.data = TRUE) |>
      JoinLayers()
  }

  export(merged_sc, cache_path)
  log_info("Merged Seurat cached ({ncol(merged_sc)} cells): {cache_path}")
  merged_sc
}


fn_plot_deg_result_celltype <- function(
  result,
  outdir,
  variant_label
) {
  if (is.null(result)) {
    return(invisible(NULL))
  }

  fs::dir_create(outdir)

  ct <- result$celltype
  safe_ct <- fn_safe_name(ct)
  cutoff <- result$cutoff
  n_low <- result$n_cells["low_af"]
  n_high <- result$n_cells["high_af"]

  p_vol <- fn_de_plot(result$markers)
  p_vol$p <- p_vol$p +
    labs(
      title = glue("{variant_label} [{ct}] AF > {cutoff} vs AF < {cutoff}"),
      subtitle = glue("high_af n={n_high}, low_af n={n_low}"),
      x = glue("avg log2FC (AF > {cutoff} vs AF < {cutoff})"),
      y = "-log10(FDR)"
    )

  saveplot(
    p_vol$p,
    filename = fs::path(outdir, glue("{safe_ct}-volcano.pdf")),
    width = 10,
    height = 6,
    device = "pdf"
  )
  fn_export_deg_excel(
    deg_dt = p_vol$markers,
    outpath = fs::path(outdir, glue("{safe_ct}-volcano.xlsx"))
  )

  purrr::walk(c("pos", "neg"), function(.dir) {
    go_list <- result[[glue("go_{.dir}")]]
    dir_label <- if (.dir == "pos") {
      glue("UP in AF > {cutoff}")
    } else {
      glue("UP in AF < {cutoff}")
    }
    dir_file <- if (.dir == "pos") "up-high-af" else "up-low-af"

    purrr::walk(c("BP", "CC", "MF"), function(.ont) {
      p_go <- fn_plot_go(
        go_list[[.ont]],
        .topn = 20,
        .ont = .ont,
        .title = glue("{variant_label} [{ct}] {dir_label} - {.ont}")
      )
      if (!is.null(p_go)) {
        saveplot(
          p_go,
          filename = fs::path(
            outdir,
            glue("{safe_ct}-go-{dir_file}-{tolower(.ont)}.pdf")
          ),
          width = 10,
          height = 8,
          device = "pdf"
        )
      }
    })
  })
  fn_export_go_excel(
    result = result,
    outpath = fs::path(outdir, glue("{safe_ct}-go.xlsx"))
  )

  invisible(NULL)
}


fn_deg_by_af_cutoff_and_celltype_level <- function(
  merged_sc,
  cutoff,
  celltype_col,
  outdir,
  variant_label = "8362T>G",
  min_cells = 15,
  celltypes_keep = NULL
) {
  load_pkg(Seurat, clusterProfiler, org.Hs.eg.db)
  conflicted::conflicts_prefer(fs::path)
  conflicted::conflict_prefer("filter", "dplyr")

  if (!celltype_col %in% colnames(merged_sc@meta.data)) {
    log_warn("Column '{celltype_col}' not in Seurat metadata - skipping")
    return(data.table())
  }

  level_dir <- if (celltype_col == "celltype_l1") "by-L1" else "by-L2"
  result_dir <- fs::path(outdir, fn_cutoff_dir(cutoff), level_dir)
  fs::dir_create(result_dir)

  meta <- merged_sc@meta.data |>
    as.data.table(keep.rownames = "barcode")
  meta[, celltype_value := as.character(get(celltype_col))]
  meta[, af_group := fn_assign_af_group(af_cell, cutoff)]

  celltypes <- sort(unique(
    meta[
      !is.na(celltype_value) & celltype_value != "",
      celltype_value
    ]
  ))
  if (!is.null(celltypes_keep)) {
    celltypes <- base::intersect(celltypes, celltypes_keep)
  }

  log_info(
    "Running {level_dir} comparisons for cutoff {cutoff} across {length(celltypes)} celltypes"
  )

  summary_list <- lapply(celltypes, function(ct) {
    safe_ct <- fn_safe_name(ct)
    cache_file <- fs::path(result_dir, glue("cache-{safe_ct}.qs"))

    keep <- which(
      meta$celltype_value == ct &
        !is.na(meta$af_group)
    )
    n_low <- sum(meta$af_group[keep] == "low_af", na.rm = TRUE)
    n_high <- sum(meta$af_group[keep] == "high_af", na.rm = TRUE)

    base_summary <- data.table(
      cutoff = cutoff,
      level = level_dir,
      celltype = ct,
      n_low_af = n_low,
      n_high_af = n_high,
      n_sig = NA_integer_,
      status = "skipped"
    )

    if (n_low < min_cells || n_high < min_cells) {
      log_info(
        "  [{ct}] skip for cutoff {cutoff} (min={min_cells}): low_af={n_low}, high_af={n_high}"
      )
      base_summary$status <- "too_few_cells"
      return(base_summary)
    }

    if (file.exists(cache_file)) {
      log_info("  [{ct}] loading cache for cutoff {cutoff}")
      ct_result <- import(cache_file)
      fn_plot_deg_result_celltype(
        ct_result,
        outdir = result_dir,
        variant_label = variant_label
      )

      sig_n <- ct_result$markers |>
        tibble::rownames_to_column("gene") |>
        dplyr::filter(
          p_val_adj < 0.05,
          abs(avg_log2FC) >= 0.25,
          pct.1 >= 0.05 | pct.2 >= 0.05
        ) |>
        nrow()
      base_summary$n_sig <- sig_n
      base_summary$status <- "cached"
      return(base_summary)
    }

    sc_ct <- merged_sc[, keep]
    sc_ct$af_group <- meta$af_group[keep]
    Idents(sc_ct) <- sc_ct$af_group

    markers <- tryCatch(
      Seurat::FindMarkers(
        sc_ct,
        ident.1 = "high_af",
        ident.2 = "low_af",
        assay = "RNA",
        test.use = "wilcox",
        min.pct = 0.1,
        logfc.threshold = 0.1
      ),
      error = function(e) {
        log_warn("  [{ct}] FindMarkers failed at cutoff {cutoff}: {e$message}")
        NULL
      }
    )

    if (is.null(markers) || nrow(markers) == 0) {
      base_summary$status <- "no_markers"
      return(base_summary)
    }

    universe <- rownames(markers)
    sig <- markers |>
      tibble::rownames_to_column("gene") |>
      dplyr::filter(
        p_val_adj < 0.05,
        abs(avg_log2FC) >= 0.25,
        pct.1 >= 0.05 | pct.2 >= 0.05
      )
    pos_genes <- sig |>
      dplyr::filter(avg_log2FC > 0) |>
      dplyr::pull(gene)
    neg_genes <- sig |>
      dplyr::filter(avg_log2FC < 0) |>
      dplyr::pull(gene)

    log_info(
      "  [{ct}] cutoff {cutoff}: {nrow(sig)} DEGs (high_af={n_high}, low_af={n_low})"
    )

    ct_result <- list(
      celltype = ct,
      cutoff = cutoff,
      markers = markers,
      n_cells = c(low_af = n_low, high_af = n_high),
      go_pos = fn_enrichGO_symbols(pos_genes, universe),
      go_neg = fn_enrichGO_symbols(neg_genes, universe)
    )
    export(ct_result, cache_file)

    fn_plot_deg_result_celltype(
      ct_result,
      outdir = result_dir,
      variant_label = variant_label
    )

    base_summary$n_sig <- nrow(sig)
    base_summary$status <- "ok"
    base_summary
  })

  summary_dt <- rbindlist(summary_list, use.names = TRUE, fill = TRUE)
  export(
    summary_dt,
    fs::path(result_dir, glue("summary-{level_dir}.qs"))
  )
  fwrite(
    summary_dt,
    fs::path(result_dir, glue("summary-{level_dir}.tsv")),
    sep = "\t"
  )
  summary_dt
}


# Main ---------------------------------------------------------------------
main <- function() {
  load_pkg(GetoptLong, logger)

  VERSION <- "v0.0.1"
  GetoptLong.options(help_style = "two-column")

  nthread <- 8
  min_cells <- 15
  variant <- "8362T>G"
  cutoffs <- "0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9"
  levels <- "L1,L2"
  celltypes <- ""

  GetoptLong(
    "variant=s",
    "Variant to compare by AF cutoff",
    "cutoffs=s",
    "Comma-separated AF cutoffs; default keeps all requested comparisons",
    "levels=s",
    "Comma-separated celltype levels to run: L1, L2, or both",
    "celltypes=s",
    "Optional comma-separated whitelist of celltypes to run within the selected level(s)",
    "nthread=i",
    "Number of threads to use",
    "min_cells=i",
    "Minimum cells required per AF group within a celltype",
    "verbose",
    "Enable verbose logging"
  )

  log_layout(layout_glue_colors)
  if (isTRUE(verbose)) {
    log_threshold(TRACE)
    log_info("Verbose mode enabled")
  } else {
    log_threshold(INFO)
  }

  dotenv()
  suppressMessages({
    conflicted::conflict_prefer("filter", "dplyr")
  })

  outdir <- fs::path(Sys.getenv("OUTDIR"))
  outdirnotuse <- fs::path(Sys.getenv("OUTDIRNOTUSE"))
  af_cutoffs <- as.numeric(trimws(strsplit(cutoffs, ",")[[1]]))
  level_specs <- trimws(strsplit(levels, ",")[[1]])
  celltypes_keep <- trimws(strsplit(celltypes, ",")[[1]])
  celltypes_keep <- celltypes_keep[celltypes_keep != ""]
  if (length(celltypes_keep) == 0) {
    celltypes_keep <- NULL
  }

  if (any(is.na(af_cutoffs))) {
    stop("All cutoffs must be numeric")
  }

  the_variant <- variant

  safe_variant <- gsub(">", "_", the_variant)
  root_outdir <- fs::path(
    outdirnotuse,
    "allvariants-prioritize",
    "17.04-af-cutoff-comparison",
    safe_variant
  )
  fs::dir_create(root_outdir)
  merged_sc_cache <- fs::path(root_outdir, glue("merged-sc-{safe_variant}.qs"))

  if (fs::file_exists(merged_sc_cache)) {
    merged_sc <- import(merged_sc_cache)
  } else {
    allvariants <- import(
      fs::path(outdir, "SAMPLE-VARIANT-CLASSIFICATION-CLUSTER-BULK-AF.xlsx")
    ) |>
      dplyr::filter(variant_type %in% c("hete", "homo"))
    metafull <- import(fs::path(outdir, "SAMPLES-METADATA-FULL.xlsx"))

    candidate_dt <- allvariants |>
      dplyr::filter(variant == the_variant) |>
      dplyr::left_join(
        metafull |>
          dplyr::select(srrid, disease),
        by = "srrid"
      ) |>
      dplyr::mutate(
        variant_type = factor(variant_type, c("hete", "homo"))
      ) |>
      dplyr::arrange(variant_type, disease)

    candidate_dt$af_cell <- if (nthread <= 1) {
      mapply(
        FUN = fn_load_vaf,
        thevariant = candidate_dt$variant,
        thegseid = candidate_dt$gseid,
        thesrrid = candidate_dt$srrid,
        SIMPLIFY = FALSE,
        USE.NAMES = FALSE
      )
    } else {
      pbmcmapply(
        FUN = fn_load_vaf,
        thevariant = candidate_dt$variant,
        thegseid = candidate_dt$gseid,
        thesrrid = candidate_dt$srrid,
        mc.cores = nthread,
        SIMPLIFY = FALSE
      )
    }

    candidate_dt <- candidate_dt |>
      as.data.table()

    if (nrow(candidate_dt) == 0) {
      stop(glue("No sample records found for variant {the_variant}"))
    }

    merged_sc <- fn_build_variant_sc(
      candidate_dt = candidate_dt,
      cache_path = merged_sc_cache
    )
  }

  if (is.null(merged_sc)) {
    stop("Merged Seurat object could not be built")
  }

  level_jobs <- list()
  if ("L1" %in% level_specs) {
    level_jobs[["L1"]] <- "celltype_l1"
  }
  if ("L2" %in% level_specs) {
    level_jobs[["L2"]] <- "celltype_l2"
  }
  if (length(level_jobs) == 0) {
    stop("levels must include at least one of: L1, L2")
  }

  summary_all <- lapply(af_cutoffs, function(cutoff) {
    rbindlist(lapply(names(level_jobs), function(level_name) {
      fn_deg_by_af_cutoff_and_celltype_level(
          merged_sc = merged_sc,
          cutoff = cutoff,
          celltype_col = level_jobs[[level_name]],
          outdir = root_outdir,
          variant_label = the_variant,
          min_cells = min_cells,
          celltypes_keep = celltypes_keep
      )
    }), use.names = TRUE, fill = TRUE)
  })

  summary_all_dt <- rbindlist(summary_all, use.names = TRUE, fill = TRUE)
  is_full_run <- identical(af_cutoffs, c(0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9)) &&
    identical(sort(level_specs), c("L1", "L2")) &&
    is.null(celltypes_keep)
  celltype_tag <- if (!is.null(celltypes_keep)) {
    paste0(
      "-celltypes-",
      gsub("[^A-Za-z0-9]+", "_", paste(celltypes_keep, collapse = "_"))
    )
  } else {
    ""
  }
  run_tag <- if (is_full_run) {
    "summary-all-cutoffs"
  } else {
    paste0(
      "summary-cutoffs-",
      gsub(",", "_", cutoffs),
      "-levels-",
      gsub(",", "_", levels),
      celltype_tag
    )
  }
  export(summary_all_dt, fs::path(root_outdir, glue("{run_tag}.qs")))
  fwrite(
    summary_all_dt,
    fs::path(root_outdir, glue("{run_tag}.tsv")),
    sep = "\t"
  )

  log_info("AF cutoff comparison outputs saved to {root_outdir}")

  if (isTRUE(verbose)) {
    sessionInfo()
  }
}


if (sys.nframe() == 0) {
  main()
}
