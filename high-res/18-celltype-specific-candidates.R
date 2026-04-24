#!/usr/bin/env Rscript
# Author: Chunjie Liu
# Contact: chunjie.sam.liu.at.gmail.com
# Date: 2026-04-22
# Description: Run single-individual AF cutoff DEG and GO analysis for celltype-specific candidate variants.
# Version: 0.1

# Reproducibility ----------------------------------------------------------
set.seed(1)

# Library ------------------------------------------------------------------
suppressMessages({
  load_pkg(jutils)
})


# Helpers ------------------------------------------------------------------
fn_format_cutoff <- function(cutoff) {
  format(cutoff, nsmall = 1, trim = TRUE)
}


fn_cutoff_dir <- function(cutoff) {
  fs::path(glue("cutoff-{fn_format_cutoff(cutoff)}"))
}


fn_cutoff_high_label <- function(cutoff) {
  glue("AF >= {fn_format_cutoff(cutoff)}")
}


fn_cutoff_low_label <- function(cutoff) {
  glue("AF < {fn_format_cutoff(cutoff)}")
}


fn_cutoff_comparison_label <- function(cutoff) {
  glue("{fn_cutoff_high_label(cutoff)} vs {fn_cutoff_low_label(cutoff)}")
}


fn_parse_cutoffs <- function(cutoffs_raw) {
  cutoff_tokens <- trimws(strsplit(cutoffs_raw, ",")[[1]])
  cutoff_tokens <- cutoff_tokens[cutoff_tokens != ""]

  if (length(cutoff_tokens) == 0) {
    stop("At least one AF cutoff must be provided")
  }

  cutoff_values <- as.numeric(cutoff_tokens)
  if (any(is.na(cutoff_values))) {
    stop("All cutoffs must be numeric")
  }

  cutoff_values[!duplicated(cutoff_values)]
}


fn_assign_af_group <- function(af_values, cutoff) {
  dplyr::case_when(
    af_values >= cutoff ~ "high_af",
    af_values < cutoff ~ "low_af",
    TRUE ~ NA_character_
  )
}


fn_safe_name <- function(x) {
  gsub("[^A-Za-z0-9]+", "_", x)
}


fn_pick_disease_label <- function(candidate_dt) {
  disease_candidates <- c(
    candidate_dt$disease[[1]],
    candidate_dt$Disease[[1]]
  )
  disease_candidates <- disease_candidates[
    !is.na(disease_candidates) &
      disease_candidates != ""
  ]

  if (length(disease_candidates) == 0) {
    return(NA_character_)
  }

  disease_candidates[[1]]
}


fn_parse_candidate_pdf <- function(filepath) {
  filename <- fs::path_file(filepath)
  parts <- strsplit(filename, "-", fixed = TRUE)[[1]]

  if (length(parts) < 4) {
    stop(glue("Could not parse variant metadata from {filename}"))
  }

  data.table(
    panel = paste(parts[seq_len(length(parts) - 3)], collapse = "-"),
    variant = parts[[length(parts) - 2]],
    gseid = parts[[length(parts) - 1]],
    srrid = sub("\\.pdf$", "", parts[[length(parts)]])
  )
}


fn_parse_variant_dir <- function(variant_dir) {
  pdfs <- fs::dir_ls(variant_dir, glob = "*.pdf", type = "file")
  if (length(pdfs) == 0) {
    stop(glue("No PDFs found in {variant_dir}"))
  }

  parsed_dt <- rbindlist(
    lapply(pdfs, fn_parse_candidate_pdf),
    use.names = TRUE,
    fill = TRUE
  )
  manifest_dt <- unique(parsed_dt[, .(variant, gseid, srrid)])

  if (nrow(manifest_dt) != 1) {
    stop(glue("Expected exactly one variant/GSE/sample tuple in {variant_dir}"))
  }

  folder_variant <- fs::path_file(variant_dir)
  if (!identical(folder_variant, manifest_dt$variant[[1]])) {
    stop(
      glue(
        "Folder name {folder_variant} does not match parsed variant {manifest_dt$variant[[1]]}"
      )
    )
  }

  data.table(
    variant_dir = as.character(fs::path_abs(variant_dir)),
    variant = manifest_dt$variant[[1]],
    gseid = manifest_dt$gseid[[1]],
    srrid = manifest_dt$srrid[[1]],
    n_panels = length(pdfs)
  )
}


fn_variant_manifest <- function(input_dir, metadata_dt) {
  variant_dirs <- fs::dir_ls(input_dir, type = "directory")
  if (length(variant_dirs) == 0) {
    stop(glue("No variant folders found in {input_dir}"))
  }

  manifest_dt <- rbindlist(
    lapply(variant_dirs, fn_parse_variant_dir),
    use.names = TRUE,
    fill = TRUE
  )
  setorder(manifest_dt, variant)

  merge(
    manifest_dt,
    metadata_dt,
    by = c("gseid", "srrid"),
    all.x = TRUE,
    sort = FALSE
  )
}


fn_load_vaf <- function(thevariant, thegseid, thesrrid) {
  conn <- db_conn(
    Sys.getenv("DUCKDB_PATH"),
    readonly = TRUE
  )
  on.exit(DBI::dbDisconnect(conn, shutdown = TRUE), add = TRUE)

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

  if (file.exists(cache_path)) {
    log_info("Loading cached Seurat from {cache_path}")
    return(import(cache_path))
  }

  gseid_i <- candidate_dt$gseid[[1]]
  srrid_i <- candidate_dt$srrid[[1]]
  disease_i <- fn_pick_disease_label(candidate_dt)
  cell_dt <- candidate_dt$af_cell[[1]]

  if (is.null(cell_dt) || nrow(cell_dt) == 0) {
    log_warn("No AF-positive cells available for {candidate_dt$variant[[1]]}")
    return(NULL)
  }

  sc_path <- fs::path(
    "/mnt/isilon/u01_project/large-scale/liuc9/raw",
    gseid_i,
    "final",
    srrid_i,
    "for_integration",
    "sc_azimuth.qs"
  )
  if (!file.exists(sc_path)) {
    stop(glue("Not found: {sc_path}"))
  }

  sc <- tryCatch(
    import(sc_path),
    error = function(e) {
      stop(glue("Failed loading {sc_path}: {e$message}"))
    }
  )

  sc <- NormalizeData(
    sc,
    assay = "RNA",
    normalization.method = "LogNormalize",
    scale.factor = 10000
  )

  valid_bc <- base::intersect(cell_dt$barcode, colnames(sc))
  if (length(valid_bc) < min_variant_cells) {
    stop(
      glue(
        "Too few variant-carrying cells for {candidate_dt$variant[[1]]}: {length(valid_bc)}"
      )
    )
  }

  sc_sub <- sc[, valid_bc]
  # match_idx <- match(valid_bc, cell_dt$barcode)
  match_idx <- match(colnames(sc_sub), cell_dt$barcode)
  sc_sub$af_cell <- cell_dt$af_cell[match_idx]
  sc_sub$depth_cell <- cell_dt$depth[match_idx]
  sc_sub$gseid_orig <- gseid_i
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

  export(sc_sub, cache_path)
  log_info("Cached {ncol(sc_sub)} variant cells to {cache_path}")
  sc_sub
}


fn_plot_deg_result_celltype <- function(
  result,
  outdir,
  variant_label,
  sample_label = NULL,
  disease_label = NULL
) {
  if (is.null(result)) {
    return(invisible(NULL))
  }

  fs::dir_create(outdir)

  ct <- result$celltype
  safe_ct <- fn_safe_name(ct)
  cutoff <- result$cutoff
  cutoff_label <- fn_cutoff_comparison_label(cutoff)
  n_low <- result$n_cells["low_af"]
  n_high <- result$n_cells["high_af"]
  subtitle_parts <- c(
    if (!is.null(sample_label) && !is.na(sample_label) && sample_label != "") {
      glue("sample={sample_label}")
    },
    if (
      !is.null(disease_label) && !is.na(disease_label) && disease_label != ""
    ) {
      glue("disease={disease_label}")
    },
    glue("high_af n={n_high}, low_af n={n_low}")
  )

  p_vol <- fn_de_plot(result$markers)
  p_vol$p <- p_vol$p +
    labs(
      title = glue("{variant_label} [{ct}] {cutoff_label}"),
      subtitle = paste(subtitle_parts, collapse = " | "),
      x = glue("avg log2FC ({cutoff_label})"),
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
      glue("UP in {fn_cutoff_high_label(cutoff)}")
    } else {
      glue("UP in {fn_cutoff_low_label(cutoff)}")
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


fn_deg_by_af_cutoff_l1 <- function(
  merged_sc,
  cutoff,
  outdir,
  variant_label,
  sample_label = NULL,
  disease_label = NULL,
  min_cells = 15,
  celltypes_keep = NULL
) {
  load_pkg(Seurat)
  conflicted::conflicts_prefer(fs::path)
  conflicted::conflict_prefer("filter", "dplyr")

  if (!"celltype_l1" %in% colnames(merged_sc@meta.data)) {
    log_warn("Column 'celltype_l1' not in Seurat metadata - skipping")
    return(data.table())
  }

  level_dir <- "by-L1"
  result_dir <- fs::path(outdir, fn_cutoff_dir(cutoff), level_dir)
  fs::dir_create(result_dir)

  meta <- merged_sc@meta.data |>
    dplyr::rename(barcode_original = barcode) |>
    as.data.table(keep.rownames = "barcode")
  meta[, celltype_value := as.character(celltype_l1)]
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

    ct_meta <- meta[
      celltype_value == ct &
        !is.na(af_group)
    ]
    n_low <- sum(ct_meta$af_group == "low_af", na.rm = TRUE)
    n_high <- sum(ct_meta$af_group == "high_af", na.rm = TRUE)

    base_summary <- data.table(
      cutoff = cutoff,
      cutoff_label = fn_cutoff_comparison_label(cutoff),
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
        variant_label = variant_label,
        sample_label = sample_label,
        disease_label = disease_label
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

    keep_idx <- match(ct_meta$barcode, colnames(merged_sc))
    keep_idx <- keep_idx[!is.na(keep_idx)]
    if (length(keep_idx) == 0) {
      base_summary$status <- "no_overlap"
      return(base_summary)
    }

    sc_ct <- merged_sc[, keep_idx]
    DefaultAssay(sc_ct) <- "RNA"
    sc_ct$af_group <- ct_meta$af_group[match(colnames(sc_ct), ct_meta$barcode)]
    Idents(sc_ct) <- sc_ct$af_group

    markers <- tryCatch(
      Seurat::FindMarkers(
        sc_ct,
        ident.1 = "high_af",
        ident.2 = "low_af",
        assay = "RNA",
        test.use = "wilcox",
        # min.pct = 0.1,
        # logfc.threshold = 0.1
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
      variant_label = variant_label,
      sample_label = sample_label,
      disease_label = disease_label
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


fn_write_variant_metadata <- function(
  candidate_dt,
  af_cell_dt,
  merged_sc,
  result_root
) {
  metadata_out <- copy(candidate_dt)
  metadata_out[, disease_label := fn_pick_disease_label(candidate_dt)]
  metadata_out[, n_variant_cells := nrow(af_cell_dt)]
  metadata_out[, n_seurat_cells := ncol(merged_sc)]
  metadata_out[, result_dir := as.character(result_root)]
  result_cols <- base::setdiff(colnames(metadata_out), "af_cell")

  fwrite(
    metadata_out[, ..result_cols],
    fs::path(result_root, "variant-metadata.tsv"),
    sep = "\t"
  )
}


fn_run_variant <- function(
  candidate_dt,
  cutoffs,
  min_cells = 15,
  celltypes_keep = NULL
) {
  variant <- candidate_dt$variant[[1]]
  gseid <- candidate_dt$gseid[[1]]
  srrid <- candidate_dt$srrid[[1]]
  disease_label <- fn_pick_disease_label(candidate_dt)
  sample_label <- glue("{gseid}/{srrid}")
  result_root <- fs::path(
    candidate_dt$variant_dir[[1]],
    "18-celltype-specific-candidates"
  )
  fs::dir_create(result_root)

  log_info("Loading AF-positive cells for {variant} ({sample_label})")
  af_cell_dt <- fn_load_vaf(
    thevariant = variant,
    thegseid = gseid,
    thesrrid = srrid
  )
  if (nrow(af_cell_dt) == 0) {
    stop(glue("No AF-positive cells found for {variant} ({sample_label})"))
  }

  candidate_dt$af_cell <- list(af_cell_dt)
  merged_sc_cache <- fs::path(result_root, "merged-sc.qs")
  merged_sc <- fn_build_variant_sc(
    candidate_dt = candidate_dt,
    cache_path = merged_sc_cache
  )
  if (is.null(merged_sc)) {
    stop(glue("Failed to build Seurat object for {variant}"))
  }

  fn_write_variant_metadata(
    candidate_dt = candidate_dt,
    af_cell_dt = af_cell_dt,
    merged_sc = merged_sc,
    result_root = result_root
  )

  summary_all <- lapply(cutoffs, function(cutoff) {
    fn_deg_by_af_cutoff_l1(
      merged_sc = merged_sc,
      cutoff = cutoff,
      outdir = result_root,
      variant_label = variant,
      sample_label = sample_label,
      disease_label = disease_label,
      min_cells = min_cells,
      celltypes_keep = celltypes_keep
    )
  })

  summary_all_dt <- rbindlist(summary_all, use.names = TRUE, fill = TRUE)
  summary_all_dt[, variant := variant]
  summary_all_dt[, gseid := gseid]
  summary_all_dt[, srrid := srrid]
  summary_all_dt[, disease := disease_label]

  default_cutoffs <- c(0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9)
  is_full_run <- identical(cutoffs, default_cutoffs) && is.null(celltypes_keep)
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
      paste(vapply(cutoffs, fn_format_cutoff, character(1)), collapse = "_"),
      celltype_tag
    )
  }

  export(summary_all_dt, fs::path(result_root, glue("{run_tag}.qs")))
  fwrite(
    summary_all_dt,
    fs::path(result_root, glue("{run_tag}.tsv")),
    sep = "\t"
  )

  log_success(
    "Finished {variant} ({sample_label}); outputs saved to {result_root}"
  )
  summary_all_dt
}


# Main ---------------------------------------------------------------------
main <- function() {
  load_pkg(GetoptLong, logger)

  VERSION <- "0.1"
  GetoptLong.options(help_style = "two-column")

  variant <- ""
  sample_id <- ""
  input_dir <- "high-res-MANUSCRIPTFIGURES-notuse/celltype-specific-candidate"
  cutoffs <- "0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9"
  celltypes <- ""
  min_cells <- 15

  GetoptLong(
    "variant=s",
    "Optional variant or comma-separated variants to run; values containing '>' are safer via sample_id",
    "sample_id=s",
    "Optional GSM/SRR sample ID or comma-separated sample IDs to run",
    "input_dir=s",
    "Directory containing per-variant folders with manuscript candidate PDFs",
    "cutoffs=s",
    "Comma-separated AF cutoffs for AF > x vs AF < x comparisons",
    "celltypes=s",
    "Optional comma-separated L1 celltype whitelist",
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
  metadata_path <- fs::path(outdir, "SAMPLES-METADATA-FULL.xlsx")
  if (!file.exists(metadata_path)) {
    stop(glue("Metadata workbook not found: {metadata_path}"))
  }

  metadata_dt <- import(metadata_path) |>
    dplyr::select(gseid, srrid, disease, Disease, status, samples) |>
    as.data.table()

  input_dir <- fs::path_abs(input_dir)
  manifest_dt <- fn_variant_manifest(
    input_dir = input_dir,
    metadata_dt = metadata_dt
  )

  variants_keep <- trimws(strsplit(variant, ",")[[1]])
  variants_keep <- variants_keep[variants_keep != ""]
  sample_ids_keep <- trimws(strsplit(sample_id, ",")[[1]])
  sample_ids_keep <- sample_ids_keep[sample_ids_keep != ""]
  if (length(variants_keep) > 0) {
    manifest_dt <- manifest_dt[variant %in% variants_keep]
  }
  if (length(sample_ids_keep) > 0) {
    manifest_dt <- manifest_dt[srrid %in% sample_ids_keep]
  }
  if (nrow(manifest_dt) == 0) {
    stop("No variants left to run after applying filters")
  }

  cutoffs_vec <- fn_parse_cutoffs(cutoffs)
  celltypes_keep <- trimws(strsplit(celltypes, ",")[[1]])
  celltypes_keep <- celltypes_keep[celltypes_keep != ""]
  if (length(celltypes_keep) == 0) {
    celltypes_keep <- NULL
  }

  failures <- character()
  for (i in seq_len(nrow(manifest_dt))) {
    candidate_dt <- copy(manifest_dt[i])
    log_info(
      "Running {candidate_dt$variant[[1]]} from {candidate_dt$gseid[[1]]}/{candidate_dt$srrid[[1]]}"
    )

    tryCatch(
      fn_run_variant(
        candidate_dt = candidate_dt,
        cutoffs = cutoffs_vec,
        min_cells = min_cells,
        celltypes_keep = celltypes_keep
      ),
      error = function(e) {
        log_error(
          "Variant {candidate_dt$variant[[1]]} failed: {e$message}"
        )
        failures <<- c(failures, candidate_dt$variant[[1]])
        NULL
      }
    )
  }

  if (length(failures) > 0) {
    stop(glue("Failed variants: {paste(unique(failures), collapse = ', ')}"))
  }

  log_success(
    "Completed {nrow(manifest_dt)} variant analyses from {input_dir}"
  )

  if (isTRUE(verbose)) {
    sessionInfo()
  }
}


if (sys.nframe() == 0) {
  main()
}
