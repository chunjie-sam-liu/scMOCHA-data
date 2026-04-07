#!/usr/bin/env Rscript
# Metainfo ----------------------------------------------------------------
# @AUTHOR: Chun-Jie Liu
# @CONTACT: chunjie.sam.liu.at.gmail.com
# @DATE: 2026-04-07 19:05:00
# @DESCRIPTION: Compare AD vs Healthy and COVID-19 vs Healthy within each celltype using all variant-carrying cells.
# @VERSION: v0.0.1

# Reproducibility ----------------------------------------------------------
set.seed(1)

# Library ------------------------------------------------------------------
suppressMessages({
  load_pkg(jutils)
})


# Shared helpers -----------------------------------------------------------
fn_script_dir <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)

  if (length(file_arg) > 0) {
    return(fs::path_dir(sub("^--file=", "", file_arg[[1]])))
  }

  fs::path_abs("high-res")
}


source(
  fs::path(
    fn_script_dir(),
    "17.05-af-cutoff-comparison-by-disease.R"
  ),
  local = TRUE
)


# Helpers ------------------------------------------------------------------
fn_disease_display_label <- function(disease) {
  dplyr::case_when(
    disease == "Alzheimer's Disease" ~ "AD",
    disease == "COVID-19" ~ "COVID-19",
    disease == "Healthy" ~ "Healthy",
    TRUE ~ as.character(disease)
  )
}


fn_disease_sheet_label <- function(disease) {
  dplyr::case_when(
    disease == "Alzheimer's Disease" ~ "AD",
    disease == "COVID-19" ~ "COVID19",
    disease == "Healthy" ~ "Healthy",
    TRUE ~ fn_safe_name(disease)
  )
}


fn_parse_comparisons <- function(comparisons_raw) {
  comparison_lookup <- list(
    ad = list(
      key = "ad",
      label = "AD vs Healthy",
      dir = "ad-vs-healthy",
      disease1 = "Healthy",
      disease2 = "Alzheimer's Disease"
    ),
    covid = list(
      key = "covid",
      label = "COVID-19 vs Healthy",
      dir = "covid-19-vs-healthy",
      disease1 = "Healthy",
      disease2 = "COVID-19"
    )
  )

  comparison_alias <- c(
    ad = "ad",
    alz = "ad",
    alzheimer = "ad",
    covid = "covid",
    covid19 = "covid",
    "covid-19" = "covid"
  )

  comparison_keys <- trimws(strsplit(comparisons_raw, ",")[[1]])
  comparison_keys <- comparison_keys[comparison_keys != ""]
  comparison_keys <- tolower(comparison_keys)
  comparison_keys <- unname(comparison_alias[comparison_keys])
  comparison_keys <- unique(comparison_keys[!is.na(comparison_keys)])

  if (length(comparison_keys) == 0) {
    stop("comparisons must include at least one of: AD, COVID")
  }

  comparison_lookup[comparison_keys]
}


fn_export_deg_excel_pairwise <- function(
  deg_dt,
  outpath,
  reference_sheet,
  target_sheet
) {
  load_pkg(writexl)

  deg_export <- deg_dt |>
    dplyr::rename(neg_log10_fdr = fdr, volcano_color = color) |>
    dplyr::mutate(
      deg_direction = dplyr::case_when(
        volcano_color == "red" ~ glue("up_{target_sheet}"),
        volcano_color == "blue" ~ glue("up_{reference_sheet}"),
        TRUE ~ "not_significant"
      )
    ) |>
    as.data.frame()

  writexl::write_xlsx(
    x = list(DEG = deg_export),
    path = outpath
  )
}


fn_export_go_excel_pairwise <- function(result, outpath) {
  load_pkg(writexl)

  go_sheets <- setNames(
    object = list(
      fn_go_to_table(result$go_pos$BP),
      fn_go_to_table(result$go_pos$CC),
      fn_go_to_table(result$go_pos$MF),
      fn_go_to_table(result$go_neg$BP),
      fn_go_to_table(result$go_neg$CC),
      fn_go_to_table(result$go_neg$MF)
    ),
    nm = c(
      glue("up_{result$target_sheet}_BP"),
      glue("up_{result$target_sheet}_CC"),
      glue("up_{result$target_sheet}_MF"),
      glue("up_{result$reference_sheet}_BP"),
      glue("up_{result$reference_sheet}_CC"),
      glue("up_{result$reference_sheet}_MF")
    )
  )

  writexl::write_xlsx(
    x = lapply(go_sheets, as.data.frame),
    path = outpath
  )
}


fn_plot_deg_result_pairwise <- function(
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
  n_reference <- unname(result$n_cells["reference"])
  n_target <- unname(result$n_cells["target"])

  p_vol <- fn_de_plot(result$markers)
  p_vol$p <- p_vol$p +
    labs(
      title = glue("{variant_label} [{result$comparison_label}] [{ct}]"),
      subtitle = glue(
        "{result$target_label} n={n_target}, {result$reference_label} n={n_reference}"
      ),
      x = glue(
        "avg log2FC ({result$target_label} vs {result$reference_label})"
      ),
      y = "-log10(FDR)"
    )

  saveplot(
    p_vol$p,
    filename = fs::path(outdir, glue("{safe_ct}-volcano.pdf")),
    width = 10,
    height = 6,
    device = "pdf"
  )
  fn_export_deg_excel_pairwise(
    deg_dt = p_vol$markers,
    outpath = fs::path(outdir, glue("{safe_ct}-volcano.xlsx")),
    reference_sheet = result$reference_sheet,
    target_sheet = result$target_sheet
  )

  purrr::walk(c("pos", "neg"), function(.dir) {
    go_list <- result[[glue("go_{.dir}")]]
    dir_label <- if (.dir == "pos") {
      glue("UP in {result$target_label}")
    } else {
      glue("UP in {result$reference_label}")
    }
    dir_file <- if (.dir == "pos") {
      glue("up-{tolower(result$target_sheet)}")
    } else {
      glue("up-{tolower(result$reference_sheet)}")
    }

    purrr::walk(c("BP", "CC", "MF"), function(.ont) {
      p_go <- fn_plot_go(
        go_list[[.ont]],
        .topn = 20,
        .ont = .ont,
        .title = glue(
          "{variant_label} [{result$comparison_label}] [{ct}] {dir_label} - {.ont}"
        )
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

  fn_export_go_excel_pairwise(
    result = result,
    outpath = fs::path(outdir, glue("{safe_ct}-go.xlsx"))
  )

  invisible(NULL)
}


fn_deg_by_disease_comparison_and_celltype_level <- function(
  merged_sc,
  disease1,
  disease2,
  comparison_label,
  comparison_dir,
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
  result_dir <- fs::path(outdir, comparison_dir, level_dir)
  fs::dir_create(result_dir)

  meta <- merged_sc@meta.data |>
    as.data.table(keep.rownames = "barcode")
  meta[, disease := as.character(disease)]
  meta[, celltype_value := as.character(get(celltype_col))]

  meta <- meta[
    !is.na(disease) &
      disease %in% c(disease1, disease2) &
      !is.na(celltype_value) &
      celltype_value != ""
  ]

  if (nrow(meta) == 0) {
    log_warn("No cells available for {comparison_label}, {level_dir}")
    return(data.table())
  }

  celltypes <- sort(unique(meta$celltype_value))
  if (!is.null(celltypes_keep)) {
    celltypes <- base::intersect(celltypes, celltypes_keep)
  }

  log_info(
    "Running {comparison_label}, {level_dir} across {length(celltypes)} celltypes"
  )

  reference_label <- fn_disease_display_label(disease1)
  target_label <- fn_disease_display_label(disease2)
  reference_sheet <- fn_disease_sheet_label(disease1)
  target_sheet <- fn_disease_sheet_label(disease2)

  summary_list <- lapply(celltypes, function(ct) {
    safe_ct <- fn_safe_name(ct)
    cache_file <- fs::path(result_dir, glue("cache-{safe_ct}.qs"))

    ct_meta <- meta[celltype_value == ct]
    n_reference <- sum(ct_meta$disease == disease1, na.rm = TRUE)
    n_target <- sum(ct_meta$disease == disease2, na.rm = TRUE)

    base_summary <- data.table(
      comparison = comparison_label,
      comparison_dir = comparison_dir,
      disease_reference = disease1,
      disease_target = disease2,
      level = level_dir,
      celltype = ct,
      n_reference = n_reference,
      n_target = n_target,
      n_sig = NA_integer_,
      status = "skipped"
    )

    if (n_reference < min_cells || n_target < min_cells) {
      log_info(
        "  [{ct}] skip for {comparison_label} (min={min_cells}): {reference_label}={n_reference}, {target_label}={n_target}"
      )
      base_summary$status <- "too_few_cells"
      return(base_summary)
    }

    if (file.exists(cache_file)) {
      log_info("  [{ct}] loading cache for {comparison_label}")
      ct_result <- import(cache_file)
      fn_plot_deg_result_pairwise(
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

    keep_idx <- match(ct_meta$barcode, colnames(merged_sc))
    keep_idx <- keep_idx[!is.na(keep_idx)]
    if (length(keep_idx) == 0) {
      base_summary$status <- "no_overlap"
      return(base_summary)
    }

    sc_ct <- merged_sc[, keep_idx]
    sc_ct$disease_compare <- ct_meta$disease[match(colnames(sc_ct), ct_meta$barcode)]
    Idents(sc_ct) <- sc_ct$disease_compare

    markers <- tryCatch(
      Seurat::FindMarkers(
        sc_ct,
        ident.1 = disease2,
        ident.2 = disease1,
        assay = "RNA",
        test.use = "wilcox",
        min.pct = 0.1,
        logfc.threshold = 0.1
      ),
      error = function(e) {
        log_warn("  [{ct}] FindMarkers failed for {comparison_label}: {e$message}")
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
      "  [{ct}] {comparison_label}: {nrow(sig)} DEGs ({target_label}={n_target}, {reference_label}={n_reference})"
    )

    ct_result <- list(
      comparison_label = comparison_label,
      comparison_dir = comparison_dir,
      celltype = ct,
      reference_label = reference_label,
      target_label = target_label,
      reference_sheet = reference_sheet,
      target_sheet = target_sheet,
      markers = markers,
      n_cells = c(reference = n_reference, target = n_target),
      go_pos = fn_enrichGO_symbols(pos_genes, universe),
      go_neg = fn_enrichGO_symbols(neg_genes, universe)
    )
    export(ct_result, cache_file)

    fn_plot_deg_result_pairwise(
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
  levels <- "L1,L2"
  celltypes <- ""
  comparisons <- "AD,COVID"

  GetoptLong(
    "levels=s",
    "Comma-separated celltype levels to run: L1, L2, or both",
    "celltypes=s",
    "Optional comma-separated whitelist of celltypes to run within the selected level(s)",
    "comparisons=s",
    "Comparisons to run: AD, COVID, or both",
    "nthread=i",
    "Number of threads to use if a merged cache must be built",
    "min_cells=i",
    "Minimum cells required per disease group within a celltype",
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
  level_specs <- trimws(strsplit(levels, ",")[[1]])
  celltypes_keep <- trimws(strsplit(celltypes, ",")[[1]])
  celltypes_keep <- celltypes_keep[celltypes_keep != ""]

  if (length(celltypes_keep) == 0) {
    celltypes_keep <- NULL
  }

  comparison_specs <- fn_parse_comparisons(comparisons)
  the_variant <- "8362T>G"
  safe_variant <- gsub(">", "_", the_variant)
  root_outdir <- fs::path(
    outdirnotuse,
    "allvariants-prioritize",
    "17.06-af-cutoff-comparison-by-covid-19-ad-vs-health",
    safe_variant
  )
  fs::dir_create(root_outdir)

  merged_sc_cache_candidates <- c(
    fs::path(root_outdir, glue("merged-sc-{safe_variant}.qs")),
    fs::path(
      outdirnotuse,
      "allvariants-prioritize",
      "17.05-af-cutoff-comparison-by-disease",
      safe_variant,
      glue("merged-sc-{safe_variant}.qs")
    ),
    fs::path(
      outdirnotuse,
      "allvariants-prioritize",
      "17.04-af-cutoff-comparison",
      safe_variant,
      glue("merged-sc-{safe_variant}.qs")
    )
  )
  merged_sc_cache_existing <- merged_sc_cache_candidates[
    fs::file_exists(merged_sc_cache_candidates)
  ]

  if (length(merged_sc_cache_existing) > 0) {
    log_info("Loading merged Seurat cache from {merged_sc_cache_existing[[1]]}")
    merged_sc <- import(merged_sc_cache_existing[[1]])
  } else {
    allvariants <- import(
      fs::path(outdir, "SAMPLE-VARIANT-CLASSIFICATION-CLUSTER-BULK-AF.xlsx")
    ) |>
      dplyr::filter(
        variant == the_variant,
        variant_type %in% c("hete", "homo")
      )
    metafull <- import(fs::path(outdir, "SAMPLES-METADATA-FULL.xlsx"))

    candidate_dt <- allvariants |>
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
      cache_path = fs::path(root_outdir, glue("merged-sc-{safe_variant}.qs"))
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

  summary_all <- lapply(comparison_specs, function(.comparison) {
    rbindlist(lapply(names(level_jobs), function(level_name) {
      fn_deg_by_disease_comparison_and_celltype_level(
        merged_sc = merged_sc,
        disease1 = .comparison$disease1,
        disease2 = .comparison$disease2,
        comparison_label = .comparison$label,
        comparison_dir = .comparison$dir,
        celltype_col = level_jobs[[level_name]],
        outdir = root_outdir,
        variant_label = the_variant,
        min_cells = min_cells,
        celltypes_keep = celltypes_keep
      )
    }), use.names = TRUE, fill = TRUE)
  })

  summary_all_dt <- rbindlist(summary_all, use.names = TRUE, fill = TRUE)

  celltype_tag <- if (!is.null(celltypes_keep)) {
    paste0(
      "-celltypes-",
      gsub("[^A-Za-z0-9]+", "_", paste(celltypes_keep, collapse = "_"))
    )
  } else {
    ""
  }
  comparison_tag <- paste(vapply(comparison_specs, `[[`, character(1), "key"), collapse = "_")
  run_tag <- paste0(
    "summary-levels-",
    gsub(",", "_", levels),
    "-comparisons-",
    comparison_tag,
    celltype_tag
  )

  export(summary_all_dt, fs::path(root_outdir, glue("{run_tag}.qs")))
  fwrite(
    summary_all_dt,
    fs::path(root_outdir, glue("{run_tag}.tsv")),
    sep = "\t"
  )

  log_info("Disease comparison outputs saved to {root_outdir}")

  if (isTRUE(verbose)) {
    sessionInfo()
  }
}


if (sys.nframe() == 0) {
  main()
}
