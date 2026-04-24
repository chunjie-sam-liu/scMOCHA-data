#!/usr/bin/env Rscript
# Author: Chunjie Liu
# Contact: chunjie.sam.liu.at.gmail.com
# Date: 2026-04-24
# Description: Run DEG and GO analyses for GSE166992 disease x AF groups.
# Version: 0.1

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
    "17.07-af-cutoff-comparison-by-disease-gse166992.R"
  ),
  local = TRUE
)


# Constants ----------------------------------------------------------------
GROUP_LABELS <- c(
  "COVID19_highAF" = "COVID-19 AF >= 0.5",
  "COVID19_lowAF" = "COVID-19 AF < 0.5",
  "Healthy_highAF" = "Health AF >= 0.5",
  "Healthy_lowAF" = "Health AF < 0.5"
)

GROUP_SHEET_LABELS <- c(
  "COVID19_highAF" = "COVID19_highAF",
  "COVID19_lowAF" = "COVID19_lowAF",
  "Healthy_highAF" = "Health_highAF",
  "Healthy_lowAF" = "Health_lowAF"
)


# Helpers ------------------------------------------------------------------
fn_format_af_cutoff <- function(af_cutoff) {
  format(as.numeric(af_cutoff), nsmall = 1, trim = TRUE)
}


fn_is_equal_cutoff <- function(af_values, af_cutoff) {
  !is.na(af_values) &
    abs(af_values - af_cutoff) <= sqrt(.Machine$double.eps)
}


fn_assign_disease_af_group <- function(disease, af_cell, af_cutoff = 0.5) {
  dplyr::case_when(
    disease == "COVID-19" & af_cell >= af_cutoff ~ "COVID19_highAF",
    disease == "COVID-19" & af_cell < af_cutoff ~ "COVID19_lowAF",
    disease == "Healthy" & af_cell >= af_cutoff ~ "Healthy_highAF",
    disease == "Healthy" & af_cell < af_cutoff ~ "Healthy_lowAF",
    TRUE ~ NA_character_
  )
}


fn_group_label <- function(group_key, af_cutoff = 0.5) {
  label <- GROUP_LABELS[[group_key]]
  if (is.null(label) || is.na(label)) {
    return(group_key)
  }

  gsub("0\\.5", fn_format_af_cutoff(af_cutoff), label)
}


fn_group_sheet_label <- function(group_key) {
  label <- GROUP_SHEET_LABELS[[group_key]]
  if (is.null(label) || is.na(label)) {
    return(fn_safe_name(group_key))
  }

  label
}


fn_make_comparison_specs <- function(af_cutoff = 0.5) {
  fmt_cutoff <- fn_format_af_cutoff(af_cutoff)

  specs <- list(
    list(
      key = "COVID19_highAF_vs_COVID19_lowAF",
      dir = "covid-19-high-af-vs-covid-19-low-af",
      target_group = "COVID19_highAF",
      reference_group = "COVID19_lowAF"
    ),
    list(
      key = "Healthy_highAF_vs_Healthy_lowAF",
      dir = "health-high-af-vs-health-low-af",
      target_group = "Healthy_highAF",
      reference_group = "Healthy_lowAF"
    ),
    list(
      key = "COVID19_highAF_vs_Healthy_highAF",
      dir = "covid-19-high-af-vs-health-high-af",
      target_group = "COVID19_highAF",
      reference_group = "Healthy_highAF"
    ),
    list(
      key = "COVID19_lowAF_vs_Healthy_lowAF",
      dir = "covid-19-low-af-vs-health-low-af",
      target_group = "COVID19_lowAF",
      reference_group = "Healthy_lowAF"
    )
  )

  lapply(specs, function(spec) {
    spec$label <- glue(
      "{fn_group_label(spec$target_group, af_cutoff)} vs ",
      "{fn_group_label(spec$reference_group, af_cutoff)}"
    )
    spec$af_cutoff <- af_cutoff
    spec$cutoff_label <- glue("AF >= {fmt_cutoff} vs AF < {fmt_cutoff}")
    spec$target_label <- fn_group_label(spec$target_group, af_cutoff)
    spec$reference_label <- fn_group_label(spec$reference_group, af_cutoff)
    spec$target_sheet <- fn_group_sheet_label(spec$target_group)
    spec$reference_sheet <- fn_group_sheet_label(spec$reference_group)
    spec
  })
}


fn_validate_comparison_specs <- function(comparison_specs) {
  observed <- vapply(
    comparison_specs,
    function(spec) paste(spec$target_group, spec$reference_group, sep = " vs "),
    character(1)
  )
  expected <- c(
    "COVID19_highAF vs COVID19_lowAF",
    "Healthy_highAF vs Healthy_lowAF",
    "COVID19_highAF vs Healthy_highAF",
    "COVID19_lowAF vs Healthy_lowAF"
  )

  if (!identical(observed, expected)) {
    stop(glue(
      "Comparison order is incorrect.\n",
      "Observed: {paste(observed, collapse = '; ')}\n",
      "Expected: {paste(expected, collapse = '; ')}"
    ))
  }

  invisible(comparison_specs)
}


fn_export_deg_excel_group <- function(
  deg_dt,
  outpath,
  target_sheet,
  reference_sheet
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


fn_export_go_excel_group <- function(result, outpath) {
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


fn_count_sig_markers <- function(markers) {
  markers |>
    tibble::rownames_to_column("gene") |>
    dplyr::filter(
      p_val_adj < 0.05,
      abs(avg_log2FC) >= 0.25,
      pct.1 >= 0.05 | pct.2 >= 0.05
    ) |>
    nrow()
}


fn_plot_deg_result_group <- function(
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
  n_target <- unname(result$n_cells["target"])
  n_reference <- unname(result$n_cells["reference"])

  p_vol <- fn_de_plot(result$markers)
  p_vol$p <- p_vol$p +
    labs(
      title = glue("{variant_label} [{result$comparison_label}] [{ct}]"),
      subtitle = glue(
        "{result$target_label} n={n_target}, ",
        "{result$reference_label} n={n_reference}"
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
  fn_export_deg_excel_group(
    deg_dt = p_vol$markers,
    outpath = fs::path(outdir, glue("{safe_ct}-volcano.xlsx")),
    target_sheet = result$target_sheet,
    reference_sheet = result$reference_sheet
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
          "{variant_label} [{result$comparison_label}] [{ct}] ",
          "{dir_label} - {.ont}"
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

  fn_export_go_excel_group(
    result = result,
    outpath = fs::path(outdir, glue("{safe_ct}-go.xlsx"))
  )

  invisible(NULL)
}


fn_load_or_build_merged_sc <- function(
  root_outdir,
  outdir,
  outdirnotuse,
  the_variant,
  safe_variant,
  nthread,
  af_cutoff
) {
  merged_sc_cache_candidates <- c(
    fs::path(root_outdir, glue("merged-sc-{safe_variant}.qs")),
    fs::path(
      outdirnotuse,
      "allvariants-prioritize",
      "17.09-af-disease-interaction-module-score-gse166992",
      safe_variant,
      glue("merged-sc-{safe_variant}.qs")
    ),
    fs::path(
      outdirnotuse,
      "allvariants-prioritize",
      "17.08-af-cutoff-comparison-by-disease-gse166992-covid19-vs-health",
      safe_variant,
      glue("merged-sc-{safe_variant}.qs")
    ),
    fs::path(
      outdirnotuse,
      "allvariants-prioritize",
      "17.07-af-cutoff-comparison-by-disease-gse166992",
      safe_variant,
      glue("merged-sc-{safe_variant}.qs")
    )
  )
  merged_sc_cache_existing <- merged_sc_cache_candidates[
    fs::file_exists(merged_sc_cache_candidates)
  ]

  if (length(merged_sc_cache_existing) > 0) {
    log_info("Loading merged Seurat cache from {merged_sc_cache_existing[[1]]}")
    return(import(merged_sc_cache_existing[[1]]))
  }

  log_info("Building merged Seurat from scratch")
  allvariants <- import(
    fs::path(outdir, "SAMPLE-VARIANT-CLASSIFICATION-CLUSTER-BULK-AF.xlsx")
  ) |>
    dplyr::filter(
      variant == the_variant,
      variant_type %in% c("hete", "homo"),
      gseid == "GSE166992"
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

  if (nrow(candidate_dt) == 0) {
    stop(glue("No GSE166992 sample records found for variant {the_variant}"))
  }

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
  export(candidate_dt, fs::path(root_outdir, "candidate_dt.qs"))

  candidate_dt |>
    dplyr::select(
      gseid,
      srrid,
      Haplogroup,
      Verbose_haplogroup,
      disease,
      af_cell
    ) |>
    tidyr::unnest(af_cell) |>
    dplyr::mutate(
      disease_af_group = fn_assign_disease_af_group(
        disease = disease,
        af_cell = af_cell,
        af_cutoff = af_cutoff
      ),
      equal_cutoff = fn_is_equal_cutoff(af_cell, af_cutoff)
    ) -> candidate_dt_af

  candidate_dt_af |>
    dplyr::count(
      gseid,
      srrid,
      Haplogroup,
      Verbose_haplogroup,
      disease,
      disease_af_group,
      equal_cutoff,
      celltype
    ) |>
    tidyr::pivot_wider(
      names_from = celltype,
      values_from = n,
      values_fill = 0
    ) |>
    dplyr::arrange(disease, disease_af_group) -> candidate_dt_af_count1

  candidate_dt_af |>
    dplyr::count(disease, disease_af_group, equal_cutoff, celltype) |>
    tidyr::pivot_wider(
      names_from = celltype,
      values_from = n,
      values_fill = 0
    ) |>
    dplyr::arrange(disease, disease_af_group) -> candidate_dt_af_count2

  export(
    list(
      by_group = candidate_dt_af_count2,
      by_sample = candidate_dt_af_count1
    ),
    fs::path(root_outdir, "candidate_dt_af_counts.xlsx")
  )

  fn_build_variant_sc(
    candidate_dt = candidate_dt,
    cache_path = fs::path(root_outdir, glue("merged-sc-{safe_variant}.qs"))
  )
}


fn_deg_by_group_comparison_and_celltype_level <- function(
  merged_sc,
  comparison_spec,
  celltype_col,
  outdir,
  variant_label = "8362T>G",
  min_cells = 15,
  celltypes_keep = NULL,
  af_cutoff = 0.5
) {
  load_pkg(Seurat, clusterProfiler, org.Hs.eg.db)
  conflicted::conflicts_prefer(fs::path)
  conflicted::conflict_prefer("filter", "dplyr")

  if (!celltype_col %in% colnames(merged_sc@meta.data)) {
    log_warn("Column '{celltype_col}' not in Seurat metadata - skipping")
    return(data.table())
  }

  level_dir <- if (celltype_col == "celltype_l1") "by-L1" else "by-L2"
  result_dir <- fs::path(outdir, comparison_spec$dir, level_dir)
  fs::dir_create(result_dir)

  meta <- merged_sc@meta.data |>
    as.data.table(keep.rownames = "barcode")
  meta[, disease := as.character(disease)]
  meta[, celltype_value := as.character(get(celltype_col))]
  meta[, disease_af_group := fn_assign_disease_af_group(
    disease = disease,
    af_cell = af_cell,
    af_cutoff = af_cutoff
  )]
  meta[, equal_cutoff := fn_is_equal_cutoff(af_cell, af_cutoff)]
  meta <- meta[
    !is.na(disease) &
      disease %in% c("COVID-19", "Healthy") &
      !is.na(celltype_value) &
      celltype_value != ""
  ]

  if (nrow(meta) == 0) {
    log_warn("No cells available for {comparison_spec$label}, {level_dir}")
    return(data.table())
  }

  celltypes <- sort(unique(meta$celltype_value))
  if (!is.null(celltypes_keep)) {
    celltypes <- base::intersect(celltypes, celltypes_keep)
  }

  log_info(
    "Running {comparison_spec$label}, {level_dir} across {length(celltypes)} celltypes"
  )

  summary_list <- lapply(celltypes, function(ct) {
    safe_ct <- fn_safe_name(ct)
    cache_file <- fs::path(result_dir, glue("cache-{safe_ct}.qs"))

    ct_meta_all <- meta[celltype_value == ct]
    ct_meta <- ct_meta_all[
      disease_af_group %in% c(
        comparison_spec$target_group,
        comparison_spec$reference_group
      )
    ]
    n_target <- sum(
      ct_meta$disease_af_group == comparison_spec$target_group,
      na.rm = TRUE
    )
    n_reference <- sum(
      ct_meta$disease_af_group == comparison_spec$reference_group,
      na.rm = TRUE
    )
    n_equal_af_excluded <- sum(ct_meta_all$equal_cutoff, na.rm = TRUE)

    base_summary <- data.table(
      variant = variant_label,
      af_cutoff = af_cutoff,
      comparison_key = comparison_spec$key,
      comparison = comparison_spec$label,
      comparison_dir = comparison_spec$dir,
      target_group = comparison_spec$target_group,
      reference_group = comparison_spec$reference_group,
      target_label = comparison_spec$target_label,
      reference_label = comparison_spec$reference_label,
      level = level_dir,
      celltype = ct,
      n_target = n_target,
      n_reference = n_reference,
      n_equal_af_excluded = n_equal_af_excluded,
      n_sig = NA_integer_,
      status = "skipped"
    )

    if (n_target < min_cells || n_reference < min_cells) {
      log_info(glue(
        "  [{ct}] skip for {comparison_spec$label} (min={min_cells}): ",
        "{comparison_spec$target_group}={n_target}, ",
        "{comparison_spec$reference_group}={n_reference}"
      ))
      base_summary$status <- "too_few_cells"
      return(base_summary)
    }

    if (file.exists(cache_file)) {
      log_info("  [{ct}] loading cache for {comparison_spec$label}")
      ct_result <- import(cache_file)
      fn_plot_deg_result_group(
        ct_result,
        outdir = result_dir,
        variant_label = variant_label
      )

      base_summary$n_sig <- fn_count_sig_markers(ct_result$markers)
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
    sc_ct$disease_af_group <- ct_meta$disease_af_group[
      match(colnames(sc_ct), ct_meta$barcode)
    ]
    Idents(sc_ct) <- sc_ct$disease_af_group

    log_info(glue(
      "  [{ct}] FindMarkers ident.1={comparison_spec$target_group}, ",
      "ident.2={comparison_spec$reference_group}"
    ))
    markers <- tryCatch(
      Seurat::FindMarkers(
        sc_ct,
        ident.1 = comparison_spec$target_group,
        ident.2 = comparison_spec$reference_group,
        assay = "RNA",
        test.use = "wilcox",
        min.pct = 0.1,
        logfc.threshold = 0.1
      ),
      error = function(e) {
        log_warn(glue(
          "  [{ct}] FindMarkers failed for {comparison_spec$label}: ",
          "{e$message}"
        ))
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

    log_info(glue(
      "  [{ct}] {comparison_spec$label}: {nrow(sig)} DEGs ",
      "({comparison_spec$target_group}={n_target}, ",
      "{comparison_spec$reference_group}={n_reference})"
    ))

    ct_result <- list(
      variant = variant_label,
      af_cutoff = af_cutoff,
      comparison_key = comparison_spec$key,
      comparison_label = comparison_spec$label,
      comparison_dir = comparison_spec$dir,
      celltype = ct,
      level = level_dir,
      target_group = comparison_spec$target_group,
      reference_group = comparison_spec$reference_group,
      target_label = comparison_spec$target_label,
      reference_label = comparison_spec$reference_label,
      target_sheet = comparison_spec$target_sheet,
      reference_sheet = comparison_spec$reference_sheet,
      markers = markers,
      n_cells = c(reference = n_reference, target = n_target),
      n_equal_af_excluded = n_equal_af_excluded,
      go_pos = fn_enrichGO_symbols(pos_genes, universe),
      go_neg = fn_enrichGO_symbols(neg_genes, universe)
    )
    export(ct_result, cache_file)

    fn_plot_deg_result_group(
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

  GetoptLong.options(help_style = "two-column")

  nthread <- 8
  min_cells <- 15
  af_cutoff <- 0.5
  levels <- "L1,L2"
  celltypes <- ""

  GetoptLong(
    "levels=s",
    "Comma-separated celltype levels to run: L1, L2, or both",
    "celltypes=s",
    "Optional comma-separated whitelist of celltypes to run",
    "nthread=i",
    "Number of threads to use if a merged cache must be built",
    "min_cells=i",
    "Minimum cells required per comparison group within a celltype",
    "af_cutoff=f",
    "AF cutoff for high/low groups; high is >= cutoff and low is < cutoff",
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

  the_variant <- "8362T>G"
  safe_variant <- gsub(">", "_", the_variant)
  root_outdir <- fs::path(
    outdirnotuse,
    "allvariants-prioritize",
    "17.21-deg-go-disease-af-gse166992",
    safe_variant
  )
  fs::dir_create(root_outdir)

  comparison_specs <- fn_make_comparison_specs(af_cutoff)
  fn_validate_comparison_specs(comparison_specs)

  log_info("Using AF groups: high is AF >= {af_cutoff}; low is AF < {af_cutoff}")
  log_info(glue(
    "Requested comparisons: ",
    "{paste(vapply(comparison_specs, `[[`, character(1), 'label'), collapse = '; ')}"
  ))

  merged_sc <- fn_load_or_build_merged_sc(
    root_outdir = root_outdir,
    outdir = outdir,
    outdirnotuse = outdirnotuse,
    the_variant = the_variant,
    safe_variant = safe_variant,
    nthread = nthread,
    af_cutoff = af_cutoff
  )

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

  summary_all <- rbindlist(
    lapply(comparison_specs, function(comparison_spec) {
      rbindlist(
        lapply(names(level_jobs), function(level_name) {
          fn_deg_by_group_comparison_and_celltype_level(
            merged_sc = merged_sc,
            comparison_spec = comparison_spec,
            celltype_col = level_jobs[[level_name]],
            outdir = root_outdir,
            variant_label = the_variant,
            min_cells = min_cells,
            celltypes_keep = celltypes_keep,
            af_cutoff = af_cutoff
          )
        }),
        use.names = TRUE,
        fill = TRUE
      )
    }),
    use.names = TRUE,
    fill = TRUE
  )

  celltype_tag <- if (!is.null(celltypes_keep)) {
    paste0(
      "-celltypes-",
      gsub("[^A-Za-z0-9]+", "_", paste(celltypes_keep, collapse = "_"))
    )
  } else {
    ""
  }
  run_tag <- paste0(
    "summary-disease-af-deg-go-af-",
    gsub("\\.", "p", fn_format_af_cutoff(af_cutoff)),
    "-levels-",
    gsub(",", "_", levels),
    celltype_tag
  )

  export(summary_all, fs::path(root_outdir, glue("{run_tag}.qs")))
  fwrite(
    summary_all,
    fs::path(root_outdir, glue("{run_tag}.tsv")),
    sep = "\t"
  )

  log_info("DEG and GO outputs saved to {root_outdir}")

  if (isTRUE(verbose)) {
    sessionInfo()
  }
}


if (sys.nframe() == 0) {
  main()
}
