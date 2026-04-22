#!/usr/bin/env Rscript
# Metainfo ----------------------------------------------------------------
# @AUTHOR: Chun-Jie Liu
# @CONTACT: chunjie.sam.liu.at.gmail.com
# @DATE: 2026-04-14 10:00:00
# @DESCRIPTION: Prioritize gene sets and genes from 17.09 module score
#   Wilcoxon test results. Reads intermediate Excel outputs (module-score
#   and per-celltype gene-expression tests), computes composite priority
#   scores, and generates diagnostic plots: bubble plot, ranked bar chart,
#   top-gene heatmap, and volcano-style dot plot.
# @VERSION: v0.0.1

# Library ------------------------------------------------------------------
load_pkg(jutils)

# args ---------------------------------------------------------------------
GetoptLong.options(help_style = "two-column")

levels = "L1,L2"
verbose = FALSE

GetoptLong(
  "levels=s",
  "Comma-separated celltype levels: L1, L2, or both (default L1,L2)",
  "verbose!",
  "Enable verbose logging"
)

# header -------------------------------------------------------------------
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
  conflicted::conflict_prefer("lag", "dplyr")
})

# constants ----------------------------------------------------------------
VARIANT <- "8362T>G"
SAFE_VARIANT <- gsub(">", "_", VARIANT)

OUTDIR_NOTUSE <- fs::path(Sys.getenv("OUTDIRNOTUSE"))
INPUT_ROOT <- fs::path(
  OUTDIR_NOTUSE,
  "allvariants-prioritize",
  "17.09-af-disease-interaction-module-score-gse166992",
  SAFE_VARIANT
)
OUTPUT_ROOT <- fs::path(
  OUTDIR_NOTUSE,
  "allvariants-prioritize",
  "17.20-prioritize-module-score-gene-sets-gse166992",
  SAFE_VARIANT
)
fs::dir_create(OUTPUT_ROOT)

COMPARISON_LABELS <- c(
  "COVID19_highAF vs COVID19_lowAF" = "COVID-19 AF >= 0.5 vs COVID-19 AF < 0.5",
  "Healthy_highAF vs Healthy_lowAF" = "Health AF >= 0.5 vs Health AF < 0.5",
  "COVID19_highAF vs Healthy_highAF" = "COVID-19 AF >= 0.5 vs Health AF >= 0.5",
  "COVID19_lowAF vs Healthy_lowAF" = "COVID-19 AF < 0.5 vs Health AF < 0.5"
)

GROUP_LABELS <- c(
  "Healthy_lowAF" = "Health AF < 0.5",
  "Healthy_highAF" = "Health AF >= 0.5",
  "COVID19_lowAF" = "COVID-19 AF < 0.5",
  "COVID19_highAF" = "COVID-19 AF >= 0.5"
)

# function -----------------------------------------------------------------

fn_apply_comparison_order <- function(dt, dt_name) {
  observed <- unique(dt$comparison)
  expected <- names(COMPARISON_LABELS)
  missing_expected <- setdiff(expected, observed)
  unexpected <- setdiff(observed, expected)

  if (length(unexpected) > 0 || length(missing_expected) > 0) {
    stop(glue(
      "{dt_name}: comparison keys do not match expected direction.\n",
      "Unexpected: {paste(unexpected, collapse = ', ')}\n",
      "Missing: {paste(missing_expected, collapse = ', ')}\n",
      "Re-run 17.09-af-disease-interaction-module-score-gse166992.R first."
    ))
  }

  dt[, comparison := factor(comparison, levels = expected)]
  dt
}

#' Read module-score Wilcoxon tests from 17.09 Excel output
fn_read_module_tests <- function(input_dir) {
  xlsx_path <- fs::path(input_dir, "module-score-wilcox-tests.xlsx")
  if (!fs::file_exists(xlsx_path)) {
    log_warn("Module-score tests not found: {xlsx_path}")
    return(data.table())
  }

  sheets <- readxl::excel_sheets(xlsx_path)
  dt_list <- lapply(sheets, function(sh) {
    d <- readxl::read_xlsx(xlsx_path, sheet = sh) |> as.data.table()
    d[, sheet := sh]
    d
  })
  rbindlist(dt_list, use.names = TRUE, fill = TRUE)
}


#' Read per-celltype gene-expression Wilcoxon tests from 17.09
fn_read_gene_expr_tests <- function(input_dir) {
  xlsx_files <- fs::dir_ls(input_dir, glob = "*-gene-expr-wilcox-tests.xlsx")
  if (length(xlsx_files) == 0) {
    log_warn("No gene-expr test Excel files in: {input_dir}")
    return(data.table())
  }

  dt_list <- lapply(xlsx_files, function(f) {
    ct_name <- gsub("-gene-expr-wilcox-tests\\.xlsx$", "", basename(f))
    sheets <- readxl::excel_sheets(f)
    sh_list <- lapply(sheets, function(sh) {
      d <- readxl::read_xlsx(f, sheet = sh) |> as.data.table()
      d[, `:=`(celltype = ct_name, gene_set_label = sh)]
      d
    })
    rbindlist(sh_list, use.names = TRUE, fill = TRUE)
  })
  rbindlist(dt_list, use.names = TRUE, fill = TRUE)
}


#' Compute gene-set-level priority scores
fn_gene_set_priority <- function(module_dt) {
  module_dt[, neg_log10p := -log10(pmax(p_value, 1e-300))]
  module_dt[, effect_size := abs(mean1 - mean2)]
  module_dt[, is_sig := p_adjusted < 0.05]

  priority <- module_dt[,
    .(
      n_tests = .N,
      n_sig = sum(is_sig, na.rm = TRUE),
      frac_sig = mean(is_sig, na.rm = TRUE),
      mean_neg_log10p = mean(neg_log10p, na.rm = TRUE),
      max_neg_log10p = max(neg_log10p, na.rm = TRUE),
      mean_effect_size = mean(effect_size, na.rm = TRUE),
      max_effect_size = max(effect_size, na.rm = TRUE),
      n_celltypes_sig = uniqueN(celltype[is_sig])
    ),
    by = .(gene_set_label)
  ]

  # Composite score: weighted z-score combination
  scale01 <- function(x) {
    rng <- range(x, na.rm = TRUE)
    if (rng[2] == rng[1]) {
      return(rep(0.5, length(x)))
    }
    (x - rng[1]) / (rng[2] - rng[1])
  }

  priority[,
    composite_score := 0.35 *
      scale01(frac_sig) +
      0.30 * scale01(mean_neg_log10p) +
      0.20 * scale01(n_celltypes_sig) +
      0.15 * scale01(mean_effect_size)
  ]
  priority[, rank := frank(-composite_score, ties.method = "min")]
  setorder(priority, rank)

  priority
}


#' Compute gene-level priority scores within each gene set
fn_gene_priority <- function(gene_dt) {
  gene_dt[, neg_log10p := -log10(pmax(p_value, 1e-300))]
  gene_dt[, effect_size := abs(mean1 - mean2)]
  gene_dt[, is_sig := p_adjusted < 0.05]

  priority <- gene_dt[,
    .(
      n_tests = .N,
      n_sig = sum(is_sig, na.rm = TRUE),
      frac_sig = mean(is_sig, na.rm = TRUE),
      mean_neg_log10p = mean(neg_log10p, na.rm = TRUE),
      max_neg_log10p = max(neg_log10p, na.rm = TRUE),
      mean_effect_size = mean(effect_size, na.rm = TRUE),
      max_effect_size = max(effect_size, na.rm = TRUE),
      n_celltypes_sig = uniqueN(celltype[is_sig])
    ),
    by = .(gene_set_label, gene)
  ]

  scale01 <- function(x) {
    rng <- range(x, na.rm = TRUE)
    if (rng[2] == rng[1]) {
      return(rep(0.5, length(x)))
    }
    (x - rng[1]) / (rng[2] - rng[1])
  }

  priority[,
    composite_score := 0.35 *
      scale01(frac_sig) +
      0.30 * scale01(mean_neg_log10p) +
      0.20 * scale01(n_celltypes_sig) +
      0.15 * scale01(mean_effect_size),
    by = gene_set_label
  ]
  priority[,
    rank := frank(-composite_score, ties.method = "min"),
    by = gene_set_label
  ]
  setorder(priority, gene_set_label, rank)

  priority
}


#' Plot A: Gene-set × celltype bubble plot
fn_plot_geneset_bubble <- function(module_dt, level_label) {
  plot_dt <- copy(module_dt)
  plot_dt[, neg_log10p := -log10(pmax(p_adjusted, 1e-300))]
  plot_dt[, mean_diff := mean1 - mean2]
  plot_dt[, comparison_label := COMPARISON_LABELS[as.character(comparison)]]
  plot_dt[, comparison_label := factor(
    comparison_label,
    levels = unname(COMPARISON_LABELS)
  )]
  plot_dt[, group1_label := GROUP_LABELS[group1]]
  plot_dt[, group2_label := GROUP_LABELS[group2]]
  plot_dt[, diff_label := paste0(group1_label, " - ", group2_label)]

  # Cap extreme values for visual clarity
  plot_dt[, neg_log10p_cap := pmin(neg_log10p, 50)]

  p <- ggplot(
    plot_dt,
    aes(
      x = celltype,
      y = gene_set_label,
      size = neg_log10p_cap,
      color = mean_diff
    )
  ) +
    geom_point(alpha = 0.8) +
    scale_size_continuous(
      name = expression(-log[10](p[adj])),
      range = c(1, 8),
      breaks = c(2, 5, 10, 20, 50)
    ) +
    scale_color_gradient2(
      name = "Mean score diff",
      low = "#2166AC",
      mid = "grey90",
      high = "#B2182B",
      midpoint = 0
    ) +
    facet_wrap(~comparison_label, ncol = 2) +
    labs(
      title = glue("Gene Set Significance - {level_label}"),
      subtitle = glue("Variant: {VARIANT} | GSE166992"),
      x = "Cell type",
      y = "Gene set"
    ) +
    theme_minimal(base_size = 11) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 9),
      axis.text.y = element_text(size = 9),
      strip.text = element_text(face = "bold", size = 10),
      plot.title = element_text(face = "bold", size = 13),
      plot.subtitle = element_text(size = 10, color = "grey40"),
      legend.position = "right",
      panel.grid.minor = element_blank()
    )

  p
}


#' Plot B: Per-comparison significance heatmap for gene sets
fn_plot_geneset_barrank <- function(module_dt, priority_dt, level_label) {
  # Summarize per gene_set × comparison
  plot_dt <- module_dt[,
    .(
      n_sig = sum(p_adjusted < 0.05, na.rm = TRUE),
      n_total = .N,
      mean_neg_log10p = mean(-log10(pmax(p_adjusted, 1e-300)), na.rm = TRUE),
      mean_effect = mean(mean1 - mean2, na.rm = TRUE)
    ),
    by = .(gene_set_label, comparison)
  ]

  plot_dt[, comparison_label := COMPARISON_LABELS[as.character(comparison)]]
  plot_dt[, comparison_label := factor(
    comparison_label,
    levels = unname(COMPARISON_LABELS)
  )]
  plot_dt[, sig_label := paste0(n_sig, "/", n_total)]
  plot_dt[, neg_log10p_cap := pmin(mean_neg_log10p, 50)]

  # Order gene sets by priority rank
  gs_order <- priority_dt$gene_set_label
  plot_dt[, gene_set_label := factor(gene_set_label, levels = rev(gs_order))]

  p <- ggplot(
    plot_dt,
    aes(
      x = comparison_label,
      y = gene_set_label,
      fill = neg_log10p_cap
    )
  ) +
    geom_tile(color = "white", linewidth = 0.5) +
    geom_text(
      aes(label = sig_label),
      size = 3.2,
      color = "black"
    ) +
    scale_fill_gradient2(
      name = expression(mean ~ -log[10](p)),
      low = "grey95",
      mid = "#FDDBC7",
      high = "#B2182B",
      midpoint = 5,
      limits = c(0, 50)
    ) +
    labs(
      title = glue("Gene Set Significance per Comparison - {level_label}"),
      subtitle = glue(
        "Variant: {VARIANT} | GSE166992 | Numbers = sig celltypes / total celltypes"
      ),
      x = "Comparison",
      y = "Gene set"
    ) +
    theme_minimal(base_size = 11) +
    theme(
      plot.title = element_text(face = "bold", size = 13),
      plot.subtitle = element_text(size = 9, color = "grey50"),
      axis.text.x = element_text(angle = 30, hjust = 1, size = 9),
      axis.text.y = element_text(size = 10),
      panel.grid = element_blank(),
      legend.position = "right"
    )

  p
}


#' Plot C: Separate heatmap per gene set (returns named list of plots)
fn_plot_gene_heatmaps <- function(gene_dt, top_n_sets = 6, level_label) {
  # Get top gene sets by mean significance
  set_order <- gene_dt[,
    .(mean_nlp = mean(-log10(pmax(p_adjusted, 1e-300)), na.rm = TRUE)),
    by = gene_set_label
  ][order(-mean_nlp)]
  top_sets <- head(set_order$gene_set_label, top_n_sets)

  base_dt <- gene_dt[gene_set_label %in% top_sets]
  base_dt[, neg_log10p_adj := -log10(pmax(p_adjusted, 1e-300))]
  base_dt[, neg_log10p_cap := pmin(neg_log10p_adj, 30)]
  base_dt[, comparison_label := COMPARISON_LABELS[as.character(comparison)]]
  base_dt[, comparison_label := factor(
    comparison_label,
    levels = unname(COMPARISON_LABELS)
  )]

  plots <- lapply(top_sets, function(gs) {
    plot_dt <- base_dt[gene_set_label == gs]

    p <- ggplot(
      plot_dt,
      aes(x = celltype, y = gene, fill = neg_log10p_cap)
    ) +
      geom_tile(color = "white", linewidth = 0.3) +
      scale_fill_gradient2(
        name = expression(-log[10](p[adj])),
        low = "grey95",
        mid = "#FDDBC7",
        high = "#B2182B",
        midpoint = 5,
        limits = c(0, 30)
      ) +
      facet_wrap(~comparison_label, ncol = 2) +
      labs(
        title = glue("Gene-Level Significance — {gs} — {level_label}"),
        subtitle = glue("Variant: {VARIANT} | GSE166992"),
        x = "Cell type",
        y = "Gene"
      ) +
      theme_minimal(base_size = 10) +
      theme(
        axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
        axis.text.y = element_text(size = 7),
        strip.text = element_text(face = "bold", size = 9),
        plot.title = element_text(face = "bold", size = 12),
        plot.subtitle = element_text(size = 9, color = "grey40"),
        legend.position = "right",
        panel.spacing = unit(0.3, "lines")
      )

    p
  })
  names(plots) <- top_sets
  plots
}


#' Plot D: Volcano-style dot plot per gene set (effect size vs -log10p)
fn_plot_gene_volcano <- function(gene_dt, top_n_sets = 6, level_label) {
  set_order <- gene_dt[,
    .(mean_nlp = mean(-log10(pmax(p_adjusted, 1e-300)), na.rm = TRUE)),
    by = gene_set_label
  ][order(-mean_nlp)]
  top_sets <- head(set_order$gene_set_label, top_n_sets)

  # Aggregate per gene across celltypes (mean effect, best p)
  agg_dt <- gene_dt[
    gene_set_label %in% top_sets,
    .(
      mean_effect = mean(mean1 - mean2, na.rm = TRUE),
      best_p_adj = min(p_adjusted, na.rm = TRUE),
      n_sig = sum(p_adjusted < 0.05, na.rm = TRUE),
      n_tests = .N
    ),
    by = .(gene_set_label, gene)
  ]

  agg_dt[, neg_log10p := -log10(pmax(best_p_adj, 1e-300))]
  agg_dt[, is_sig := best_p_adj < 0.05]
  agg_dt[, gene_set_label := factor(gene_set_label, levels = top_sets)]

  # Label top genes per set
  agg_dt[, top_gene := FALSE]
  agg_dt[,
    top_gene := frank(-neg_log10p, ties.method = "first") <= 5,
    by = gene_set_label
  ]

  p <- ggplot(
    agg_dt,
    aes(x = mean_effect, y = neg_log10p, color = is_sig)
  ) +
    geom_point(aes(size = n_sig), alpha = 0.7) +
    ggrepel::geom_text_repel(
      data = agg_dt[top_gene == TRUE],
      aes(label = gene),
      size = 2.8,
      max.overlaps = 15,
      segment.color = "grey60",
      show.legend = FALSE
    ) +
    geom_hline(
      yintercept = -log10(0.05),
      linetype = "dashed",
      color = "grey50"
    ) +
    scale_color_manual(
      name = "Significant\n(p_adj < 0.05)",
      values = c("TRUE" = "#E18727FF", "FALSE" = "grey70")
    ) +
    scale_size_continuous(
      name = "N sig\ncelltypes",
      range = c(1.5, 5)
    ) +
    facet_wrap(~gene_set_label, scales = "free", ncol = 3) +
    labs(
      title = glue(
        "Gene Volcano — Top {length(top_sets)} Gene Sets — {level_label}"
      ),
      subtitle = glue(
        "Variant: {VARIANT} | x = mean score difference, y = best -log10(p_adj)"
      ),
      x = "Mean Score Difference",
      y = expression(-log[10](p[adj]))
    ) +
    theme_minimal(base_size = 10) +
    theme(
      strip.text = element_text(face = "bold", size = 10),
      plot.title = element_text(face = "bold", size = 12),
      plot.subtitle = element_text(size = 9, color = "grey40"),
      legend.position = "right",
      panel.grid.minor = element_blank()
    )

  p
}


#' Process one level (L1 or L2)
fn_process_level <- function(level_name, input_root, output_root) {
  level_dir <- if (level_name == "L1") "by-L1" else "by-L2"
  level_label <- if (level_name == "L1") "L1 (major)" else "L2 (detailed)"
  input_dir <- fs::path(input_root, level_dir)
  result_dir <- fs::path(output_root, level_dir)
  fs::dir_create(result_dir)

  if (!fs::dir_exists(input_dir)) {
    log_warn("{level_name}: input directory not found: {input_dir}")
    return(invisible(NULL))
  }

  log_info("{level_name}: Reading module-score Wilcoxon tests")
  module_dt <- fn_read_module_tests(input_dir)
  if (nrow(module_dt) == 0) {
    log_warn("{level_name}: No module-score tests found — skipping")
    return(invisible(NULL))
  }
  module_dt <- fn_apply_comparison_order(module_dt, dt_name = "module_dt")

  log_info("{level_name}: Reading gene-expression Wilcoxon tests")
  gene_dt <- fn_read_gene_expr_tests(input_dir)
  if (nrow(gene_dt) > 0) {
    gene_dt <- fn_apply_comparison_order(gene_dt, dt_name = "gene_dt")
  }
  log_info(
    "{level_name}: Loaded {nrow(module_dt)} module tests, {nrow(gene_dt)} gene tests"
  )

  # Compute priorities -------------------------------------------------------
  log_info("{level_name}: Computing gene-set priority scores")
  gs_priority <- fn_gene_set_priority(module_dt)
  log_info(
    "{level_name}: Top gene sets: {paste(head(gs_priority$gene_set_label, 5), collapse = ', ')}"
  )

  gene_priority <- data.table()
  if (nrow(gene_dt) > 0) {
    log_info("{level_name}: Computing gene-level priority scores")
    gene_priority <- fn_gene_priority(gene_dt)
  }

  # Export priority tables ---------------------------------------------------
  log_info("{level_name}: Exporting priority tables")
  export(gs_priority, fs::path(result_dir, "gene-set-priority.qs2"))
  fwrite(gs_priority, fs::path(result_dir, "gene-set-priority.tsv"), sep = "\t")

  # Gene-set priority Excel
  export(gs_priority, fs::path(result_dir, "gene-set-priority.xlsx"))

  # Gene priority by set — one Excel, sheets per gene set
  if (nrow(gene_priority) > 0) {
    export(gene_priority, fs::path(result_dir, "gene-priority-by-set.qs2"))

    gene_set_names <- unique(gene_priority$gene_set_label)
    gene_sheets <- lapply(gene_set_names, function(gs) {
      gene_priority[gene_set_label == gs]
    })
    names(gene_sheets) <- substr(
      gsub("[^A-Za-z0-9_ -]", "", gene_set_names),
      1,
      31
    )
    writexl::write_xlsx(
      gene_sheets,
      fs::path(result_dir, "gene-priority-by-set.xlsx")
    )
  }

  # Diagnostic plots ---------------------------------------------------------
  log_info("{level_name}: Creating diagnostic plots")

  # Plot A — Gene Set Bubble Plot
  p_bubble <- fn_plot_geneset_bubble(module_dt, level_label)
  n_gs <- uniqueN(module_dt$gene_set_label)
  n_ct <- uniqueN(module_dt$celltype)
  saveplot(
    fs::path(result_dir, "A-geneset-bubble-plot.pdf"),
    p_bubble,
    width = max(8, n_ct * 0.8 + 4),
    height = max(6, n_gs * 0.5 + 4),
    device = "pdf"
  )

  # Plot B — Per-comparison significance heatmap
  p_barrank <- fn_plot_geneset_barrank(module_dt, gs_priority, level_label)
  saveplot(
    fs::path(result_dir, "B-geneset-priority-barrank.pdf"),
    p_barrank,
    width = 10,
    height = max(4, n_gs * 0.45 + 1),
    device = "pdf"
  )

  # Plot C — Per-gene-set heatmaps (separate pages)
  if (nrow(gene_dt) > 0) {
    top_n <- min(6, uniqueN(gene_dt$gene_set_label))
    heatmap_list <- fn_plot_gene_heatmaps(
      gene_dt,
      top_n_sets = top_n,
      level_label = level_label
    )

    # Save each gene set heatmap as a separate PDF page
    pdf_path <- fs::path(result_dir, "C-gene-heatmap-top-sets.pdf")
    n_ct <- uniqueN(gene_dt$celltype)
    pdf(pdf_path, width = max(10, n_ct * 1.2 + 4), height = 8)
    for (gs_name in names(heatmap_list)) {
      n_genes <- gene_dt[gene_set_label == gs_name, uniqueN(gene)]
      print(heatmap_list[[gs_name]])
    }
    dev.off()
    log_info(
      "{level_label}: Saved {length(heatmap_list)} gene-set heatmaps to C-gene-heatmap-top-sets.pdf"
    )

    # Plot D — Gene Volcano
    p_volcano <- fn_plot_gene_volcano(
      gene_dt,
      top_n_sets = top_n,
      level_label = level_label
    )
    saveplot(
      fs::path(result_dir, "D-gene-volcano-top-sets.pdf"),
      p_volcano,
      width = 14,
      height = max(6, ceiling(top_n / 3) * 5),
      device = "pdf"
    )
  }

  log_info("{level_name}: Done — outputs in {result_dir}")

  list(
    level = level_name,
    gs_priority = gs_priority,
    gene_priority = gene_priority
  )
}


# body ---------------------------------------------------------------------
level_specs <- trimws(strsplit(levels, ",")[[1]])
log_info("Levels to process: {paste(level_specs, collapse = ', ')}")
log_info("Input root: {INPUT_ROOT}")
log_info("Output root: {OUTPUT_ROOT}")

if (!fs::dir_exists(INPUT_ROOT)) {
  stop(glue("Input root does not exist: {INPUT_ROOT}"))
}

results <- lapply(level_specs, function(lv) {
  fn_process_level(lv, INPUT_ROOT, OUTPUT_ROOT)
})
results <- results[!vapply(results, is.null, logical(1))]

# Combined summary across levels -------------------------------------------
if (length(results) > 0) {
  all_gs <- rbindlist(
    lapply(results, function(r) {
      r$gs_priority[, level := r$level]
    }),
    use.names = TRUE,
    fill = TRUE
  )
  export(all_gs, fs::path(OUTPUT_ROOT, "gene-set-priority-all-levels.qs2"))
  fwrite(
    all_gs,
    fs::path(OUTPUT_ROOT, "gene-set-priority-all-levels.tsv"),
    sep = "\t"
  )
  export(all_gs, fs::path(OUTPUT_ROOT, "gene-set-priority-all-levels.xlsx"))

  all_genes <- rbindlist(
    lapply(results, function(r) {
      if (nrow(r$gene_priority) > 0) {
        r$gene_priority[, level := r$level]
      } else {
        data.table()
      }
    }),
    use.names = TRUE,
    fill = TRUE
  )
  if (nrow(all_genes) > 0) {
    export(all_genes, fs::path(OUTPUT_ROOT, "gene-priority-all-levels.qs2"))

    # Cross-level comparison plot: gene-set rank per level
    gs_rank_compare <- all_gs[, .(gene_set_label, level, rank, composite_score)]
    if (uniqueN(gs_rank_compare$level) > 1) {
      gs_wide <- dcast(
        gs_rank_compare,
        gene_set_label ~ level,
        value.var = "rank"
      )
      p_rank_compare <- ggplot(
        gs_wide,
        aes(x = L1, y = L2, label = gene_set_label)
      ) +
        geom_point(size = 3, color = "#E18727FF") +
        ggrepel::geom_text_repel(size = 3.2, max.overlaps = 20) +
        geom_abline(
          slope = 1,
          intercept = 0,
          linetype = "dashed",
          color = "grey50"
        ) +
        scale_x_reverse() +
        scale_y_reverse() +
        labs(
          title = "Gene Set Rank Concordance: L1 vs L2",
          subtitle = glue("Variant: {VARIANT} | 1 = highest priority"),
          x = "Rank (L1 — major celltypes)",
          y = "Rank (L2 — detailed celltypes)"
        ) +
        theme_minimal(base_size = 11) +
        theme(
          plot.title = element_text(face = "bold", size = 13),
          plot.subtitle = element_text(size = 10, color = "grey40")
        )

      saveplot(
        fs::path(OUTPUT_ROOT, "E-rank-concordance-L1-vs-L2.pdf"),
        p_rank_compare,
        width = 8,
        height = 7,
        device = "pdf"
      )
    }
  }
}

log_info("All outputs saved to {OUTPUT_ROOT}")

if (isTRUE(verbose)) {
  sessionInfo()
}
