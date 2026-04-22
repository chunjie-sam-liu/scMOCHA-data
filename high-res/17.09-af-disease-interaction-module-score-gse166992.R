#!/usr/bin/env Rscript
# Metainfo ----------------------------------------------------------------
# @AUTHOR: Chun-Jie Liu
# @CONTACT: chunjie.sam.liu.at.gmail.com
# @DATE: 2026-04-13 14:00:00
# @DESCRIPTION: Disease x AF interaction module score analysis for GSE166992.
#   Tests whether COVID-19 + high 8362T>G heteroplasmy amplifies inflammatory
#   gene expression and reduces mitochondrial gene expression.
#   2x2 design: {Healthy, COVID-19} x {AF < 0.5, AF >= 0.5}
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


# Gene sets ----------------------------------------------------------------
INFLAMMATORY_GENES <- c(
  # Cytokines / Chemokines
  "IL6",
  "IL1B",
  "TNF",
  "CXCL8",
  "CCL2",
  "CCL3",
  "CCL4",
  "CXCL1",
  "CXCL2",
  "CXCL3",
  "CXCL10",
  # NF-kB targets
  "NFKBIA",
  "NFKB1",
  "REL",
  "RELA",
  # Interferon response
  "ISG15",
  "ISG20",
  "MX1",
  "MX2",
  "OAS1",
  "IFIT1",
  "IFIT2",
  "IFIT3",
  "IRF1",
  "STAT1",
  # Stress / acute phase
  "FOS",
  "JUN",
  "EGR1",
  "ATF3",
  "SOCS3",
  "TNFAIP3"
)

MITO_GENES <- c(
  # 13 MT-encoded protein genes
  "MT-ND1",
  "MT-ND2",
  "MT-ND3",
  "MT-ND4",
  "MT-ND4L",
  "MT-ND5",
  "MT-ND6",
  "MT-CO1",
  "MT-CO2",
  "MT-CO3",
  "MT-ATP6",
  "MT-ATP8",
  "MT-CYB",
  # MT rRNA
  "MT-RNR1",
  "MT-RNR2"
)

OXPHOS_GENES <- c(
  # Complex I (NADH dehydrogenase)
  "NDUFA1",
  "NDUFA2",
  "NDUFA4",
  "NDUFB8",
  "NDUFB10",
  "NDUFS1",
  "NDUFS2",
  "NDUFS3",
  "NDUFS7",
  "NDUFS8",
  # Complex II (Succinate dehydrogenase)
  "SDHA",
  "SDHB",
  "SDHC",
  "SDHD",
  # Complex III (Cytochrome bc1)
  "UQCRC1",
  "UQCRC2",
  "UQCRB",
  "UQCRFS1",
  # Complex IV (Cytochrome c oxidase)
  "COX4I1",
  "COX5A",
  "COX5B",
  "COX6A1",
  "COX6B1",
  "COX7A2",
  "COX7C",
  # Complex V (ATP synthase)
  "ATP5F1A",
  "ATP5F1B",
  "ATP5F1C",
  "ATP5F1D",
  "ATP5MC1",
  "ATP5MC2",
  "ATP5MC3",
  "ATP5MG"
)

APOPTOSIS_GENES <- c(
  # Pro-apoptotic
  "BAX",
  "BAK1",
  "BID",
  "BIM",
  "BAD",
  "PUMA",
  "NOXA1",
  "CASP3",
  "CASP7",
  "CASP8",
  "CASP9",
  "CYCS",
  "APAF1",
  "DIABLO",
  # Anti-apoptotic
  "BCL2",
  "BCL2L1",
  "MCL1",
  # Death receptors
  "FAS",
  "FASLG",
  "TNFRSF10A",
  "TNFRSF10B",
  # Executioner
  "DFFA",
  "DFFB"
)

IFN_TYPE1_GENES <- c(
  # Interferon-stimulated genes
  "ISG15",
  "ISG20",
  "MX1",
  "MX2",
  "OAS1",
  "OAS2",
  "OAS3",
  "IFIT1",
  "IFIT2",
  "IFIT3",
  "IFIT5",
  "IFITM1",
  "IFITM2",
  "IFITM3",
  "IFI6",
  "IFI27",
  "IFI35",
  "IFI44",
  "IFI44L",
  # Type I IFN signaling
  "IRF1",
  "IRF7",
  "IRF9",
  "STAT1",
  "STAT2",
  # Effectors
  "RSAD2",
  "BST2",
  "HERC5",
  "USP18",
  "TRIM22"
)

ANTIGEN_PRESENTATION_GENES <- c(
  # MHC class I
  "HLA-A",
  "HLA-B",
  "HLA-C",
  "HLA-E",
  "HLA-F",
  "B2M",
  "TAP1",
  "TAP2",
  "TAPBP",
  # MHC class II
  "HLA-DRA",
  "HLA-DRB1",
  "HLA-DPA1",
  "HLA-DPB1",
  "HLA-DQA1",
  "HLA-DQB1",
  # Antigen processing
  "PSMB8",
  "PSMB9",
  "PSMB10",
  "PSME1",
  "PSME2",
  "CD74",
  "CIITA"
)

MITO_TRANSLATION_GENES <- c(
  # Elongation factors
  "TUFM",
  "GFM1",
  "GFM2",
  "MTIF2",
  # Mitoribosomal proteins
  "MRPS12",
  "MRPS28",
  "MRPL11",
  "MRPL24",
  "MRPL44",
  # Mitochondrial aminoacyl-tRNA synthetases
  "MARS2",
  "KARS1",
  "YARS2",
  "AARS2",
  "EARS2",
  "DARS2",
  "TARS2"
)

MTUPR_GENES <- c(
  # Mitochondrial chaperones
  "HSPA9",
  "HSPD1",
  "HSPE1",
  "DNAJA3",
  # Mitochondrial proteases
  "CLPP",
  "LONP1",
  "AFG3L2",
  "SPG7",
  # Cytosolic stress chaperones
  "HSPB1",
  "HSP90AA1",
  # Integrated stress response
  "ATF4",
  "ATF5",
  "DDIT3"
)

GLYCOLYSIS_GENES <- c(
  "HK1",
  "HK2",
  "GPI",
  "PFKL",
  "PFKM",
  "ALDOA",
  "TPI1",
  "GAPDH",
  "PGK1",
  "PGAM1",
  "ENO1",
  "PKM",
  "LDHA",
  "LDHB",
  # Glucose transporters
  "SLC2A1",
  "SLC2A3"
)

OXIDATIVE_STRESS_GENES <- c(
  # Superoxide dismutases
  "SOD1",
  "SOD2",
  # Catalase / peroxidases
  "CAT",
  "GPX1",
  "GPX4",
  # Peroxiredoxins / thioredoxins
  "PRDX1",
  "PRDX2",
  "TXN",
  "TXNRD1",
  # NRF2 pathway
  "NFE2L2",
  "NQO1",
  "HMOX1",
  # Glutathione synthesis
  "GSR",
  "GCLC",
  "GCLM"
)

NFKB_GENES <- c(
  # NF-kB subunits
  "NFKB1",
  "NFKB2",
  "RELA",
  "RELB",
  # IkB / IKK complex
  "NFKBIA",
  "IKBKG",
  "CHUK",
  "IKBKB",
  # Target / feedback genes
  "TNFAIP3",
  "BIRC3",
  "TRAF1",
  "TRAF2"
)

# All module gene sets as a named list
MODULE_GENE_SETS <- list(
  inflammatory = list(
    genes = INFLAMMATORY_GENES,
    label = "Inflammatory",
    score_col = "inflammatory_score1"
  ),
  mito = list(
    genes = MITO_GENES,
    label = "Mitochondrial",
    score_col = "mito_score1"
  ),
  oxphos = list(
    genes = OXPHOS_GENES,
    label = "OXPHOS",
    score_col = "oxphos_score1"
  ),
  apoptosis = list(
    genes = APOPTOSIS_GENES,
    label = "Apoptosis",
    score_col = "apoptosis_score1"
  ),
  ifn_type1 = list(
    genes = IFN_TYPE1_GENES,
    label = "Type I IFN Response",
    score_col = "ifn_type1_score1"
  ),
  antigen_presentation = list(
    genes = ANTIGEN_PRESENTATION_GENES,
    label = "Antigen Presentation",
    score_col = "antigen_presentation_score1"
  ),
  mito_translation = list(
    genes = MITO_TRANSLATION_GENES,
    label = "Mito Translation",
    score_col = "mito_translation_score1"
  ),
  mtupr = list(
    genes = MTUPR_GENES,
    label = "mtUPR",
    score_col = "mtupr_score1"
  ),
  glycolysis = list(
    genes = GLYCOLYSIS_GENES,
    label = "Glycolysis",
    score_col = "glycolysis_score1"
  ),
  oxidative_stress = list(
    genes = OXIDATIVE_STRESS_GENES,
    label = "Oxidative Stress",
    score_col = "oxidative_stress_score1"
  ),
  nfkb = list(
    genes = NFKB_GENES,
    label = "NF-kB Signaling",
    score_col = "nfkb_score1"
  )
)


GROUP_ORDER <- c(
  "Healthy_lowAF",
  "Healthy_highAF",
  "COVID19_lowAF",
  "COVID19_highAF"
)

GROUP_COLORS <- c(
  "Healthy_lowAF" = "#7FBCE3",
  "Healthy_highAF" = "#0072B5FF",
  "COVID19_lowAF" = "#F5C28C",
  "COVID19_highAF" = "#E18727FF"
)

GROUP_LABELS <- c(
  "Healthy_lowAF" = "Healthy\nAF<0.5",
  "Healthy_highAF" = "Healthy\nAF\u22650.5",
  "COVID19_lowAF" = "COVID-19\nAF<0.5",
  "COVID19_highAF" = "COVID-19\nAF\u22650.5"
)

PAIRWISE_COMPARISONS <- list(
  c("COVID19_highAF", "COVID19_lowAF"),
  c("Healthy_highAF", "Healthy_lowAF"),
  c("COVID19_highAF", "Healthy_highAF"),
  c("COVID19_lowAF", "Healthy_lowAF")
)


# Helpers ------------------------------------------------------------------
fn_export_gene_sets_excel <- function(gene_sets, available_genes, outpath) {
  load_pkg(writexl)

  sheets <- lapply(names(gene_sets), function(nm) {
    gs <- gene_sets[[nm]]
    genes <- gs$genes
    data.frame(
      gene = genes,
      present_in_data = genes %in% available_genes,
      stringsAsFactors = FALSE
    )
  })
  names(sheets) <- vapply(
    gene_sets,
    function(gs) gs$label,
    character(1)
  )

  # Add a summary sheet
  summary_df <- data.frame(
    gene_set = vapply(gene_sets, function(gs) gs$label, character(1)),
    n_total = vapply(gene_sets, function(gs) length(gs$genes), integer(1)),
    n_present = vapply(
      gene_sets,
      function(gs) sum(gs$genes %in% available_genes),
      integer(1)
    ),
    n_missing = vapply(
      gene_sets,
      function(gs) sum(!gs$genes %in% available_genes),
      integer(1)
    ),
    score_column = vapply(gene_sets, function(gs) gs$score_col, character(1)),
    stringsAsFactors = FALSE
  )
  sheets <- c(list(Summary = summary_df), sheets)

  writexl::write_xlsx(x = sheets, path = outpath)
  log_info("Gene sets exported to {outpath}")
}


fn_filter_gene_set <- function(gene_set, available_genes, set_name) {
  present <- gene_set[gene_set %in% available_genes]
  missing <- gene_set[!gene_set %in% available_genes]

  log_info(
    "{set_name}: {length(present)}/{length(gene_set)} genes found"
  )
  if (length(missing) > 0) {
    log_info(
      "  Missing: {paste(missing, collapse = ', ')}"
    )
  }

  present
}


fn_assign_disease_af_group <- function(disease, af_cell, af_cutoff = 0.5) {
  dplyr::case_when(
    disease == "Healthy" & af_cell < af_cutoff ~ "Healthy_lowAF",
    disease == "Healthy" & af_cell >= af_cutoff ~ "Healthy_highAF",
    disease == "COVID-19" & af_cell < af_cutoff ~ "COVID19_lowAF",
    disease == "COVID-19" & af_cell >= af_cutoff ~ "COVID19_highAF",
    TRUE ~ NA_character_
  )
}


fn_plot_module_violin <- function(
  plot_dt,
  score_col,
  score_label,
  celltype_label,
  variant_label,
  comparisons = PAIRWISE_COMPARISONS
) {
  load_pkg(ggpubr)

  plot_dt <- plot_dt |>
    dplyr::filter(!is.na(disease_af_group)) |>
    dplyr::mutate(
      disease_af_group = factor(disease_af_group, levels = GROUP_ORDER)
    )

  group_counts <- plot_dt |>
    dplyr::count(disease_af_group) |>
    dplyr::mutate(
      label = glue("{GROUP_LABELS[as.character(disease_af_group)]}\nn={n}")
    )

  x_labels <- setNames(group_counts$label, group_counts$disease_af_group)

  p <- ggplot(
    plot_dt,
    aes(x = disease_af_group, y = .data[[score_col]], fill = disease_af_group)
  ) +
    geom_violin(trim = FALSE, alpha = 0.7, scale = "width") +
    geom_boxplot(width = 0.15, outlier.size = 0.3, alpha = 0.9) +
    ggpubr::stat_compare_means(
      comparisons = comparisons,
      method = "wilcox.test",
      label = "p.signif",
      size = 3.5,
      step.increase = 0.08,
      tip.length = 0.01
    ) +
    scale_fill_manual(values = GROUP_COLORS, guide = "none") +
    scale_x_discrete(labels = x_labels) +
    labs(
      title = glue("{score_label} — {celltype_label}"),
      subtitle = glue("Variant: {variant_label} | GSE166992"),
      x = NULL,
      y = score_label
    ) +
    theme_classic(base_size = 12) +
    theme(
      plot.title = element_text(face = "bold", size = 13),
      plot.subtitle = element_text(size = 10, color = "grey40"),
      axis.text.x = element_text(size = 9)
    )

  p
}


fn_plot_gene_expr_violin <- function(
  plot_dt,
  gene,
  celltype_label,
  variant_label,
  comparisons = PAIRWISE_COMPARISONS
) {
  load_pkg(ggpubr)

  plot_dt <- plot_dt |>
    dplyr::filter(!is.na(disease_af_group)) |>
    dplyr::mutate(
      disease_af_group = factor(disease_af_group, levels = GROUP_ORDER)
    )

  group_counts <- plot_dt |>
    dplyr::count(disease_af_group) |>
    dplyr::mutate(
      label = glue("{GROUP_LABELS[as.character(disease_af_group)]}\nn={n}")
    )

  x_labels <- setNames(group_counts$label, group_counts$disease_af_group)

  p <- ggplot(
    plot_dt,
    aes(x = disease_af_group, y = .data[[gene]], fill = disease_af_group)
  ) +
    geom_violin(trim = FALSE, alpha = 0.7, scale = "width") +
    geom_boxplot(width = 0.15, outlier.size = 0.3, alpha = 0.9) +
    ggpubr::stat_compare_means(
      comparisons = comparisons,
      method = "wilcox.test",
      label = "p.signif",
      size = 3.5,
      step.increase = 0.08,
      tip.length = 0.01
    ) +
    scale_fill_manual(values = GROUP_COLORS, guide = "none") +
    scale_x_discrete(labels = x_labels) +
    labs(
      title = glue("{gene} - {celltype_label}"),
      subtitle = glue("Variant: {variant_label} | GSE166992"),
      x = NULL,
      y = "Expression"
    ) +
    theme_classic(base_size = 12) +
    theme(
      plot.title = element_text(face = "bold", size = 13),
      plot.subtitle = element_text(size = 10, color = "grey40"),
      axis.text.x = element_text(size = 9)
    )

  p
}


fn_plot_geneset_expr_pages <- function(
  merged_sc,
  cell_barcodes,
  genes,
  gene_set_label,
  celltype_label,
  variant_label,
  outpath
) {
  # Extract expression for genes from RNA assay
  expr_mat <- Seurat::GetAssayData(
    merged_sc,
    assay = "RNA",
    layer = "data"
  )[genes, cell_barcodes, drop = FALSE]
  expr_dt <- as.data.table(
    as.matrix(t(expr_mat)),
    keep.rownames = "barcode"
  )

  # Merge with metadata for group info
  meta_dt <- merged_sc@meta.data[
    cell_barcodes,
    "disease_af_group",
    drop = FALSE
  ] |>
    as.data.table(keep.rownames = "barcode")
  expr_dt <- merge(expr_dt, meta_dt, by = "barcode")

  # One violin plot per gene
  plot_list <- lapply(genes, function(g) {
    fn_plot_gene_expr_violin(
      plot_dt = expr_dt,
      gene = g,
      celltype_label = celltype_label,
      variant_label = variant_label
    )
  })

  # Save list of plots as multi-page PDF
  saveplot(
    plot_list,
    filename = outpath,
    width = 7,
    height = 6,
    device = "pdf"
  )
}


fn_wilcox_pairwise <- function(
  dt,
  value_col,
  group_col = "disease_af_group",
  comparisons = PAIRWISE_COMPARISONS
) {
  results <- lapply(comparisons, function(pair) {
    g1 <- pair[1]
    g2 <- pair[2]
    vals1 <- dt[get(group_col) == g1][[value_col]]
    vals2 <- dt[get(group_col) == g2][[value_col]]

    base_row <- data.table(
      group1 = g1,
      group2 = g2,
      comparison = paste(g1, "vs", g2),
      n1 = length(vals1),
      n2 = length(vals2),
      mean1 = mean(vals1, na.rm = TRUE),
      mean2 = mean(vals2, na.rm = TRUE),
      median1 = median(vals1, na.rm = TRUE),
      median2 = median(vals2, na.rm = TRUE)
    )

    if (length(vals1) < 3 || length(vals2) < 3) {
      base_row[, `:=`(statistic = NA_real_, p_value = NA_real_)]
      return(base_row)
    }

    wt <- wilcox.test(vals1, vals2, exact = FALSE)
    base_row[, `:=`(
      statistic = unname(wt$statistic),
      p_value = wt$p.value
    )]
    base_row
  })
  rbindlist(results)
}


fn_export_module_tests_excel <- function(
  module_test_list,
  outpath
) {
  load_pkg(writexl)

  all_dt <- rbindlist(module_test_list, use.names = TRUE, fill = TRUE)
  if (nrow(all_dt) == 0) {
    return(invisible(NULL))
  }

  gene_set_names <- unique(all_dt$gene_set)
  sheets <- lapply(gene_set_names, function(gs) {
    gs_dt <- all_dt[gene_set == gs]
    gs_dt[, p_adjusted := p.adjust(p_value, method = "BH")]
    gs_dt[,
      significance := dplyr::case_when(
        p_adjusted < 0.001 ~ "***",
        p_adjusted < 0.01 ~ "**",
        p_adjusted < 0.05 ~ "*",
        p_adjusted < 0.1 ~ ".",
        TRUE ~ "ns"
      )
    ]
    gs_dt
  })
  labels <- vapply(
    gene_set_names,
    function(gs) {
      MODULE_GENE_SETS[[gs]]$label
    },
    character(1)
  )
  # Sanitize sheet names (max 31 chars, no special chars)
  labels <- substr(gsub("[^A-Za-z0-9_ -]", "", labels), 1, 31)
  names(sheets) <- labels

  writexl::write_xlsx(x = sheets, path = outpath)
  log_info("Module score tests exported to {outpath}")
}


fn_export_gene_expr_tests_excel <- function(
  merged_sc,
  cell_barcodes,
  filtered_gene_sets,
  celltype,
  outpath,
  comparisons = PAIRWISE_COMPARISONS
) {
  load_pkg(writexl)

  all_genes <- unique(unlist(filtered_gene_sets))
  expr_mat <- Seurat::GetAssayData(
    merged_sc,
    assay = "RNA",
    layer = "data"
  )[all_genes, cell_barcodes, drop = FALSE]
  expr_dt <- as.data.table(
    as.matrix(t(expr_mat)),
    keep.rownames = "barcode"
  )

  meta_dt <- merged_sc@meta.data[
    cell_barcodes,
    "disease_af_group",
    drop = FALSE
  ] |>
    as.data.table(keep.rownames = "barcode")
  expr_dt <- merge(expr_dt, meta_dt, by = "barcode")

  sheets <- lapply(names(filtered_gene_sets), function(gs_name) {
    gs <- MODULE_GENE_SETS[[gs_name]]
    gs_genes <- filtered_gene_sets[[gs_name]]
    gene_results <- lapply(gs_genes, function(gene) {
      test_dt <- fn_wilcox_pairwise(expr_dt, gene, comparisons = comparisons)
      test_dt[, gene := gene]
      test_dt
    })
    result <- rbindlist(gene_results, use.names = TRUE)
    result[, p_adjusted := p.adjust(p_value, method = "BH")]
    result[,
      significance := dplyr::case_when(
        p_adjusted < 0.001 ~ "***",
        p_adjusted < 0.01 ~ "**",
        p_adjusted < 0.05 ~ "*",
        p_adjusted < 0.1 ~ ".",
        TRUE ~ "ns"
      )
    ]
    # Reorder columns
    setcolorder(
      result,
      c(
        "gene",
        "comparison",
        "group1",
        "group2",
        "n1",
        "n2",
        "mean1",
        "mean2",
        "median1",
        "median2",
        "statistic",
        "p_value",
        "p_adjusted",
        "significance"
      )
    )
    result
  })

  labels <- vapply(
    names(filtered_gene_sets),
    function(gs) {
      MODULE_GENE_SETS[[gs]]$label
    },
    character(1)
  )
  labels <- substr(gsub("[^A-Za-z0-9_ -]", "", labels), 1, 31)
  names(sheets) <- labels

  writexl::write_xlsx(x = sheets, path = outpath)
  log_info("Gene expression tests for {celltype} exported to {outpath}")
}


fn_module_score_by_celltype <- function(
  merged_sc,
  celltype_col,
  outdir,
  filtered_gene_sets = list(),
  variant_label = "8362T>G",
  min_cells = 15,
  celltypes_keep = NULL,
  af_cutoff = 0.5
) {
  load_pkg(Seurat)
  conflicted::conflicts_prefer(fs::path)

  if (!celltype_col %in% colnames(merged_sc@meta.data)) {
    log_warn("Column '{celltype_col}' not in Seurat metadata - skipping")
    return(data.table())
  }

  level_dir <- if (celltype_col == "celltype_l1") "by-L1" else "by-L2"
  result_dir <- fs::path(outdir, level_dir)
  fs::dir_create(result_dir)

  meta <- merged_sc@meta.data |>
    as.data.table(keep.rownames = "cell_id")
  # Remove original barcode col if present to avoid duplication
  if ("barcode" %in% colnames(meta)) {
    meta[, barcode := NULL]
  }
  data.table::setnames(meta, "cell_id", "barcode")
  meta[, disease := as.character(disease)]
  meta[, celltype_value := as.character(get(celltype_col))]
  meta[,
    disease_af_group := fn_assign_disease_af_group(
      disease,
      af_cell,
      af_cutoff
    )
  ]

  meta <- meta[
    !is.na(disease_af_group) &
      !is.na(celltype_value) &
      celltype_value != ""
  ]

  if (nrow(meta) == 0) {
    log_warn("No cells available for {level_dir}")
    return(data.table())
  }

  celltypes <- sort(unique(meta$celltype_value))
  if (!is.null(celltypes_keep)) {
    celltypes <- base::intersect(celltypes, celltypes_keep)
  }

  log_info(
    "Running module score analysis, {level_dir}: {length(celltypes)} celltypes"
  )

  module_test_collector <- list()

  summary_list <- lapply(celltypes, function(ct) {
    safe_ct <- fn_safe_name(ct)
    ct_meta <- meta[celltype_value == ct]

    group_counts <- ct_meta[, .N, by = disease_af_group]
    n_sufficient <- sum(group_counts$N >= min_cells)

    base_summary <- data.table(
      level = level_dir,
      celltype = ct,
      n_Healthy_lowAF = sum(ct_meta$disease_af_group == "Healthy_lowAF"),
      n_Healthy_highAF = sum(ct_meta$disease_af_group == "Healthy_highAF"),
      n_COVID19_lowAF = sum(ct_meta$disease_af_group == "COVID19_lowAF"),
      n_COVID19_highAF = sum(ct_meta$disease_af_group == "COVID19_highAF"),
      n_total = nrow(ct_meta),
      status = "skipped"
    )

    if (n_sufficient < 3) {
      log_info(
        "  [{ct}] skip: only {n_sufficient}/4 groups with >= {min_cells} cells"
      )
      base_summary$status <- "too_few_groups"
      return(base_summary)
    }

    # Add module score columns from merged Seurat metadata
    all_score_cols <- vapply(
      MODULE_GENE_SETS,
      function(gs) gs$score_col,
      character(1)
    )
    for (sc in all_score_cols) {
      if (sc %in% colnames(merged_sc@meta.data)) {
        ct_meta[[sc]] <- merged_sc@meta.data[ct_meta$barcode, sc]
      }
    }

    # Compute per-group summary stats for all scores
    stat_list <- lapply(names(MODULE_GENE_SETS), function(gs_name) {
      gs <- MODULE_GENE_SETS[[gs_name]]
      sc <- gs$score_col
      if (!sc %in% colnames(ct_meta)) {
        return(NULL)
      }
      ct_meta[
        !is.na(disease_af_group),
        .(
          score_name = gs_name,
          score_label = gs$label,
          n = .N,
          mean_score = mean(get(sc), na.rm = TRUE),
          median_score = median(get(sc), na.rm = TRUE),
          sd_score = sd(get(sc), na.rm = TRUE)
        ),
        by = disease_af_group
      ]
    })
    stat_dt <- rbindlist(
      stat_list[!vapply(stat_list, is.null, logical(1))],
      use.names = TRUE
    )

    # Violin plots for each gene set
    plot_list <- list()
    for (gs_name in names(MODULE_GENE_SETS)) {
      gs <- MODULE_GENE_SETS[[gs_name]]
      if (!gs$score_col %in% colnames(ct_meta)) {
        next
      }
      plot_list[[gs_name]] <- fn_plot_module_violin(
        plot_dt = ct_meta,
        score_col = gs$score_col,
        score_label = glue("{gs$label} Score"),
        celltype_label = ct,
        variant_label = variant_label
      )
    }

    n_plots <- length(plot_list)
    ncol_p <- min(n_plots, 3)
    nrow_p <- ceiling(n_plots / ncol_p)

    p_combined <- patchwork::wrap_plots(plot_list, ncol = ncol_p) +
      patchwork::plot_annotation(
        title = glue("{variant_label} — {ct}"),
        subtitle = "Disease x AF interaction module scores",
        theme = theme(
          plot.title = element_text(face = "bold", size = 14),
          plot.subtitle = element_text(size = 11, color = "grey40")
        )
      )

    saveplot(
      p_combined,
      filename = fs::path(result_dir, glue("{safe_ct}-module-violin.pdf")),
      width = 7 * ncol_p,
      height = 6 * nrow_p,
      device = "pdf"
    )

    # Gene expression plots per gene set (multi-page PDF)
    for (gs_name in names(filtered_gene_sets)) {
      gs <- MODULE_GENE_SETS[[gs_name]]
      gs_genes <- filtered_gene_sets[[gs_name]]
      fn_plot_geneset_expr_pages(
        merged_sc = merged_sc,
        cell_barcodes = ct_meta$barcode,
        genes = gs_genes,
        gene_set_label = gs$label,
        celltype_label = ct,
        variant_label = variant_label,
        outpath = fs::path(
          result_dir,
          glue("{safe_ct}-{gs_name}-gene-expr.pdf")
        )
      )
    }

    # Wilcoxon tests for module scores — collect for level-wide Excel
    module_test_dt <- lapply(names(MODULE_GENE_SETS), function(gs_name) {
      gs <- MODULE_GENE_SETS[[gs_name]]
      sc <- gs$score_col
      if (!sc %in% colnames(ct_meta)) {
        return(NULL)
      }
      test_dt <- fn_wilcox_pairwise(ct_meta, sc)
      test_dt[, `:=`(
        celltype = ct,
        gene_set = gs_name,
        gene_set_label = gs$label
      )]
      test_dt
    })
    module_test_dt <- rbindlist(
      module_test_dt[!vapply(module_test_dt, is.null, logical(1))],
      use.names = TRUE
    )
    module_test_collector[[ct]] <<- module_test_dt

    # Gene expression Wilcoxon tests — one Excel per celltype
    fn_export_gene_expr_tests_excel(
      merged_sc = merged_sc,
      cell_barcodes = ct_meta$barcode,
      filtered_gene_sets = filtered_gene_sets,
      celltype = ct,
      outpath = fs::path(
        result_dir,
        glue("{safe_ct}-gene-expr-wilcox-tests.xlsx")
      )
    )

    # Export per-celltype stats
    export(stat_dt, fs::path(result_dir, glue("{safe_ct}-stats.qs")))
    fwrite(
      stat_dt,
      fs::path(result_dir, glue("{safe_ct}-stats.tsv")),
      sep = "\t"
    )

    log_info(
      "  [{ct}] done: {nrow(ct_meta)} cells, {nrow(stat_dt)/4} groups, {n_plots} scores"
    )

    base_summary$status <- "ok"
    # Add mean difference for each score
    for (gs_name in names(MODULE_GENE_SETS)) {
      gs_stat <- stat_dt[score_name == gs_name]
      if (nrow(gs_stat) == 0) {
        next
      }
      base_summary[[paste0(gs_name, "_mean_diff_lowAF")]] <- tryCatch(
        {
          gs_stat[disease_af_group == "COVID19_lowAF"]$mean_score -
            gs_stat[disease_af_group == "Healthy_lowAF"]$mean_score
        },
        error = function(e) NA_real_
      )
      base_summary[[paste0(gs_name, "_mean_diff_highAF")]] <- tryCatch(
        {
          gs_stat[disease_af_group == "COVID19_highAF"]$mean_score -
            gs_stat[disease_af_group == "Healthy_highAF"]$mean_score
        },
        error = function(e) NA_real_
      )
    }

    base_summary
  })

  summary_dt <- rbindlist(summary_list, use.names = TRUE, fill = TRUE)
  export(summary_dt, fs::path(result_dir, glue("summary-{level_dir}.qs")))
  fwrite(
    summary_dt,
    fs::path(result_dir, glue("summary-{level_dir}.tsv")),
    sep = "\t"
  )

  # Export module score Wilcoxon tests — one Excel per level
  fn_export_module_tests_excel(
    module_test_list = module_test_collector,
    outpath = fs::path(result_dir, "module-score-wilcox-tests.xlsx")
  )

  # Combined multi-celltype faceted figure for celltypes with status=ok
  ok_cts <- if (nrow(summary_dt) > 0 && "status" %in% names(summary_dt)) {
    summary_dt[status == "ok"]$celltype
  } else {
    character(0)
  }
  if (length(ok_cts) > 0) {
    # Build per-score, per-celltype plot lists
    score_ct_plots <- list()
    for (gs_name in names(MODULE_GENE_SETS)) {
      score_ct_plots[[gs_name]] <- list()
    }

    for (ct in ok_cts) {
      ct_meta <- meta[celltype_value == ct]
      all_score_cols <- vapply(
        MODULE_GENE_SETS,
        function(gs) gs$score_col,
        character(1)
      )
      for (sc in all_score_cols) {
        if (sc %in% colnames(merged_sc@meta.data)) {
          ct_meta[[sc]] <- merged_sc@meta.data[ct_meta$barcode, sc]
        }
      }

      for (gs_name in names(MODULE_GENE_SETS)) {
        gs <- MODULE_GENE_SETS[[gs_name]]
        if (!gs$score_col %in% colnames(ct_meta)) {
          next
        }
        score_ct_plots[[gs_name]][[ct]] <- fn_plot_module_violin(
          plot_dt = ct_meta,
          score_col = gs$score_col,
          score_label = gs$label,
          celltype_label = ct,
          variant_label = variant_label
        ) +
          theme(plot.subtitle = element_blank())
      }
    }

    n_cts <- length(ok_cts)
    ncol_grid <- min(n_cts, 3)
    nrow_grid <- ceiling(n_cts / ncol_grid)

    for (gs_name in names(score_ct_plots)) {
      plot_list <- score_ct_plots[[gs_name]]
      if (length(plot_list) == 0) {
        next
      }
      gs <- MODULE_GENE_SETS[[gs_name]]
      p_grid <- patchwork::wrap_plots(plot_list, ncol = ncol_grid) +
        patchwork::plot_annotation(
          title = glue(
            "{gs$label} Score — {variant_label} ({level_dir})"
          ),
          theme = theme(
            plot.title = element_text(face = "bold", size = 14)
          )
        )
      saveplot(
        p_grid,
        filename = fs::path(
          result_dir,
          glue("combined-{gs_name}-violin.pdf")
        ),
        width = 6 * ncol_grid,
        height = 5 * nrow_grid,
        device = "pdf"
      )
    }

    log_info(
      "Combined figures saved for {level_dir}: {n_cts} celltypes, {length(MODULE_GENE_SETS)} scores"
    )
  }

  summary_dt
}


# Main ---------------------------------------------------------------------
main <- function() {
  load_pkg(GetoptLong, logger, Seurat)

  VERSION <- "v0.0.1"
  GetoptLong.options(help_style = "two-column")

  nthread = 8
  min_cells = 15
  af_cutoff = 0.5
  levels = "L1,L2"
  celltypes = ""

  GetoptLong(
    "levels=s",
    "Comma-separated celltype levels to run: L1, L2, or both",
    "celltypes=s",
    "Optional comma-separated whitelist of celltypes",
    "nthread=i",
    "Number of threads for building merged Seurat if needed",
    "min_cells=i",
    "Minimum cells per group within a celltype (default 15)",
    "af_cutoff=f",
    "AF cutoff for high/low groups (default 0.5)",
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
    "17.09-af-disease-interaction-module-score-gse166992",
    safe_variant
  )
  fs::dir_create(root_outdir)

  # Reuse merged Seurat from 17.07 or 17.08
  merged_sc_cache_candidates <- c(
    fs::path(root_outdir, glue("merged-sc-{safe_variant}.qs")),
    fs::path(
      outdirnotuse,
      "allvariants-prioritize",
      "17.07-af-cutoff-comparison-by-disease-gse166992",
      safe_variant,
      glue("merged-sc-{safe_variant}.qs")
    ),
    fs::path(
      outdirnotuse,
      "allvariants-prioritize",
      "17.08-af-cutoff-comparison-by-disease-gse166992-covid19-vs-health",
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
        metafull |> dplyr::select(srrid, disease),
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

    candidate_dt <- candidate_dt |> as.data.table()
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
      unnest(af_cell) |>
      mutate(
        af_group = if_else(
          af_cell < af_cutoff,
          "AF<0.5",
          "AF>=0.5"
        )
      ) -> candidate_dt_af
    candidate_dt_af |>
      count(
        gseid,
        srrid,
        Haplogroup,
        Verbose_haplogroup,
        disease,
        af_group,
        celltype
      ) |>
      pivot_wider(
        names_from = celltype,
        values_from = n,
        values_fill = 0
      ) |>
      arrange(disease, af_group) -> candidate_dt_af_count1

    candidate_dt_af |>
      count(disease, af_group, celltype) |>
      pivot_wider(
        names_from = celltype,
        values_from = n,
        values_fill = 0
      ) |>
      arrange(disease, af_group) -> candidate_dt_af_count2

    export(
      list(
        by_group = candidate_dt_af_count2,
        by_sample = candidate_dt_af_count1
      ),
      fs::path(root_outdir, "candidate_dt_af_counts.xlsx")
    )

    if (nrow(candidate_dt) == 0) {
      stop(glue("No GSE166992 sample records found for variant {the_variant}"))
    }

    merged_sc <- fn_build_variant_sc(
      candidate_dt = candidate_dt,
      cache_path = fs::path(root_outdir, glue("merged-sc-{safe_variant}.qs"))
    )
  }

  if (is.null(merged_sc)) {
    stop("Merged Seurat object could not be built")
  }

  # Filter gene sets to available genes and export to Excel
  available_genes <- rownames(merged_sc)

  filtered_gene_sets <- list()
  for (gs_name in names(MODULE_GENE_SETS)) {
    gs <- MODULE_GENE_SETS[[gs_name]]
    filtered <- fn_filter_gene_set(
      gs$genes,
      available_genes,
      gs$label
    )
    if (length(filtered) < 3) {
      log_warn("{gs$label}: only {length(filtered)} genes found, skipping")
      next
    }
    filtered_gene_sets[[gs_name]] <- filtered
  }

  if (length(filtered_gene_sets) == 0) {
    stop("No gene sets have enough genes in the Seurat object")
  }

  # Export gene sets to Excel
  fn_export_gene_sets_excel(
    gene_sets = MODULE_GENE_SETS,
    available_genes = available_genes,
    outpath = fs::path(root_outdir, "module-gene-sets.xlsx")
  )

  # Compute module scores
  for (gs_name in names(filtered_gene_sets)) {
    gs <- MODULE_GENE_SETS[[gs_name]]
    gs_genes <- filtered_gene_sets[[gs_name]]
    score_name <- sub("1$", "", gs$score_col)
    log_info(
      "Computing AddModuleScore for {gs$label} ({length(gs_genes)} genes)"
    )
    merged_sc <- Seurat::AddModuleScore(
      merged_sc,
      features = list(gs_genes),
      name = score_name,
      ctrl = min(100, length(gs_genes) * 5)
    )
  }

  # Assign disease x AF groups
  merged_sc$disease_af_group <- fn_assign_disease_af_group(
    disease = merged_sc$disease,
    af_cell = merged_sc$af_cell,
    af_cutoff = af_cutoff
  )

  group_table <- table(merged_sc$disease_af_group)
  log_info(
    "Group sizes: {paste(names(group_table), group_table, sep = '=', collapse = ', ')}"
  )

  # Level jobs
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
    lapply(names(level_jobs), function(level_name) {
      fn_module_score_by_celltype(
        merged_sc = merged_sc,
        celltype_col = level_jobs[[level_name]],
        outdir = root_outdir,
        filtered_gene_sets = filtered_gene_sets,
        variant_label = the_variant,
        min_cells = min_cells,
        celltypes_keep = celltypes_keep,
        af_cutoff = af_cutoff
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
    "summary-levels-",
    gsub(",", "_", levels),
    celltype_tag
  )

  export(summary_all, fs::path(root_outdir, glue("{run_tag}.qs")))
  fwrite(
    summary_all,
    fs::path(root_outdir, glue("{run_tag}.tsv")),
    sep = "\t"
  )

  log_info("Module score analysis outputs saved to {root_outdir}")

  if (isTRUE(verbose)) {
    sessionInfo()
  }
}


if (sys.nframe() == 0) {
  main()
}
