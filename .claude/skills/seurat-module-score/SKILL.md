---
name: seurat-module-score
description: "Workflow for Seurat-based module score and gene expression analysis with 2×2 factorial design (Disease × AF group). Use when asked to: create module score analysis scripts, add/modify gene sets for AddModuleScore, plot module score or gene expression violins across groups, export Wilcoxon test results to Excel, build Disease×AF interaction analyses, or work with scMOCHA variant-level Seurat objects. Triggers: 'module score', 'gene set', 'AddModuleScore', 'Disease AF interaction', 'violin by group', 'gene expression comparison', 'Wilcoxon test Excel'."
---

# Seurat Module Score & Gene Expression Analysis

Workflow for analyzing gene set module scores and individual gene expression across a 2×2 factorial design (Disease × AF allele frequency groups) using Seurat objects.

## Architecture Overview

```
main()
├── Load merged Seurat (cache from 17.07/17.08)
├── Filter gene sets → available genes
├── Export gene set summary → Excel
├── AddModuleScore for each gene set
├── Assign Disease×AF groups
└── Per-level (L1/L2):
    └── fn_module_score_by_celltype()
        └── Per-celltype:
            ├── Module score violin plots (combined PDF)
            ├── Gene expression violin plots (multi-page PDF per gene set)
            ├── Module Wilcoxon tests → collect for level Excel
            ├── Gene expression Wilcoxon tests → per-celltype Excel
            └── Per-group summary stats → TSV/QS
```

## Key Patterns

### Gene Set Definition

Define gene sets as named character vectors, register in `MODULE_GENE_SETS` list:

```r
MY_GENES <- c("GENE1", "GENE2", "GENE3")

MODULE_GENE_SETS <- list(
  my_set = list(
    genes = MY_GENES,
    label = "My Gene Set",          # Human-readable label for plots/sheets
    score_col = "my_set_score1"     # AddModuleScore appends "1" to name
  )
)
```

- `score_col` must end with `1` — Seurat's `AddModuleScore` appends index suffix
- Filter genes against `rownames(merged_sc)` before scoring; skip sets with <3 genes
- The `label` is used for plot titles, Excel sheet names, and display

### 2×2 Group Assignment

```r
fn_assign_disease_af_group <- function(disease, af_cell, af_cutoff = 0.5) {
  dplyr::case_when(
    disease == "Healthy" & af_cell < af_cutoff ~ "Healthy_lowAF",
    disease == "Healthy" & af_cell >= af_cutoff ~ "Healthy_highAF",
    disease == "COVID-19" & af_cell < af_cutoff ~ "COVID19_lowAF",
    disease == "COVID-19" & af_cell >= af_cutoff ~ "COVID19_highAF",
    TRUE ~ NA_character_
  )
}
```

Standard pairwise comparisons (4 pairs):
1. Healthy_lowAF vs COVID19_lowAF (disease effect at low AF)
2. Healthy_highAF vs COVID19_highAF (disease effect at high AF)
3. Healthy_lowAF vs Healthy_highAF (AF effect in healthy)
4. COVID19_lowAF vs COVID19_highAF (AF effect in disease)

### Module Score Violin Plot

Use `ggplot2` + `ggpubr::stat_compare_means()` with Wilcoxon test:

```r
fn_plot_module_violin <- function(plot_dt, score_col, score_label,
                                   celltype_label, variant_label,
                                   comparisons = PAIRWISE_COMPARISONS) {
  ggplot(plot_dt, aes(x = disease_af_group, y = .data[[score_col]],
                      fill = disease_af_group)) +
    geom_violin(trim = FALSE, alpha = 0.7, scale = "width") +
    geom_boxplot(width = 0.15, outlier.size = 0.3, alpha = 0.9) +
    ggpubr::stat_compare_means(
      comparisons = comparisons, method = "wilcox.test",
      label = "p.signif", size = 3.5,
      step.increase = 0.08, tip.length = 0.01
    ) +
    scale_fill_manual(values = GROUP_COLORS, guide = "none")
}
```

### Gene Expression Multi-Page PDF

Extract expression from RNA assay, create one violin per gene, save as multi-page PDF:

```r
fn_plot_geneset_expr_pages <- function(merged_sc, cell_barcodes, genes,
                                       gene_set_label, celltype_label,
                                       variant_label, outpath) {
  expr_mat <- Seurat::GetAssayData(merged_sc, assay = "RNA", layer = "data")[
    genes, cell_barcodes, drop = FALSE
  ]
  expr_dt <- as.data.table(as.matrix(t(expr_mat)), keep.rownames = "barcode")

  meta_dt <- merged_sc@meta.data[cell_barcodes, "disease_af_group", drop = FALSE] |>
    as.data.table(keep.rownames = "barcode")
  expr_dt <- merge(expr_dt, meta_dt, by = "barcode")

  plot_list <- lapply(genes, function(g) {
    fn_plot_gene_expr_violin(expr_dt, gene = g, ...)
  })

  saveplot(plot_list, filename = outpath, width = 7, height = 6, device = "pdf")
}
```

- `saveplot()` from `jutils` accepts a list of plots → multi-page PDF
- Use `GetAssayData(merged_sc, assay = "RNA", layer = "data")` for normalized expression

### Wilcoxon Test Helper

Generic pairwise Wilcoxon test function returning a `data.table`:

```r
fn_wilcox_pairwise <- function(dt, value_col, group_col = "disease_af_group",
                                comparisons = PAIRWISE_COMPARISONS) {
  results <- lapply(comparisons, function(pair) {
    vals1 <- dt[get(group_col) == pair[1]][[value_col]]
    vals2 <- dt[get(group_col) == pair[2]][[value_col]]
    # ... wilcox.test(vals1, vals2, exact = FALSE)
    # Returns: group1, group2, comparison, n1, n2, mean1, mean2,
    #          median1, median2, statistic, p_value
  })
  rbindlist(results)
}
```

### Excel Export — Module Score Tests

One Excel per level (by-L1, by-L2), sheets by gene set. Rows = celltype × comparison.

```r
fn_export_module_tests_excel <- function(module_test_list, outpath) {
  all_dt <- rbindlist(module_test_list, use.names = TRUE, fill = TRUE)
  # Split by gene_set, add BH-adjusted p-values per sheet
  # Sheet names from MODULE_GENE_SETS[[gs]]$label (sanitized, max 31 chars)
  writexl::write_xlsx(x = sheets, path = outpath)
}
```

Collect tests inside the per-celltype loop using `<<-` into a collector list:

```r
module_test_collector <- list()
# ... inside lapply(celltypes, function(ct) { ... })
module_test_collector[[ct]] <<- module_test_dt
```

### Excel Export — Gene Expression Tests

One Excel per celltype, sheets by gene set. Rows = gene × comparison.

```r
fn_export_gene_expr_tests_excel <- function(merged_sc, cell_barcodes,
                                             filtered_gene_sets, celltype, outpath) {
  # Extract all genes at once, run fn_wilcox_pairwise per gene
  # Add p_adjusted (BH) and significance columns per sheet
  # Column order: gene, comparison, group1, group2, n1, n2,
  #               mean1, mean2, median1, median2, statistic,
  #               p_value, p_adjusted, significance
}
```

## Seurat Gotchas

- **Duplicate barcode column**: When converting `meta.data` to `data.table`, use `as.data.table(keep.rownames = "cell_id")` then delete old `barcode` col, rename `cell_id` → `barcode`
- **Celltype names**: L1 uses short names ("Mono", not "Monocyte"). Always check `unique(meta$celltype_l1)` first
- **AddModuleScore naming**: Appends `1` to the `name` argument. Set `name = sub("1$", "", score_col)` so the final column is `score_col`
- **GetAssayData layer**: Seurat v5 uses `layer = "data"` for normalized, `layer = "counts"` for raw
- **saveplot device**: Always pass `device = "pdf"` explicitly — see [saveplot convention](references/saveplot-convention.md)

## Output File Naming

Per celltype in `by-{level}/`:
- `{safe_ct}-module-violin.pdf` — Combined module score violins
- `{safe_ct}-{gs_name}-gene-expr.pdf` — Multi-page gene expression violins
- `{safe_ct}-gene-expr-wilcox-tests.xlsx` — Gene expression tests
- `{safe_ct}-stats.qs` / `.tsv` — Per-group summary statistics

Per level in `by-{level}/`:
- `module-score-wilcox-tests.xlsx` — Module score tests across celltypes
- `combined-{gs_name}-violin.pdf` — Multi-celltype grid figure
- `summary-{level}.qs` / `.tsv` — Level summary

Root:
- `module-gene-sets.xlsx` — Gene set definitions with presence check

## Tmux Launch Pattern

Use a bash launcher script for long-running jobs:

```bash
SESSION_PREFIX="af8362_1709"
TIMESTAMP="$(date '+%Y%m%d_%H%M%S')"
SESSION_NAME="${SESSION_PREFIX}_${TIMESTAMP}"
LOG_DIR="${OUTPUT_DIR}/logs-launch_${TIMESTAMP}"

# Write job script, then:
tmux new-session -d -s "${SESSION_NAME}" "${JOB_SCRIPT}"
```

Workflow: test with single celltype first (`--celltypes=Mono`), then launch full run.

## Bioinformatics Context

### Variant: 8362T>G in MT-TK (tRNA-Lys)

- Affects mitochondrial **translation** globally (not just one protein)
- Relevant gene sets: Mito Translation, mtUPR, OXPHOS, Glycolysis, Oxidative Stress
- NOT in MT-ATP6 — common misconception about position 8362

### Standard Gene Sets (11)

See [references/gene-sets.md](references/gene-sets.md) for full gene lists.

| Name                 | Label                | N genes | Relevance                           |
| -------------------- | -------------------- | ------- | ----------------------------------- |
| inflammatory         | Inflammatory         | 31      | Cytokines, NF-κB, IFN, stress       |
| mito                 | Mitochondrial        | 15      | 13 MT-encoded proteins + rRNA       |
| oxphos               | OXPHOS               | 33      | Complex I-V nuclear-encoded         |
| apoptosis            | Apoptosis            | 23      | Pro/anti-apoptotic, death receptors |
| ifn_type1            | Type I IFN Response  | 29      | ISGs, IFN signaling                 |
| antigen_presentation | Antigen Presentation | 22      | MHC I/II, processing                |
| mito_translation     | Mito Translation     | 16      | MT elongation, ribosomes, aaRS      |
| mtupr                | mtUPR                | 13      | MT chaperones, proteases, ISR       |
| glycolysis           | Glycolysis           | 16      | Full pathway + transporters         |
| oxidative_stress     | Oxidative Stress     | 15      | SOD, catalase, NRF2, glutathione    |
| nfkb                 | NF-κB Signaling      | 12      | Subunits, IKK, feedback             |

## Checklist

Before running:
- [ ] Gene sets defined and registered in `MODULE_GENE_SETS`
- [ ] `score_col` ends with `1` for each gene set
- [ ] Group constants defined: `GROUP_ORDER`, `GROUP_COLORS`, `GROUP_LABELS`, `PAIRWISE_COMPARISONS`
- [ ] `load_pkg(jutils)` first, then `Seurat`, `ggpubr`, `writexl`
- [ ] Output directories created with `fs::dir_create()`
- [ ] `saveplot` called with `device = "pdf"` explicitly
- [ ] Test with single celltype before full run
- [ ] Clean old results before re-run (`rm -rf by-L1 by-L2 summary-*`)
