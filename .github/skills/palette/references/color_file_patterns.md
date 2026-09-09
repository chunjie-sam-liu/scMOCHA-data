# Color File Patterns

Use these patterns for project-level R color files:

- `src/color.R`, and `src_<variant>/color.R` in a repo with more than one
  analysis track
- `src/colors.R`
- `src/plot_colors.R`

If one exists, use it. If multiple exist, prefer the one already sourced by plotting scripts. If none exist, ask before creating one.

In a multi-track repo, load the file belonging to the track the script lives
in. The files are independent; they may hold the same object names today and
diverge later, and nothing warns you when the wrong one is loaded.

```r
source("src/color.R")            # a src/ script
source("src_<variant>/color.R")  # a src_<variant>/ script
```

## Canonical Structure

Every object carries the `color_` prefix when that is the convention the
existing file uses. Check the file before assuming it.

```r
color_celltype <- c(
  "B cell" = "#4E79A7",
  "T cell" = "#F28E2B",
  "Monocyte" = "#E15759"
)

color_celltype_light <- prismatic::clr_lighten(color_celltype, 0.25)
color_celltype_alpha <- prismatic::clr_alpha(color_celltype, 0.75)

color_expression <- function(n = 256) {
  paletteer::paletteer_c("scico::batlow", n = n)
}

color_signed_score <- function(n = 256) {
  paletteer::paletteer_c("scico::vik", n = n)
}
```

## ggplot2

```r
source("src/color.R")

p + scale_color_manual(values = color_celltype)
p + scale_fill_manual(values = color_celltype, na.value = "grey80")
p + scale_color_gradientn(colors = color_expression())
```

## ComplexHeatmap

```r
source("src/color.R")

ha <- ComplexHeatmap::HeatmapAnnotation(
  celltype = celltype,
  col = list(celltype = color_celltype)
)
```

For continuous heatmap values, create a separate color function in the color file:

```r
color_expr_fun <- circlize::colorRamp2(
  c(0, 1, 2),
  color_expression(3)
)
```

## openxlsx2

Workbook chrome and block fills are colors too, and come from the same file.

```r
source("src/color.R")

wb$add_fill(
  sheet = sheet, dims = dims,
  color = openxlsx2::wb_color(hex = color_xlsx_hdr)
)
wb$add_font(
  sheet = sheet, dims = dims,
  color = openxlsx2::wb_color(hex = color_xlsx_white), bold = TRUE
)

# A cohort column block reuses the cohort's own color, never a new hex.
wb$add_fill(
  sheet = sheet, dims = dims,
  color = openxlsx2::wb_color(hex = color_newgroup[["AA-Primary"]])
)
```

`wb_color(hex = )` needs a 6-digit hex. Strip an alpha suffix first:
`substr(color_group$nhw, 1, 7)`.

## pheatmap

```r
source("src/color.R")

pheatmap::pheatmap(
  mat,
  annotation_col = annotation_col,
  annotation_colors = list(celltype = color_celltype),
  color = color_expression(100)
)
```

## Seurat

```r
source("src/color.R")

DimPlot(obj, group.by = "celltype", cols = color_celltype)
FeaturePlot(obj, features = "GENE", cols = color_expression(100))
```

Before using a named vector, check that names match observed levels:

```r
setdiff(levels(obj$celltype), names(color_celltype))
setdiff(names(color_celltype), levels(obj$celltype))
```

