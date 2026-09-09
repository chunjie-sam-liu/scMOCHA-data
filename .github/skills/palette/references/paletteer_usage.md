# Paletteer Usage

Use `paletteer` to select palettes and retrieve stable hex colors for R scientific plots.

## Choosing Palette Type

- Discrete categories: use `paletteer_d()` or `paletteer_dynamic()`.
- Many categories with variable `n`: prefer `paletteer_dynamic()` when the palette supports it.
- Continuous nonnegative values: use `paletteer_c()` with a sequential palette.
- Signed values centered at zero: use `paletteer_c()` with a diverging palette and make the midpoint explicit in the plot scale.

## Exploratory ggplot Scales

Inline scales are acceptable for quick exploration:

```r
p + paletteer::scale_color_paletteer_d("ggsci::default_jama")
p + paletteer::scale_fill_paletteer_c("scico::batlow")
```

Do not use these inline scales as the default for final shared figure code when stable category mappings matter.

## Fixed Discrete Colors

For final/shared figures, retrieve colors once and store a named vector in the project color file:

```r
celltype_colors <- paletteer::paletteer_d("ggsci::default_jama")
celltype_colors <- as.character(celltype_colors[seq_along(celltype_levels)])
names(celltype_colors) <- celltype_levels
```

After fixing the mapping, plots should use:

```r
p + scale_color_manual(values = celltype_colors)
p + scale_fill_manual(values = celltype_colors)
```

## Continuous Palette Functions

Store continuous palettes as functions:

```r
expression_colors <- function(n = 256) {
  paletteer::paletteer_c("scico::batlow", n = n)
}
```

Use with ggplot:

```r
p + scale_color_gradientn(colors = expression_colors())
p + scale_fill_gradientn(colors = expression_colors())
```

## Diverging Palette Functions

For signed statistics, document the intended midpoint:

```r
signed_score_colors <- function(n = 256) {
  paletteer::paletteer_c("scico::vik", n = n)
}
```

Use with ggplot:

```r
p + scale_color_gradient2(
  low = signed_score_colors(3)[1],
  mid = signed_score_colors(3)[2],
  high = signed_score_colors(3)[3],
  midpoint = 0
)
```

For smoother diverging scales, use `scale_color_gradientn()` with explicit limits and midpoint-aware rescaling when needed.

