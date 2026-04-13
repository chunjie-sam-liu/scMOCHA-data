# saveplot Convention

In high-res analysis scripts, `saveplot` (from `jutils`) should be called with:

```r
saveplot(
  plot_or_list,
  filename = "path/to/output.pdf",   # Use filename= (NOT file=)
  width = 7,
  height = 6,
  device = "pdf"                      # Always set explicitly
)
```

## Key Rules

1. **Use `filename=`** not `file=` — the parameter name is `filename`
2. **Set `device = "pdf"` explicitly** — avoids "Unsupported device NA" errors when filename/device inference fails
3. **List of plots → multi-page PDF**: When passed a list of ggplot objects, each plot becomes a separate page
4. **Single plot**: Pass the plot object directly

## Multi-page PDF Example

```r
plot_list <- lapply(genes, function(g) {
  ggplot(...) + ...
})

saveplot(
  plot_list,
  filename = fs::path(outdir, "gene-expression.pdf"),
  width = 7,
  height = 6,
  device = "pdf"
)
```

## Combined Panel Example

For a single page with multiple panels, use `patchwork::wrap_plots()` first:

```r
p_combined <- patchwork::wrap_plots(plot_list, ncol = 3)
saveplot(
  p_combined,
  filename = fs::path(outdir, "combined.pdf"),
  width = 7 * ncol,
  height = 6 * nrow,
  device = "pdf"
)
```
