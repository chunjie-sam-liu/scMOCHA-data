# Plotting helpers

## fn_xy_breaks_limits() — pretty axis breaks

```r
fn_xy_breaks_limits(vec, step = NULL, n_breaks = 5, max = TRUE)
```

Returns `list(breaks, limits, labels, step)` for use with ggplot2 scales.

```r
y <- c(0.5, 2.3, 4.7, 8.1)
info <- fn_xy_breaks_limits(y)

ggplot(data, aes(x, y)) +
  geom_point() +
  scale_y_continuous(breaks = info$breaks, limits = info$limits)

# Custom step
fn_xy_breaks_limits(1:100, step = 20)

# Custom number of breaks
fn_xy_breaks_limits(1:100, n_breaks = 10)

# Without max value in breaks
fn_xy_breaks_limits(1:100, max = FALSE)
```

---

## human_read() — format numbers for readability

```r
human_read(x)
```

Formats numbers with adaptive significant digits. Vectorized.

```r
human_read(123.456)       # "120"
human_read(0.0456)        # "0.046"
human_read(0.0000123)     # "1.23e-05"
human_read(-42.5)         # "-42"
human_read(0)             # "0"
human_read(c(0.5, 0.05))  # c("0.5", "0.05")
```

---

## human_read_latex_pval() — format p-values for LaTeX/plots

```r
human_read_latex_pval(x, s = NA, tex = TRUE)
```

Returns a `latex2exp::TeX()` expression object when `tex = TRUE`.
Pass directly to ggplot2 `label` — do NOT use `parse = TRUE`.

```r
# Basic
label <- human_read_latex_pval("0.05")
label <- human_read_latex_pval("1e-5")

# With statistic prefix
label <- human_read_latex_pval("0.01", s = "R = 0.85")

# As character string
str <- human_read_latex_pval("1e-5", tex = FALSE)

# Typical usage with ggplot2
cor_result <- cor.test(data$x, data$y)
pval_label <- human_read_latex_pval(
  human_read(cor_result$p.value),
  s = paste0("R = ", round(cor_result$estimate, 2))
)
ggplot(data, aes(x, y)) +
  geom_point() +
  annotate("text", x = Inf, y = Inf,
           label = pval_label, hjust = 1.1, vjust = 1.5)
```

---

## saveplot() — save plots to file

```r
saveplot(filename, plot = NULL, device = NULL, path = NULL,
         scale = 1, width = NA, height = NA, units = "in",
         dpi = 300, bg = "white", create.dir = TRUE, ...)
```

Saves single plots or lists of plots. Multi-page formats (PDF, TIFF, PS/EPS)
support multiple plots per file. Single-page formats (PNG, JPEG, BMP, SVG)
produce numbered files for lists.

```r
p <- ggplot(mtcars, aes(wt, mpg)) + geom_point()

# Single plot
saveplot("figure.pdf", p)
saveplot("figure.png", p, width = 8, height = 6)
saveplot("figure.tiff", p, dpi = 600, compression = "lzw")

# Multiple plots → multi-page PDF
plots <- list(
  ggplot(mtcars, aes(wt, mpg)) + geom_point(),
  ggplot(mtcars, aes(hp, mpg)) + geom_point()
)
saveplot("figures.pdf", plots)

# Multiple plots → numbered PNGs (figure_1.png, figure_2.png)
saveplot("figure.png", plots)

# Last plot (default)
saveplot("last.pdf")
```

**Supported formats:** pdf, png, tiff/tif, jpeg/jpg, bmp, svg, eps, ps.
