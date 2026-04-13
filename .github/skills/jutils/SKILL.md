---
name: jutils
description: Guide for using the jutils R package - a personal utility toolkit for data import/export, DuckDB database workflows, parallel processing, and plotting helpers. Use when writing R code that needs to read/write files (CSV, TSV, Parquet, Excel, JSON, YAML, RDS, QS, FST), work with DuckDB databases, run parallel computations with progress bars, format numbers for plots, save plots, or load multiple packages efficiently.
---

# jutils R Package

jutils is a personal R utility toolkit. `library(jutils)` auto-loads 20+
common packages (data.table, dplyr, ggplot2, arrow, etc.) and resolves
conflicts (dplyr wins for filter/select/mutate/etc.; data.table wins for
transpose; purrr wins for set_names).

## Exported functions

Detailed signatures, parameters, and examples are in `references/`.

### Data I/O — [references/io.md](references/io.md)

| Function    | Purpose                      | Key example                               |
| ----------- | ---------------------------- | ----------------------------------------- |
| `import()`  | Read any file (10+ formats)  | `import("data.csv")` → lazy arrow Dataset |
| `export()`  | Write any file (10+ formats) | `export(df, "out.parquet")`               |
| `convert()` | Convert between formats      | `convert("data.csv", "data.parquet")`     |

- Lazy by default for csv/tsv/parquet; call `collect()` to materialize
- Use `.qs2` not `.qs` — qs format is deprecated
- Excel auto-styles; errors if > 1M rows or > 16K columns

### DuckDB database — [references/db.md](references/db.md)

| Function               | Purpose                        | Key example                                    |
| ---------------------- | ------------------------------ | ---------------------------------------------- |
| `db_conn()`            | Get/create connection (pooled) | `conn <- db_conn("my.duckdb")`                 |
| `db_disconn()`         | Close all connections          | `db_disconn()`                                 |
| `conn_ls()`            | List active connections        | `conn_ls()`                                    |
| `db_attach()`          | Attach external database       | `db_attach(conn, "other.duckdb", alias = "o")` |
| `db_detach()`          | Detach database                | `db_detach(conn, "o")`                         |
| `tbl_import()`         | Create lazy view from file     | `tbl_import(conn, "sales", "sales.csv")`       |
| `tbl_export()`         | Export dbplyr tbl to file      | `tbl_export(conn, tbl, "out.parquet")`         |
| `tbl_ls()`             | List tables/views              | `tbl_ls(conn, show_colnames = TRUE)`           |
| `tbl_drop()`           | Drop table/view                | `tbl_drop(conn, "x", confirm = TRUE)`          |
| `tbl_analyze()`        | Update table statistics        | `tbl_analyze(conn, "big_table")`               |
| `tbl_register_arrow()` | Register Arrow dataset as view | `tbl_register_arrow(conn, "v", ds)`            |

- `tbl_drop()` requires `confirm = TRUE` in non-interactive sessions
- Use `types =` (not `columns =`) in `tbl_import()` for partial type override
- `tbl_export()` uses DuckDB COPY — more efficient than `collect()` + `export()`

### Parallel processing — [references/parallel.md](references/parallel.md)

| Function       | Purpose                        | Key example                           |
| -------------- | ------------------------------ | ------------------------------------- |
| `pbmclapply()` | Parallel lapply + progress bar | `pbmclapply(1:100, fn, mc.cores = 8)` |
| `pbmcmapply()` | Parallel mapply + progress bar | `pbmcmapply(fn, x, y, mc.cores = 8)`  |

- Falls back to sequential on Windows or `mc.cores = 1`
- `mc.preschedule = TRUE` (default): faster but failure affects whole batch
- `mc.preschedule = FALSE`: isolates errors per job

### Plotting helpers — [references/plot.md](references/plot.md)

| Function                  | Purpose                        | Key example                                  |
| ------------------------- | ------------------------------ | -------------------------------------------- |
| `fn_xy_breaks_limits()`   | Pretty axis breaks for ggplot2 | `fn_xy_breaks_limits(vec, n_breaks = 5)`     |
| `human_read()`            | Format numbers readably        | `human_read(0.0456)` → `"0.046"`             |
| `human_read_latex_pval()` | P-values for LaTeX/plots       | `human_read_latex_pval("1e-5", s = "R=0.9")` |
| `saveplot()`              | Save plots (single/multi-page) | `saveplot("fig.pdf", p, width = 8)`          |

- Pass `human_read_latex_pval()` directly to `label` (no `parse = TRUE`)
- `saveplot()` multi-page for PDF/TIFF; numbered files for PNG/JPEG

### Utilities — [references/utils.md](references/utils.md)

| Function     | Purpose                       | Key example                             |
| ------------ | ----------------------------- | --------------------------------------- |
| `dotenv()`   | Load .env files               | `dotenv(".env.prod", override = FALSE)` |
| `load_pkg()` | Load packages with CLI output | `load_pkg(ggplot2, dplyr, tidyr)`       |

---

## Deprecated — never use in new code

| Deprecated        | Replacement                           |
| ----------------- | ------------------------------------- |
| `view_create()`   | `tbl_import()`                        |
| `view_drop()`     | `tbl_drop(conn, name, type = "view")` |
| `table_ls()`      | `tbl_ls()`                            |
| `table_analyze()` | `tbl_analyze()`                       |

---

## Common patterns

### Large file processing

```r
library(jutils)
ds <- import("large.csv")
result <- ds |> filter(category == "A") |> collect()
convert("large.csv", "large.parquet")
```

### DuckDB analytics pipeline

```r
library(jutils)
conn <- db_conn("analytics.duckdb", read_only = FALSE)
tbl_import(conn, "sales", "sales.parquet")
tbl_import(conn, "products", "products.csv",
           types = list(product_id = "VARCHAR", price = "DOUBLE"))
result <- tbl(conn, "sales") |>
  left_join(tbl(conn, "products"), by = "product_id") |>
  group_by(category) |>
  summarise(revenue = sum(price * quantity)) |>
  collect()
export(result, "revenue.xlsx")
db_disconn()
```

### Parallel file processing

```r
library(jutils)
files <- fs::dir_ls("data/", glob = "*.csv")
results <- pbmclapply(files, function(f) {
  df <- import(f, lazy = FALSE)
  nrow(df)
}, mc.cores = 8)
```

### Publication-ready plot

```r
library(jutils)
data <- import("results.csv", lazy = FALSE)
cor_result <- cor.test(data$x, data$y)
p <- ggplot(data, aes(x, y)) +
  geom_point() +
  geom_smooth(method = "lm") +
  annotate("text", x = Inf, y = Inf,
           label = human_read_latex_pval(
             human_read(cor_result$p.value),
             s = paste0("R = ", round(cor_result$estimate, 2))
           ),
           hjust = 1.1, vjust = 1.5)
saveplot("correlation.pdf", p, width = 8, height = 6)
```
