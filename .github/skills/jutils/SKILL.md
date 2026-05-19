---
name: jutils
description: Guide for using the jutils R package - a personal utility toolkit for data import/export, DuckDB database workflows, parallel processing, and plotting helpers. Use when writing R code that needs to read/write files (CSV, TSV, Parquet, Excel, JSON, YAML, RDS, QS, FST), work with DuckDB databases, run parallel computations with progress bars, format numbers for plots, save plots, or load multiple packages efficiently.
---

# jutils R Package

jutils is a personal R utility toolkit. Calling `library(jutils)` auto-loads
a curated set of common packages and resolves namespace conflicts via
`conflicted`, so a single `library(jutils)` replaces dozens of
`library()` calls.

---

## Auto-loaded packages

When you call `library(jutils)`, the following packages are automatically
loaded and attached (via `load_pkg()` in `.onAttach`):

### Plotting / visualization

- **ggplot2** — grammar of graphics
- **patchwork** — composing multiple ggplots

### Colors

- **prismatic** — color manipulation utilities
- **paletteer** — comprehensive color palette collection

### Data core

- **data.table** — fast data manipulation (fread, fwrite, `:=`, .SD, etc.)
- **dplyr** — data manipulation verbs (filter, select, mutate, etc.)
- **dtplyr** — data.table backend for dplyr (lazy translation)
- **dbplyr** — database backend for dplyr
- **tidyr** — tidy data reshaping (pivot_longer, pivot_wider, etc.)
- **purrr** — functional programming (map, walk, etc.)
- **tibble** — modern data frames
- **rlang** — tidy evaluation and metaprogramming
- **arrow** — Apache Arrow for columnar data (Parquet, CSV, etc.)

### Infrastructure

- **here** — project-relative file paths
- **glue** — string interpolation
- **fs** — cross-platform filesystem operations
- **parallel** — base R parallel computing
- **GetoptLong** — command-line argument parsing

### Logging / UX

- **logger** — structured logging
- **cli** — rich CLI output (progress bars, colors, etc.)

### Summary

After `library(jutils)`, you have 20 packages available without extra
`library()` calls: ggplot2, patchwork, prismatic, paletteer, data.table,
dplyr, dtplyr, dbplyr, tidyr, purrr, tibble, rlang, arrow, here, glue, fs,
parallel, GetoptLong, logger, cli.

---

## Conflict resolution

jutils uses `conflicted` to set **deterministic** conflict winners. These
are the resolved preferences:

| Function      | Winner         | Losers                    |
| ------------- | -------------- | ------------------------- |
| `filter()`    | **dplyr**      | stats, data.table         |
| `select()`    | **dplyr**      | MASS, data.table          |
| `mutate()`    | **dplyr**      | data.table (hypothetical) |
| `summarise()` | **dplyr**      | data.table (hypothetical) |
| `arrange()`   | **dplyr**      | data.table (hypothetical) |
| `lag()`       | **dplyr**      | stats                     |
| `first()`     | **dplyr**      | data.table                |
| `last()`      | **dplyr**      | data.table                |
| `between()`   | **dplyr**      | data.table                |
| `sql()`       | **dbplyr**     | DBI                       |
| `transpose()` | **data.table** | purrr                     |
| `set_names()` | **purrr**      | stats                     |

**Rule of thumb:** dplyr semantics are the default for common verbs.
data.table wins for `transpose()`. purrr wins for `set_names()`.

When writing scripts with `library(jutils)`, you do NOT need to prefix
`dplyr::filter()` etc. — the conflicts are already resolved.

---

## Exported functions

Detailed signatures, parameters, and examples are in `references/`.

### Data I/O — [references/io.md](references/io.md)

| Function    | Purpose                      | Key example                                                                  |
| ----------- | ---------------------------- | ---------------------------------------------------------------------------- |
| `import()`  | Read any file (10+ formats)  | `import("data.csv")` → lazy arrow Dataset                                    |
| `export()`  | Write any file (10+ formats) | `export(df, "out.parquet")` or `export(df, "out", format = c("csv", "fst"))` |
| `convert()` | Convert between formats      | `convert("data.csv", "data.parquet")`                                        |

- Lazy by default for csv/tsv/parquet; call `collect()` to materialize
- Use `.qs2` not `.qs` — qs format is deprecated
- Excel auto-styles; errors if > 1M rows or > 16K columns
- Compressed files auto-detected (.gz, .zst, .bz2)
- URLs supported for import (forces eager mode)
- `export(format = c("csv", "fst"))` exports to multiple formats at once
- `export(format = "both")` is legacy shorthand for `c("csv", "fst")`
- Named list of data.frames → multi-sheet Excel

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
- Connection pool is per-process; parallel workers get their own connections

### Parallel processing — [references/parallel.md](references/parallel.md)

| Function       | Purpose                        | Key example                           |
| -------------- | ------------------------------ | ------------------------------------- |
| `pbmclapply()` | Parallel lapply + progress bar | `pbmclapply(1:100, fn, mc.cores = 8)` |
| `pbmcmapply()` | Parallel mapply + progress bar | `pbmcmapply(fn, x, y, mc.cores = 8)`  |

- Falls back to sequential on Windows or `mc.cores = 1`
- `mc.preschedule = TRUE` (default): faster but failure affects whole batch
- `mc.preschedule = FALSE`: isolates errors per job
- Default `mc.cores` is `getOption("mc.cores", 8L)`

### Plotting helpers — [references/plot.md](references/plot.md)

| Function                  | Purpose                        | Key example                                  |
| ------------------------- | ------------------------------ | -------------------------------------------- |
| `fn_xy_breaks_limits()`   | Pretty axis breaks for ggplot2 | `fn_xy_breaks_limits(vec, n_breaks = 5)`     |
| `human_read()`            | Format numbers readably        | `human_read(0.0456)` → `"0.046"`             |
| `human_read_latex_pval()` | P-values for LaTeX/plots       | `human_read_latex_pval("1e-5", s = "R=0.9")` |
| `saveplot()`              | Save plots (single/multi-page) | `saveplot("fig.pdf", p, width = 8)`          |

- Pass `human_read_latex_pval()` directly to `label` (no `parse = TRUE`)
- `saveplot()` multi-page for PDF/TIFF; numbered files for PNG/JPEG
- `saveplot()` auto-creates output directories

### Utilities — [references/utils.md](references/utils.md)

| Function     | Purpose                       | Key example                             |
| ------------ | ----------------------------- | --------------------------------------- |
| `dotenv()`   | Load .env files               | `dotenv(".env.prod", override = FALSE)` |
| `load_pkg()` | Load packages with CLI output | `load_pkg(ggplot2, dplyr, tidyr)`       |

- `dotenv()` supports comments, quoted values, multiline (`"""`), variable
  expansion (`${VAR}`), escape sequences, `export` prefix
- `load_pkg()` accepts unquoted names, strings, character vectors, or any mix

---

## Deprecated — never use in new code

| Deprecated        | Replacement                           |
| ----------------- | ------------------------------------- |
| `view_create()`   | `tbl_import()`                        |
| `view_drop()`     | `tbl_drop(conn, name, type = "view")` |
| `table_ls()`      | `tbl_ls()`                            |
| `table_analyze()` | `tbl_analyze()`                       |

---

## Coding style guide for jutils scripts

When writing R scripts that use `library(jutils)`, follow these conventions:

### Script setup

- Start with `library(jutils)` — this is the ONLY `library()` call needed
  for the 20 auto-loaded packages
- Only add extra `library()` calls for packages NOT auto-loaded by jutils
- Use `dotenv()` to load environment variables from `.env` files
- Use `here::here()` for project-relative paths (loaded by jutils)

### Pipe and anonymous functions

- Use base pipe `|>` (NOT magrittr pipe `%>%`)
- Use `\() ...` for single-line anonymous functions
- Use `function() { ... }` for multi-line anonymous functions
- Do NOT use `_$x` or `_$[["x"]]` with pipe placeholder

### Data I/O patterns

- Use `import()` / `export()` — never call `fread()`, `fwrite()`,
  `read_parquet()`, etc. directly
- Use `convert()` for format conversion (memory-efficient)
- `import()` returns lazy arrow Datasets for csv/tsv/parquet by default —
  always `collect()` before using as data.table/data.frame
- For DuckDB workflows, prefer `tbl_import()` + dplyr over `import()` to
  keep data out of R memory

### Naming and style

- Use `snake_case` for all names
- 2-space indentation, 80-character line width
- Format all code with `air format .`

### Error and messaging

- Use `cli::cli_abort()` for user-facing errors
- Use `cli::cli_alert_info()` / `cli::cli_alert_success()` for messages
- Use `rlang::check_installed()` before using optional packages
- Use explicit namespace prefix (e.g., `DBI::dbGetQuery()`) for packages
  NOT auto-loaded by jutils

### Parallel processing

- Use `pbmclapply()` / `pbmcmapply()` instead of `parallel::mclapply()`
- Wrap worker functions in `tryCatch()` for error isolation
- Always clean up DuckDB connections in parallel workers:
  `on.exit(db_disconn())`

### Plotting

- Use ggplot2 + patchwork (both auto-loaded)
- Use `saveplot()` instead of `ggsave()` — supports multi-page output
- Use `fn_xy_breaks_limits()` for axis scaling
- Use `human_read_latex_pval()` for p-value annotations (pass directly to
  `label`, NOT with `parse = TRUE`)

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

### Parallel DuckDB processing

```r
library(jutils)
files <- fs::dir_ls("data/", glob = "*.parquet")
results <- pbmclapply(files, function(f) {
  conn <- db_conn("analysis.duckdb", read_only = TRUE)
  on.exit(db_disconn())
  tbl_import(conn, "tmp", f, overwrite = TRUE)
  dplyr::tbl(conn, "tmp") |>
    filter(pval < 5e-8) |>
    collect()
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

### Multi-page figure export

```r
library(jutils)
plots <- lapply(unique(mtcars$cyl), function(c) {
  ggplot(mtcars |> filter(cyl == c), aes(wt, mpg)) +
    geom_point() +
    labs(title = glue("Cylinders: {c}"))
})
saveplot("all_figures.pdf", plots, width = 8, height = 6)
```

### Environment variables

```r
library(jutils)
dotenv(".env")
db_host <- Sys.getenv("DB_HOST")
api_key <- Sys.getenv("API_KEY")
```

### Typical data analysis script

```r
library(jutils)
# jutils auto-loads: data.table, dplyr, ggplot2, arrow, patchwork,
# purrr, tidyr, tibble, fs, glue, here, cli, logger, etc.

# Load data
df <- import(here("data", "raw", "experiment.csv"), lazy = FALSE)

# Process
result <- df |>
  filter(!is.na(value)) |>
  group_by(group) |>
  summarise(
    mean_val = mean(value),
    sd_val = sd(value),
    n = n()
  )

# Plot
p <- ggplot(result, aes(group, mean_val)) +
  geom_col() +
  geom_errorbar(aes(ymin = mean_val - sd_val, ymax = mean_val + sd_val))

# Save
export(result, here("output", "summary.csv"))
saveplot(here("output", "figure.pdf"), p, width = 8, height = 6)
```
