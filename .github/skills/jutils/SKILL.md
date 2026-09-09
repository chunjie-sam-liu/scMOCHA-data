---
name: jutils
description: Guide for using the jutils R package - a personal utility toolkit for data import/export, DuckDB database workflows, parallel processing, and plotting helpers. Use when writing R code that needs to read/write files (CSV, TSV, Parquet, Excel, JSON, YAML, RDS, QS, FST), work with DuckDB databases, run parallel computations with progress bars, format numbers for plots, save plots, or load multiple packages efficiently.
---

# jutils R Package

jutils is a personal R utility toolkit. Attaching it auto-loads a curated set
of common packages and resolves namespace conflicts via `conflicted`, so one
entry line replaces dozens of `library()` calls. That entry line is
`suppressMessages({library(jutils)})`, placed in the `# Library` section at the
top of the script. `jutils` is the only package ever loaded with `library()`;
everything else goes through `load_pkg(...)`.

---

## Auto-loaded packages

Attaching jutils loads these 20 packages via `load_pkg()` in `.onAttach`.
No extra `library()` or `load_pkg()` call is needed for any of them.

| Area           | Packages                                                              |
| -------------- | --------------------------------------------------------------------- |
| Data core      | data.table, dplyr, dtplyr, dbplyr, tidyr, purrr, tibble, rlang, arrow |
| Plotting       | ggplot2, patchwork                                                    |
| Colors         | prismatic, paletteer                                                  |
| Infrastructure | here, glue, fs, parallel, GetoptLong                                  |
| Logging / UX   | logger, cli                                                           |

13 of these are only in `Suggests`. If one is missing, attaching jutils does
not fail and, under `Rscript`, does not even print a message — the error
surfaces later at the first call. See [references/gotchas.md](references/gotchas.md).

---

## Conflict resolution

jutils uses `conflicted` to set **deterministic** conflict winners. There are
12 rules, registered with `conflicted::conflict_prefer(name, pkg)`. No losers
are specified, so the winner takes precedence over **every** other package:

| Function                                                   | Winner         |
| ---------------------------------------------------------- | -------------- |
| `filter()` `select()` `mutate()` `summarise()` `arrange()` | **dplyr**      |
| `lag()` `first()` `last()` `between()`                     | **dplyr**      |
| `sql()`                                                    | **dbplyr**     |
| `transpose()`                                              | **data.table** |
| `set_names()`                                              | **purrr**      |

**Rule of thumb:** dplyr semantics are the default for common verbs.
data.table wins for `transpose()`. purrr wins for `set_names()`.

These rules are enforced by a `.conflicts` shim environment that `conflicted`
puts at the top of the search path — note that `package:conflicted` itself is
**not** attached. Conflicts outside this table remain ambiguous and raise an
error at call time; use an explicit namespace prefix for those.

When writing scripts that attach jutils, you do NOT need to prefix
`dplyr::filter()` etc. — the conflicts are already resolved.

---

## Exported functions

Detailed signatures, parameters, and examples are in `references/`.

**Read [references/gotchas.md](references/gotchas.md) before writing jutils
code.** Several functions silently produce wrong files or wrong types.

### Data I/O — [references/io.md](references/io.md)

| Function    | Purpose                      | Key example                                                                  |
| ----------- | ---------------------------- | ---------------------------------------------------------------------------- |
| `import()`  | Read any file (10+ formats)  | `import("data.csv")` → lazy arrow Dataset                                    |
| `export()`  | Write any file (10+ formats) | `export(df, "out.parquet")` or `export(df, "out", format = c("csv", "fst"))` |
| `convert()` | Convert between formats      | `convert("data.csv", "data.parquet")`                                        |

- Lazy by default for csv/tsv/parquet; call `collect()` to materialize
- Eager types differ: csv/tsv → data.table, parquet/xlsx → tibble, fst →
  data.frame. Call `setDT()` when you need data.table semantics
- **Name qs files `.qs`, never `.qs2`** — `export(df, "out.qs2")` writes
  `out.qs`, so the round trip breaks
- **`.gz` is the only compression suffix that works** with `import()`/
  `export()`, and only on the data.frame path. `.zst`/`.bz2`/`.xz`/`.zip` and
  Arrow inputs all error. Use `tbl_export()` or parquet for zstd
- Excel auto-styles; errors if > 1M rows or > 16K columns
- URLs supported for import (silently forces eager mode)
- `export(format = c("csv", "fst"))` exports to multiple formats at once
- `export(format = "both")` is legacy shorthand for `c("csv", "fst")`
- Named list of data.frames → multi-sheet Excel
- `convert()` forwards `...` to **both** `import()` and `export()`

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

- `db_conn()` defaults to `read_only = TRUE`
- **`read_only` is not part of the connection pool key** — reusing a path
  silently returns the existing connection in its original mode. Call
  `db_disconn()` before switching modes
- `tbl_drop()` requires `confirm = TRUE` in non-interactive sessions
- Use `types =` (not `columns =`) in `tbl_import()` for partial type override
- `tbl_import()`'s parquet branch ignores `header`/`delim`/`columns`/`types`
- `format = "auto"` falls back to `csv` for unknown extensions, it does not error
- `tbl_ls()` enforces `check_dots_empty()`; use named arguments only
- `tbl_export()` uses DuckDB COPY — more efficient than `collect()` + `export()`,
  and the only writer that supports `.zst` for CSV/TSV
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
- **Failures never abort** — failed elements return as `try-error` objects in
  the result list; check them explicitly
- Progress bar requires a TTY. Redirected to a log file you get only 50%/90%
  milestones — that is normal, not a hang

### Plotting helpers — [references/plot.md](references/plot.md)

| Function                  | Purpose                        | Key example                                  |
| ------------------------- | ------------------------------ | -------------------------------------------- |
| `fn_xy_breaks_limits()`   | Pretty axis breaks for ggplot2 | `fn_xy_breaks_limits(vec, n_breaks = 5)`     |
| `human_read()`            | Format numbers readably        | `human_read(0.0456)` → `"0.046"`             |
| `human_read_latex_pval()` | P-values for LaTeX/plots       | `human_read_latex_pval("1e-5", s = "R=0.9")` |
| `saveplot()`              | Save plots (single/multi-page) | `saveplot("fig.pdf", p, width = 8)`          |

- Pass `human_read_latex_pval()` directly to `label` (no `parse = TRUE`)
- `human_read()` **errors on `NA`** — filter before calling
- `saveplot()` multi-page for PDF/TIFF; numbered files for PNG/JPEG,
  zero-padded to the plot count (`fig_01.png` … `fig_10.png`)
- `saveplot()` writes a blank page and only warns when a plot fails to render
- `saveplot()` auto-creates output directories
- `fn_xy_breaks_limits()` returns `limits, breaks, labels, step` — index by name

### Utilities — [references/utils.md](references/utils.md)

| Function     | Purpose                             | Key example                             |
| ------------ | ----------------------------------- | --------------------------------------- |
| `dotenv()`   | Load .env files                     | `dotenv(".env.prod", override = FALSE)` |
| `load_pkg()` | Load extra packages with CLI output | `load_pkg(ComplexHeatmap, circlize)`    |

- `dotenv()` supports comments, quoted values, multiline (`"""`), variable
  expansion (`${VAR}`), escape sequences, `export` prefix
- `dotenv()` defaults to `override = TRUE` — it replaces existing environment
  variables unless you pass `override = FALSE`
- `dotenv()` skips malformed lines **silently**; a typo yields no error
- `load_pkg()` accepts unquoted names, strings, character vectors, or any mix
- `load_pkg()` is silent under `Rscript` and returns `FALSE` for a missing
  package instead of erroring — use `rlang::check_installed()` for hard
  dependencies in scripts

### Developing jutils itself — [references/dev.md](references/dev.md)

Use **pixi** for every command in this repository (`pixi run test`,
`pixi run doc`, `pixi run check`, `pixi run format`). Never use conda or mamba.

---

## Deprecated — never use in new code

All four are still exported and emit a `lifecycle::deprecate_warn()`.

| Deprecated        | Replacement                           |
| ----------------- | ------------------------------------- |
| `view_create()`   | `tbl_import()`                        |
| `view_drop()`     | `tbl_drop(conn, name, type = "view")` |
| `table_ls()`      | `tbl_ls()`                            |
| `table_analyze()` | `tbl_analyze()`                       |

---

## Coding style guide for jutils scripts

When writing R scripts that attach jutils, follow these conventions:

### Script setup

- Start every script with `suppressMessages({library(jutils)})` in the
  `# Library` section. `jutils` is the only package ever loaded with
  `library()`.
- The `rmeta` snippet repeats a bare `library(jutils)` in the `# Load data`
  section. That repeat is intentional and harmless — attaching is idempotent.
  Keep it; do not delete it as a duplicate.
- Load anything jutils does not auto-attach with `load_pkg(...)`, which takes
  unquoted names, strings, or a character vector:
  `load_pkg(ComplexHeatmap, circlize)`.
- Use `dotenv()` to load environment variables from `.env` files
- Source the stage's `config.R` in the `Source` section, right after
  `dotenv()`, then call `stage_paths()`. Stage constants, output paths, and
  helpers shared by two or more steps live there, never inline in a step. The
  full contract is in the `analysis-pipeline` skill,
  `references/script-templates.md` section G.
- Build every path from the env-file variable that owns that root, resolved
  through `stage_paths()`. **Do not rebuild a path from the repository root**,
  with `here::here()` or otherwise: the roots are configured per project and
  are frequently symlinks pointing outside the repository, so
  `here("data", ...)` and `${datadir}` are not the same directory. `here()` is
  loaded by jutils and is fine for a throwaway interactive lookup, never in a
  script that ships. -> `data-result-layout`
- The `rmeta` editor snippet emits the full section skeleton: Metainfo /
  Reproducibility / Library / Args / Logger / Load data / Source / Conn /
  Function / Main / Save / Session info.

### Pipe and anonymous functions

- Use the base pipe `|>` (never the magrittr pipe `%>%`)
- Use `\() ...` for single-line anonymous functions
- Use `function() { ... }` for multi-line anonymous functions
- Do NOT use `_$x` or `_[["x"]]` with the `_` pipe placeholder

### Data I/O patterns

- Default to `import()` / `export()`. Do not reach for `read_parquet()`,
  `read.csv()`, or `readxl` directly.
- `fread()` / `fwrite()` are the one documented exception, for cases where the
  exact delimiter, quoting, `na` string, or header form matters. PLINK2
  `#FID IID ...` files are the standard case; see the `analysis-pipeline`
  script templates.
- Use `convert()` for format conversion (memory-efficient)
- `import()` returns lazy arrow Datasets for csv/tsv/parquet by default —
  always `collect()` before using as data.table/data.frame
- After an eager `import()`, only csv/tsv are data.tables. Call
  `data.table::setDT()` for parquet, fst, and xlsx
- Name qs files `.qs`, never `.qs2`
- `.gz` is the only compression suffix `import()`/`export()` handle, and only
  for data.frame input
- For DuckDB workflows, prefer `tbl_import()` + dplyr over `import()` to
  keep data out of R memory

### Naming and style

- Use `snake_case` for all names
- 2-space indentation, 80-character line width
- Format with the repository's own format task when it has one. In jutils that
  is `pixi run format`, which runs `air format .` across the package — correct
  here because the whole tree is the package. In an analysis repository
  without such a task, run `air format <file>` on the files you touched and
  never `air format .` on the whole tree
- In this repository, run every command through **pixi**. Never use
  `conda activate`, `mamba activate`, or `conda run -n renv`

### Error and messaging

- Use `cli::cli_abort()` for user-facing errors
- Use `cli::cli_alert_info()` / `cli::cli_alert_success()` for messages
- Use `rlang::check_installed()` before using optional packages
- Use explicit namespace prefix (e.g., `DBI::dbGetQuery()`) for packages
  NOT auto-loaded by jutils

### Parallel processing

- Use `pbmclapply()` / `pbmcmapply()` instead of `parallel::mclapply()`
- Wrap worker functions in `tryCatch()` for error isolation
- Check the result for `try-error` elements — these functions never abort
- Decide `read_only` once per process; `db_conn()` ignores it for an
  already-pooled path
- Always clean up DuckDB connections in parallel workers:
  `on.exit(db_disconn())`

### Plotting

- Use ggplot2 + patchwork (both auto-loaded)
- Use `saveplot()` instead of `ggsave()` — supports multi-page output
- Use `fn_xy_breaks_limits()` for axis scaling
- Use `human_read_latex_pval()` for p-value annotations (pass directly to
  `label`, NOT with `parse = TRUE`)
- Guard `NA` before calling `human_read()`
- Do not ignore warnings from `saveplot()` — a failed plot becomes a blank page

---

## Canonical patterns

More worked examples (large-file streaming, parallel DuckDB, publication
plots, multi-page export, `dotenv()`) are in `references/`. Known traps are
collected in [references/gotchas.md](references/gotchas.md).

### Typical analysis script

```r
suppressMessages({library(jutils)})

df <- import(paths$raw_experiment_csv, lazy = FALSE)

result <- df |>
  filter(!is.na(value)) |>
  group_by(group) |>
  summarise(
    mean_val = mean(value),
    sd_val = sd(value),
    n = n()
  )

p <- ggplot(result, aes(group, mean_val)) +
  geom_col() +
  geom_errorbar(aes(ymin = mean_val - sd_val, ymax = mean_val + sd_val))

export(result, fs::path(paths$tabledir, "summary.csv"))
saveplot(fs::path(paths$figuredir, "figure.pdf"), p, width = 8, height = 6)
```

### DuckDB pipeline (keeps data out of R memory)

```r
suppressMessages({library(jutils)})
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
