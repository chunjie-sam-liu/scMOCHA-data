---
name: jutils
description: Guide for using the jutils R package - a personal utility toolkit for data import/export, DuckDB database workflows, parallel processing, and plotting helpers. Use when writing R code that needs to read/write files (CSV, TSV, Parquet, Excel, JSON, YAML, RDS, QS, FST), work with DuckDB databases, run parallel computations with progress bars, format numbers for plots, or load multiple packages efficiently.
---

# jutils R Package

jutils is a personal R utility toolkit that provides unified interfaces for common data science tasks. Use this skill when writing R code that involves file I/O, database operations, parallel processing, or plotting utilities.

## Overview

The package provides:
1. **Unified Data Import/Export** - Single `import()`/`export()` functions for 10+ file formats
2. **DuckDB Database Workflow** - Connection management, lazy views, and efficient exports
3. **Parallel Processing** - `pbmclapply`/`pbmcmapply` with progress bars
4. **Plotting Helpers** - Axis breaks, human-readable numbers, LaTeX p-values
5. **Utility Functions** - Package loading, environment variables

## Package Loading

When `library(jutils)` or `load_pkg(jutils)` is called, it automatically loads common packages:
- **Data**: data.table, dplyr, dtplyr, dbplyr, tidyr, purrr, tibble, rlang, arrow
- **Visualization**: ggplot2, patchwork, prismatic, paletteer
- **Infrastructure**: here, glue, fs, parallel, GetoptLong
- **Logging**: logger, cli

It also resolves package conflicts using the `conflicted` package with these
explicit preferences:
- **dplyr** wins: `filter`, `select`, `mutate`, `summarise`, `arrange`, `lag`,
  `first`, `last`, `between`
- **dbplyr** wins: `sql`
- **data.table** wins: `transpose`
- **purrr** wins: `set_names`

---

## Data Import/Export

### import() - Read files with automatic format detection

```r
# Basic usage - format auto-detected from extension
df <- import("data.csv")              # CSV (lazy arrow Dataset by default)
df <- import("data.tsv")              # TSV (lazy arrow Dataset by default)
df <- import("data.parquet")          # Parquet (lazy arrow Dataset by default)
df <- import("data.xlsx")             # Excel (always eager, returns data.table)
df <- import("data.json")             # JSON (always eager, returns list/data.frame)
df <- import("data.yaml")             # YAML (always eager, returns list)
df <- import("data.rds")              # RDS (always eager, returns R object)
df <- import("data.qs")               # QS2 format (always eager)
df <- import("data.fst")              # FST format (always eager, returns data.table)

# Eager loading (load all data into memory)
df <- import("data.csv", lazy = FALSE)      # Returns data.table (fread, nThread=8)
df <- import("data.parquet", lazy = FALSE)  # Returns arrow Table

# Lazy loading (default for csv/tsv/parquet)
ds <- import("large_data.csv", lazy = TRUE)  # Returns arrow Dataset
# Work with dplyr verbs, then collect
result <- ds |> filter(col > 100) |> collect()

# URLs are supported
df <- import("https://example.com/data.csv")

# Compressed files auto-detected
df <- import("data.csv.gz")
df <- import("data.tsv.zst")
```

**Lazy loading only applies to csv, tsv, and parquet.** All other formats
(qs, qs2, rds, fst, xlsx, json, yaml) ignore the `lazy` parameter and always
load eagerly.

**Supported formats:**
| Extension | Lazy Support | Eager Return Type |
| --------- | ------------ | ----------------- |
| .csv      | Yes          | data.table        |
| .tsv      | Yes          | data.table        |
| .parquet  | Yes          | arrow Table       |
| .xlsx     | No           | data.table        |
| .json     | No           | list/data.frame   |
| .yaml     | No           | list              |
| .rds      | No           | R object          |
| .qs       | No           | R object          |
| .fst      | No           | data.table        |
| .env      | No           | named list        |

### export() - Write files with automatic format detection

```r
# Basic usage - format auto-detected from extension
export(df, "output.csv")
export(df, "output.tsv")
export(df, "output.parquet")
export(df, "output.xlsx")
export(df, "output.json")
export(df, "output.yaml")
export(df, "output.rds")
export(df, "output.qs2")   # use qs2, not qs (qs is deprecated)
export(df, "output.fst")

# Compressed output
export(df, "output.csv.gz")
export(df, "output.tsv.zst")

# Export to both CSV and FST simultaneously
export(df, "output", format = "both")  # Creates output.csv and output.fst

# Parquet options
export(df, "output.parquet", lazy = FALSE)  # Single file
export(df, "output_dir", format = "parquet", lazy = TRUE)  # Dataset folder
```

**IMPORTANT: `format = "qs"` is deprecated.** Always use `format = "qs2"` (or
`.qs2` extension) for new files. Using `format = "qs"` emits a warning and
saves as qs2 anyway.

**Excel export auto-styling:** When exporting to `.xlsx`, the output is
automatically styled:
- Headers: bold, centered, light gray background (`#EDEDED`)
- Numeric columns: `#,##0.00` number format, right-aligned, conditional
  color-scale formatting (light-to-dark gradient, or diverging if column has
  both negative and positive values)
- Categorical columns (factor or character with ≤ 20 unique values): each
  level gets a distinct HCL color fill
- Text columns (character with > 20 unique values): left-aligned

**Excel size limits:** Export errors if data exceeds Excel's hard limits:
- Rows: 1,048,575 maximum
- Columns: 16,384 maximum

Use parquet or CSV for larger datasets.

### convert() - Convert between formats

```r
# Memory-efficient conversion (streaming when possible)
convert("input.csv", "output.parquet")
convert("input.parquet", "output.csv")

# Lazy-to-lazy conversions avoid loading all data
convert("large.csv", "large.parquet", lazy_in = TRUE, lazy_out = TRUE)
```

---

## DuckDB Database Functions

### Connection Management

```r
# Get or create a connection (reuses existing valid connection)
conn <- db_conn("mydb.duckdb")  # File-based database
conn <- db_conn()               # In-memory database

# Read-write connection
conn <- db_conn("mydb.duckdb", read_only = FALSE)

# List all active connections
conn_ls()

# Close all connections
db_disconn()
```

### Working with Tables and Views

```r
# Import a file as a lazy view (no data loaded into memory)
lazy_tbl <- tbl_import(conn, "my_view", "data.csv")
lazy_tbl <- tbl_import(conn, "my_view", "data.parquet")

# Partial type override (recommended — only named columns are forced)
lazy_tbl <- tbl_import(conn, "sales", "sales.csv",
                       types = list(date = "DATE", amount = "DOUBLE"))

# Full schema (advanced — ALL columns must be specified)
lazy_tbl <- tbl_import(conn, "sales", "sales.csv",
                       columns = list(id = "INTEGER", name = "VARCHAR",
                                      amount = "DOUBLE"))

# Overwrite existing view (default errors if view already exists)
lazy_tbl <- tbl_import(conn, "sales", "new_sales.csv", overwrite = TRUE)

# Work with dplyr verbs (lazy evaluation)
result <- lazy_tbl |>
  filter(col > 100) |>
  group_by(category) |>
  summarise(total = sum(value)) |>
  collect()

# Export directly from DuckDB (more efficient than collect + export)
tbl_export(conn, lazy_tbl, "output.parquet")
tbl_export(conn, lazy_tbl, "output.csv")

# Overwrite existing output file
tbl_export(conn, lazy_tbl, "output.parquet", overwrite = TRUE)

# With partitioning
tbl_export(conn, lazy_tbl, "output_dir",
           format = "parquet",
           partition_by = c("year", "month"))

# List all tables/views
tbl_ls(conn)
tbl_ls(conn, show_colnames = TRUE)  # Include column names

# Drop a table or view
# In non-interactive sessions, confirm = TRUE is required or the call errors
tbl_drop(conn, "my_view", type = "view", confirm = TRUE)
tbl_drop(conn, "my_table", type = "table", confirm = TRUE)

# Analyze table for query optimization
tbl_analyze(conn, "my_table")
```

### Multi-Database Queries

```r
# Attach another database
db_attach(conn, "other.duckdb", alias = "other_db")

# Query across databases
DBI::dbGetQuery(conn, "SELECT * FROM other_db.main.table_name")

# Detach when done
db_detach(conn, "other_db")
```

### Arrow Dataset Integration

```r
# Register an Arrow dataset as a DuckDB view
ds <- import("data.parquet", lazy = TRUE)
tbl_register_arrow(conn, "arrow_view", ds)

# Now query it with DuckDB
result <- tbl(conn, "arrow_view") |>
  filter(col > 100) |>
  collect()
```

---

## Deprecated Functions

The following functions are deprecated and emit warnings. Use the replacements:

| Deprecated        | Replacement                           |
| ----------------- | ------------------------------------- |
| `view_create()`   | `tbl_import()`                        |
| `view_drop()`     | `tbl_drop(conn, name, type = "view")` |
| `table_ls()`      | `tbl_ls()`                            |
| `table_analyze()` | `tbl_analyze()`                       |

Never use deprecated functions in new code.

---

## Parallel Processing

### pbmclapply() - Parallel lapply with progress bar

```r
# Basic usage (8 cores by default)
results <- pbmclapply(1:100, function(x) {
  Sys.sleep(0.1)
  x^2
})

# Custom number of cores
results <- pbmclapply(items, process_fn, mc.cores = 4)

# Disable progress bar
results <- pbmclapply(items, process_fn, .progress = FALSE)
```

**`mc.preschedule = FALSE` is the default** (unlike base `mclapply` which
defaults to `TRUE`). This means each job runs as an independent process,
isolating errors so a failure in one job does not affect others.
Trade-off: slightly higher process-spawn overhead vs `mc.preschedule = TRUE`
where jobs are pre-assigned to cores in round-robin batches (faster for many
small jobs, but a failing job corrupts the whole batch on its core).

**Features:**
- Real-time progress bar with percentage and ETA
- Milestone alerts at 50% and 90%
- Error isolation per job (default `mc.preschedule = FALSE`)
- Automatic fallback to sequential on Windows
- Completion summary with success/failure count
- In non-interactive sessions (scripts), progress bar is skipped automatically

### pbmcmapply() - Parallel mapply with progress bar

```r
# Multiple vectorized arguments
results <- pbmcmapply(
  function(x, y) x + y,
  1:100,
  101:200,
  mc.cores = 8
)

# With additional fixed arguments
results <- pbmcmapply(
  process_fn,
  files,
  params,
  MoreArgs = list(verbose = TRUE),
  mc.cores = 4
)
```

---

## Plotting Helpers

### fn_xy_breaks_limits() - Generate axis breaks and limits

```r
# Generate pretty breaks for a numeric vector
y_data <- c(0.5, 2.3, 4.7, 8.1)
breaks_info <- fn_xy_breaks_limits(y_data)
# Returns: list(breaks, limits, labels, step)

# Use with ggplot2
ggplot(data, aes(x, y)) +
  geom_point() +
  scale_y_continuous(
    breaks = breaks_info$breaks,
    limits = breaks_info$limits
  )

# Custom step size
breaks_info <- fn_xy_breaks_limits(y_data, step = 1)

# Custom number of breaks
breaks_info <- fn_xy_breaks_limits(y_data, n_breaks = 10)
```

### human_read() - Format numbers for readability

```r
human_read(123.456)     # "120"
human_read(0.0456)      # "0.046"
human_read(0.0000123)   # "1.23e-05"
human_read(-5.67)       # "-5.7"
```

### human_read_latex_pval() - Format p-values for LaTeX

Returns a `latex2exp::TeX()` expression object (an R `expression` built from
a LaTeX string like `$\textit{P}=0.05$`). Pass it directly to `label`
**without** `parse = TRUE` — the expression is already parsed.

```r
# Returns latex2exp expression object (pass directly, no parse = TRUE needed)
label <- human_read_latex_pval("0.05")

# With statistic prefix
label <- human_read_latex_pval("1e-5", s = "R = 0.95")

# As character string (for non-plot contexts)
label <- human_read_latex_pval("1e-5", tex = FALSE)
# Returns: "R = 0.95, $\textit{P}=1 \times 10^{-5}$"

# Use in ggplot2 - pass expression directly (no parse = TRUE)
ggplot(data, aes(x, y)) +
  geom_point() +
  annotate("text", x = 1, y = 1,
           label = human_read_latex_pval("0.001"))

# Combine with human_read for correlation test
cor_result <- cor.test(data$x, data$y)
pval_label <- human_read_latex_pval(
  human_read(cor_result$p.value),
  s = paste0("R = ", round(cor_result$estimate, 2))
)
```

---

## Utility Functions

### load_pkg() - Load packages with colorful output

```r
# Single package
load_pkg(ggplot2)

# Multiple packages
load_pkg(ggplot2, dplyr, tidyr)

# As character vector
load_pkg(c("ggplot2", "dplyr"))

# Mixed input
load_pkg(ggplot2, c("dplyr", "tidyr"), "patchwork")

# Silent loading
load_pkg(ggplot2, verbose = FALSE)
```

### dotenv() - Load environment variables from .env file

```r
# Load from default .env file
dotenv()

# Load from custom file
dotenv("config/.env.production")

# Don't override existing variables
dotenv(override = FALSE)

# Verbose mode
dotenv(verbose = TRUE)
```

**Supported .env features:**
- Comments (# full line and inline)
- Quoted values (single and double quotes)
- Multiline values (""" or ''')
- Variable expansion (${VAR} and $VAR)
- Escape sequences (\n, \t, \\, \", \')

---

## Common Patterns

### Large File Processing

```r
library(jutils)

# Lazy load large CSV
ds <- import("large_file.csv", lazy = TRUE)

# Process with dplyr (lazy)
result <- ds |>
  filter(category == "A") |>
  group_by(year) |>
  summarise(total = sum(value)) |>
  collect()

# Or convert to Parquet for faster future reads
convert("large_file.csv", "large_file.parquet")
```

### DuckDB Analytics Pipeline

```r
library(jutils)

# Create connection
conn <- db_conn("analytics.duckdb", read_only = FALSE)

# Import data as views (use types= to control column types)
tbl_import(conn, "sales", "sales.parquet")
tbl_import(conn, "products", "products.csv",
           types = list(product_id = "VARCHAR", price = "DOUBLE"))

# Query with dplyr
result <- tbl(conn, "sales") |>
  left_join(tbl(conn, "products"), by = "product_id") |>
  group_by(category) |>
  summarise(revenue = sum(price * quantity)) |>
  collect()

# Export result
export(result, "revenue_by_category.xlsx")

# Clean up
db_disconn()
```

### Parallel Processing with Progress

```r
library(jutils)

# Process many files in parallel
files <- fs::dir_ls("data/", glob = "*.csv")

results <- pbmclapply(files, function(f) {
  df <- import(f, lazy = FALSE)
  # Process...
  nrow(df)
}, mc.cores = 8)
```

### Publication-Ready Plots

```r
library(jutils)

# Prepare data
data <- import("results.csv", lazy = FALSE)

# Calculate correlation
cor_result <- cor.test(data$x, data$y)

# Create plot with formatted p-value
pval_expr <- human_read_latex_pval(
  human_read(cor_result$p.value),
  s = paste0("R = ", round(cor_result$estimate, 2))
)

ggplot(data, aes(x, y)) +
  geom_point() +
  geom_smooth(method = "lm") +
  annotate("text", x = Inf, y = Inf,
           label = pval_expr,
           hjust = 1.1, vjust = 1.5)
```

---

## Checklist

When using jutils in R code:

- [ ] Use `import()` instead of `read.csv()`, `fread()`, `read_parquet()`, etc.
- [ ] Use `export()` instead of `write.csv()`, `fwrite()`, `write_parquet()`, etc.
- [ ] Use `lazy = TRUE` for large CSV/TSV/Parquet (default); other formats ignore it
- [ ] Use `qs2` format (not `qs`) — `qs` is deprecated
- [ ] Use `db_conn()` for DuckDB connections (manages global connections)
- [ ] Use `tbl_import()` to create lazy views from files
- [ ] Use `types =` in `tbl_import()` to control column types (partial override)
- [ ] Use `overwrite = TRUE` in `tbl_import()` / `tbl_export()` when replacing
- [ ] Use `confirm = TRUE` in `tbl_drop()` for non-interactive sessions
- [ ] Never use deprecated functions: `view_create`, `view_drop`, `table_ls`,
  `table_analyze`
- [ ] Use `pbmclapply()` instead of `mclapply()` for progress bars
- [ ] Pass `human_read_latex_pval()` expression directly to `label` (no
  `parse = TRUE`)
- [ ] Call `db_disconn()` when done with database work
- [ ] Remember that `library(jutils)` auto-loads common packages
- [ ] Excel export errors if data > 1,048,575 rows or > 16,384 columns
