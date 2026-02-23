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
3. **Parallel Processing** - `mclapply`/`mcmapply` with progress bars
4. **Plotting Helpers** - Axis breaks, human-readable numbers, LaTeX p-values
5. **Utility Functions** - Package loading, environment variables

## Package Loading

When `library(jutils)` is called, it automatically loads common packages:
- **Data**: data.table, dplyr, dtplyr, dbplyr, tidyr, purrr, tibble, rlang, arrow
- **Visualization**: ggplot2, patchwork, prismatic, paletteer
- **Infrastructure**: here, glue, fs, parallel, GetoptLong
- **Logging**: logger, cli

It also resolves package conflicts, preferring dplyr functions (filter, select, mutate, etc.).

---

## Data Import/Export

### import() - Read files with automatic format detection

```r
# Basic usage - format auto-detected from extension
df <- import("data.csv")              # CSV (lazy arrow Dataset)
df <- import("data.tsv")              # TSV (lazy arrow Dataset)
df <- import("data.parquet")          # Parquet (lazy arrow Dataset)
df <- import("data.xlsx")             # Excel
df <- import("data.json")             # JSON
df <- import("data.yaml")             # YAML
df <- import("data.rds")              # RDS
df <- import("data.qs")               # QS2 format
df <- import("data.fst")              # FST format

# Eager loading (load all data into memory)
df <- import("data.csv", lazy = FALSE)  # Returns data.table
df <- import("data.parquet", lazy = FALSE)  # Returns arrow Table

# Lazy loading for large files (memory-efficient)
ds <- import("large_data.csv", lazy = TRUE)  # Returns arrow Dataset
# Work with dplyr verbs, then collect
result <- ds |> filter(col > 100) |> collect()

# URLs are supported
df <- import("https://example.com/data.csv")

# Compressed files auto-detected
df <- import("data.csv.gz")
df <- import("data.tsv.zst")
```

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
export(df, "output.qs")
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

# Work with dplyr verbs (lazy evaluation)
result <- lazy_tbl |>
  filter(col > 100) |>
  group_by(category) |>
  summarise(total = sum(value)) |>
  collect()

# Export directly from DuckDB (more efficient)
tbl_export(conn, lazy_tbl, "output.parquet")
tbl_export(conn, lazy_tbl, "output.csv")

# With partitioning
tbl_export(conn, lazy_tbl, "output_dir",
           format = "parquet",
           partition_by = c("year", "month"))

# List all tables/views
tbl_ls(conn)
tbl_ls(conn, show_colnames = TRUE)  # Include column names

# Drop a table or view
tbl_drop(conn, "my_view", type = "view")
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

**Features:**
- Real-time progress bar with percentage and ETA
- Milestone alerts at 50% and 90%
- Error isolation (one failure doesn't stop others)
- Automatic fallback to sequential on Windows
- Completion summary with success/failure count

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

```r
# For ggplot2 labels (returns expression)
label <- human_read_latex_pval("0.05")
# Returns: expression(italic(P) == 0.05)

# With statistic prefix
label <- human_read_latex_pval("1e-5", s = "R = 0.95")
# Returns: expression("R = 0.95, " * italic(P) == 1 %*% 10^{-5})

# As character string
label <- human_read_latex_pval("1e-5", tex = FALSE)
# Returns: "P=1 × 10^{-5}"

# Use in ggplot2
ggplot(data, aes(x, y)) +
  geom_point() +
  annotate("text", x = 1, y = 1,
           label = human_read_latex_pval("0.001"),
           parse = TRUE)
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

# Import data as views
tbl_import(conn, "sales", "sales.parquet")
tbl_import(conn, "products", "products.csv")

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
ggplot(data, aes(x, y)) +
  geom_point() +
  geom_smooth(method = "lm") +
  labs(
    title = "Correlation Analysis",
    subtitle = human_read_latex_pval(
      format(cor_result$p.value, scientific = TRUE),
      s = paste0("R = ", round(cor_result$estimate, 2))
    )
  )
```

---

## Checklist

When using jutils in R code:

- [ ] Use `import()` instead of `read.csv()`, `fread()`, `read_parquet()`, etc.
- [ ] Use `export()` instead of `write.csv()`, `fwrite()`, `write_parquet()`, etc.
- [ ] Use `lazy = TRUE` for large files to avoid loading all data into memory
- [ ] Use `db_conn()` for DuckDB connections (manages global connections)
- [ ] Use `tbl_import()` to create lazy views from files
- [ ] Use `pbmclapply()` instead of `mclapply()` for progress bars
- [ ] Use `human_read_latex_pval()` for publication-quality p-value formatting
- [ ] Call `db_disconn()` when done with database work
- [ ] Remember that `library(jutils)` auto-loads common packages
