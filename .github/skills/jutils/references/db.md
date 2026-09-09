# DuckDB database functions

## db_conn() — get/create DuckDB connection

```r
db_conn(dbpath = NULL, read_only = TRUE, ...)
```

Manages a global connection pool in a package-level environment. Reuses
existing valid connections. Each parallel worker gets its own pool.
The default is **read-only**.

```r
conn <- db_conn()                                  # in-memory
conn <- db_conn("my.duckdb")                       # file-based, read-only
conn <- db_conn("my.duckdb", read_only = FALSE)    # read-write
```

**Trap: `read_only` is not part of the pool key and is not stored.** Once a
path is pooled, later calls return that same connection regardless of the
`read_only` you ask for — with no warning.

```r
c1 <- db_conn("x.duckdb", read_only = FALSE)
c2 <- db_conn("x.duckdb", read_only = TRUE)   # returns c1, still read-write
c3 <- db_conn("x.duckdb", read_only = FALSE)
identical(c1, c3)                              # TRUE
```

Decide the mode once per process. To switch modes, close first:

```r
db_disconn()
conn <- db_conn("x.duckdb", read_only = FALSE)
```

---

## db_disconn() — close all DuckDB connections

```r
db_disconn(shutdown = TRUE)
```

Closes all connections managed by `db_conn()` in the current process.

```r
db_disconn()
```

---

## conn_ls() — list active connections

```r
conn_ls()
```

Shows all active DuckDB connections with validity status, access mode, and
variable names pointing to each connection. Returns the keys **invisibly** —
assign the result if you need the value.

```r
conn_ls()
# ℹ Active connections: 2
# ✔ :memory: (valid) [read-write] [conn]
# ✔ /path/to/my.duckdb (valid) [read-only] [conn2]

keys <- conn_ls()
```

---

## db_attach() / db_detach() — multi-database queries

```r
db_attach(conn, dbpath, alias = NULL, read_only = TRUE)
db_detach(conn, alias)
```

Attach external DuckDB files for cross-database queries.

```r
conn <- db_conn()
db_attach(conn, "other.duckdb", alias = "other")
DBI::dbGetQuery(conn, "SELECT * FROM other.main.my_table LIMIT 10")
db_detach(conn, "other")
db_disconn()
```

---

## tbl_import() — create lazy view from file

```r
tbl_import(conn, view_name, file,
           format = c("auto", "csv", "tsv", "parquet"),  # match.arg -> "auto"
           header = TRUE,
           delim = NULL,        # auto: "," for csv, "\t" for tsv
           columns = NULL,      # full schema (all columns, advanced)
           types = NULL,        # partial type override (recommended)
           overwrite = FALSE,
           ...)
```

Creates a DuckDB view backed by a file and returns a lazy dbplyr tbl.
Data is NOT loaded into memory until `collect()`.

```r
conn <- db_conn("my.duckdb", read_only = FALSE)

# Auto-detect format from extension
tbl_sales <- tbl_import(conn, "sales", "sales.csv")
tbl_data  <- tbl_import(conn, "data", "data.parquet")

# Partial type override (recommended)
tbl_sales <- tbl_import(conn, "sales", "sales.csv",
                        types = list(date = "DATE", amount = "DOUBLE"))

# Full schema (advanced — must specify ALL columns)
tbl_sales <- tbl_import(conn, "sales", "sales.csv",
                        columns = list(id = "INTEGER", name = "VARCHAR",
                                       amount = "DOUBLE"))

# Overwrite existing view
tbl_sales <- tbl_import(conn, "sales", "new_sales.csv", overwrite = TRUE)

# Query lazily with dplyr
result <- tbl_sales |>
  filter(amount > 100) |>
  group_by(category) |>
  summarise(total = sum(amount)) |>
  collect()
```

**The parquet branch ignores `header`, `delim`, `columns`, `types` and `...`.**
Those only apply to csv/tsv.

**`format = "auto"` falls back to `csv` for any unrecognised extension** — it
does not error. Pass `format` explicitly for files like `.txt` or `.dat`.

---

## tbl_export() — export dbplyr table to file

```r
tbl_export(conn, tbl, file,
           format = c("auto", "csv", "tsv", "parquet"),  # match.arg -> "auto"
           header = TRUE,
           partition_by = NULL,  # parquet only
           overwrite = FALSE)
```

Uses DuckDB COPY command — more efficient than `collect()` + `export()`.
Supports partitioned parquet output. Returns the path **invisibly** as an
`fs_path` object, not a plain string.

```r
tbl <- dplyr::tbl(conn, "sales") |> filter(year == 2024)

# Auto-detect format
tbl_export(conn, tbl, "filtered.csv")
tbl_export(conn, tbl, "filtered.parquet")
tbl_export(conn, tbl, "filtered.tsv")

# Compressed output — .gz and .zst only; .bz2 and .zip abort
tbl_export(conn, tbl, "filtered.csv.gz")
tbl_export(conn, tbl, "filtered.csv.zst")

# Overwrite — physically deletes the existing file/directory first
tbl_export(conn, tbl, "data.csv", overwrite = TRUE)

# Partitioned parquet
tbl_all <- dplyr::tbl(conn, "sales")
tbl_export(conn, tbl_all, "by_year",
           format = "parquet", partition_by = "year")
tbl_export(conn, tbl_all, "partitioned",
           format = "parquet", partition_by = c("year", "month"))
```

This is the only jutils writer that supports zstd for CSV/TSV; `export()`
supports gzip only. Parquet output is always zstd. Like `tbl_import()`,
`format = "auto"` falls back to `csv` for unrecognised extensions.

---

## tbl_ls() — list tables and views

```r
tbl_ls(conn, show_colnames = FALSE, ...)
```

Returns a data.frame with columns: database_alias, schema_name, table_name,
table_type, is_insertable_into, nrows, ncols (and colnames if requested,
semicolon-separated). Includes tables from attached databases.

`...` is enforced empty by `rlang::check_dots_empty()` — passing any extra
argument errors. Use named arguments only.

```r
tbl_ls(conn)
tbl_ls(conn, show_colnames = TRUE)
# tbl_ls(conn, TRUE, 1)   # ERROR: `...` must be empty
```

`nrows` is computed with a per-table `SELECT COUNT(*)`, not read from
metadata. This is slow on databases with many or large tables.

---

## tbl_drop() — drop a table or view

```r
tbl_drop(conn, name, type = c("table", "view"),
         if_exists = FALSE, confirm = FALSE)
```

Requires interactive confirmation by default. In non-interactive sessions
(scripts, agents), you MUST pass `confirm = TRUE`. The check is `isTRUE()`,
so `confirm = 1` or `confirm = "yes"` do not count as confirmed.

```r
tbl_drop(conn, "my_view", type = "view", confirm = TRUE)
tbl_drop(conn, "my_table", type = "table", confirm = TRUE)
tbl_drop(conn, "maybe", type = "table", if_exists = TRUE, confirm = TRUE)
```

---

## tbl_analyze() — update table statistics

```r
tbl_analyze(conn, table_name)
```

Runs ANALYZE on a table to improve DuckDB query optimizer decisions.

```r
tbl_analyze(conn, "large_table")
```

---

## tbl_register_arrow() — register Arrow dataset as DuckDB view

```r
tbl_register_arrow(conn, name, ds, overwrite = FALSE)
```

Registers an Arrow object (Dataset, Table, RecordBatchReader) as a DuckDB
view and returns a lazy dplyr tbl.

```r
ds <- import("data.parquet", lazy = TRUE)
tbl_data <- tbl_register_arrow(conn, "data", ds)

result <- tbl_data |>
  filter(pval < 5e-8) |>
  collect()

# Overwrite existing
tbl_data <- tbl_register_arrow(conn, "data", ds, overwrite = TRUE)
```
