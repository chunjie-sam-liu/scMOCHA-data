# Data import/export

## import() — read any file

```r
import(file, format = NULL, lazy = TRUE, ...)
```

Auto-detects format from extension. Lazy loading (arrow Dataset) is default
for csv, tsv, parquet. All other formats always load eagerly.

Eager return types are **not** uniform — only csv/tsv give a data.table.

| Extension | Lazy? | Eager return type      | Backend                          |
| --------- | ----- | ---------------------- | -------------------------------- |
| .csv      | Yes   | `data.table`           | fread (nThread=8) / open_dataset |
| .tsv      | Yes   | `data.table`           | fread (nThread=8) / open_dataset |
| .parquet  | Yes   | `tbl_df` (tibble)      | read_parquet / open_dataset      |
| .xlsx     | No    | `tbl_df` (tibble)      | readxl::read_xlsx                |
| .json     | No    | list/data.frame        | jsonlite::fromJSON               |
| .yaml     | No    | list                   | yaml::read_yaml                  |
| .rds      | No    | R object               | readr::read_rds                  |
| .qs       | No    | R object               | qs2::qs_read (nthreads=8)        |
| .fst      | No    | `data.frame`           | fst::read_fst                    |
| .env      | No    | named list (sets env!) | dotenv()                         |

Call `data.table::setDT()` when you need data.table semantics from parquet,
fst, or xlsx.

```r
# Lazy (default for csv/tsv/parquet) — returns arrow Dataset
ds <- import("data.csv")
ds <- import("data.parquet")
result <- ds |> filter(x > 10) |> collect()

# Eager
df <- import("data.csv", lazy = FALSE)     # data.table
df <- import("data.parquet", lazy = FALSE) # tibble, NOT a data.table

# Other formats (always eager, lazy param ignored)
df <- import("data.xlsx")
df <- import("data.rds")
df <- import("data.qs")
df <- import("data.fst")
df <- import("data.json")
df <- import("data.yaml")

# Compression: ONLY .gz, .zip, .tar.gz are stripped during format detection
df <- import("data.csv.gz")   # OK, lazy and eager both work
# import("data.tsv.zst")      # ERROR: Unsupported format: zst

# URLs supported — silently forces lazy = FALSE, without warning
df <- import("https://example.com/data.csv")

# Override format explicitly
df <- import("data.txt", format = "csv")

# Extra args passed to backend
df <- import("data.csv", lazy = FALSE, select = c("id", "value"))
```

---

## export() — write any file

```r
export(x, file, format = NULL, lazy = TRUE, create.dir = TRUE, ...)
```

Auto-detects format from extension. `format` accepts a single format string
or a character vector of multiple formats. Creates output directory if needed.

```r
# Basic — format from extension
export(df, "out.csv")
export(df, "out.tsv")
export(df, "out.parquet")              # dataset folder (lazy=TRUE default)
export(df, "out.parquet", lazy = FALSE) # single file
export(df, "out.xlsx")
export(df, "out.json")
export(df, "out.yaml")
export(df, "out.rds")
export(df, "out.qs")                   # write .qs — see warning below
export(df, "out.fst")

# Compressed output — .gz ONLY, and only for csv/tsv
export(df, "out.csv.gz")
# export(df, "out.tsv.zst")            # ERROR: 'arg' should be one of
                                       # "auto", "none", "gzip"

# Multiple formats at once (uses same base filename)
export(df, "out", format = c("csv", "fst"))
export(df, "out", format = c("qs2", "tsv", "xlsx"))

# Legacy shorthand for c("csv", "fst")
export(df, "out", format = "both")

# Named list of data.frames → multi-sheet Excel
export(list(cars = mtcars, flowers = iris), "out.xlsx")
```

**Excel auto-styling:** bold headers, `#,##0.00` number format with color
gradient, categorical columns get HCL fill, text left-aligned. Errors if
data > 1,048,575 rows or > 16,384 columns.

**Always name qs files `.qs`, never `.qs2`.** The `qs2` engine is used either
way, but `export(df, "out.qs2")` writes a file called **`out.qs`**, so a
later `import("out.qs2")` fails with "File does not exist".

**Compression support is narrow and asymmetric:**

| Suffix | `export()` csv/tsv | `import()`           | `tbl_export()` csv/tsv |
| ------ | ------------------ | -------------------- | ---------------------- |
| `.gz`  | yes                | yes                  | yes                    |
| `.zst` | error              | error                | yes                    |
| `.bz2` | error              | error                | error                  |
| `.xz`  | error              | error                | error                  |
| `.zip` | error              | stripped, unreadable | error                  |

`export()` compression only applies on the data.frame path (`fwrite`). With an
Arrow input, even `.csv.gz` errors with "`compression` is not a valid argument
for your chosen `format`" — `collect()` first, or write parquet.

For zstd output use `tbl_export()`, or write `.parquet` — parquet is written
with zstd compression by default. Compression suffixes on non-csv/tsv formats
(`out.rds.gz`) are accepted but **not** applied.

**`create.dir = FALSE` is not validated.** Unlike `saveplot()`, `export()`
does not check that the directory exists; it fails inside the writer instead.

---

## convert() — convert between formats

```r
convert(filein, fileout, lazy_in = TRUE, lazy_out = TRUE, ...)
```

Streams data lazily when both formats support it (csv, tsv, parquet).
Supports every format `import()` and `export()` support, including qs and env
as input.

```r
convert("data.csv", "data.parquet")                     # lazy → lazy
convert("data.parquet", "data.csv")                     # lazy → lazy
convert("data.csv", "data.rds")                         # lazy → eager (materializes)
convert("large.csv", "large.parquet", lazy_in = TRUE, lazy_out = TRUE)
```

**`...` is forwarded to both `import()` and `export()`.** An argument that is
valid for only one side makes the other side error. Pass no `...` unless the
argument is meaningful to both.
