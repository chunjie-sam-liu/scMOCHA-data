# Data import/export

## import() — read any file

```r
import(file, format = NULL, lazy = TRUE, ...)
```

Auto-detects format from extension. Lazy loading (arrow Dataset) is default
for csv, tsv, parquet. All other formats always load eagerly.

| Extension  | Lazy? | Eager return type | Backend                          |
| ---------- | ----- | ----------------- | -------------------------------- |
| .csv       | Yes   | data.table        | fread (nThread=8) / open_dataset |
| .tsv       | Yes   | data.table        | fread (nThread=8) / open_dataset |
| .parquet   | Yes   | arrow Table       | read_parquet / open_dataset      |
| .xlsx      | No    | data.table        | readxl::read_xlsx                |
| .json      | No    | list/data.frame   | jsonlite::fromJSON               |
| .yaml      | No    | list              | yaml::read_yaml                  |
| .rds       | No    | R object          | readr::read_rds                  |
| .qs / .qs2 | No    | R object          | qs2::qs_read (nthreads=8)        |
| .fst       | No    | data.table        | fst::read_fst                    |
| .env       | No    | named list        | dotenv()                         |

```r
# Lazy (default for csv/tsv/parquet) — returns arrow Dataset
ds <- import("data.csv")
ds <- import("data.parquet")
result <- ds |> filter(x > 10) |> collect()

# Eager — returns data.table / arrow Table
df <- import("data.csv", lazy = FALSE)
df <- import("data.parquet", lazy = FALSE)

# Other formats (always eager, lazy param ignored)
df <- import("data.xlsx")
df <- import("data.rds")
df <- import("data.qs")
df <- import("data.fst")
df <- import("data.json")
df <- import("data.yaml")

# Compressed files auto-detected
df <- import("data.csv.gz")
df <- import("data.tsv.zst")

# URLs supported (forces eager)
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
export(df, "out.qs")                   # saves as qs2 (qs is deprecated)
export(df, "out.fst")

# Compressed output
export(df, "out.csv.gz")
export(df, "out.tsv.zst")

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

**qs format is deprecated.** Always use `.qs2` extension for new files.

---

## convert() — convert between formats

```r
convert(filein, fileout, lazy_in = TRUE, lazy_out = TRUE, ...)
```

Streams data lazily when both formats support it (csv, tsv, parquet).

```r
convert("data.csv", "data.parquet")                     # lazy → lazy
convert("data.parquet", "data.csv")                     # lazy → lazy
convert("data.csv", "data.rds")                         # lazy → eager (materializes)
convert("large.csv", "large.parquet", lazy_in = TRUE, lazy_out = TRUE)
```
