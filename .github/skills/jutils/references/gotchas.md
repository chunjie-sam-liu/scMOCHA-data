# Gotchas — behaviour that looks correct but is not

Every item here was verified by running the code against the package in this
repository. These are the failure modes that produce **wrong files, wrong
types, or silently wrong results** rather than an obvious error.

---

## 1. `export()` to `.qs2` writes a `.qs` file

```r
export(df, "out.qs2")   # actually creates out.qs
import("out.qs2")       # ERROR: File does not exist
```

Both the `qs2` and `qs` branches call `build_path("qs")`, so the extension is
always `.qs`. The _content_ is written by `qs2::qs_save()`.

**Rule:** always name the file `.qs`. Never write `.qs2` in a path.

```r
export(df, "out.qs")    # correct: file name and content agree
import("out.qs")
```

---

## 2. CSV/TSV compression: only `.gz` works

`export()` maps `gz|bz2|xz|zst|zstd` to a `compress=` argument, but
`data.table::fwrite()` only accepts `"auto"`, `"none"`, `"gzip"`. Everything
except gzip aborts. `.zip` fails even earlier, at format detection.

```r
export(df, "out.csv.gz")    # OK
export(df, "out.csv.zst")   # ERROR: 'arg' should be one of "auto","none","gzip"
export(df, "out.tsv.bz2")   # ERROR (same)
export(df, "out.csv.zip")   # ERROR: Cannot detect format from file extension
```

With an **Arrow** input even `.gz` fails, because `compression` is not a valid
argument for Arrow's CSV writer:

```r
export(arrow::arrow_table(df), "out.csv.gz")
# ERROR: `compression` is not a valid argument for your chosen `format`.
```

`import()` is asymmetric — it only strips `.gz`, `.zip`, `.tar.gz` when
detecting format, and of those only gzip is actually readable:

```r
import("data.csv.gz")       # OK (lazy and eager both work)
import("data.tsv.zst")      # ERROR: Unsupported format: zst
```

**Rule:** with `import()`/`export()`, `.gz` on a data.frame is the only
working combination. For zstd, use `tbl_export()` (DuckDB COPY supports `.gz`
and `.zst`) or write plain `.parquet`, which is zstd-compressed internally by
default.

---

## 3. Eager `import()` return types are not all data.table

| Call                                | Actual class                |
| ----------------------------------- | --------------------------- |
| `import("x.csv", lazy = FALSE)`     | `data.table` `data.frame`   |
| `import("x.tsv", lazy = FALSE)`     | `data.table` `data.frame`   |
| `import("x.parquet", lazy = FALSE)` | `tbl_df` `tbl` `data.frame` |
| `import("x.fst")`                   | `data.frame`                |
| `import("x.xlsx")`                  | `tbl_df` `tbl` `data.frame` |

Only csv/tsv give you a data.table. Applying `dt[, .(x)]` or `:=` to a
parquet/fst/xlsx import fails.

**Rule:** call `data.table::setDT()` explicitly when you need data.table
semantics from parquet, fst, or xlsx.

---

## 4. `db_conn()` ignores `read_only` for an already-pooled path

The pool key is the normalised path only. `read_only` is **not** part of the
key and is **not** stored.

```r
c1 <- db_conn("x.duckdb", read_only = FALSE)
c2 <- db_conn("x.duckdb", read_only = TRUE)   # returns c1 — still read-write
c3 <- db_conn("x.duckdb", read_only = FALSE)
identical(c1, c3)                              # TRUE
```

Requesting a different mode silently returns the existing connection. The
usual symptom is a write failing with a read-only error long after the
`db_conn()` call.

**Rule:** decide the mode once per process. To switch, close first:

```r
db_disconn()
conn <- db_conn("x.duckdb", read_only = FALSE)
```

Also note the default is `read_only = TRUE`.

---

## 5. `human_read()` errors on `NA`

```r
human_read(c(1, NA))   # ERROR: missing value where TRUE/FALSE needed
```

P-value and correlation vectors routinely contain `NA`.

**Rule:** filter or guard first.

```r
out <- rep(NA_character_, length(x))
ok <- !is.na(x)
out[ok] <- human_read(x[ok])
```

---

## 6. `pbmclapply()` never throws — failures hide in the result

Failed elements come back as `try-error` objects inside the result list. The
call itself succeeds. A completion summary counts them, but nothing aborts.

**Rule:** check explicitly, or return a sentinel from inside the worker.

```r
res <- pbmclapply(items, process_fn, mc.cores = 8)
bad <- vapply(res, \(x) inherits(x, "try-error"), logical(1))
if (any(bad)) cli::cli_abort("{sum(bad)} job{?s} failed")
```

---

## 7. Progress bars disappear when output is redirected

The progress bar is chosen by `isatty(stdout())`. Under the tmux + log
redirection workflow (`... > logs/job.log 2>&1`) there is no TTY, so you get
only two milestone lines at 50% and 90% plus the final summary.

That is expected, not a hang. Use `CPU_USED` vs `RUN_TIME` to tell a slow job
from a stalled one.

---

## 8. `saveplot()` swallows render errors and writes blank pages

A plot that fails to `print()` is replaced by a blank page and a warning; the
file is still written and `saveplot()` returns normally. A `NULL` element in a
plot list does the same.

**Rule:** treat warnings from `saveplot()` as errors when the figure matters.
Numbered output is zero-padded to the width of the count:

```r
saveplot("fig.png", plots)   # 2 plots  -> fig_1.png, fig_2.png
                             # 10 plots -> fig_01.png ... fig_10.png
```

`saveplot()` returns the original `filename`, not the numbered file names.

---

## 9. `load_pkg()` is completely silent under `Rscript`

Every message is wrapped in `if (interactive())`. In a script, a missing
package produces no output at all — `load_pkg()` just returns `FALSE` for it
and execution continues until the first call to the missing function.

This also applies to `library(jutils)` itself: 13 of the 20 auto-attached
packages are in `Suggests`, so a missing `ggplot2` or `dbplyr` is not
reported when a script starts.

**Rule:** for a hard dependency in a script, use
`rlang::check_installed("pkg")` instead of relying on `load_pkg()`.

---

## 10. `dotenv()` overwrites existing environment variables by default

`override = TRUE` is the default. Additionally `override = FALSE` tests
`nzchar(Sys.getenv(key))`, so a variable that exists but is empty counts as
absent and gets overwritten anyway.

**Rule:** pass `override = FALSE` when the shell environment should win.

---

## 11. Small API traps

| Call                                    | Trap                                                                   |
| --------------------------------------- | ---------------------------------------------------------------------- |
| `tbl_ls(conn, TRUE, 1)`                 | `...` must be empty (`check_dots_empty()`); errors on any extra arg    |
| `tbl_drop(...)` in a script             | aborts unless `confirm = TRUE`; `confirm` must be exactly `TRUE`       |
| `tbl_import(conn, n, "x.parquet", ...)` | parquet branch ignores `header`, `delim`, `columns`, `types`, `...`    |
| `tbl_import()` on `.txt`                | unknown extensions fall back to `csv`, they do not error               |
| `human_read_latex_pval(x, tex = FALSE)` | returns a `glue` object, not plain character                           |
| `human_read_latex_pval(x, s = vec)`     | `s` must be length 1; longer vectors error                             |
| `fn_xy_breaks_limits()`                 | list order is `limits, breaks, labels, step`; always index by name     |
| `conn_ls()`                             | returns **invisibly**; assign it if you need the value                 |
| `tbl_export()`                          | returns an `fs_path`, not a plain string                               |
| `tbl_export(..., "x.csv.bz2")`          | aborts; CSV/TSV accept `.gz` and `.zst` only                           |
| `export(..., create.dir = FALSE)`       | no validation — fails later inside the writer, unlike `saveplot()`     |
| `convert(a, b, ...)`                    | `...` goes to **both** `import()` and `export()`; one-sided args error |
| `import()` on a URL                     | silently forces `lazy = FALSE`, without warning                        |
