# Parallel processing

## pbmclapply() — parallel lapply with progress bar

```r
pbmclapply(X, FUN, ...,
           mc.cores = getOption("mc.cores", 8L),
           mc.preschedule = TRUE,
           mc.set.seed = TRUE,
           mc.cleanup = TRUE,
           mc.allow.recursive = TRUE,
           .progress = TRUE)
```

Wraps `parallel::mclapply()` with a cli progress bar. Falls back to
sequential `lapply()` on Windows or when `mc.cores = 1`.

```r
# Basic (8 cores default)
results <- pbmclapply(1:100, function(x) {
  Sys.sleep(0.1)
  x^2
})

# Custom cores
results <- pbmclapply(items, process_fn, mc.cores = 4)

# No progress bar
results <- pbmclapply(items, process_fn, .progress = FALSE)

# Error handling with tryCatch
results <- pbmclapply(1:100, function(i) {
  tryCatch(
    {
      if (i == 50) stop("oops")
      i^2
    },
    error = function(e) NA
  )
}, mc.cores = 4)
```

**mc.preschedule:** `TRUE` (default) pre-assigns jobs to cores (faster for
many small jobs, but a failure affects the whole batch on that core). `FALSE`
isolates errors to individual jobs but has higher fork overhead.

**Failures never abort.** Failed elements come back as `try-error` objects
inside the result list and the call returns normally. Check explicitly:

```r
res <- pbmclapply(items, process_fn, mc.cores = 8)
bad <- vapply(res, \(x) inherits(x, "try-error"), logical(1))
if (any(bad)) cli::cli_abort("{sum(bad)} job{?s} failed")
```

**The progress bar needs a TTY.** The display is selected by
`isatty(stdout())`. Under the tmux + log redirection workflow
(`... > logs/job.log 2>&1`) there is no TTY, so you get only milestone lines
at 50% and 90% plus a final summary. That is expected, not a hang.

**`.progress = FALSE` bypasses the wrapper entirely** and calls
`parallel::mclapply()` directly — no progress file, no monitor process, no
completion summary.

---

## pbmcmapply() — parallel mapply with progress bar

```r
pbmcmapply(FUN, ...,
           MoreArgs = NULL,
           mc.cores = getOption("mc.cores", 8L),
           mc.preschedule = TRUE,
           mc.set.seed = TRUE,
           mc.cleanup = TRUE,
           SIMPLIFY = TRUE,
           USE.NAMES = TRUE,
           .progress = TRUE)
```

Wraps `parallel::mcmapply()` with a cli progress bar. Note there is no
`mc.allow.recursive` argument here, unlike `pbmclapply()`. The same
`try-error`, TTY, and `.progress = FALSE` notes above apply.

```r
# Multiple vectorized arguments
results <- pbmcmapply(function(x, y) x + y, 1:100, 101:200, mc.cores = 8)

# With fixed arguments
results <- pbmcmapply(
  process_fn, files, params,
  MoreArgs = list(verbose = TRUE),
  mc.cores = 4
)
```
