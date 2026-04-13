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

Wraps `parallel::mcmapply()` with a cli progress bar.

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
