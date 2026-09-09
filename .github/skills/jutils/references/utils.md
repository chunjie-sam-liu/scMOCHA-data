# Utility functions

## dotenv() — load .env files

```r
dotenv(file = ".env", override = TRUE, verbose = FALSE, encoding = "UTF-8")
```

Parses .env files and sets system environment variables via `Sys.setenv()`.
**`override = TRUE` is the default**, so existing environment variables are
replaced unless you opt out.

```r
# Load default .env (overwrites existing vars)
dotenv()

# Custom file, don't override existing vars
dotenv("config/.env.production", override = FALSE)

# Verbose mode
dotenv(verbose = TRUE)

# Access loaded variables
db_host <- Sys.getenv("DB_HOST")
```

**Supported .env features:** comments (`#`), quoted values (single/double),
multiline values (`"""`/`'''`), variable expansion (`${VAR}`, `$VAR`),
escape sequences (`\n`, `\t`, `\\`), `export` prefix.

Notes:

- Malformed lines and invalid keys are skipped **silently** — nothing is
  logged, so a typo produces no error and no variable.
- A missing file emits a cli alert and returns an empty list. It is a plain
  message, not an R condition, so `tryCatch(warning = )` will not catch it.
- `override = FALSE` tests `nzchar(Sys.getenv(key))`, so a variable that
  exists but is empty is treated as absent and overwritten anyway.
- Inline comments require whitespace before `#`, so `KEY=secret#123` keeps
  `#123` as part of the value.
- Variable expansion runs on every value, including single-quoted and
  multiline ones.

---

## load_pkg() — load packages with CLI output

```r
load_pkg(..., verbose = TRUE)
```

Accepts unquoted names, strings, character vectors, or any mix. Returns an
invisible named logical vector (TRUE = loaded, FALSE = not installed).
Packages are attached to the search path, sequentially. `verbose` is accepted
but currently has no effect.

```r
load_pkg(ggplot2)
load_pkg(ggplot2, dplyr, tidyr)
load_pkg(c("ggplot2", "dplyr"))
load_pkg(ggplot2, c("dplyr", "tidyr"), "patchwork")
```

**Silent under `Rscript`.** Every message is wrapped in `if (interactive())`,
and a missing package is not an error — it just returns `FALSE`. In a script
you get no indication until the first call to a missing function.

For a hard dependency in a script, check explicitly:

```r
rlang::check_installed("ComplexHeatmap")
load_pkg(ComplexHeatmap)
```
