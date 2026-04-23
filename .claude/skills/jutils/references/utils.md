# Utility functions

## dotenv() — load .env files

```r
dotenv(file = ".env", override = TRUE, verbose = FALSE, encoding = "UTF-8")
```

Parses .env files and sets system environment variables via `Sys.setenv()`.

```r
# Load default .env
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

---

## load_pkg() — load packages with CLI output

```r
load_pkg(..., verbose = TRUE)
```

Accepts unquoted names, strings, character vectors, or any mix. Returns an
invisible named logical vector (TRUE = loaded, FALSE = not installed).

```r
load_pkg(ggplot2)
load_pkg(ggplot2, dplyr, tidyr)
load_pkg(c("ggplot2", "dplyr"))
load_pkg(ggplot2, c("dplyr", "tidyr"), "patchwork")
```
