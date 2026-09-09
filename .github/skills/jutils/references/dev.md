# Developing the jutils package itself

This applies when editing the jutils source in this repository, not when
using jutils from an analysis script.

## Environment

**Always use pixi. Never use conda or mamba.**

```bash
pixi run <task>            # run a defined task
pixi run Rscript -e "..."  # ad-hoc R command
pixi shell                 # interactive shell in the environment
```

Do not run `conda activate renv`, `mamba activate renv`, or
`conda run -n renv ...`. The `.vscode/tasks.json` entries that wrap
`conda run -n renv task ...` are stale; use the pixi tasks below instead.

## Tasks

Defined in `pixi.toml`. Run from the repository root.

| Task                  | Command                                                  |
| --------------------- | -------------------------------------------------------- |
| `pixi run doc`        | `devtools::document()` — regenerate NAMESPACE and `man/` |
| `pixi run test`       | `devtools::test()`                                       |
| `pixi run dev`        | `doc` then `test`                                        |
| `pixi run check`      | `dev` then `devtools::check(error_on = 'error')`         |
| `pixi run check-cran` | `devtools::check(clean = TRUE, args = '--as-cran')`      |
| `pixi run build`      | `R CMD build . --no-build-vignettes`                     |
| `pixi run install`    | `build` then `R CMD INSTALL`                             |
| `pixi run format`     | `air format .`                                           |
| `pixi run clean`      | remove `*.tar.gz`, `*.Rcheck`                            |

`air` is **not** a pixi dependency. `pixi run format` only works because the
Posit Air VS Code extension puts `air` on `PATH`. Outside that environment,
install `air` separately or the task fails.

Ad-hoc commands:

```bash
pixi run Rscript -e "devtools::load_all()"
pixi run Rscript -e "devtools::test(filter = 'dotenv', reporter = 'llm')"
pixi run Rscript -e "testthat::test_file('tests/testthat/test-dotenv.R')"
pixi run Rscript -e "lintr::lint_package()"
```

## Source layout

| File                 | Exports                                                               |
| -------------------- | --------------------------------------------------------------------- |
| `R/io.R`             | `import` `export` `convert`                                           |
| `R/db.R`             | `db_*` `tbl_*` `conn_ls` + 4 deprecated                               |
| `R/parallel.R`       | `pbmclapply` `pbmcmapply` + 11 internal helpers                       |
| `R/plot.R`           | `fn_xy_breaks_limits` `human_read` `human_read_latex_pval` `saveplot` |
| `R/dotenv.R`         | `dotenv`                                                              |
| `R/load_pkg.R`       | `load_pkg`                                                            |
| `R/utils.R`          | none — internal Excel styling helpers only                            |
| `R/jutils-package.R` | `.onAttach` (20 auto-loaded packages, 12 conflict rules)              |

Tests mirror source names: `R/name.R` -> `tests/testthat/test-name.R`.
testthat edition 3.

## Conventions

- Base pipe `|>`, never `%>%`.
- `\() ...` for single-line anonymous functions, `function() {}` for
  multi-line.
- `snake_case`, 2-space indent. `air.toml` sets `line-width = 80`; `.lintr`
  allows up to 120 via `line_length_linter(120)`. Write to 80.
- Run `pixi run format` after any code change. Do not hand-format.
- Explicit namespace prefixes in package code (`cli::cli_abort()`), never
  `library()`.
- `cli::cli_abort()` for user-facing errors, `rlang::abort()` only when an
  error class is needed. Note `rlang::abort()` does **not** interpret cli
  inline markup — `{.file {x}}` is emitted literally.
- Roxygen with markdown for every export; keep `_pkgdown.yml` in sync.
- Run `pixi run doc` after changing roxygen. Never edit `NAMESPACE` by hand.
- Do not edit `NEWS.md`; maintainers manage releases.

## Known source-level defects

The user-facing symptoms are documented in [gotchas.md](gotchas.md). As bugs
in this repository, the root causes are:

- `export()` writes `.qs` when asked for `.qs2`; its `compression_map` offers
  suffixes `fwrite()` rejects; `import()` strips a different suffix set, so
  the two are not round-trip compatible.
- `db_conn()` omits `read_only` from the pool key.
- `human_read()` errors on `NA`.
- `tbl_export()` passes cli markup to `rlang::abort()`, which prints it
  literally.
- Roxygen drift: `load_pkg()` claims `mclapply` but uses `sapply()`;
  `dotenv()` documents `file = "~/.Renviron"` but defaults to `".env"`.
- Undeclared dependencies used at runtime: `qs`, `grid`, `tools`, `stats`,
  `grDevices`, `utils`. `lifecycle` is called by the deprecated `db.R`
  functions but is only in `Suggests`.
- `title_case_names()` in `R/utils.R` is dead code.
