# 01-load-harmonize.R

## Purpose

Join the two `variant_stats` tables on `variant`, tag every variant with the
stage it survives under each rule set, and extract the two cell-level tables
the later steps need: per-cell depth at candidate variant positions, and the
cell-by-variant AF submatrix restricted to the S1 union.

Nothing downstream reads the raw inputs again.

## Inputs

All under `${ISILON_BASE}/compare-with-mgatk`:

- `cell.variant_stats.tsv.gz` - scMOCHA per-variant statistics
- `cell.variant_stats_mgatk_original.tsv.gz` - original mgatk per-variant statistics
- `cell.cell_heteroplasmic_df.tsv.gz` - scMOCHA cell-by-variant AF, post-C3
- `cell.cell_heteroplasmic_df_mgatk_original.tsv.gz` - mgatk cell-by-variant AF, post-C3
- `cell.cell_heteroplasmic_df_raw.tsv.gz` - scMOCHA cell-by-variant AF, all candidates
- `cell.depthTable.txt` - per-cell mean MT coverage
- `cell.coverage.txt.gz` - per-cell per-position coverage, `pos,barcode,cov`

## Outputs

All under `${ISILON_BASE}/compare-with-mgatk/derived`:

- `01-variant-joined.qs` - one row per candidate variant, both callers' statistics side by side, plus `s1_*`, `s2_mgatk`, `candidate_*`, `blacklisted`
- `01-cell-depth.qs` - per-cell mean coverage and `kept_by_mgatk`
- `01-cell-pos-coverage.qs` - per-cell depth at every S1-union variant position
- `01-cell-af-s1.qs` - cell-by-variant AF for the S1-union variants

## Run

```bash
cd "$(git rev-parse --show-toplevel)"
pixi run Rscript newplots/compare-with-mgatk/01-load-harmonize.R
```

Takes roughly half a minute, most of it streaming `cell.coverage.txt.gz`.

## Verify

```bash
ls -la ~/project/scmocha/compare-with-mgatk/derived/
pixi run Rscript -e 'suppressMessages(library(jutils)); d <- import("~/project/scmocha/compare-with-mgatk/derived/01-variant-joined.qs"); cat(nrow(d), sum(d$s1_scmocha), sum(d$s1_mgatk), sum(d$s2_mgatk), "\n")'
pixi run Rscript -e 'suppressMessages(library(jutils)); a <- import("~/project/scmocha/compare-with-mgatk/derived/01-cell-af-s1.qs"); cat(nrow(a), ncol(a) - 1L, "\n")'
```

Expected: 25746 variants, 736 scMOCHA S1, 1247 mgatk S1, 230 mgatk S2; the AF
submatrix is 7210 cells by 1249 variants.

## Notes

Two schema traps in these files, both already handled here.

The AF matrices have an **empty first header field** - the barcode column has
no name. Reading the barcode column with `header = TRUE` makes `fread` consume
the first barcode as the column name and silently return one cell too few. The
script reads with `header = FALSE`, and for the wide submatrix uses
`skip = 1` with explicit `col.names`.

Barcode **order differs** between `cell.depthTable.txt` and the AF matrices, so
the assertions compare sets with `setequal()`, never row counts.

`cut` emits fields in file order regardless of the order requested, so the
column names for the AF submatrix are taken from the header in that same order
rather than from the requested variant vector.
