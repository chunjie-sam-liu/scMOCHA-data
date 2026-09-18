# 01-load-harmonize.R

## Purpose

Join the two `variant_stats` tables on `variant`, tag every variant with the
stage it survives under each rule set, and extract the two cell-level tables
the later steps need: per-cell depth at candidate variant positions, and the
cell-by-variant AF submatrix restricted to the S1 union.

Nothing downstream reads the raw inputs again.

One sample per run, selected with `--sample`. Every path below is per sample,
so the five samples never share a cache file.

## Inputs

All under `${ISILON_BASE}/compare-with-mgatk/samples/<sample_id>`, written by
`00-extract-archives.sh`:

- `cell.variant_stats.tsv.gz` - scMOCHA per-variant statistics
- `cell.variant_stats_mgatk_original.tsv.gz` - original mgatk per-variant statistics
- `cell.cell_heteroplasmic_df.tsv.gz` - scMOCHA cell-by-variant AF, post-C3
- `cell.cell_heteroplasmic_df_mgatk_original.tsv.gz` - mgatk cell-by-variant AF, post-C3
- `cell.cell_heteroplasmic_df_raw.tsv.gz` - scMOCHA cell-by-variant AF, all candidates
- `cell.depthTable.txt` - per-cell mean MT coverage
- `cell.coverage.txt.gz` - per-cell per-position coverage, `pos,barcode,cov`

## Outputs

All under `${ISILON_BASE}/compare-with-mgatk/derived/<sample_id>`:

- `01-variant-joined.qs` - one row per candidate variant, both callers' statistics side by side, plus `s1_*`, `s2_mgatk`, `candidate_*`, `blacklisted`
- `01-cell-depth.qs` - per-cell mean coverage, `kept_by_mgatk` and `meancov_over_cutoff`
- `01-cell-pos-coverage.qs` - per-cell depth at every S1-union variant position
- `01-cell-af-s1.qs` - cell-by-variant AF for the S1-union variants

## Run

```bash
cd "$(git rev-parse --show-toplevel)"
pixi run Rscript newplots/compare-with-mgatk/01-load-harmonize.R --sample=GSE181279_GSM5494116_5PPE
```

`--sample` is required and must be one of the `sample_id` values in `SAMPLES`
in `config.R`; anything else stops with the valid list. Takes roughly half a
minute per sample, most of it streaming `cell.coverage.txt.gz`.

## Verify

```bash
sample=GSE181279_GSM5494116_5PPE
ls -la ~/project/scmocha/compare-with-mgatk/derived/"${sample}"/
pixi run Rscript -e 'suppressMessages(library(jutils)); d <- import("~/project/scmocha/compare-with-mgatk/derived/GSE181279_GSM5494116_5PPE/01-variant-joined.qs"); cat(nrow(d), sum(d$s1_scmocha), sum(d$s1_mgatk), sum(d$s2_mgatk), "\n")'
pixi run Rscript -e 'suppressMessages(library(jutils)); a <- import("~/project/scmocha/compare-with-mgatk/derived/GSE181279_GSM5494116_5PPE/01-cell-af-s1.qs"); cat(nrow(a), ncol(a) - 1L, "\n")'
```

Expected **for GSE181279_GSM5494116_5PPE**: 25746 variants, 736 scMOCHA S1,
1247 mgatk S1, 230 mgatk S2; the AF submatrix is 7210 cells by 1249 variants.

Across the five samples the scMOCHA S0 counts are 1897, 7077, 2484, 25746 and
7540, the scMOCHA S1 counts 20, 9, 20, 736 and 2, and the mgatk S1 counts 19,
4, 21, 1247 and 6, in registry order.

## Notes

**`cell.depthTable.txt` rounds mean coverage to two decimals, so mgatk's cell
filter cannot be reproduced from it.** A cell written as `10.0` may sit on
either side of `mean coverage > 10`, and in GSE163314 two such cells do.
Recomputing the rule there keeps 274 cells against the 276 mgatk itself used.
The kept set is therefore taken as membership in the barcode set of
`cell.cell_heteroplasmic_df_mgatk_original.tsv.gz`, which is the only
unambiguous record of which cells mgatk used, and stored as `kept_by_mgatk`.
The recomputed rule is kept alongside it as `meancov_over_cutoff` so the
criterion is still reported and still checkable. Any disagreement between the
two is logged with the coverage values involved and asserted to be confined to
cells within 0.01 of the cutoff; the other four samples agree exactly.

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

Two conditions **warn rather than fail**, because the shallow samples meet
them legitimately. A candidate position with no row in `cell.coverage.txt.gz`
has depth 0 in every cell, which the detection helper already handles. An S1
variant absent from the raw AF matrix is one caller's candidate that the other
never wrote a column for; it is carried with no per-cell AF. Both are logged
with a count so a large number is still visible.
