# 02-cell-inclusion.R

## Purpose

Criteria C1 and C2. Quantify what original mgatk's cell filter removes, and
what scMOCHA's stricter per-cell read requirement costs.

- **02a** per-cell MT coverage, cells split by the mgatk `mean coverage > 10` filter
- **02b** cells entering variant calling under each caller
- **02c** lowest heteroplasmy a single cell can support under each confident-detection rule
- **02d** cell-level detections and variants lost with the discarded cells

Panel 02c is where the cost of scMOCHA's extra requirement is stated: mgatk
needs 4 alt reads in a cell, scMOCHA needs 10, so scMOCHA's per-cell floor is
`10 / depth` against mgatk's `4 / depth`.

## Inputs

All under `${ISILON_BASE}/compare-with-mgatk/derived/<sample_id>`:

- `01-variant-joined.qs`
- `01-cell-depth.qs`
- `01-cell-af-s1.qs`
- `01-cell-pos-coverage.qs`

## Outputs

- `newplots/compare-with-mgatk/figures/<sample_id>/02a-cell-coverage.pdf`
- `newplots/compare-with-mgatk/figures/<sample_id>/02b-cells-retained.pdf`
- `newplots/compare-with-mgatk/figures/<sample_id>/02c-min-detectable-af.pdf`
- `newplots/compare-with-mgatk/figures/<sample_id>/02d-detections-lost.pdf`
- `newplots/compare-with-mgatk/tables/<sample_id>/02-cell-inclusion.tsv`

The table is `metric` / `value` in long form with `sample` as its first
column. Step 07 casts the five copies of it into the cross-sample cell-filter
figure and table.

## Run

```bash
cd "$(git rev-parse --show-toplevel)"
pixi run Rscript newplots/compare-with-mgatk/02-cell-inclusion.R --sample=GSE181279_GSM5494116_5PPE
```

`--sample` is required and must be one of the `sample_id` values in `SAMPLES`
in `config.R`.

## Verify

```bash
sample=GSE181279_GSM5494116_5PPE
cat newplots/compare-with-mgatk/tables/"${sample}"/02-cell-inclusion.tsv
ls -la newplots/compare-with-mgatk/figures/"${sample}"/02*.pdf
```

Expected **for GSE181279_GSM5494116_5PPE**: 7210 cells total, 6724 kept by
mgatk, 486 dropped; 127,998 cell-level detections of which 4,130 sit in the
discarded cells; 14 variants fall below the 10-cell reliability threshold when
those cells are removed.

Cells total and cells dropped by mgatk across the eight samples, in registry
order: 721 / 389, 5805 / 4381, 7949 / 7673, 191 / 21, 7210 / 486,
17919 / 16099, 5639 / 5018, 8645 / 7155. Most of the shallow 3' samples lose
the bulk of their cells to the coverage filter, which is why this step is per
sample and its eight tables are compared in step 07.

## Notes

`detection` here uses the scMOCHA gate definition (cell AF >= 0.05 and depth
>= 10) because the question is what the cell filter costs the scMOCHA gate.
Any panel about heteroplasmy *level* uses the uncensored measure instead; see
`config.R.md`.

`kept_by_mgatk` comes from mgatk's own AF matrix, not from recomputing
`mean coverage > 10` on `cell.depthTable.txt`; the recomputed rule is carried
as `meancov_over_cutoff` and the two can disagree by a couple of cells. See
`01-load-harmonize.md`.

Every panel in this step is defined on cells rather than on arms, so nothing
here depends on an arm being non-empty. The empty-arm handling starts at step
03.
