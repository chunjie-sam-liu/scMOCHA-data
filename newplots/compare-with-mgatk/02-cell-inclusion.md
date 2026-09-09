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

- `${ISILON_BASE}/compare-with-mgatk/derived/01-variant-joined.qs`
- `${ISILON_BASE}/compare-with-mgatk/derived/01-cell-depth.qs`
- `${ISILON_BASE}/compare-with-mgatk/derived/01-cell-af-s1.qs`
- `${ISILON_BASE}/compare-with-mgatk/derived/01-cell-pos-coverage.qs`

## Outputs

- `newplots/compare-with-mgatk/figures/02a-cell-coverage.pdf`
- `newplots/compare-with-mgatk/figures/02b-cells-retained.pdf`
- `newplots/compare-with-mgatk/figures/02c-min-detectable-af.pdf`
- `newplots/compare-with-mgatk/figures/02d-detections-lost.pdf`
- `newplots/compare-with-mgatk/tables/02-cell-inclusion.tsv`

## Run

```bash
cd "$(git rev-parse --show-toplevel)"
pixi run Rscript newplots/compare-with-mgatk/02-cell-inclusion.R
```

## Verify

```bash
cat newplots/compare-with-mgatk/tables/02-cell-inclusion.tsv
ls -la newplots/compare-with-mgatk/figures/02*.pdf
```

Expected: 7210 cells total, 6724 kept by mgatk, 486 dropped; 127,998 cell-level
detections of which 4,130 sit in the discarded cells; 14 variants fall below
the 10-cell reliability threshold when those cells are removed.

## Notes

`detection` here uses the scMOCHA gate definition (cell AF >= 0.05 and depth
>= 10) because the question is what the cell filter costs the scMOCHA gate.
Any panel about heteroplasmy *level* uses the uncensored measure instead; see
`config.R.md`.
