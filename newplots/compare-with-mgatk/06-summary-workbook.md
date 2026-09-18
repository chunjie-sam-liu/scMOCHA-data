# 06-summary-workbook.R

## Purpose

Collect every table this stage produced, **for every sample**, into one styled
workbook for the response letter, and write the criterion comparison as a
machine-readable table. Adds nothing new: it reads the `.tsv` files the
earlier steps wrote and builds the samples, criteria and definitions sheets
from the registry and constants in `config.R`.

One workbook with a `sample` column on every topic sheet, not five workbooks
of the same sheets. The five per-sample copies of each table are bound into a
single sheet, and the step 07 tables come in as they are.

A table that is missing is skipped with a warning rather than failing, so the
workbook can be rebuilt after re-running only part of the stage.

## Inputs

Per sample, from `newplots/compare-with-mgatk/tables/<sample_id>/`, for each
of the five `sample_id` values in `SAMPLES`:

- `02-cell-inclusion.tsv`
- `03-funnel-counts.tsv`
- `03-gate-crossapplied.tsv`
- `03-exclusion-reasons.tsv`
- `03-arm-combinations.tsv`
- `03-variant-membership.tsv`
- `04-af-bins.tsv`
- `04-af-bins-own-definition.tsv`
- `04-tests.tsv`
- `04-measure-sensitivity.tsv`
- `05-vmr-strand-rejected.tsv`
- `05-arm-in-mgatk-plane.tsv`
- `05-read-support.tsv`

Cross-sample, from `newplots/compare-with-mgatk/tables/cross-sample/`, all
written by step 07:

- `07-sample-overview.tsv`
- `07-arm-yield.tsv`
- `07-cell-filter.tsv`
- `07-strand-support.tsv`
- `07-exclusion-reasons.tsv`
- `07-gate-test.tsv`

## Outputs

- `newplots/compare-with-mgatk/tables/cross-sample/06-compare-with-mgatk.xlsx` - 22 sheets
- `newplots/compare-with-mgatk/tables/cross-sample/06-criteria-comparison.tsv`

Sheets, in order: `00_Samples`, `01_Criteria`, `02_Definitions`,
`03_Overview`, `04_Arm_yield`, `05_Cell_filter`, `06_Strand_correlation`,
`07_Exclusion_pooled`, `08_Gate_test_per_sample`, `09_Cell_inclusion`,
`10_Funnel_counts`, `11_Gate_crossapplied`, `12_Exclusion_reasons`,
`13_Arm_combinations`, `14_AF_bins`, `15_AF_bins_own_definition`,
`16_Gate_test`, `17_Measure_sensitivity`, `18_VMR_strand_rejected`,
`19_Arm_in_mgatk_plane`, `20_Read_support`, `21_Variant_membership`.

Sheets `00` to `02` describe the stage, `03` to `08` are the cross-sample
layer from step 07, and `09` to `21` are the per-sample tables bound across
samples.

## Run

```bash
cd "$(git rev-parse --show-toplevel)"
pixi run Rscript newplots/compare-with-mgatk/06-summary-workbook.R
```

No `--sample`: it spans them. **Run it last, and after step 07**, because it
binds the tables steps 01 to 05 wrote for every sample plus the six step 07
wrote. It has no dependency check, and a stale workbook looks identical to a
fresh one.

## Verify

```bash
pixi run Rscript -e 'cat(openxlsx2::wb_load("newplots/compare-with-mgatk/tables/cross-sample/06-compare-with-mgatk.xlsx")$get_sheet_names(), sep = "\n")'
cat newplots/compare-with-mgatk/tables/cross-sample/06-criteria-comparison.tsv
```

Expected: 22 sheet names in the order above, and 7 criterion rows C1 to C7.
The log line at the end names every sheet written, so a sheet dropped for a
missing table is visible there.

## Notes

Every criterion string in the `01_Criteria` sheet is built with `glue()` from
the `config.R` constants, so a threshold cannot be described in the workbook
differently from how the code applied it.

**A header-only table is dropped, not bound.** A `.tsv` with no rows reads
back with every column typed logical, which collides with the other samples on
`rbindlist`. Those samples contribute no rows anyway and their zero is already
recorded in the count tables, so they are logged and skipped.

**The bind uses `ignore.attr = TRUE`.** A column that happens to be
whole-numbered in a shallow sample and fractional in a deep one comes back
with a different class from each file, which would otherwise stop the bind on
a difference that carries no meaning.
