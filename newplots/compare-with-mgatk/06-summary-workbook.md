# 06-summary-workbook.R

## Purpose

Collect every table this stage produced into one styled workbook for the
response letter, and write the criterion comparison as a machine-readable
table. Adds nothing new: it reads the `.tsv` files the earlier steps wrote and
builds the criteria and definitions sheets from the constants in `config.R`.

A step whose table is missing is skipped with a warning rather than failing, so
the workbook can be rebuilt after re-running only part of the stage.

## Inputs

From `newplots/compare-with-mgatk/tables/`:

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

## Outputs

- `newplots/compare-with-mgatk/tables/06-compare-with-mgatk.xlsx` - 15 sheets
- `newplots/compare-with-mgatk/tables/06-criteria-comparison.tsv`

Sheets: `00_Criteria`, `01_Definitions`, `02_Cell_inclusion`,
`03_Funnel_counts`, `04_Gate_crossapplied`, `05_Exclusion_reasons`,
`06_Arm_combinations`, `07_AF_bins`, `08_AF_bins_own_definition`,
`09_Gate_test`, `10_Measure_sensitivity`, `11_VMR_strand_rejected`,
`12_Arm_in_mgatk_plane`, `13_Read_support`, `14_Variant_membership`.

## Run

```bash
cd "$(git rev-parse --show-toplevel)"
pixi run Rscript newplots/compare-with-mgatk/06-summary-workbook.R
```

Run it last; it reads what the other steps wrote.

## Verify

```bash
pixi run Rscript -e 'cat(openxlsx2::wb_load("newplots/compare-with-mgatk/tables/06-compare-with-mgatk.xlsx")$get_sheet_names(), sep = "\n")'
cat newplots/compare-with-mgatk/tables/06-criteria-comparison.tsv
```

Expected: 15 sheet names in the order above, and 7 criterion rows C1 to C7.

## Notes

Every criterion string in the `00_Criteria` sheet is built with `glue()` from
the `config.R` constants, so a threshold cannot be described in the workbook
differently from how the code applied it.
