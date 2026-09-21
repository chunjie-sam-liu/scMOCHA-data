# config.R

## Purpose

Stage constants, path construction and shared helpers for the comparison
between scMOCHA-updated mgatk and original mgatk variant calling. Sourced by
every step in this stage. Defines no colours; those come from the track colour
file at `high-res/00-colors.R`, reached through `stage_paths()$colorfile`.

## Inputs

- `.env` at the repository root, read by `jutils::dotenv()` before this file is
  sourced. `REPODIR` and `HIGHRESDIR` carry a leading `~`, so every env-derived
  path is passed through `fs::path_expand()`.

## Outputs

None. Sourcing this file defines objects only.

## Contents

### Identity

- `STAGE` - directory name of this stage.
- `SAMPLES` - the sample registry, and the single source of truth for which
  samples exist: `sample_id`, `archive`, `gse`, `gsm`, `chemistry`. `sample_id`
  normalises the `-` that two archive names carry between GSE and GSM, so one
  token is safe as a path, a factor level and a file name.
- `SAMPLE_IDS` - `SAMPLES$sample_id`, the accepted values of `--sample`.
- `CROSS_SAMPLE` - the output leaf used by anything that spans samples.
- `SAMPLE_LABEL` - set per run by each step from `fn_sample_label()`. Every
  figure subtitle reads it through `fn_sample_note()`.
- `ARCHIVE_MEMBERS` - the seven files extracted from each archive. The same
  list appears in `00-extract-archives.sh`; keep them in step.

### Paths

`stage_paths(sample_id)` returns `repodir`, `root`, `indir`, `cachedir`,
`stagedir`, `figroot`, `tabroot`, `figdir`, `tabdir`, `colorfile`. Passing
`NULL`, the default, selects the `cross-sample` leaf instead of a sample.
`stage_inputs()` returns the seven input files under a given `paths$indir`.

### Criterion constants

Original mgatk, from `newplots/thecode/original-mgatk-variant-calling.py` and
the canonical mgatk/Signac downstream gate: `CUTOFF_CELL_MEANCOV`,
`CUTOFF_ALT_STRAND`, `CUTOFF_NCELLS_CONF`, `CUTOFF_VMR_MGATK`,
`CUTOFF_STRAND_MGATK`.

scMOCHA, from `newplots/thecode/scmocha-mgatk-variant-calling.py`:
`CUTOFF_ALT_READS`.

scMOCHA reliability gate, copied verbatim from
`src/06.1-collect-variants-new.R`: `CUTOFF_HETEROPLASMIC`,
`CUTOFF_HOMOPLASMIC`, `CUTOFF_MIN_READS`, `CUTOFF_NOTRELIABLE`,
`POS_RNA_EDITING`, `POS_MISSALIGNMENT_ERROR`. That script is the production
classification; if it changes, these must change with it.

`CUTOFF_MIN_ALT_OBS` belongs to this stage alone. It is the read-support floor
for the uncensored AF measure and appears in neither caller.

### Helpers

- `fn_theme()` - the shared ggplot2 theme.
- `fn_af_bin()` - AF binning, with an explicit `No detected cell` level so a
  variant with no carrier cell is shown rather than silently dropped.
- `fn_sample_note()` - the sample clause used in every figure subtitle.
- `fn_check_sample()` - validates `--sample`, stopping with the valid list.
- `fn_sample_label()` - the `GSE GSM (chemistry)` caption for one sample. It
  looks the row up with `match()` rather than a data.table `i` expression,
  because `i` cannot see a function argument through `..name`.
- `fn_empty_panel()` / `fn_or_empty()` - the placeholder drawn when an arm or a
  group is empty. `fn_or_empty()` takes the real plot as a promise, so a panel
  that would fail on an empty group is never evaluated. Seven of the ten
  samples retain no original-mgatk variant at all, and a zero there is a result
  that has to be drawn rather than skipped.
- `fn_testable()` - TRUE when a two-group comparison has at least
  `CUTOFF_MIN_GROUP` variants on both sides. Below that the medians are still
  reported and the p-value is `NA`.
- `fn_export_tab()` - writes a per-sample table with `sample` as its first
  column, so the cross-sample step can bind the five files without re-deriving
  provenance.
- `fn_is_blacklisted()` - position blacklist test.
- `fn_detection_long()` - melts the cell-by-variant AF matrix, attaches the
  cell's depth at the variant position, and flags scMOCHA detections. A
  `(position, barcode)` pair absent from the coverage file means no reads, so
  its depth is set to 0 rather than left `NA`.
- `fn_scmocha_gate()` - the scMOCHA reliability gate.
- `fn_background_rate()` - the background alt-read rate. A shallow sample can
  carry no background alt read, or no qualifying background cell at all, which
  makes the estimate 0 or undefined. Either value breaks the carrier test in
  opposite directions - rate 0 calls every cell holding one alt read a carrier,
  rate 1 calls nothing a carrier - so it is floored at one alt read over all
  observed depth, with a warning. Shallow samples with no qualifying background
  cell hit that floor - `GSE155673_GSM4712895_3PV3` and
  `GSE220189_GSM6793474_3PV3` did on the 2026-09-18 run; the deep GSE181279
  does not. Check the step 03 log for the `not estimable` warning rather than
  assuming.
- `fn_carrier_stats()` - carrier-level heteroplasmy on three definitions.

## The three heteroplasmy measures

This is the one thing to understand before reading any figure in this stage.
Every one of them is a "median AF across the cells carrying the variant"; they
differ only in which cells count as carriers, and that choice changes the
answer.

| Column | Carrier definition | Use for |
| --- | --- | --- |
| `mean` (from either caller) | all cells | nothing here |
| `af_gate_*` | cell AF >= 0.05 and depth >= 10 | describing the scMOCHA gate itself |
| `af_loose_*` | depth >= 10 and >= 2 alt reads | the sensitivity table only |
| `af_carrier_*` | alt count inconsistent with the background error rate | **everything else** |

**`mean`** is total alt reads over total coverage across every cell. A variant
at 60% heteroplasmy in 12 of 7,210 cells scores below 0.002, so the column
ranks variants by prevalence, not by heteroplasmy level.

**`af_gate_*`** is censored at 0.05, because 0.05 is the AF floor inside its own
carrier definition. Using it for an AF distribution deletes the region this
stage exists to examine, and it does so silently: the analysis runs cleanly and
returns a null result.

**`af_loose_*`** removes the AF floor but replaces it with a floor so weak that
background reads dominate. At the measured background rate of 0.026%, two alt
reads in a cell of depth 1,000 is expected by chance, so a variant with 19
genuine carriers above 5% and 876 background cells gets a median of 0.0017.
That is the background rate, not the variant's heteroplasmy.

**`af_carrier_*`** keeps a cell when its alt count is inconsistent with the
background error rate under a binomial test, Bonferroni-corrected across all
9,005,290 cell-by-variant observations. The threshold scales with depth, so it
excludes background cells without imposing an AF floor. This is the adopted
measure.

`fn_background_rate()` estimates the background from cells that are nowhere
near carrying the variant, **including the zero-alt cells**; dropping them
inflates the estimate by more than twofold (0.026% to 0.059%).

Panel `04f` reports the gate test under all four definitions so the choice is
auditable. The direction is the same in every one; only the censored measure
fails to resolve it. See `D12` and `D17` in `DECISION.md` for the two wrong
turns that produced this section.
