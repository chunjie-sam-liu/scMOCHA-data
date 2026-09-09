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
- `SAMPLE_LABEL` - the GSE/GSM/SRR identifier of the compared sample, currently
  `<pending>`. Every figure subtitle reads it through `fn_sample_note()`, so
  one edit here relabels the whole stage.

### Paths

`stage_paths()` returns `repodir`, `indir`, `cachedir`, `stagedir`, `figdir`,
`tabdir`, `colorfile`. `stage_inputs()` returns the eight input files for a
given sample prefix, defaulting to `cell`.

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
- `fn_is_blacklisted()` - position blacklist test.
- `fn_detection_long()` - melts the cell-by-variant AF matrix, attaches the
  cell's depth at the variant position, and flags scMOCHA detections. A
  `(position, barcode)` pair absent from the coverage file means no reads, so
  its depth is set to 0 rather than left `NA`.
- `fn_scmocha_gate()` - the scMOCHA reliability gate.
- `fn_carrier_stats()` - carrier-level heteroplasmy on two definitions.

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
