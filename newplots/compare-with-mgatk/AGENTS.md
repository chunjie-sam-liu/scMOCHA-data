# AGENTS.md - newplots/compare-with-mgatk

Orientation for working in this stage. Run state lives in `PROGRESS.md`,
design in `PLAN.md`, rationale in `DECISION.md`, and the reader's guide to the
figures and results in `README.md`. When this file disagrees with
`PROGRESS.md`, the progress file wins.

Active set: `PLAN.md` + `PROGRESS.md` + `DECISION.md` at this directory root.
No dated campaign files yet.

## What this stage is

A criterion-by-criterion comparison of **three calling arms** on one sample
where the callers ran on the same allele counts. It exists to answer the Cell
Metabolism editorial concern recorded in `EDITOR.md`.

| Arm | What it is | Variants |
| --- | --- | --- |
| Original mgatk | cell filter + VMR/strand gate | 230 |
| scMOCHA variant call | what the caller emits, no AF filter | 736 |
| scMOCHA AF>5% | plus the downstream gate in `src/06.1-collect-variants-new.R` | 216 |

`DIAGRAM.md` holds the Mermaid version of all three arms with their cutoffs and
counts.

## Run

Everything runs through pixi from the repository root. Conda is not installed
on this host.

```bash
cd "$(git rev-parse --show-toplevel)"
pixi run Rscript newplots/compare-with-mgatk/01-load-harmonize.R
pixi run Rscript newplots/compare-with-mgatk/02-cell-inclusion.R
pixi run Rscript newplots/compare-with-mgatk/03-variant-funnel-overlap.R
pixi run Rscript newplots/compare-with-mgatk/04-heteroplasmy-spectrum.R
pixi run Rscript newplots/compare-with-mgatk/05-vmr-strand.R
pixi run Rscript newplots/compare-with-mgatk/06-summary-workbook.R
```

Order matters. Steps 02 to 05 read step 01's cache, step 04 and 05 read step
03's cache, and step 06 reads the `.tsv` files the others wrote. No script
takes arguments. The whole stage runs in about a minute.

## Layout

- `config.R` - constants, paths, shared helpers. No colours.
- `NN-*.R` + `NN-*.md` - one step per file, each with its contract.
- `README.md` - reader's guide: how to read every figure and table.
- `DIAGRAM.md` - Mermaid comparison of the three arms.
- `figures/` - 22 PDFs, `02a` through `05e`. Rewritten on every run.
- `tables/` - 15 files including `06-compare-with-mgatk.xlsx`.

## Paths

| What | Where | Variable |
| --- | --- | --- |
| Inputs, read only | `${ISILON_BASE}/compare-with-mgatk` | `stage_paths()$indir` |
| Cache | `${ISILON_BASE}/compare-with-mgatk/derived` | `stage_paths()$cachedir` |
| Figures | `newplots/compare-with-mgatk/figures` | `stage_paths()$figdir` |
| Tables | `newplots/compare-with-mgatk/tables` | `stage_paths()$tabdir` |
| Colours | `high-res/00-colors.R` | `stage_paths()$colorfile` |

Nothing is written into the input root except the `derived/` subdirectory.

## Pitfalls

- **The AF matrices have an empty first header field.** Reading the barcode
  column with `fread(header = TRUE)` consumes the first barcode as a column
  name and silently loses one cell. Use `header = FALSE`, or `skip = 1` with
  explicit `col.names` for the wide read.
- **Barcode order differs** between `cell.depthTable.txt` and the AF matrices.
  Compare barcode sets with `setequal()`, never row counts.
- **Never use the mgatk `mean` column as heteroplasmy.** It is total alt reads
  over total coverage across every cell, so it measures prevalence. A variant
  at 60% heteroplasmy in 12 of 7,210 cells scores below 0.002.
- **Never define carriers by a fixed AF floor or a fixed read floor.** An AF
  floor censors the region under study (`af_gate_*`, null result); a fixed read
  floor lets background cells dominate the median (`af_loose_*`, background
  rate reported as the variant's heteroplasmy). Use `af_carrier_*`, the
  binomial test against the measured background rate. This cost three analysis
  passes; see `D12` and `D17` in `DECISION.md`.
- **Estimate the background rate including the zero-alt cells.** Dropping them
  inflates it from 0.026% to 0.059%.
- **`04g` bars do not share a measure.** Each arm is summarised on its own AF
  definition by design (`D19`), so the AF>5% bar cannot fall below 0.05 while
  the other two can. Use `04c` for any cross-arm statement.
- **Save with `device = cairo_pdf`.** Subtitles contain `\u2265` and `\u00b7`;
  the default `pdf()` device substitutes ASCII and only warns.
- **Wrap subtitles at save time** with `fn_wrap_labs(p)`. Subtitles carry N and
  every threshold, so they run past the panel width and ggplot2 does not wrap
  them.
- **Step 06 must run after any step whose table it reads.** It has no
  dependency check; a stale workbook looks identical to a fresh one.
- **Sources are ASCII.** Write `\u00b7`, not the character. Check with
  `grep -nP '[^\x00-\x7F]' newplots/compare-with-mgatk/*.R`.
- **`export()` writes qs2 but wants a `.qs` name.** `export(x, "f.qs2")` writes
  `f.qs` and breaks the round trip.
- `cut` emits fields in file order regardless of the order requested.

## Criteria under comparison

Audited from `newplots/thecode/original-mgatk-variant-calling.py`,
`newplots/thecode/scmocha-mgatk-variant-calling.py`, and
`src/06.1-collect-variants-new.R`.

| | Original mgatk | scMOCHA |
| --- | --- | --- |
| Cell inclusion | mean MT cov > 10 | none |
| Per-cell confident call | fwd >= 2 and rev >= 2 | plus fwd + rev >= 10 |
| Variant retention | n_cells_conf >= 3 | identical |
| Reliability gate | vmr > 0.01 and strand r > 0.65 | none in the call; the AF>5% arm adds blacklist plus >= 10 cells at AF >= 0.05 and depth >= 10 |

`src/06.1-collect-variants-new.R` is the production classification and applies
**no** VMR or strand-correlation filter. The 0.3 cutoff in
`preprocessing/high_quality_variant/02-update-high-quality-variant.py` is a
side branch feeding a summary flextable, not the calling path; do not present
it as a scMOCHA criterion.

## Known drift

- `SAMPLE_LABEL` in `config.R` is `<pending>`. Figure subtitles say so.
- The four scMOCHA gate constants and the two position blacklists in `config.R`
  are copied from `src/06.1-collect-variants-new.R`. If that script changes,
  they must be updated here by hand.
