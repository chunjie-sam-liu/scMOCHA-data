# AGENTS.md - newplots/compare-with-mgatk

Orientation for working in this stage. Run state lives in `PROGRESS.md`,
design in `PLAN.md`, rationale in `DECISION.md`, and the reader's guide to the
figures and results in `README.md`. When this file disagrees with
`PROGRESS.md`, the progress file wins.

Active set: `2026-09-17-multi-sample-comparison.PLAN.md` + `.PROGRESS.md` +
`.DECISION.md`. The stage-level `PLAN.md` / `PROGRESS.md` / `DECISION.md`
remain the charter and are still read first on resume.

## What this stage is

A criterion-by-criterion comparison of **three calling arms** on **five
samples** where the callers ran on the same allele counts. It exists to answer
the Cell Metabolism editorial concern recorded in `EDITOR.md`.

| Arm | What it is |
| --- | --- |
| Original mgatk | cell filter + VMR/strand gate |
| scMOCHA variant call | what the caller emits, no AF filter |
| scMOCHA AF>5% | plus the downstream gate in `src/06.1-collect-variants-new.R` |

Variants retained, per sample:

| Sample | Chemistry | Original mgatk | scMOCHA call | scMOCHA AF>5% |
| --- | --- | ---: | ---: | ---: |
| GSE149689_GSM4509019_3PV3 | SC3Pv3 | **0** | 20 | 20 |
| GSE163314_GSM4976997_3PV2 | SC3Pv2 | **0** | 9 | 9 |
| GSE163668_GSM4995445_5PR2 | SC5P-R2 | **0** | 20 | 18 |
| GSE181279_GSM5494116_5PPE | SC5P-PE | 230 | 736 | 216 |
| GSE271107_GSM8369876_3PV3 | SC3Pv3 | **0** | 2 | 0 |

`DIAGRAM.md` holds the Mermaid version of all three arms with their cutoffs.

## Run

Everything runs through pixi from the repository root. Conda is not installed
on this host.

```bash
cd "$(git rev-parse --show-toplevel)"
bash newplots/compare-with-mgatk/run-all.sh
```

That is the whole stage: extract, five samples x five steps, then the
cross-sample step and the workbook. About four minutes after extraction.
`run-all.sh GSE181279_GSM5494116_5PPE` restricts it to named samples but still
rebuilds the cross-sample layer, which needs every sample's cache.

One sample, one step at a time:

```bash
bash newplots/compare-with-mgatk/00-extract-archives.sh
pixi run Rscript newplots/compare-with-mgatk/01-load-harmonize.R --sample=<id>
pixi run Rscript newplots/compare-with-mgatk/02-cell-inclusion.R --sample=<id>
pixi run Rscript newplots/compare-with-mgatk/03-variant-funnel-overlap.R --sample=<id>
pixi run Rscript newplots/compare-with-mgatk/04-heteroplasmy-spectrum.R --sample=<id>
pixi run Rscript newplots/compare-with-mgatk/05-vmr-strand.R --sample=<id>
pixi run Rscript newplots/compare-with-mgatk/07-cross-sample.R
pixi run Rscript newplots/compare-with-mgatk/06-summary-workbook.R
```

Order matters. Steps 02 to 05 read step 01's cache, 04 and 05 read step 03's
cache, **07 reads every sample's step 03 cache and step 02 table**, and **06
runs last** because it binds the tables 01 to 05 and 07 wrote. `<id>` is one of
the `sample_id` values in `SAMPLES` in `config.R`; anything else stops with the
valid list.

## Layout

- `config.R` - `SAMPLES` registry, constants, paths, shared helpers. No colours.
- `00-extract-archives.sh` - unpack the archives, idempotent.
- `NN-*.R` + `NN-*.md` - one step per file, each with its contract.
- `run-all.sh` - the driver.
- `README.md` - reader's guide: how to read every figure and table.
- `figures/<sample_id>/` - 24 PDFs each; `figures/cross-sample/` - 5 PDFs.
- `tables/<sample_id>/` - 14 TSVs each; `tables/cross-sample/` - the workbook
  plus 7 files.
- `figures_stale-2026-09-17/`, `tables_stale-2026-09-17/` - the single-sample
  outputs this campaign replaced. History; do not write there.

## Paths

| What | Where | Variable |
| --- | --- | --- |
| Archives | `${ISILON_BASE}/compare-with-mgatk/*.zip` | `SAMPLES$archive` |
| Inputs, read only | `${ISILON_BASE}/compare-with-mgatk/samples/<id>` | `stage_paths(<id>)$indir` |
| Cache | `${ISILON_BASE}/compare-with-mgatk/derived/<id>` | `stage_paths(<id>)$cachedir` |
| Figures | `newplots/compare-with-mgatk/figures/<id>` | `stage_paths(<id>)$figdir` |
| Tables | `newplots/compare-with-mgatk/tables/<id>` | `stage_paths(<id>)$tabdir` |
| Cross-sample | same roots, leaf `cross-sample` | `stage_paths()` |
| Colours | `high-res/00-colors.R` | `stage_paths()$colorfile` |

`stage_paths(NULL)` is the cross-sample leaf. Nothing is written into the input
root except `samples/` and `derived/`.

## Pitfalls

Schema and input traps:

- **The AF matrices have an empty first header field.** Reading the barcode
  column with `fread(header = TRUE)` consumes the first barcode as a column
  name and silently loses one cell. Use `header = FALSE`, or `skip = 1` with
  explicit `col.names` for the wide read.
- **Barcode order differs** between `cell.depthTable.txt` and the AF matrices.
  Compare barcode sets with `setequal()`, never row counts.
- **`cell.depthTable.txt` rounds mean coverage to two decimals, so mgatk's cell
  filter cannot be reproduced from it.** GSE163314 has four cells written as
  `10.0`, two of which mgatk kept and two it did not; recomputing
  `mean_coverage > 10` gives 274 kept against mgatk's own 276. `kept_by_mgatk`
  is therefore membership in the mgatk AF matrix barcode set, and the coverage
  rule is only cross-checked. See `M3`.
- **A mgatk S1 variant can be absent from the scMOCHA raw AF matrix.** Step 01
  warns and carries on with depth 0 rather than asserting.

Degenerate arms, which four of the five samples have:

- **Original mgatk retains zero variants in four samples.** That is the result,
  not a failure. Every panel that groups or tests must tolerate an empty group;
  `fn_or_empty()` and `fn_testable()` in `config.R` are how. See `M1`.
- **Never put a count on a log scale here.** A zero bar and its label both
  vanish silently, and the zero is the headline. Use
  `scales::transform_pseudo_log(base = 10)`.
- **A `by = sample` aggregation drops samples with no rows.** A missing row
  looks identical to a failed run, so cross-sample tables are merged back onto
  `SAMPLES` and carry an explicit zero.
- **The background alt rate can be 0 or undefined in a shallow sample.** Rate 0
  makes every cell with one alt read a carrier; the degenerate fallback rate 1
  makes nothing a carrier. `fn_background_rate()` floors it at one alt read
  over all observed depth and warns. See `M5`.
- **A header-only TSV reads back with every column typed logical**, which
  collides with the other samples on `rbindlist`. Step 06 drops empty parts and
  binds with `ignore.attr = TRUE`.

Measurement and plotting:

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
  inflates it from 0.026% to 0.059% in GSE181279.
- **`04g` bars do not share a measure.** Each arm is summarised on its own AF
  definition by design (`D19`), so the AF>5% bar cannot fall below 0.05 while
  the other two can. Use `04c` for any cross-arm statement.
- **Strand correlation does not track per-variant coverage within a sample.**
  Spearman rho is -0.013 (P = 0.66) in GSE181279. `07c` is descriptive - where
  the variants sit against the floor - and must not be captioned as a
  mechanism. See `M4`.
- **Save with `device = cairo_pdf`.** Subtitles contain `\u2265` and `\u00b7`;
  the default `pdf()` device substitutes ASCII and only warns.
- **Wrap subtitles at save time** with `fn_wrap_labs(p)`. Subtitles carry N and
  every threshold, so they run past the panel width and ggplot2 does not wrap
  them.
- **Sources are ASCII.** Write `\u00b7`, not the character. Check with
  `grep -nP '[^\x00-\x7F]' newplots/compare-with-mgatk/*.R`.

Mechanics:

- **`--sample` binds a global named `sample`, which shadows `base::sample`.**
  Every step calls `rm(sample)` right after validating it. Keep that line.
- **Step 06 must run after every step whose table it reads, 07 included.** It
  has no dependency check; a stale workbook looks identical to a fresh one.
- **`export()` writes qs2 but wants a `.qs` name.** `export(x, "f.qs2")` writes
  `f.qs` and breaks the round trip.
- `cut` emits fields in file order regardless of the order requested.
- A data.table `i` expression cannot see a function argument through `..name`;
  look the value up with `match()` outside the `[` instead.

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

- The four scMOCHA gate constants and the two position blacklists in `config.R`
  are copied from `src/06.1-collect-variants-new.R`. If that script changes,
  they must be updated here by hand.
- `ARCHIVE_MEMBERS` in `config.R` and the `members` array in
  `00-extract-archives.sh` are the same list in two languages.
- The provenance of the four archives added on 2026-09-11 is unknown; `DATA.md`
  at the repository root records them as such.
