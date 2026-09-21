# PROGRESS: extend the mgatk comparison from one sample to five

Stage: `newplots/compare-with-mgatk/`
Campaign: `2026-09-17-multi-sample-comparison`
Last updated: 2026-09-18

Files: `2026-09-17-multi-sample-comparison.PLAN.md`,
`.DECISION.md`, this file.

---

## State

**Done.** Plan approved 2026-09-17 with all recommendations accepted;
implemented, run to completion, and verified.

| Item | Status |
| --- | --- |
| Input verification | done, 2026-09-17 |
| `PLAN.md`, `Q1`..`Q8` | approved; `Q8` still open, see below |
| `DECISION.md` | `M1`..`M6` |
| `high-res/00-colors.R`, `color_sample` | done |
| `config.R` registry, paths, helpers | done |
| `00-extract-archives.sh`, `run-all.sh` | done |
| `01`..`05` per sample | done, all ten samples |
| `07-cross-sample.R` | done |
| `06-summary-workbook.R` | done, 23 sheets |
| Paired `.md`, `AGENTS.md` | done |
| `README.md`, `DIAGRAM.md` | done, 2026-09-18 |
| `05f` / `05g`, scMOCHA decision plane | done, 2026-09-18 |
| Presentation pass: GSM labels, linear 07a, red retained, arm-specific planes | done, 2026-09-18, see `M8` |
| Three samples added (GSE155673, GSE188632, GSE220189) | done, 2026-09-18, see `M9` |
| Cross-sample figures restricted to the six 3' samples | done, 2026-09-18, see `M9` |
| Registry-driven extraction; 04f guard fixed | done, 2026-09-18, see `M10` |

Outputs: 22 figures and 13 tables per sample, 5 figures and 7 tables
cross-sample, one 22-sheet workbook. No zero-byte file. All sources ASCII.

## Results

### Yield, the headline

**Original mgatk retains zero variants in seven of the ten samples, and one in
two more.**

| Sample | Chemistry | Cells | Dropped by mgatk | S1 scMOCHA | S1 mgatk | mgatk final | scMOCHA call | scMOCHA AF>5% |
| --- | --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| GSE149689_GSM4509019_3PV3 | SC3Pv3 | 721 | 389 (54%) | 20 | 19 | **0** | 20 | 20 |
| GSE163314_GSM4976997_3PV2 | SC3Pv2 | 7,949 | 7,673 (97%) | 9 | 4 | **0** | 9 | 9 |
| GSE163668_GSM4995445_5PR2 | SC5P-R2 | 191 | 21 (11%) | 20 | 21 | **0** | 20 | 18 |
| GSE181279_GSM5494116_5PPE | SC5P-PE | 7,210 | 486 (6.7%) | 736 | 1,247 | 230 | 736 | 216 |
| GSE271107_GSM8369876_3PV3 | SC3Pv3 | 8,645 | 7,155 (83%) | 2 | 6 | **0** | 2 | 0 |

The cause is the strand-correlation floor, not missing data. Of the mgatk S1
variants, the number reaching `strand r > 0.65` is 0, 0, 0, 300 and 0; the
highest strand correlation observed anywhere in the four shallow samples is
**0.590**, and many values are negative.

### The editor's question

Computable in one sample only, because the other four have no comparison group:
mgatk's gate passes nothing there, so the split has a single level.

| Sample | Passes mgatk gate | Rejected | Median AF passed | Median AF rejected | P |
| --- | ---: | ---: | ---: | ---: | ---: |
| GSE149689 | 0 | 20 | - | 1.000 | not testable |
| GSE163314 | 0 | 11 | - | 1.000 | not testable |
| GSE163668 | 0 | 18 | - | 1.000 | not testable |
| GSE181279 | 54 | 215 | 0.205 | 0.069 | **0.0029** |
| GSE271107 | 0 | 0 | - | - | not testable |

GSE181279 reproduces its previously published values exactly, which is the
check that the refactor changed no behaviour: 7,210 cells, 486 dropped, S1
736 / 1,247, S2 230 / 216, 127,998 detections, 4,130 in dropped cells, 14
variants lost, 0.069 vs 0.205 at P = 0.0029, background rate 0.000259.

### What is not claimed

No pooled or combined test across samples (`M2`). No mechanism linking the
strand floor to coverage: Spearman rho is **-0.013, P = 0.66** in the deep
sample, so panel `07c` is descriptive only (`M4`).

## Open

- **`Q8`, provenance: answered in part.** The archives come from Ting (stated
  2026-09-18) and `DATA.md` records that. Still outstanding: the scMOCHA and
  mgatk versions used, and written confirmation that both arms ran on identical
  allele counts in all ten samples. Stage decision `D1` still rests on that
  assumption; schema and candidate counts are consistent with it, but
  consistency is not confirmation.
- **The four shallow samples cannot answer the editor's question directly.**
  They answer a different and simpler one: original mgatk returns nothing at
  all on them. Do not present the two as the same result.
- **Still no orthogonal truth set** (`D4`). Nothing here is a sensitivity or
  recall measurement.

## Error log

| Date | Symptom | Cause | Fix |
| --- | --- | --- | --- |
| 2026-09-17 | Step 01 died with `object '..sample_id' not found` | A data.table `i` expression cannot reach a function argument through `..name`; that prefix works only for column selection | Look the row up with `match()` outside the `[` |
| 2026-09-17 | Step 06 died on `rbindlist`: class attribute mismatch | A column whole-numbered in a shallow sample and fractional in a deep one reads back with a different class; a header-only file reads back all-logical | Drop empty parts, bind with `ignore.attr = TRUE` |
| 2026-09-17 | GSE163314 reported `background_alt_rate = 1` and vanished from the gate table | It has no qualifying background cell at all, so `0 / 0` floored to `1 / 1`, and a rate of 1 makes no cell a carrier | Floor at one alt read over all observed depth. See `M5` |
| 2026-09-17 | Planned panel 07c asserted a mechanism the data reject | rho = -0.013, P = 0.66 in the only sample large enough to test it | Redrew 07c as descriptive and renamed the file. See `M4` |

## Resume

Everything is current. To rebuild from scratch:

```bash
cd "$(git rev-parse --show-toplevel)"
bash newplots/compare-with-mgatk/run-all.sh
```

Read this file, then `.DECISION.md`, then `.PLAN.md`. Read the stage
`DECISION.md` before changing any AF measure.
