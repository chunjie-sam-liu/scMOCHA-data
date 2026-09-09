# PROGRESS: scMOCHA-updated mgatk vs original mgatk variant calling

Stage: `newplots/compare-with-mgatk/`
Last updated: 2026-09-09

---

## State

All seven steps implemented, syntax-checked, run to completion, and verified.
22 figures and 15 tables produced, plus `README.md` and `DIAGRAM.md`. One
sample, three arms.

| Step | Status | Last run |
| --- | --- | --- |
| `config.R` | done | - |
| `01-load-harmonize.R` | done | 2026-09-09 16:36 |
| `02-cell-inclusion.R` | done | 2026-09-09 16:37 |
| `03-variant-funnel-overlap.R` | done | 2026-09-09 16:55 |
| `04-heteroplasmy-spectrum.R` | done | 2026-09-09 16:55 |
| `05-vmr-strand.R` | done | 2026-09-09 16:55 |
| `06-summary-workbook.R` | done | 2026-09-09 16:55 |

Last full pass ran clean: no warnings, no errors, all sources ASCII.

## Results

Sample: 7,210 cells, median per-cell MT coverage 38.9. Identity pending (Q1).

**Cell inclusion, C1.** Original mgatk discards 486 of 7,210 cells (6.74%) for
mean coverage <= 10. Those cells hold 4,130 of 127,998 cell-level detections,
and removing them drops 14 variants below the 10-cell reliability threshold.

**Three-arm funnel.**

| Stage | Original mgatk | scMOCHA variant call | scMOCHA AF>5% |
| --- | --- | --- | --- |
| S0 candidates | 25,580 | 25,746 | 25,746 |
| S1, n_cells_conf >= 3 | 1,247 | 736 | 736 |
| S2, reliability gate | 230 | none applied | 216 |

Arm 1 against arm 3: 45 shared, 185 mgatk only, 171 AF>5% only.

**Why each arm misses the other's variants.**

| Direction | Cutoff responsible | n |
| --- | --- | --- |
| in AF>5%, not in mgatk (171) | strand r <= 0.65 alone | **149 (87%)** |
| | vmr <= 0.01 and strand r <= 0.65 | 17 |
| | vmr <= 0.01 alone | 5 |
| in mgatk, not in AF>5% (185) | < 10 cells at AF >= 0.05 | 117 (63%) |
| | < 3 cells with fwd+rev alt >= 10 | 68 (37%) |

Across all 1,247 mgatk S1 variants the gate rejects 628 on strand correlation
alone, 319 on both, and 70 on VMR alone.

**Gates cross-applied.** Holding the S1 set fixed, the scMOCHA gate retains
more variants than mgatk's in both directions: 269 vs 230 on mgatk's S1, and
216 vs 164 on scMOCHA's S1.

**The editor's question.** Among variants passing the scMOCHA reliability
gate, mgatk's `vmr > 0.01 & strand r > 0.65` gate rejects 215 and keeps 54.
The rejected variants sit at median carrier heteroplasmy **0.069** against
**0.205** for those it keeps (Wilcoxon P = 0.0029). Carriers are cells whose
alt count is inconsistent with the measured 0.026% background error rate,
Bonferroni-corrected over 9,005,290 cell-by-variant observations.

The direction holds under every carrier definition tried
(`04-measure-sensitivity.tsv`): P = 0.59 with an AF >= 0.05 floor, which cannot
resolve it because it has no values below 0.05; 0.043 with a 2-alt-read floor,
which is diluted by background; 0.0029 with the error model; 0.0036 on the
carrier maximum.

**The counter-direction, C2.** scMOCHA's extra 10-alt-read requirement is why
its S1 set is smaller. The variants only mgatk reports carry a median of
**9** alt reads per carrier cell, against **48** for shared variants and
**38** for scMOCHA-only variants.

## Open

- **Q1, sample identity.** `SAMPLE_LABEL` in `config.R` is `<pending>`. Every
  figure subtitle reads it; one edit relabels the stage. Needed before the
  figures go into a response letter.
- **One sample only.** The editor's criticism is about the whole resource.
  Extending to more samples needs paired `*_mgatk_original.tsv.gz` outputs,
  which currently exist for this sample alone.
- **No orthogonal truth set.** Nothing here is a sensitivity or recall
  measurement, and the figures must not be described as one.
- **P = 0.044 on one sample.** The gate test is nominally significant but not
  robust to a single sample. Treat it as a demonstration of direction, not as
  an established effect size, until more samples are added.

## Error log

| Date | Symptom | Cause | Fix |
| --- | --- | --- | --- |
| 2026-09-09 | `stopifnot` failed on barcode count in step 01 | The AF matrices have an empty first header field, so `fread(header = TRUE)` consumed the first barcode as a column name and returned 7,209 of 7,210 cells | Read with `header = FALSE`; assert with `setequal()` on barcode sets rather than on row counts |
| 2026-09-09 | Gate test returned P = 0.82, then P = 0.59, both null | The AF measure was censored: first the population `mean` column (prevalence, not heteroplasmy), then a carrier definition with a 0.05 AF floor, which deletes exactly the region under test | Added an uncensored measure with a read-support floor. Result became 0.019 vs 0.053, P = 0.043. See D12 |
| 2026-09-09 | scMOCHA AF>5% appeared to have 29% of variants at carrier AF < 0.01, for a set gated on 10 cells above 0.05 | The D12 fix over-corrected: a fixed 2-alt-read floor does not scale with depth, so at a 0.026% background rate two alt reads at depth 1,000 counts as a carrier and the median reports the background | Carriers now defined by a binomial test against the measured background, Bonferroni-corrected. Result 0.069 vs 0.205, P = 0.0029; the sub-0.01 fraction fell to 6 of 216. See D17 |
| 2026-09-09 | `log-10 transformation introduced infinite values` on every run | Zero-count bins on a log-scaled count axis in panel 04d | Plot density instead of counts. See D14 |

## Resume

Everything is current. To rebuild from scratch:

```bash
cd "$(git rev-parse --show-toplevel)"
for s in 01-load-harmonize 02-cell-inclusion 03-variant-funnel-overlap \
         04-heteroplasmy-spectrum 05-vmr-strand 06-summary-workbook; do
  pixi run Rscript "newplots/compare-with-mgatk/${s}.R" || break
done
```

Steps 02 to 06 depend on step 01's cache; step 06 must run last. Read
`DECISION.md` before changing any AF measure.
