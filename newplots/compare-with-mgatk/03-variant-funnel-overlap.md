# 03-variant-funnel-overlap.R

## Purpose

Criteria C3, C4 and C6 across **three calling arms**. Count variants at each
filtering stage, show how the three sets overlap, and attribute every
arm-specific variant to the single cutoff that excluded it from the other arm.

- **03a** variants surviving S0 -> S1 -> S2 in each arm
- **03b** three-way overlap of the arms
- **03c** each reliability gate applied to each S1 set
- **03d** why original mgatk misses the scMOCHA AF>5% variants
- **03e** why scMOCHA AF>5% misses the original mgatk variants
- **03f** every arm membership combination

This step also computes the carrier-level heteroplasmy statistics and the
exclusion attributions, and writes them into the cached variant table, so
steps 04 and 05 inherit both.

## The three arms

| | Arm 1, original mgatk | Arm 2, scMOCHA variant call | Arm 3, scMOCHA AF>5% |
| --- | --- | --- | --- |
| Cells | mean MT cov > 10 | all | all |
| Confident cell | fwd >= 2 and rev >= 2 | plus fwd + rev >= 10 | same as arm 2 |
| S1 | `n_cells_conf_detected >= 3` | same | same |
| S2 | `vmr > 0.01` and `strand r > 0.65` | **none** | blacklist, then >= 10 cells at AF >= 0.05 and depth >= 10 |

Arm 3 is a strict subset of arm 2. The AF>5% rule lives in
`src/06.1-collect-variants-new.R`, downstream of variant calling, which is why
it is a separate arm rather than part of the call.

## Exclusion attribution

Each excluded variant is assigned to the **first** criterion it fails, in that
arm's own order, so the bars in 03d and 03e sum to the set size and no variant
is double-counted. The orders are fixed in `fn_exclusion_mgatk()`,
`fn_exclusion_scmocha_af5()` and `fn_exclusion_scmocha_call()` in `config.R`.
A variant failing both VMR and strand correlation is reported as failing both,
not arbitrarily assigned to one.

## Inputs

- `${ISILON_BASE}/compare-with-mgatk/derived/01-variant-joined.qs`
- `${ISILON_BASE}/compare-with-mgatk/derived/01-cell-af-s1.qs`
- `${ISILON_BASE}/compare-with-mgatk/derived/01-cell-pos-coverage.qs`

## Outputs

- `newplots/compare-with-mgatk/figures/03a-variant-funnel.pdf`
- `newplots/compare-with-mgatk/figures/03b-variant-overlap.pdf`
- `newplots/compare-with-mgatk/figures/03c-gate-crossapplied.pdf`
- `newplots/compare-with-mgatk/figures/03d-excluded-from-mgatk.pdf`
- `newplots/compare-with-mgatk/figures/03e-excluded-from-scmocha-af5.pdf`
- `newplots/compare-with-mgatk/figures/03f-arm-combinations.pdf`
- `newplots/compare-with-mgatk/tables/03-funnel-counts.tsv`
- `newplots/compare-with-mgatk/tables/03-gate-crossapplied.tsv`
- `newplots/compare-with-mgatk/tables/03-exclusion-reasons.tsv`
- `newplots/compare-with-mgatk/tables/03-arm-combinations.tsv`
- `newplots/compare-with-mgatk/tables/03-variant-membership.tsv`
- `${ISILON_BASE}/compare-with-mgatk/derived/03-variant-gated.qs`

## Run

```bash
cd "$(git rev-parse --show-toplevel)"
pixi run Rscript newplots/compare-with-mgatk/03-variant-funnel-overlap.R
```

## Verify

```bash
cat newplots/compare-with-mgatk/tables/03-funnel-counts.tsv
cat newplots/compare-with-mgatk/tables/03-exclusion-reasons.tsv
cat newplots/compare-with-mgatk/tables/03-arm-combinations.tsv
```

Expected arms: mgatk 25,580 -> 1,247 -> 230; scMOCHA call 25,746 -> 736;
scMOCHA AF>5% 25,746 -> 736 -> 216. Arm 1 against arm 3: 45 shared, 185 mgatk
only, 171 AF>5% only. Of those 171, **149 (87%) fail strand correlation
alone**; of the 185, 117 fail the 10-cell rule and 68 fail the 10-alt-read
rule.

## Notes

The scMOCHA gate is evaluated on **all 7,210 cells for every arm**, so panel
03c isolates the gate. Folding the cell filter in as well would conflate C1
with C4/C6; C1 is quantified on its own in step 02.
