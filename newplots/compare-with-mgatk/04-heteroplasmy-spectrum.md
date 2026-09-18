# 04-heteroplasmy-spectrum.R

## Purpose

The editor's question, answered directly: is a filter biased toward
higher-heteroplasmy variants?

- **04a** ECDF of carrier heteroplasmy by variant class
- **04b** density of maximum carrier heteroplasmy by variant class
- **04c** heteroplasmy composition of each of the three arms
- **04d** cell-level heteroplasmy of the observations each set carries
- **04e** heteroplasmy of the variants mgatk's VMR/strand gate discards, with a
  Wilcoxon test
- **04f** the same test under four carrier definitions
- **04g** the 04c composition with each arm on its own AF definition

Panel 04e is the load-bearing result. It is restricted to variants that pass
the scMOCHA reliability gate, so both arms are variants scMOCHA considers
reliable and the only difference is mgatk's decision about them.

**04c against 04g.** 04c uses one measure for all three arms, so the bars are
comparable. 04g measures each arm on its own definition: the AF>5% arm on
cells at AF >= 0.05, because that rule is what defines it, and the other two
arms on all carrier cells, because they apply no AF rule. The three bars in
04g therefore do **not** share a measure - the AF>5% bar cannot fall below 0.05
while the other two can, so part of its apparent shift is definitional rather
than biological. The panel says so in its subtitle and the `measured` column of
`04-af-bins-own-definition.tsv` records which measure each arm used. Use 04c
for any like-for-like statement.

## Heteroplasmy measure

Every AF axis and the test use `af_carrier_median`: the median AF over cells
whose alt count is inconsistent with the background error rate under a
binomial test, Bonferroni-corrected over all cell-by-variant observations. No
AF floor, so the sub-0.05 region stays visible; no fixed read floor, so
background cells do not dilute the estimate.

Read `config.R.md` before changing this. Two other definitions were tried
first and both are wrong in opposite directions: the gate-based one is censored
at 0.05, and the 2-alt-read one is dominated by background. Panel 04f shows
what each would have concluded.

## Inputs

All under `${ISILON_BASE}/compare-with-mgatk/derived/<sample_id>`:

- `03-variant-gated.qs`
- `01-cell-af-s1.qs`
- `01-cell-pos-coverage.qs`

## Outputs

Figures in `newplots/compare-with-mgatk/figures/<sample_id>/`:

- `04a-af-ecdf.pdf`
- `04b-max-heteroplasmy.pdf`
- `04c-af-bins.pdf`
- `04d-cell-level-af.pdf`
- `04e-gate-rejected-af.pdf`
- `04f-measure-sensitivity.pdf`
- `04g-af-bins-own-definition.pdf`

Tables in `newplots/compare-with-mgatk/tables/<sample_id>/`, each with
`sample` as its first column:

- `04-af-bins.tsv`
- `04-af-bins-own-definition.tsv`
- `04-tests.tsv`
- `04-measure-sensitivity.tsv`

## Run

```bash
cd "$(git rev-parse --show-toplevel)"
pixi run Rscript newplots/compare-with-mgatk/04-heteroplasmy-spectrum.R --sample=GSE181279_GSM5494116_5PPE
```

`--sample` is required and must be one of the `sample_id` values in `SAMPLES`
in `config.R`.

## Verify

```bash
sample=GSE181279_GSM5494116_5PPE
cat newplots/compare-with-mgatk/tables/"${sample}"/04-tests.tsv
cat newplots/compare-with-mgatk/tables/"${sample}"/04-af-bins.tsv
cat newplots/compare-with-mgatk/tables/"${sample}"/04-measure-sensitivity.tsv
```

Expected **for GSE181279_GSM5494116_5PPE**: 45 variants in both sets, 185
mgatk only, 171 scMOCHA only. The gate test compares 54 variants mgatk keeps
against 215 it rejects, median carrier AF 0.205 versus 0.069, Wilcoxon
P = 0.0029, at a background alt rate of 0.000259. Panel 04c covers all three
arms: mgatk 230, scMOCHA variant call 736, scMOCHA AF>5% 216. Sensitivity:
P = 0.59 censored, 0.043 loose, 0.0029 adopted, 0.0036 on the max. In panel
04g the AF>5% arm is 37.0% in the 0.05-0.10 bin and 29.2% at >= 0.50, with
nothing below 0.05 by construction.

GSE181279 is the only sample where this test is computable. In the other four
`04-tests.tsv` carries `testable` FALSE and `p_value` NA, because mgatk's gate
passes nothing there and the comparison has only one group.

## Notes

This is a criterion-to-criterion comparison on identical reads. There is no
orthogonal truth set in the intermediate data, so nothing here is a
sensitivity, recall, or false-negative rate, and the figures must not be
described as such.

**A two-group test needs two groups.** `fn_testable()` requires at least
`CUTOFF_MIN_GROUP` = 3 variants on each side; below that no test is run. The
`testable` column in `04-tests.tsv` records that decision and `p_value` is NA,
and every row of `04-measure-sensitivity.tsv` likewise reports `p_value` NA.
The medians, the two group sizes and the background rate are still written, so
an untestable sample is still described rather than blank. Reporting NA rather
than omitting the row matters: a missing row cannot be told apart from a step
that never ran.

The panels go through `fn_or_empty()` for the same reason. In the four samples
where mgatk retains nothing, 04e and 04f are drawn as a panel stating why they
are empty instead of being skipped.
