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

- `${ISILON_BASE}/compare-with-mgatk/derived/03-variant-gated.qs`
- `${ISILON_BASE}/compare-with-mgatk/derived/01-cell-af-s1.qs`
- `${ISILON_BASE}/compare-with-mgatk/derived/01-cell-pos-coverage.qs`

## Outputs

- `newplots/compare-with-mgatk/figures/04a-af-ecdf.pdf`
- `newplots/compare-with-mgatk/figures/04b-max-heteroplasmy.pdf`
- `newplots/compare-with-mgatk/figures/04c-af-bins.pdf`
- `newplots/compare-with-mgatk/figures/04d-cell-level-af.pdf`
- `newplots/compare-with-mgatk/figures/04e-gate-rejected-af.pdf`
- `newplots/compare-with-mgatk/figures/04f-measure-sensitivity.pdf`
- `newplots/compare-with-mgatk/figures/04g-af-bins-own-definition.pdf`
- `newplots/compare-with-mgatk/tables/04-af-bins.tsv`
- `newplots/compare-with-mgatk/tables/04-af-bins-own-definition.tsv`
- `newplots/compare-with-mgatk/tables/04-tests.tsv`
- `newplots/compare-with-mgatk/tables/04-measure-sensitivity.tsv`

## Run

```bash
cd "$(git rev-parse --show-toplevel)"
pixi run Rscript newplots/compare-with-mgatk/04-heteroplasmy-spectrum.R
```

## Verify

```bash
cat newplots/compare-with-mgatk/tables/04-tests.tsv
cat newplots/compare-with-mgatk/tables/04-af-bins.tsv
cat newplots/compare-with-mgatk/tables/04-measure-sensitivity.tsv
```

Expected: 45 variants in both sets, 185 mgatk only, 171 scMOCHA only. The gate
test compares 54 variants mgatk keeps against 215 it rejects, median carrier AF
0.205 versus 0.069, Wilcoxon P = 0.0029, at a background alt rate of 0.000259.
Panel 04c covers all three arms: mgatk 230, scMOCHA variant call 736, scMOCHA
AF>5% 216. Sensitivity: P = 0.59 censored, 0.043 loose, 0.0029 adopted, 0.0036
on the max. In panel 04g the AF>5% arm is 37.0% in the 0.05-0.10 bin and 29.2%
at >= 0.50, with nothing below 0.05 by construction.

## Notes

This is a criterion-to-criterion comparison on identical reads. There is no
orthogonal truth set in the intermediate data, so nothing here is a
sensitivity, recall, or false-negative rate, and the figures must not be
described as such.
