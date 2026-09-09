# 05-vmr-strand.R

## Purpose

Criterion C4 in detail. VMR and strand correlation are computed by both callers
but acted on only by mgatk, so plotting scMOCHA's variants in the same plane
shows exactly which of them mgatk's gate would have discarded and at what
heteroplasmy.

- **05a** the VMR / strand-correlation plane for each rule set, mgatk's cutoffs drawn on both, coloured by carrier heteroplasmy (`af_carrier_median`, see `config.R.md`)
- **05b** heteroplasmy composition of what mgatk's gate keeps versus discards
- **05c** per-cell alt-read support behind each variant class
- **05d** the same plane on mgatk's own coordinates, coloured by which arm reports the variant
- **05e** 05d split into one facet per arm, with all variants repeated in grey

Panel 05c is the honesty check on criterion C2. scMOCHA's extra 10-alt-read
requirement drops variants that mgatk reports; this panel asks whether those
are weakly supported calls or genuine low-heteroplasmy ones.

Panels 05d and 05e answer where the arm-specific variants sit in the plane
mgatk actually decides on. Everything mgatk reports is in the upper-right
quadrant by construction, so the informative points are the scMOCHA AF>5%
variants outside it.

## Inputs

- `${ISILON_BASE}/compare-with-mgatk/derived/03-variant-gated.qs`
- `${ISILON_BASE}/compare-with-mgatk/derived/01-cell-af-s1.qs`
- `${ISILON_BASE}/compare-with-mgatk/derived/01-cell-pos-coverage.qs`

## Outputs

- `newplots/compare-with-mgatk/figures/05a-vmr-strand-plane.pdf`
- `newplots/compare-with-mgatk/figures/05b-gate-rejected-af-bins.pdf`
- `newplots/compare-with-mgatk/figures/05c-read-support.pdf`
- `newplots/compare-with-mgatk/figures/05d-arm-in-mgatk-plane.pdf`
- `newplots/compare-with-mgatk/figures/05e-arm-plane-facets.pdf`
- `newplots/compare-with-mgatk/tables/05-vmr-strand-rejected.tsv`
- `newplots/compare-with-mgatk/tables/05-arm-in-mgatk-plane.tsv`
- `newplots/compare-with-mgatk/tables/05-read-support.tsv`

## Run

```bash
cd "$(git rev-parse --show-toplevel)"
pixi run Rscript newplots/compare-with-mgatk/05-vmr-strand.R
```

## Verify

```bash
cat newplots/compare-with-mgatk/tables/05-vmr-strand-rejected.tsv
cat newplots/compare-with-mgatk/tables/05-read-support.tsv
cat newplots/compare-with-mgatk/tables/05-arm-in-mgatk-plane.tsv
```

Expected: the rejection rate is highest in the lowest AF bin (76 of 84 below
0.01, 90.5%) and lower in the middle bins. Median alt reads per carrier cell is
9 for mgatk-only variants against 48 for shared and 38 for scMOCHA-only. In
mgatk's plane: 45 both arms, 185 mgatk only, 168 scMOCHA AF>5% only, 403
scMOCHA call only, 445 neither.

## Notes

Alt reads per cell are reconstructed as `AF x depth`. Neither caller stores the
alt count itself in these outputs; they store AF and total coverage. The
product is exact up to the float rounding in the AF column, and the figure is
read as a distribution rather than as exact integers.

Panels 05d and 05e drop variants with no mgatk `vmr` or a `vmr` of zero, so
"scMOCHA AF>5% only" is 168 there against 171 in the step 03 overlap. The
three missing variants are ones mgatk never proposed as candidates.
