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
- **05f** the mirror of 05d in scMOCHA's decision plane: cells at AF >= 0.05 with depth >= 10 against median alt reads per carrying cell
- **05g** 05f split into one facet per arm, with all variants repeated in grey

Panel 05c is the honesty check on criterion C2. scMOCHA's extra 10-alt-read
requirement drops variants that mgatk reports; this panel asks whether those
are weakly supported calls or genuine low-heteroplasmy ones.

Panels 05d and 05e answer where the arm-specific variants sit in the plane
mgatk actually decides on. Everything mgatk reports is in the upper-right
quadrant by construction, so the informative points are the scMOCHA AF>5%
variants outside it.

Panels 05f and 05g ask the same question the other way round, which is what
makes the pair diagnostic: given scMOCHA's criteria, where do **mgatk's**
variants land. The x axis is exactly the quantity scMOCHA's reliability gate
thresholds, so every scMOCHA AF>5% variant is right of the line by
construction and the informative points are the mgatk variants left of it.

**The y axis is not a criterion.** Median alt reads per carrying cell
summarises the support the per-cell confident-call rule acts on, but the rule
itself asks for >= 3 individual cells each at >= 10 alt reads, which a median
cannot express. The horizontal line is orientation only; first-failed-criterion
attribution lives in step 03's `03e`, not here.

## Inputs

All under `${ISILON_BASE}/compare-with-mgatk/derived/<sample_id>`:

- `03-variant-gated.qs`
- `01-cell-af-s1.qs`
- `01-cell-pos-coverage.qs`

## Outputs

Figures in `newplots/compare-with-mgatk/figures/<sample_id>/`:

- `05a-vmr-strand-plane.pdf`
- `05b-gate-rejected-af-bins.pdf`
- `05c-read-support.pdf`
- `05d-arm-in-mgatk-plane.pdf`
- `05e-arm-plane-facets.pdf`
- `05f-arm-in-scmocha-plane.pdf`
- `05g-scmocha-plane-facets.pdf`

Tables in `newplots/compare-with-mgatk/tables/<sample_id>/`, each with
`sample` as its first column:

- `05-vmr-strand-rejected.tsv`
- `05-arm-in-mgatk-plane.tsv`
- `05-arm-in-scmocha-plane.tsv`
- `05-read-support.tsv`

## Run

```bash
cd "$(git rev-parse --show-toplevel)"
pixi run Rscript newplots/compare-with-mgatk/05-vmr-strand.R --sample=GSE181279_GSM5494116_5PPE
```

`--sample` is required and must be one of the `sample_id` values in `SAMPLES`
in `config.R`.

## Verify

```bash
sample=GSE181279_GSM5494116_5PPE
cat newplots/compare-with-mgatk/tables/"${sample}"/05-vmr-strand-rejected.tsv
cat newplots/compare-with-mgatk/tables/"${sample}"/05-read-support.tsv
cat newplots/compare-with-mgatk/tables/"${sample}"/05-arm-in-mgatk-plane.tsv
cat newplots/compare-with-mgatk/tables/"${sample}"/05-arm-in-scmocha-plane.tsv
```

Expected **for GSE181279_GSM5494116_5PPE**: the rejection rate is highest in
the lowest AF bin (76 of 84 below 0.01, 90.5%) and lower in the middle bins.
Median alt reads per carrier cell is 9 for mgatk-only variants against 48 for
shared and 38 for scMOCHA-only. In mgatk's plane: 45 both arms, 185 mgatk
only, 168 scMOCHA AF>5% only, 403 scMOCHA call only, 445 neither.

In scMOCHA's plane the arm counts reconcile with step 03 exactly: 45 both arms
(30 + 15 by read support), 185 mgatk only (173 failing both axes, 9 failing
read support alone, 3 failing the cell gate alone), and 171 scMOCHA AF>5% only,
all of them right of the cell gate by construction. **176 of the 185 mgatk-only
variants fall short of scMOCHA's 10-cell gate.**

In the other four samples mgatk retains nothing, so the mgatk-only and
both-arms classes are empty and the tables are correspondingly short or
header-only.

## Notes

Alt reads per cell are reconstructed as `AF x depth`. Neither caller stores the
alt count itself in these outputs; they store AF and total coverage. The
product is exact up to the float rounding in the AF column, and the figure is
read as a distribution rather than as exact integers.

Panels 05d and 05e drop variants with no mgatk `vmr` or a `vmr` of zero, so
"scMOCHA AF>5% only" is 168 there against 171 in the step 03 overlap for
GSE181279. The three missing variants are ones mgatk never proposed as
candidates. Panels 05f and 05g do not use mgatk coordinates, so all 171 are
present there; that difference is why the two tables disagree by three.

`n_cells_gate` is 0 for a variant no cell carries above the AF floor, which a
log axis cannot show, so the x axis of 05f and 05g is pseudo-log. A median over
small integer counts puts most low-support variants on y = 1 or 2, where the
arms occlude each other in the combined panel; 05g exists because faceting is
the only honest fix, as jittering a median would move points off their own
value.

All seven panels go through `fn_or_empty()`. Four of the five samples have an
empty original-mgatk arm, so a panel with no points is drawn with a note
saying why rather than skipped; a missing figure would read as a failed run.
