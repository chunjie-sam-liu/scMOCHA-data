# 07-cross-sample.R

## Purpose

What holds across all ten samples rather than inside one. Steps 01 to 05
answer the comparison per sample; this step puts the ten answers side by side
so a claim can be made about the method instead of about GSE181279.

- **07a** variants retained by each of the three arms, in every sample
- **07b** what the mgatk cell filter costs per sample
- **07c** strand correlation of mgatk's own S1 variants against its 0.65 floor
- **07d** which mgatk cutoff rejects each of its own S1 variants
- **07e** heteroplasmy of the variants mgatk's gate discards, per sample

The headline is 07a: original mgatk retains **zero** variants in seven of the
ten samples and exactly one in two more, while the scMOCHA call retains 20, 6,
9, 20, 14, 736, 24, 1, 2 and 25 in registry order. 07c and 07d say which cutoff
produced those zeros - in every one of them it is the strand-correlation floor,
not VMR and not the cell count.

## Inputs

For every sample in `SAMPLES`:

- `${ISILON_BASE}/compare-with-mgatk/derived/<sample_id>/03-variant-gated.qs`
- `newplots/compare-with-mgatk/tables/<sample_id>/02-cell-inclusion.tsv`

Both are required for all ten samples. A missing cache stops the step with
the `--sample=` command that would rebuild it, because a silently absent
sample would be indistinguishable from a sample with no variants.

## Outputs

Figures in `newplots/compare-with-mgatk/figures/cross-sample/`:

- `07a-arm-yield.pdf`
- `07b-cell-filter.pdf`
- `07c-strand-correlation.pdf`
- `07d-exclusion-reasons.pdf`
- `07e-gate-rejected-af.pdf`

Tables in `newplots/compare-with-mgatk/tables/cross-sample/`:

- `07-sample-overview.tsv` - one row per sample: registry columns, S0 and S1
  counts for both callers, what each arm retained, cells total and dropped
- `07-arm-yield.tsv` - sample by arm variant counts, the 07a data
- `07-cell-filter.tsv` - cells total, cells dropped, detections in the dropped
  cells, variants lost by the cell filter, and the dropped fraction
- `07-strand-support.tsv` - per sample: n, median, max and above-floor count of
  mgatk S1 strand correlation, plus the Spearman rho between strand
  correlation and per-variant coverage with its P value
- `07-exclusion-reasons.tsv` - sample by first-failed-criterion counts over
  mgatk's own S1 variants
- `07-gate-test.tsv` - per sample: n passed, n rejected, the two medians,
  `testable`, and `p_value`

Step 06 reads all six into the workbook, which is why 07 runs before it.

## Run

```bash
cd "$(git rev-parse --show-toplevel)"
pixi run Rscript newplots/compare-with-mgatk/07-cross-sample.R
```

No `--sample`: it is the step that spans them. Run it after steps 01 to 05
have completed for **every** sample, and before step 06.

## Verify

```bash
cat newplots/compare-with-mgatk/tables/cross-sample/07-arm-yield.tsv
cat newplots/compare-with-mgatk/tables/cross-sample/07-strand-support.tsv
cat newplots/compare-with-mgatk/tables/cross-sample/07-gate-test.tsv
ls -la newplots/compare-with-mgatk/figures/cross-sample/
```

Expected retained counts, original mgatk / scMOCHA call / scMOCHA AF>5%:

| Sample | mgatk | scMOCHA call | scMOCHA AF>5% |
| --- | ---: | ---: | ---: |
| GSE149689_GSM4509019_3PV3 | 0 | 20 | 20 |
| GSE163314_GSM4976997_3PV2 | 0 | 9 | 9 |
| GSE163668_GSM4995445_5PR2 | 0 | 20 | 18 |
| GSE181279_GSM5494116_5PPE | 230 | 736 | 216 |
| GSE271107_GSM8369876_3PV3 | 0 | 2 | 0 |

In `07-strand-support.tsv`, `n_above_floor` is 0, 0, 0, 300, 0 in that order.
The four zero samples reach a maximum strand correlation of 0.569, 0.006, 0.59
and 0.048, all under the 0.65 floor. GSE181279 is the only sample where the
gate test is computable, so `testable` is TRUE in one row of
`07-gate-test.tsv` and FALSE in four, with `p_value` NA in those four. Five
PDFs and six TSVs.

## Notes

**07c is descriptive.** It shows where mgatk's own S1 variants sit relative to
the strand-correlation floor that rejects them, and nothing more. It must not
be captioned as showing that shallow coverage drives strand correlation down:
within GSE181279, Spearman rho between a variant's strand correlation and its
coverage is **-0.013 with P = 0.66**. The `rho_strand_vs_coverage` and
`rho_p_value` columns of `07-strand-support.tsv` report that correlation per
sample as measured, precisely so the mechanism claim cannot be made from this
panel by accident.

**No pooled or combined P value is computed anywhere in this step.** The gate
test needs both a passing and a rejected group; four samples have no passing
variant at all, so they contribute no comparison. A stratified or pooled test
across the five would be arithmetic over one informative stratum and would
collapse to GSE181279's own P value while reading as a five-sample result.
`07-gate-test.tsv` therefore reports one row per sample with an explicit
`testable` flag and `p_value = NA` where there is nothing to test.

**Every panel tolerates an empty group.** Four samples have an empty
original-mgatk arm, so a panel that would error or silently drop a group is
drawn through `fn_or_empty()` with a note stating why it is empty, rather than
skipped. A skipped figure looks like a failed run; a zero is the result.

**Counts are never read off a log axis.** 07a uses a plain linear axis with
`scales::comma`, so the four zero bars sit on the baseline and carry a printed
`0` label. A log scale drops a zero bar and its label without warning, and the
zero is the headline of this step.

**Samples are labelled by GSE and GSM, never by chemistry.** Several samples
share a chemistry, so chemistry cannot identify a sample; axis labels use
`{gse}\n{gsm}` and every cross-sample table carries `sample` as `{gse}_{gsm}`
alongside a `sample_id` column that keeps the join back to the per-sample
directories. Chemistry stays in the `SAMPLES` registry and in the workbook's
`00_Samples` sheet, which is where a reader looks it up.

**The panels draw a subset; the tables never do.** `CROSS_SAMPLE_FIG_EXCLUDE`
in `config.R` lists the samples kept out of the figures - currently the two 5'
libraries, so the panels compare within one library family. `fn_fig_subset()`
is applied at each `ggplot()` call and nowhere else; every `export()` below
uses the unfiltered object, so all ten samples keep their rows. Any subtitle
that quotes a sample count quotes the count actually drawn, not
`length(SAMPLE_IDS)`.

**That exclusion removes the only informative sample.** GSE181279 is the one
sample where mgatk retains variants and the only one where the gate test is
computable. With it out of the panels, `07a` shows mgatk at zero everywhere and
`07e` shows no computable test at all. Neither is a claim about the resource;
both are a property of the six samples drawn. `07-gate-test.tsv` and
`07-arm-yield.tsv` carry the full picture.

**07d colours the retained slice red.** `Retained` is the one outcome a reader
looks for first, so it carries the NEJM red in `color_exclusion`; every
exclusion reason is a non-red hue. The four classes that co-occur in 07d -
retained, `strand r only`, `VMR and strand r`, `VMR only` - are mutually
distinguishable.

**Aggregations are merged back onto `SAMPLES`.** A `by = sample` summary emits
no row for a sample with no qualifying variants, and a missing row cannot be
told apart from a sample that never ran, so the overview and gate-test tables
carry an explicit zero for every sample in the registry.
