# 08-call-af-depth.R

## Purpose

The allele-frequency spectrum of the scMOCHA variant calls in one sample, and
the read support behind it. These are the per-sample versions of the two
panels; step 07 draws the same two pooled over `CROSS_SAMPLE_FIG_IDS`, using
the same builders in `config.R` so the two can never drift into drawing the
same quantity two different ways.

- **08a** one dot per sample-variant call, stacked by AF on a fixed log10 grid,
  drawn twice: once on the prevalence measure and once on the carrier
  heteroplasmy measure
- **08b** one dot per cell in which the variant is callable, cell AF against
  that cell's depth

**08a deliberately draws two panels rather than picking one.** The two measures
disagree about exactly the thing the panel exists to show. Pooled over the
eight figure samples, the prevalence measure puts 5 of 101 calls below the 5%
cutoff and reaches down to 0.065%; the carrier-heteroplasmy measure puts
**none** below it and bottoms out at 11.8%. Neither is wrong - they answer
different questions - but a single-panel version would silently commit the
reader to one answer. See `D12`, `D17` and `M14`.

**Prevalence is not heteroplasmy.** `mean` is total alt reads over total
coverage across every cell, so a variant at 60% heteroplasmy in 12 of 7,210
cells scores below 0.002. Quote the prevalence panel as prevalence.

**08b is a statement about detectability, not about biology.** A cell enters it
when it clears scMOCHA's per-cell rule: `CUTOFF_ALT_STRAND` = 2 alt reads on
**each strand**, at `CUTOFF_MIN_READS` = 10 reads of **depth**. The depth floor
is why nothing appears left of 10 on the x axis.

**The strand pair is approximated by its sum.** The cache carries AF and depth,
not the forward/reverse split, so alt reads are rebuilt as `AF x depth` and the
pair condition becomes `>= 4` alt reads in total. That is permissive: 4 alt
reads all on one strand pass here and would fail in the caller.

**The caller's source disagrees with its own comment on the third condition.**
`newplots/thecode/scmocha-mgatk-variant-calling.py` comments
`# C.J. minimum total coverage >=10` but the code is
`(fwd >= 2) & (rev >= 2) & ((fwd + rev) >= low_coverage_threshold)` evaluated on
the per-base **alt** matrices, which is 10 alt reads rather than 10 reads of
depth - `total_coverage_variant_df` never enters the condition. This panel uses
the depth reading, as intended. The two are not interchangeable: 10 alt reads
is far stricter, and which one the shipped `variant_stats` actually used is
worth settling before the figure is published. See `M19`.

**Almost nothing survives below the cutoff, and that is the finding.** Pooled
over the eight figure samples the median cell depth at a called position is 4
reads. The rule leaves 108,993 pairs with a minimum AF of 3.77% and **1** below
5%. In these shallow 3' libraries there is essentially no cell deep enough to
call a low-heteroplasmy variant - which is what the panel title says.

**The strand half of the rule is not reconstructable.** The cached data carries
AF and depth, not forward and reverse alt counts, so alt reads are rebuilt as
`AF x depth`. A cell clearing both floors is therefore a permissive upper bound
on "scMOCHA would have called it here".

## Inputs

All under `${ISILON_BASE}/compare-with-mgatk/derived/<sample_id>`:

- `03-variant-gated.qs`
- `01-cell-af-s1.qs`
- `01-cell-pos-coverage.qs`

## Outputs

Figures in `newplots/compare-with-mgatk/figures/<sample_id>/`:

- `08a-call-af-spectrum.pdf`
- `08b-cell-af-depth.pdf`

Tables in `newplots/compare-with-mgatk/tables/<sample_id>/`, each with
`sample` as its first column:

- `08-call-af.tsv` - one row per call per measure
- `08-cell-af-depth.tsv` - one row per cell observation

`08-cell-af-depth.tsv` is the large one: 108,993 rows across the eight figure
samples, 225,803 for GSE181279 alone. It is **not** a workbook sheet for that
reason; it stays a TSV beside the figure it backs. Step 07 reads both tables to
build the pooled panels, so this step runs before 07.

## Run

```bash
cd "$(git rev-parse --show-toplevel)"
pixi run Rscript newplots/compare-with-mgatk/08-call-af-depth.R --sample=GSE279945_GSM8583916_3PV3
```

`--sample` is required and must be one of the `sample_id` values in `SAMPLES`
in `config.R`.

## Verify

```bash
sample=GSE279945_GSM8583916_3PV3
head -3 newplots/compare-with-mgatk/tables/"${sample}"/08-call-af.tsv
wc -l newplots/compare-with-mgatk/tables/"${sample}"/08-cell-af-depth.tsv
```

Expected: `08-call-af.tsv` holds twice the call count, one row per measure.
Call counts per sample in registry order are 20, 6, 9, 20, 14, 736, 24, 1, 2
and 25; callable cell observations 4,302, 18,261, 13,038, 2,564, 27,197,
248,994, 59,554, 5,566, 6,382 and 62,858.

## Notes

Alt reads are reconstructed as `AF x depth`; neither caller stores the alt
count. A cell enters 08b when that product reaches `CUTOFF_ALT_READS` **and**
its depth reaches `CUTOFF_MIN_READS`, so background observations and
under-covered cells are both excluded. Carrier membership is still a different
question and is settled by the binomial test in `fn_carrier_stats()`; this
panel asks only whether scMOCHA could call the variant in that cell.

Both axes are decade breaks via `scales::label_log()`. The y axis always covers
at least 0.01% to 100% so the requested decades are on the panel, but the lower
bound extends past that whenever the data go lower - a fixed 0.01% floor
silently dropped 126 observations in GSE181279, where cells reach 70,107 reads
and so can report an AF of 0.0063%.

The legend carries a key for every AF band even when a band holds no cells, via
an invisible placeholder point per level; ggplot otherwise draws a bare label
with no glyph for an absent level.

The stacking grid for 08a is fixed in `config.R` (`CALL_AF_LOG_MIN`,
`CALL_AF_STACK_STEP`) rather than derived from the data, so bin width is
identical across samples and across the two measures and the panels stay
comparable. Calls below 1e-5 are clamped into the first bin and counted in the
subtitle.

Both panels go through `fn_or_empty()`. A sample whose call arm is empty gets a
panel saying so rather than a missing file.
