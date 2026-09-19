# DECISION: extend the mgatk comparison from one sample to five

Stage: `newplots/compare-with-mgatk/`
Campaign: `2026-09-17-multi-sample-comparison`
Append-only. Each entry: what was chosen, what it beat, why, evidence, date.

Parent: `DECISION.md` (stage log, `D1`..`D19`). Entries here are numbered
`M1..Mn` so the two never collide.

---

## M1 - Original mgatk retaining zero variants in four samples is a result, not a data problem

**Date:** 2026-09-17
**Chosen:** Treat `S2 mgatk = 0` in GSE149689, GSE163314, GSE163668 and
GSE271107 as the measured outcome, report it, and make every script tolerate an
empty arm.
**Beat:** Treating it as a broken input, a parsing bug, or grounds for dropping
those samples from the comparison.
**Why:** The rejection is explicit and not marginal. Every mgatk S1 variant in
those four samples fails `strand_correlation > 0.65`, the highest observed value
anywhere in the four is 0.59, and many are negative. Both `vmr` and
`strand_correlation` are populated, so it is not an `NA` artefact. Dropping the
samples would keep only the one sample where mgatk happens to work, which is the
selection bias the editor is entitled to object to.
**Evidence:** `awk -F'\t' 'NR>1 && $7>=3 && $4>0.01 && $13>0.65'` over
`cell.variant_stats_mgatk_original.tsv.gz` returns 0 rows for each of the four,
and 230 for GSE181279. Per-variant listing of the mgatk S1 sets confirms the
strand-correlation range.
**Status:** approved 2026-09-17; implemented and verified.

---

## M2 - The cross-sample claim is yield plus mechanism, not a pooled test

**Date:** 2026-09-17
**Chosen:** State the five-sample result as (a) variants retained per arm per
sample and (b) a pooled strand-correlation-versus-read-support panel over all
five samples' mgatk S1 variants. No combined p-value.
**Beat:** A stratified Wilcoxon or a Fisher combination of the per-sample gate
tests.
**Why:** Step 04's test compares variants the mgatk gate rejects against those
it keeps, among variants passing the scMOCHA gate. Four samples have an empty
"keeps" group, so the test is undefined there. A stratified test over five
strata where four contribute nothing is the single-sample test relabelled, and
Fisher's method over four undefined p-values combines nothing. Reporting either
would overstate the evidence in a response to an editor.
**Evidence:** `S2 mgatk` column of the verified inventory in the plan: 0, 0, 0,
230, 0.
**Status:** approved 2026-09-17 as `Q6`; no combined p-value appears anywhere
in the stage.

---

## M3 - `kept_by_mgatk` comes from the mgatk AF matrix, not from recomputing the coverage rule

**Date:** 2026-09-17
**Chosen:** Define the mgatk cell set as the barcode set of
`cell.cell_heteroplasmic_df_mgatk_original.tsv.gz`, keep
`mean_coverage > CUTOFF_CELL_MEANCOV` as the reported criterion, and log the
symmetric difference instead of asserting `setequal`.
**Beat:** Keeping the recomputed rule and loosening the assertion to a
tolerance.
**Why:** `cell.depthTable.txt` stores mean coverage rounded to two decimals.
Cells whose true mean coverage straddles 10 are all written as `10.0`, so the
rule is not reproducible from that file. The caller's own output is the only
non-ambiguous record of which cells it used, and every downstream count in step
02 is a statement about what mgatk actually did.
**Evidence:** GSE163314 has four cells recorded as `10.0`. Recomputing
`> 10` gives 274 kept cells; `>= 10` gives 278; mgatk's own AF matrix carries
276. Barcode-set `comm` shows the two extra cells are among the four ties. The
other four samples agree exactly under `> 10`, so this is a no-op for them and
for the published GSE181279 numbers.
**Status:** approved 2026-09-17; implemented and verified - GSE181279 still
reports 7,210 cells with 486 dropped.

---

## M4 - Panel 07c is descriptive, because the coverage mechanism is not in the data

**Date:** 2026-09-17
**Chosen:** Draw 07c as the distribution of mgatk's own S1 strand correlations
against its 0.65 floor, per sample, and caption it as that and nothing more.
Report the Spearman correlation between strand correlation and per-variant
coverage in `07-strand-support.tsv` as measured, including where it is null.
**Beat:** The panel this campaign was planned with: a pooled scatter of strand
correlation against coverage captioned "mgatk's strand-correlation floor tracks
coverage".
**Why:** The planned caption was a mechanistic claim, and the first run
falsified it in the only sample with enough variants to test it. Shipping it
would have put an unsupported causal statement into a response to an editor,
which is exactly the kind of claim that gets checked. What the data do support
is a flat descriptive fact, and that fact is strong on its own.
**Evidence:** Spearman rho between `cov_mgatk` and `strand_mgatk` over the
1,247 mgatk S1 variants of GSE181279 is **-0.013, P = 0.66**. The smaller
samples give +0.66 (n = 19, P = 0.003), +0.36 (n = 21, P = 0.11), and no
estimate at n = 4 and n = 6. The descriptive fact is unambiguous instead:
0 of 19, 0 of 4, 0 of 21 and 0 of 6 variants reach the 0.65 floor in the four
shallow samples, with maxima 0.569, 0.006, 0.590 and 0.048, while 300 of 1,247
reach it in GSE181279.
**Consequence:** `07c-strand-vs-coverage.pdf` became
`07c-strand-correlation.pdf`. Nothing in the stage now claims a mechanism for
the bias; it reports that the floor is unreachable in shallow libraries.

---

## M5 - Floor an unestimable background rate instead of letting it degenerate

**Date:** 2026-09-17
**Chosen:** In `fn_background_rate()`, when the estimate is 0 or undefined,
fall back to one alt read over all observed depth, and warn.
**Beat:** (a) leaving the raw estimate, and (b) the first fix, which floored at
`1 / background_depth` only.
**Why:** The carrier definition is a binomial test against this rate, and both
degenerate values break it, in opposite directions. A rate of 0 makes
`pbinom(alt - 1, depth, 0, lower.tail = FALSE)` zero for any cell holding a
single alt read, so every such cell becomes a carrier and the measure reverts
to the noise-dominated definition that `D17` exists to remove. A rate of 1
makes the same quantity 1, so no cell is ever a carrier and the sample silently
drops out of every downstream table. Option (b) produced exactly that.
**Evidence:** Before the fix, `04-tests.tsv` carried `background_alt_rate = 1`
for GSE163314 - which has no qualifying background cell at all, so `0 / 0`
floored to `1 / 1` - and that sample reported 0 gate-rejected variants in
`07-gate-test.tsv`. After it the rate is 2.54e-06 and the sample reports 11.
GSE149689 moved from a floored 0.0172 to 1.47e-06 and GSE271107 to 6.58e-06.
GSE163668 (0.00277) and GSE181279 (0.000259) estimate a real rate and are
untouched, so no published number moved.

---

## M6 - Zeros are written, not omitted

**Date:** 2026-09-17
**Chosen:** Cross-sample tables are merged back onto `SAMPLES` so every sample
has a row, and an empty arm gets a drawn panel carrying a stated note.
**Beat:** Letting `by = sample` return only the groups it observed, and letting
an empty panel be skipped.
**Why:** With four samples retaining nothing in one arm, an omitted row and a
crashed run look identical from the output directory. The zero is the finding
this campaign exists to report, so it has to be written down explicitly.
**Evidence:** `07-gate-test.tsv` initially held three of five samples; it now
holds five, two of them explicit zeros. Every sample directory holds the same
22 figures and 13 tables however degenerate the sample is.
**Status:** approved as `Q5`.

---

## M7 - The scMOCHA decision plane is a mirror pair, and its y axis is not a criterion

**Date:** 2026-09-18
**Chosen:** Add `05f` (combined scatter) and `05g` (faceted) plotting every S1
variant on scMOCHA's own quantities: x = cells at AF >= 0.05 with depth >= 10,
y = median alt reads per carrying cell. Label the x axis as the reliability
gate it literally is, and label the y axis explicitly as a summary that is
**not** a criterion.
**Beat:** (a) a single combined panel only; (b) using the S1 confident-cell
rule itself as the y axis; (c) reusing `nconf_scmocha` for y.
**Why:** `05d` answers "where do scMOCHA's variants sit in mgatk's plane"; the
reverse question is what a reader actually asks when judging what scMOCHA's
criteria do to mgatk's output, and nothing in the stage answered it. The x axis
is exact - the gate thresholds that very count - so every scMOCHA AF>5% variant
lies right of the line by construction and the mgatk variants left of it are
the finding. The y axis cannot be exact: the rule asks for >= 3 individual
cells each at >= 10 alt reads, which no per-variant median can express.
Labelling it as a criterion would have invited the reader to read
first-failed-criterion attribution off a panel that cannot carry it, and that
attribution already exists in `03e`. (b) was rejected because a per-variant
plane cannot hold a per-cell rule; (c) because `nconf_scmocha` is strongly
correlated with the x axis and would have produced a diagonal band carrying one
piece of information, not two.
**Evidence:** Arm counts in `05-arm-in-scmocha-plane.tsv` reconcile exactly
with step 03 for GSE181279: 45 both arms (30 + 15), 185 mgatk only (173 short
on both axes, 9 on read support alone, 3 on the cell gate alone), 171 AF>5%
only, all right of the gate. 176 of the 185 mgatk-only variants fall short of
the 10-cell gate. Counts hold for all five samples; the four shallow ones plot
3, 11, 21 and 20 variants.
**Consequence:** `05f` alone is unreadable where it matters. A median over
small integer counts lands most low-support variants on y = 1 or 2, so the arms
occlude each other exactly in the region the figure exists to show. `05g`
facets by arm for that reason; jittering was rejected because it would move a
median off its own value.

---

## M8 - Cross-sample presentation: GSM labels, linear counts, red retained, arm-specific planes

**Date:** 2026-09-18
**Chosen:** Four presentation changes, all requested by the user on 2026-09-18.

1. Chemistry is no longer used anywhere in the cross-sample layer. Panels label
   samples `{gse}\n{gsm}`; tables carry `sample` as `{gse}_{gsm}` plus a
   `sample_id` join key.
2. 07a moves from `transform_pseudo_log` to a plain linear axis.
3. `color_exclusion["Retained"]` becomes the NEJM red `#BC3C29`; the hue it
   displaced moves to `VMR and strand r`, and `< 10 cells at AF >= 0.05` takes
   the grey that `Retained` gave up.
4. 05d, 05e, 05f and 05g plot only the three arm-specific classes; `Both arms`
   and `Neither arm` are dropped from the panels.

**Beat:** Keeping chemistry as the axis label; keeping the pseudo-log axis;
keeping grey for `Retained`; keeping all five arm classes in the planes.
**Why:** (1) Two of the five samples are `SC3Pv3`, so chemistry does not
identify a sample and a reader cannot tell those two panels apart. GSE plus GSM
is unique and is what a reader can look up. (2) The pseudo-log axis was adopted
so a zero bar would still be drawn, but with a maximum of 736 a linear axis
shows the zeros just as well and does not distort the 230-vs-736 comparison,
which is the thing the panel is for. (3) `Retained` is the outcome a reader
looks for first in 07d and red is where the eye goes; the constraint this
creates - every exclusion reason must be non-red - is satisfied and checked
against the four classes that actually co-occur in 07d. (4) `Both arms` and
`Neither arm` do not separate the callers, and `Neither arm` is the largest
class, so it rendered as a grey mass covering the arm-specific points the
panels exist to show.
**Evidence:** `grep -l "SC5P\|SC3Pv" tables/cross-sample/*.tsv` returns
nothing. 07a renders all four zero bars with printed `0` labels on a linear
axis to 800. In 07d the four co-occurring classes are red, steel blue, pink and
pale yellow. Step 05 and step 07 re-ran for all five samples at exit 0 with no
warnings; the workbook rebuilt at 23 sheets.
**Consequence:** The exported tables are unchanged in content - `tab_f` and the
mgatk-plane table are still built from the unfiltered data, so `Both arms` and
`Neither arm` keep their counts. Only the panels are filtered. `sample_id` is
retained in every cross-sample table because it is the join key back to
`figures/<sample_id>/` and `tables/<sample_id>/`; it still carries the
chemistry suffix as part of the directory name.

---

## M9 - Eight samples, and the cross-sample figures restricted to the 3' libraries

**Date:** 2026-09-18
**Chosen:** Add GSE155673_GSM4712895_3PV3, GSE188632_GSM5687372_3PV3 and
GSE220189_GSM6793474_3PV3 to the registry, taking the stage to eight samples.
Exclude `GSE163668_GSM4995445_5PR2` and `GSE181279_GSM5494116_5PPE` from the
cross-sample FIGURES only, via `CROSS_SAMPLE_FIG_EXCLUDE` in `config.R`; every
cross-sample table still carries all eight.
**Beat:** Excluding the two 5' samples from the tables as well, which was the
other reading of the request.
**Why:** The user asked for figures only and confirmed it explicitly. Keeping
the tables complete means the exclusion is a presentation choice a reader can
undo, not a deletion. Restricting the panels to the 3' libraries makes them a
comparison within one library family rather than across three chemistries.
**Consequence, and it is a real one:** GSE181279 is the only sample in which
mgatk retains any variant (230) and the only one in which the gate test is
computable, and it is one of the two excluded. The panels as drawn therefore
show mgatk retaining nothing in 6 of 6 and no computable test anywhere. That is
a property of the six samples drawn, not of the resource, and every affected
subtitle now says how many samples it is describing. `07-arm-yield.tsv` and
`07-gate-test.tsv` carry the full eight-sample picture. Anyone quoting a panel
without the table will understate what mgatk does on deep 5' data.
**Evidence:** Verified before use, per archive: all three new archives carry
the same 17-file layout and the same 14-column `variant_stats` schema. New
counts - GSE155673 s1_scmocha 6 / s1_mgatk 19, GSE188632 24 / 32, GSE220189
1 / 6; `s2_mgatk` is 0 in all three, with max strand correlation 0.196, 0.503
and 0.013 against the 0.65 floor. Across all eight, mgatk retains nothing in
seven and no variant reaches the strand floor in seven.

---

## M10 - Two defects found by adding samples

**Date:** 2026-09-18

**(a) The extraction script kept its own copy of the sample list.**
`00-extract-archives.sh` held a hard-coded `archives=()` array duplicating the
`SAMPLES` registry, so adding three registry rows extracted nothing and the
three new samples had to be unpacked by hand. A second copy of a list is a
second source of truth and this one failed silently: the script exited 0 and
reported skipping five samples. It now reads `sample_id` and `archive` from
`config.R` at run time and prints how many the registry lists. Verified: the
script reports `registry lists 8 samples` and skips all eight as present.

**(b) 04f guarded on the wrong condition.** `fn_or_empty(nrow(d_e) > 0L, ...)`
is true whenever the rejected group has variants, so in a sample where mgatk's
gate passes nothing the panel still drew: all four `median_passed` values NA,
`geom_line()` with one observation per group, and `Removed 4 rows containing
missing values` from the log scale. The rendered panel was four lone points
that read as a real comparison. The guard now requires at least one measure
with both medians present. Verified: the warning is gone for all eight samples
and the panel states why it is empty.
