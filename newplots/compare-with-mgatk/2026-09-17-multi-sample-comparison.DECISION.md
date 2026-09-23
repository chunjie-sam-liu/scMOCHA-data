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

---

## M11 - Ten samples, a chemistry-agnostic palette, and a headline that had to change

**Date:** 2026-09-20
**Chosen:** Add `GSE175499_GSM5335510_3PV3` and `GSE279945_GSM8583916_3PV3`,
taking the stage to ten samples. Restate the headline: original mgatk retains
zero variants in **seven of ten** and exactly **one** in two more, rather than
"zero in all but the deep sample". Replace the chemistry-anchored
`color_sample` with the categorical `ggthemes::Classic_10`.
**Beat:** Keeping the "mgatk retains nothing anywhere except GSE181279"
wording; keeping the chemistry-anchored palette.
**Why:** The two new samples are the first outside GSE181279 in which a variant
clears mgatk's strand floor - max strand correlation 0.673 and 0.919 against
the 0.65 cutoff, one variant each. The previous wording would have been false
for them, and this is a rebuttal to an editor, so the claim has to match the
data exactly. The weaker statement is still decisive: across ten samples mgatk
yields 0, 0, 0, 0, 1, 230, 0, 0, 0, 1 while the scMOCHA call yields 1 to 736.
On the palette: seven of the ten samples are SC3Pv3, so anchoring colour to
chemistry would have meant seven greens, and `M8` already removed chemistry
from the cross-sample layer entirely, so the anchoring no longer carried
meaning. `Classic_10` has the best worst-case separation of the ten-colour
candidates checked under `clr_deutan()` (minimum pairwise distance 26.3 against
8.1 for Tableau_10); sample is also on the x axis wherever the palette is used,
so colour is a redundant encoding.
**Evidence:** Per-archive verification before use: both new archives carry the
17-file layout and the 14-column `variant_stats` schema. GSE175499 s1_scmocha
14 / s1_mgatk 20 / s2_mgatk 1; GSE279945 25 / 29 / 1. Cross-sample run reports
`10 samples in the tables, 8 in the figures; mgatk retains nothing in 7`.
**Consequence:** `fn_testable()` still refuses both new samples - one variant
is below the three-per-side minimum - so the gate test remains computable in
GSE181279 alone and `M2` stands unchanged.

---

## M12 - A group of exactly one broke two panels quietly

**Date:** 2026-09-20
**Chosen:** In 04b build the density from classes with at least two variants and
name any dropped class in the subtitle; in 04e compute the violin layer only
over groups with at least two members, leaving the boxplot to cover the rest.
**Beat:** Leaving ggplot to drop them with a warning.
**Why:** The two new samples are the first with a variant class of size exactly
one. `geom_density()` dropped that class from 04b entirely - the class vanished
from the panel with nothing but a console warning to say so, which is the same
silent-omission failure `M10(b)` fixed in 04f. 04e was less severe because the
boxplot still drew the one-variant group, but the warning is noise that hides
real ones.
**Evidence:** Before: `Groups with fewer than two data points have been
dropped` and `Removed 1 row containing missing values` on 04b, and `Groups with
fewer than two datapoints have been dropped` on 04e, for GSE175499 and
GSE279945. After: step 04 re-run across all ten samples, every sample CLEAN.

---

## M13 - The background-rate warning was reporting depth as reads

**Date:** 2026-09-20
**Chosen:** Reword the `fn_background_rate()` fallback warning to say
"{total_depth} background depth carried 0 alt reads" instead of "not estimable
from {total_depth} background reads".
**Beat:** Leaving it.
**Why:** The variable interpolated is the summed *depth* of background cells,
not an alt-read count. On GSE175499 it printed "not estimable from 22845
background reads", which reads as a contradiction - 22,845 reads is plenty to
estimate from - and cost an investigation to establish that the logic was in
fact correct: 22,845 depth carrying zero alt reads. A log line that makes
correct code look broken is a defect in its own right.
**Evidence:** `fn_background_rate()` floors only when `rate <= 0`, and rate is
`sum(round(af * depth)) / total_depth`, so the trigger is zero alt reads, not
low depth.

---

## M14 - The call AF spectrum is drawn on both measures, never one

**Date:** 2026-09-23
**Chosen:** Add the call-level AF spectrum (panel C of the user's reference
figure) and the cell-level AF-versus-depth panel (panel D), per sample as
`08a` / `08b` and pooled over `CROSS_SAMPLE_FIG_IDS` as `07f` / `07g`. Draw the
spectrum on **both** the prevalence measure and the carrier-heteroplasmy
measure, as two facets of one figure, rather than choosing one.
**Beat:** Picking one measure, which is what the reference figure does.
**Why:** The two measures disagree about exactly the quantity the panel exists
to show. Pooled over the eight figure samples: prevalence puts 5 of 101 calls
below the 5% cutoff and reaches 0.065%; carrier heteroplasmy puts **zero**
below it and bottoms out at 11.8%. Choosing prevalence would label a prevalence
statistic "Variant AF", which `D12` and the README explicitly forbid - `mean`
is total alt reads over total coverage across all cells, so it tracks how many
cells carry a variant, not how strongly. Choosing carrier heteroplasmy would
delete the low-AF tail that is the panel's entire subject, which is the
censoring failure `D17` already fixed once. Drawing both is the only option
that neither mislabels nor censors, and it makes the measure-dependence the
reader's to judge.
**Evidence:** `07-call-af-bins.tsv`. Prevalence: 1 call `<0.1%`, 2 `0.1-1%`,
2 `1-5%`, 1 `5-20%`, 95 `20-100%`. Heteroplasmy: 1 `5-20%`, 100 `20-100%`.
**Reconciliation with the reference figure:** our pooled call count is 101 over
the same eight samples, matching it exactly, as do 89 for the AF>5% arm and 2
for mgatk. Two numbers do not match and were not forced to: the reference
reports 65 distinct variants for the call arm where we count 70 - our 65 is the
AF>5% arm's distinct count, so the reference appears to have shifted that
column by one row - and it reports 12 of 101 calls below 5%, which no measure
in this stage reproduces; prevalence gives 5. Those two numbers should be
checked against whoever produced the figure before either version is published.

---

## M15 - run-all.sh held the same duplicated sample list that M10 fixed elsewhere

**Date:** 2026-09-23
**Chosen:** Read `SAMPLE_IDS` from `config.R` at run time in `run-all.sh`, the
same way `00-extract-archives.sh` now does, and print how many samples it is
about to run.
**Beat:** Leaving it.
**Why:** `M10(a)` fixed a hard-coded `archives=()` array in the extraction
script but missed the identical `samples=()` array in `run-all.sh`. That array
still listed **five** samples while the registry held ten, so `bash run-all.sh`
with no arguments would silently have processed half the stage and exited 0.
The defect had been latent since the registry went from five to eight.
**Evidence:** The array listed GSE149689, GSE163314, GSE163668, GSE181279 and
GSE271107 only. After the change the script prints `running 10 samples`.
**Consequence:** Step 08 was added to the per-sample step list in the same
edit, and the comment above the cross-sample call now records that 07 depends
on 08.

---

## M16 - The depth panel shows callable cells only, and the sub-5% region empties

**Date:** 2026-09-23
**Chosen:** A cell enters `08b` / `07g` only when the variant is callable in it,
meaning at least `CUTOFF_CELL_ALT_CALLED` = 2 x `CUTOFF_ALT_STRAND` = 4 alt
reads. Both axes use decade breaks; the y axis always spans at least 0.01% to
100% but extends lower whenever the data do.
**Beat:** Keeping every cell with one or more alt reads, which is what the
panel did when it was first added.
**Why:** The user asked for cells that hold the variant, defined as the variant
being called there by mgatk or scMOCHA. mgatk's per-cell rule is 2 alt reads on
each strand and scMOCHA's is that plus 10 in total, so the union of the two is
mgatk's and the floor is 4. The strand half cannot be rebuilt from AF and
depth, so this is a permissive upper bound on "the caller would have called it
here", and the caption says so.
**Consequence, and it is the substantive one:** the sub-5% region nearly
vanishes. Pooled over the eight figure samples, median cell depth at a called
position is 4 reads. Below-5% observations go from 196 at >= 1 alt read, to 21
at >= 2, to **1** at >= 4, to **0** at >= 10, while n falls from 382,773 to
197,158. The earlier version's low-AF cloud was background, not calls. The
honest reading is that these shallow 3' libraries contain almost no cell deep
enough to call a low-heteroplasmy variant, which is the panel's own title read
back at the data.
**Evidence:** Per-cell alt-read sweep over `CROSS_SAMPLE_FIG_IDS`, and the
`n_low` count now printed in every subtitle.

---

## M17 - A fixed axis floor silently cropped 126 observations

**Date:** 2026-09-23
**Chosen:** Compute the y lower bound as
`min(0.01, min(af_pct)) * 0.9` rather than pinning it at 0.01%.
**Beat:** The fixed `limits = c(0.01, 100)` written when the decade axis was
first added.
**Why:** GSE181279 reaches 70,107 reads in a cell, so 4 alt reads there is an
AF of 0.0057%. A hard floor at 0.01% put 126 real observations outside the
scale, and ggplot removed them with a warning rather than an error. Same
failure mode as `M10(b)` and `M12`: the panel looked fine and quietly held less
data than it claimed.
**Evidence:** `Removed 126 rows containing missing values` on GSE181279, min AF
0.0063% against the 0.01% floor. After the change that sample runs warning-free
and all ten samples report zero warnings.
**Consequence:** The requested decades are still guaranteed - the bound only
ever moves down, never up - so a sample with nothing below 0.01% still shows
the full 0.01%-100% span and the empty band under the cutoff stays visible as
an absence.

---

## M18 - M16 dropped scMOCHA's depth floor; the panel used mgatk's rule instead

**Date:** 2026-09-23
**Chosen:** A cell enters `08b` / `07g` only when **scMOCHA** would call the
variant there: `round(AF x depth) >= CUTOFF_ALT_READS` (10) **and**
`depth >= CUTOFF_MIN_READS` (10). `CUTOFF_CELL_ALT_CALLED` is removed.
**Beat:** The `M16` rule, `>= 2 x CUTOFF_ALT_STRAND` = 4 alt reads with no
depth condition.
**Why:** `M16` reasoned that "callable by mgatk or scMOCHA" is the union of the
two per-cell rules, and that the union is mgatk's because it is the weaker one.
That is true as set algebra and wrong as a figure: mgatk imposes **no depth
floor at all**, so the union silently discards scMOCHA's `depth >= 10`
requirement. The panel then drew a large cloud between depth 4 and 10 - cells
scMOCHA would never have called - on an x axis labelled as read depth. The user
spotted it as "why are there so many dots below 10 on the x axis".
**Rejected middle option:** `alt >= 4 & depth >= 10` keeps one more sub-5%
observation, but it is mgatk's read floor bolted onto scMOCHA's depth floor and
belongs to neither caller, so it cannot be described in one sentence.
**Evidence:** Pooled over the eight figure samples, cell-variant pairs and the
observed depth minimum: `alt>=4` 197,158 pairs at min depth **4**;
`alt>=4 & depth>=10` 108,993 at min depth 10; `alt>=10 & depth>=10` 98,704 at
min depth 10, min AF 6.97%, 0 below the 5% cutoff. After the change the pooled
table reports min depth 10 exactly.
**Consequence:** The sub-5% region is now empty rather than holding a single
point, and the cell count backing the panel falls from 49,512 to 34,939 of
59,181. Both are the honest consequence of applying the caller's own floors.

---

## M19 - The scMOCHA caller's code and its comment disagree on C2's third condition

**Date:** 2026-09-23
**Chosen:** The per-cell rule for `08b` / `07g` is
`CUTOFF_ALT_STRAND` alt reads on each strand at `CUTOFF_MIN_READS` reads of
**depth**, approximated as `round(AF x depth) >= 2 x CUTOFF_ALT_STRAND` and
`depth >= CUTOFF_MIN_READS`. This supersedes `M18`, which read the third
condition as 10 **alt** reads.
**Beat:** `alt >= CUTOFF_ALT_READS & depth >= CUTOFF_MIN_READS` (`M18`).
**Why:** The user, who is the caller's author, states the criterion is
`fwd >= 2 AND rev >= 2 AND total depth >= 10`. The source supports the comment
but not the code:

```python
# C.J. This the number of reads supporting the variant, requires >= 2 on both strands
# C.J. minimum total coverage >=10
variant_n_cells_conf_detected = (
    (fwd_cell_variant_df >= 2)
    & (rev_cell_variant_df >= 2)
    & ((fwd_cell_variant_df + rev_cell_variant_df) >= low_coverage_threshold)
).sum()
```

`fwd_cell_variant_df` and `rev_cell_variant_df` are per-base **alt** matrices -
`base_coverage_dict[base][0]` and `[1]` - and `heteroplasmic_df` is
`(fwd + rev) / total_coverage_variant_df`, so the third condition is 10 alt
reads. `total_coverage_variant_df`, the depth, never enters the condition. The
comment says coverage, the code says alt reads, and they are not
interchangeable.
**Open, and it matters:** the shipped `variant_stats` tables were produced by
the code, so `n_cells_conf_detected` in the data reflects the alt reading
whatever the intent was. The panel now uses the depth reading per the author's
instruction, which means the panel and the upstream `n_cells_conf_detected`
column are built on different thirds of C2. Settle which is intended before
publishing, and note that `PLAN.md` section 2 and the README both describe C2
as `(fwd + rev) >= 10` alt reads, matching the code, not the comment.
**Evidence:** Pooled over the eight figure samples: `alt>=4 & depth>=10` gives
108,993 pairs, min depth 10, min AF 3.77%, 1 below the 5% cutoff;
`alt>=10 & depth>=10` gives 98,704, min AF 6.97%, 0 below. All ten samples
re-run at exit 0 with no warnings.

---

## M20 - "25 - 21 = 4" is not four low-AF variants

**Date:** 2026-09-23
**Finding, recorded because the inference is natural and wrong:** in
GSE279945 the scMOCHA call arm holds 25 variants and the AF>5% arm 21, but the
four dropped are not four variants below 5%. Three are **blacklisted
positions** dropped at high AF - `16192C>T` at 99.1% prevalence (mis-alignment
window 16182-16194), `2617A>G` at 37.5% and `2617A>T` at 44.3% (RNA-editing
list) - and only `385A>G` is dropped for cell count, with 7 qualifying cells
against the 10 the gate requires.
**Why it matters:** `385A>G` is also a clean illustration of `M14`: its
prevalence AF is **0.065%** and its carrier AF is **96.9%**. Whether it counts
as a low-AF variant depends entirely on which measure is quoted, which is why
`07f` draws both and neither is presented alone.
**Consequence:** the difference between the call arm and the AF>5% arm cannot
be read as an AF statement. It is `blacklist + cell count`, and
`03-exclusion-reasons.tsv` carries the split per sample.
