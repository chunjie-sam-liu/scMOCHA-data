# Reading this stage

A criterion-by-criterion comparison of three mitochondrial variant-calling
arms across **eight samples**, built to answer the Cell Metabolism editorial
concern in `EDITOR.md`:

> "...concerns regarding the reliance on mgatk, which is biased toward
> detecting higher-heteroplasmy mutations ... and missing many mutations that
> would be tolerated at lower heteroplasmy."

Start here, then open `DIAGRAM.md` for the Mermaid version of the three arms.
`AGENTS.md` says how to run the code. The active campaign is
`2026-09-17-multi-sample-comparison`: read its `.PROGRESS.md`, then
`.DECISION.md`, then `.PLAN.md`. The stage-level `PLAN.md` and `DECISION.md`
hold the original single-sample design and the criterion audit.

---

## 1. The one-paragraph answer

The manuscript does not use original mgatk. It uses a modified caller whose
criteria differ at four points, and the difference the editor is describing is
mgatk's post-hoc gate, `vmr > 0.01 AND strand correlation > 0.65`. **scMOCHA
does not apply it at all.**

Across eight samples that gate does not merely bias the output, it eliminates
it: **original mgatk retains zero variants in seven of the eight samples**,
while scMOCHA reports 1 to 736. In the one deep sample where mgatk does return
a set, the variants its gate discards sit at median carrier heteroplasmy 0.069
against 0.205 for the ones it keeps (Wilcoxon P = 0.0029), and of the 171
variants scMOCHA AF>5% reports and mgatk does not, **149 (87%) are lost to the
strand-correlation cutoff alone**.

Running the other way: scMOCHA requires 10 alt reads in a cell where mgatk
requires 4, which is why scMOCHA's call set is smaller at the S1 stage in the
deep sample (736 vs 1,247). The variants only mgatk reports there carry a
median of **9** alt reads per carrier cell against **48** for the shared ones,
so what that criterion removes is weak read support, not low heteroplasmy.

---

## 2. The eight samples

Read from `tables/cross-sample/07-sample-overview.tsv`. **Fig** marks the six
samples the cross-sample panels show; the two 5' libraries are excluded from
those panels and kept in every table.

| Sample | Chemistry | Fig | Cells | Dropped by mgatk | S1 scMOCHA | S1 mgatk | **mgatk final** | scMOCHA call | scMOCHA AF>5% |
| --- | --- | :-: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| GSE149689_GSM4509019 | SC3Pv3 | yes | 721 | 389 (54%) | 20 | 19 | **0** | 20 | 20 |
| GSE155673_GSM4712895 | SC3Pv3 | yes | 5,805 | 4,381 (75%) | 6 | 19 | **0** | 6 | 6 |
| GSE163314_GSM4976997 | SC3Pv2 | yes | 7,949 | 7,673 (97%) | 9 | 4 | **0** | 9 | 9 |
| GSE163668_GSM4995445 | SC5P-R2 | no | 191 | 21 (11%) | 20 | 21 | **0** | 20 | 18 |
| GSE181279_GSM5494116 | SC5P-PE | no | 7,210 | 486 (6.7%) | 736 | 1,247 | **230** | 736 | 216 |
| GSE188632_GSM5687372 | SC3Pv3 | yes | 17,919 | 16,099 (90%) | 24 | 32 | **0** | 24 | 22 |
| GSE220189_GSM6793474 | SC3Pv3 | yes | 5,639 | 5,018 (89%) | 1 | 6 | **0** | 1 | 1 |
| GSE271107_GSM8369876 | SC3Pv3 | yes | 8,645 | 7,155 (83%) | 2 | 6 | **0** | 2 | 0 |

GSE181279 is the only deep sample and the only one that can carry a
distribution comparison. The other seven are shallow 10x libraries where 11% to
97% of cells fail mgatk's own coverage filter. **Treat the two groups as
answering different questions** - see section 8.

### Why mgatk returns nothing in seven samples

Not missing data. Every one of the 107 mgatk S1 variants in those seven samples
fails the strand-correlation floor:

| Sample | mgatk S1 | median strand r | max strand r | reaching r > 0.65 |
| --- | ---: | ---: | ---: | ---: |
| GSE149689 | 19 | 0.078 | 0.569 | 0 |
| GSE155673 | 19 | -0.062 | 0.196 | 0 |
| GSE163314 | 4 | -0.106 | 0.006 | 0 |
| GSE163668 | 21 | 0.007 | 0.590 | 0 |
| GSE181279 | 1,247 | 0.229 | 1.000 | 300 |
| GSE188632 | 32 | 0.102 | 0.503 | 0 |
| GSE220189 | 6 | -0.089 | 0.013 | 0 |
| GSE271107 | 6 | -0.206 | 0.048 | 0 |

The highest strand correlation seen anywhere in the seven shallow samples is
**0.590**, and many values are negative. Both `vmr` and `strand_correlation`
are populated, so this is a real rejection, not an `NA` artefact.

---

## 3. The three arms

| Arm | What it is |
| --- | --- |
| **Original mgatk** | cell filter, then the VMR/strand gate |
| **scMOCHA variant call** | what the caller emits; no AF filter, no VMR, no strand filter |
| **scMOCHA AF>5%** | the call plus the downstream gate in `src/06.1-collect-variants-new.R` |

Arm 3 is a strict subset of arm 2. The 5% AF rule belongs to the downstream
analysis, not to variant calling, which is why it is a separate arm.

### Where they differ

| Step | Original mgatk | scMOCHA variant call | scMOCHA AF>5% |
| --- | --- | --- | --- |
| Cell inclusion | mean MT cov > 10 | none | none |
| Confident cell | fwd >= 2 AND rev >= 2 | plus fwd+rev >= 10 | same |
| Variant retention | n_cells_conf >= 3 | same | same |
| Reliability gate | vmr > 0.01 AND strand r > 0.65 | **none** | blacklist, then >= 10 cells at AF >= 0.05 with depth >= 10 |

---

## 4. Cross-sample figures

5 PDFs in `figures/cross-sample/`. **These carry the eight-sample claim.**

**The panels show six samples; the tables show all eight.** The two 5'
libraries, `GSE163668_GSM4995445_5PR2` and `GSE181279_GSM5494116_5PPE`, are
excluded from the cross-sample panels so those panels compare within one
library family. The list is `CROSS_SAMPLE_FIG_EXCLUDE` in `config.R`.

**This matters when reading the panels.** GSE181279 is the only sample in which
mgatk retains anything at all and the only one in which the gate test is
computable, and it is one of the two excluded. So the panels alone show mgatk
retaining nothing anywhere and no testable comparison; that is a property of
the six samples drawn, not of the resource. Read the tables alongside them.

| Figure | Shows | Read it as |
| --- | --- | --- |
| **`07a-arm-yield`** | variants retained per arm per sample, linear axis | the headline: mgatk 0 in all six panels shown, and in seven of the eight samples overall |
| `07b-cell-filter` | fraction of cells mgatk's coverage filter discards, per sample | 6.7% in the deep sample, 54% to 97% across the shallow ones |
| `07c-strand-correlation` | distribution of mgatk strand correlation per sample, floor drawn | **descriptive only** - see the caveat below |
| `07d-exclusion-reasons` | which mgatk cutoff rejects each S1 variant, per sample; **retained is red**, every exclusion reason a non-red hue | strand correlation is implicated in every rejection in all eight samples |
| **`07e-gate-rejected-af`** | carrier AF of gate-rejected vs gate-passed variants, per sample | testable in GSE181279 only; the other four have no passed group |

**Samples are labelled by GSE and GSM.** Two of the five share `SC3Pv3`, so
chemistry cannot identify a sample and is not used on any cross-sample panel or
table. Every cross-sample table carries `sample` as `{gse}_{gsm}` plus a
`sample_id` column that joins back to the per-sample directories; chemistry is
in the `SAMPLES` registry and the workbook's `00_Samples` sheet.

**`07a` is on a linear axis on purpose.** A zero bar and its label both vanish
silently on a log scale, and four zeros are the headline of this figure.

**The `07c` caveat.** The plan expected this panel to show the strand floor
tracking read depth, which would explain the bias mechanically. **The data do
not support that.** In the only sample large enough to test it, Spearman rho
between strand correlation and per-variant coverage is **-0.013, P = 0.66**
(`07-strand-support.tsv`). The two small samples where rho is positive have
n = 19 and n = 21. So `07c` documents *that* the floor rejects everything in
shallow data; it does **not** establish *why*. Do not quote it as a mechanism.
See `M4` in the campaign decision log.

**The `07e` caveat.** Seven of the eight samples have zero variants passing
mgatk's gate, so the split has one level and the Wilcoxon is undefined. Those
panels show the rejected group only and are labelled as having no comparison
group. The one testable sample, GSE181279, is excluded from the panels, so
`07e` as drawn contains no computable test at all; `07-gate-test.tsv` carries
it. There is **no pooled or combined test across samples** - with one sample
contributing both groups, a stratified test collapses to that sample and a
combined p-value would be one p-value wearing a meta-analysis label (`M2`).

| Sample | Passed | Rejected | Median AF passed | Median AF rejected | P |
| --- | ---: | ---: | ---: | ---: | ---: |
| **GSE181279** (not in panels) | 54 | 215 | **0.205** | **0.069** | **0.0029** |
| all seven others | 0 | 0 to 20 | - | - | not testable |

---

## 5. How heteroplasmy is measured, and why it matters

The single thing most likely to be misread. Every measure below is a "median
AF across the cells carrying the variant". They differ only in **which cells
count as carriers**, and that choice changes the answer.

| Column | Carrier definition | Status |
| --- | --- | --- |
| `mean` (either caller) | all cells | measures **prevalence**, not heteroplasmy |
| `af_gate_*` | cell AF >= 0.05 and depth >= 10 | **censored at 0.05** |
| `af_loose_*` | depth >= 10 and >= 2 alt reads | **background-dominated** |
| `af_carrier_*` | alt count inconsistent with the background error rate | **adopted** |

`mean` is total alt reads over total coverage across every cell, so a variant
at 60% heteroplasmy in 12 cells scores below 0.002.

`af_gate_*` cannot see below 0.05, because 0.05 is the floor inside its own
carrier definition. Using it for an AF distribution deletes the region under
study.

`af_loose_*` removes the AF floor but replaces it with one so weak that
background reads take over. In GSE181279 the measured background alt rate is
**0.026%**, so two alt reads in a cell of depth 1,000 is expected by chance.
Variant `1801A>G` has 19 genuine carriers above 5% and 876 background cells;
its `af_loose_median` is 0.0017, which is the background rate, not the variant.

`af_carrier_*` keeps a cell when its alt count is inconsistent with that
background rate under a binomial test, Bonferroni-corrected across all
cell-by-variant observations (9,005,290 of them in GSE181279). The threshold
scales with depth, so it excludes background without imposing an AF floor.

**A shallow sample can have no qualifying background cell at all**, which makes
the rate `0 / 0`. It is floored at one alt read over all observed depth rather
than allowed to reach 1, which would make every cell a non-carrier and silently
delete the sample from the gate table. This was a real defect, caught on
GSE163314; see `M5` in the campaign decision log.

**`04f-measure-sensitivity` reports the headline test under all four
definitions** (GSE181279), so the choice is auditable rather than asserted:

| Carrier definition | median rejected | median passed | P |
| --- | --- | --- | --- |
| AF >= 0.05 (censored) | 0.103 | 0.107 | 0.59 |
| >= 2 alt reads (background-dominated) | 0.019 | 0.053 | 0.043 |
| **error model (adopted)** | **0.069** | **0.205** | **0.0029** |
| max across error-model carriers | 0.82 | 1.00 | 0.0036 |

The direction is identical in all four. Only the censored measure fails to
resolve it, and it fails for a structural reason rather than a biological one.
Both wrong measures were adopted first and both returned a clean null; see
`D12` and `D17` in the stage `DECISION.md`.

### Why an AF>5% variant can still have a carrier median below 0.05

The two rules count different cells. The gate asks for **at least 10** cells
above 5%; the carrier median is taken over **every** cell that genuinely
carries the variant, including the ones below 5%. A variant with 19 cells at
7% and 40 more real carriers at 1-3% passes the gate and has a median under
0.05. That is a real pattern, not an artifact - and it is the pattern the
editor's question is about. Six of the 216 AF>5% variants sit below 0.01 on
this measure; all six have a carrier maximum between 0.09 and 0.53.

---

## 6. Per-sample figures

24 PDFs in `figures/<sample_id>/`, the same set for every sample. Panel letters
follow the step that produced them. **Every number quoted below is GSE181279**,
the only sample where these panels are well populated; the same panels exist
for the other four and are mostly near-empty by construction.

In the seven shallow samples `03b`, `03d`, `03e`, `04e` and `04f` are drawn
with an explicit "no variants in this arm" annotation rather than left blank,
so a missing comparison is distinguishable from a failed run.

### 02 - what the cell filter costs (criteria C1, C2)

| Figure | Shows | Read it as |
| --- | --- | --- |
| `02a-cell-coverage` | per-cell mean MT coverage, split at the mgatk cutoff of 10 | 486 of 7,210 cells (6.7%) fall below the line and never enter mgatk's call |
| `02b-cells-retained` | cells entering variant calling | 6,724 vs 7,210 |
| `02c-min-detectable-af` | lowest AF a single cell can support, `4/depth` for mgatk against `10/depth` for scMOCHA, over the observed coverage range | the cost of scMOCHA's stricter per-cell rule, stated plainly: its floor is higher in any single cell |
| `02d-detections-lost` | detections and variants lost with the discarded cells | 4,130 of 127,998 detections and 14 variants |

`02c` is the panel that argues against scMOCHA. It is here on purpose.

### 03 - the funnel and where the sets diverge (C3, C4, C6)

| Figure | Shows | Read it as |
| --- | --- | --- |
| `03a-variant-funnel` | S0 -> S1 -> S2 for all three arms, log y | mgatk 25,580 -> 1,247 -> 230; scMOCHA 25,746 -> 736, then -> 216 with the AF>5% gate |
| `03b-variant-overlap` | three-way Venn | arm 1 vs arm 3: 45 shared, 185 mgatk only, 171 AF>5% only |
| `03c-gate-crossapplied` | each reliability gate applied to each S1 set | holding S1 fixed, the scMOCHA gate keeps more: 269 vs 230 on mgatk's S1, 216 vs 164 on scMOCHA's |
| **`03d-excluded-from-mgatk`** | for the 171 variants mgatk misses, the cutoff responsible | **149 strand r alone, 17 both, 5 VMR alone** |
| `03e-excluded-from-scmocha-af5` | for the 185 variants scMOCHA AF>5% misses, the cutoff responsible | 117 fail the 10-cell rule, 68 fail the 10-alt-read rule |
| `03f-arm-combinations` | every membership combination | 403 call only, 171 call + AF>5%, 117 mgatk + call, 68 mgatk alone, 45 all three |

`03d` and `03e` are the pair that answers "why do the two sets differ". Each
variant is attributed to the **first** criterion it fails in that arm's order,
so the bars sum to the set size and nothing is double-counted. A variant
failing both VMR and strand correlation gets its own category rather than being
folded into either, which is what makes the 87% figure trustworthy.

### 04 - the heteroplasmy spectrum (the editor's question)

| Figure | Shows | Read it as |
| --- | --- | --- |
| `04a-af-ecdf` | ECDF of carrier heteroplasmy per variant class | the three classes' distributions, uncensored |
| `04b-max-heteroplasmy` | density of maximum carrier AF | whether any class is confined to low heteroplasmy |
| `04c-af-bins` | AF-bin composition of all three arms | the raw scMOCHA call is broadest (70% below AF 0.05), mgatk sits between (48%), and the AF>5% gate concentrates its set on higher heteroplasmy (30% below 0.05, 34% at >= 0.50) |
| `04d-cell-level-af` | density of cell-level AF, both sets | the sub-0.05 tail, visible because there is no AF floor |
| **`04e-gate-rejected-af`** | carrier AF of variants mgatk's gate rejects vs keeps, among variants scMOCHA considers reliable | **median 0.069 rejected vs 0.205 kept, Wilcoxon P = 0.0029** |
| **`04f-measure-sensitivity`** | the same test under four carrier definitions | the direction holds in all four; only the censored measure cannot resolve it |
| `04g-af-bins-own-definition` | `04c` with each arm on its own AF definition: AF>5% arm on cells at AF >= 0.05, the other two on all carrier cells | what each arm delivers on its own terms; **not** a like-for-like comparison |

`04e` is the direct test. Both arms are variants scMOCHA already calls
reliable, so the only thing differing between them is mgatk's decision.
`04f` is what makes `04e` trustworthy: it shows the result is a property of the
data, not of the measure that was picked.

**`04c` against `04g`.** They answer different questions. `04c` uses one
measure across all three arms, so the bars are comparable. `04g` measures each
arm on the rule that defines it - the 5% AF rule belongs to the AF>5% arm
alone, so only that arm is summarised on cells at AF >= 0.05. The consequence
is that the three bars in `04g` do **not** share a measure: the AF>5% bar
cannot fall below 0.05 while the other two can, so part of its apparent shift
is definitional rather than biological. The `measured` column of
`04-af-bins-own-definition.tsv` records which measure each arm used. Quote
`04c`, not `04g`, for any statement comparing arms.

### 05 - inside the VMR / strand-correlation plane (C4 in detail)

Both callers **compute** VMR and strand correlation; only mgatk **acts** on
them. That is what makes this plane readable.

| Figure | Shows | Read it as |
| --- | --- | --- |
| `05a-vmr-strand-plane` | the plane per rule set, coloured by carrier AF, mgatk's cutoffs drawn on both | where each heteroplasmy band sits relative to the gate |
| `05b-gate-rejected-af-bins` | AF-bin composition of what the gate keeps vs discards | rejection is 90.5% in the `<0.01` bin (76/84) against 58% in `0.05-0.10` (11/19) |
| **`05d-arm-in-mgatk-plane`** | mgatk's own coordinates, coloured by which arm reports the variant | orange sits in the upper-right quadrant by construction; the informative points are the **dark blue scMOCHA AF>5% variants to the left of the strand r = 0.65 line but above the VMR line** - high VMR, discarded purely for strand correlation |
| `05e-arm-plane-facets` | `05d` split one facet per arm, all variants repeated in grey | the same, easier to compare densities |
| `05c-read-support` | alt reads per carrier cell by variant class | mgatk-only 9, scMOCHA-only 38, shared 48 |
| **`05f-arm-in-scmocha-plane`** | the mirror of `05d`: **scMOCHA's** decision plane, x = cells at AF >= 0.05 with depth >= 10, y = median alt reads per carrying cell | where **mgatk's** variants land under scMOCHA's criteria: **176 of the 185 mgatk-only variants fall short of the 10-cell gate** |
| `05g-scmocha-plane-facets` | `05f` split one facet per arm, all variants repeated in grey | the readable version of `05f`; see the caveat below |

`05d` is the figure to put next to `03d`: one shows the count, the other shows
the geometry behind it.

`05d` and `05e` drop variants with no mgatk `vmr`, so "scMOCHA AF>5% only" is
168 there against 171 in `03b`. The three missing variants are ones mgatk never
proposed as candidates. `05f` and `05g` use no mgatk coordinate, so all 171 are
present there.

**`05f` and `05g` are the diagnostic pair.** `05d` asks where scMOCHA's
variants sit in mgatk's plane; `05f` asks the reverse, which is the question to
put to anyone who wants to know what scMOCHA's criteria do to mgatk's output.
Arm counts reconcile with `03b` exactly: 45 both arms, 185 mgatk only (173
short on both axes, 9 on read support alone, 3 on the cell gate alone), and 171
AF>5% only, every one of them right of the cell gate by construction.

**Read the y axis carefully.** The x axis *is* scMOCHA's reliability gate, so
that vertical line is a real cutoff. The y axis is only a **summary** of the
read support the per-cell rule acts on - that rule asks for >= 3 individual
cells each at >= 10 alt reads, which a median cannot express. The horizontal
line is orientation, not a criterion; `03e` carries the attribution. `05g`
exists because the median puts most low-support variants on y = 1 or 2, where
the arms occlude each other in `05f`, and jittering a median would move points
off their own value.

**All four plane panels show only the three arm-specific classes.** `05d`,
`05e`, `05f` and `05g` drop `Both arms` and `Neither arm`: neither separates
the callers, and `Neither arm` is the largest group, so leaving it in draws a
grey mass over the points the panels exist to show. Both classes keep their
counts in `05-arm-in-mgatk-plane.tsv` and `05-arm-in-scmocha-plane.tsv`, which
are built from the unfiltered data.

---

## 7. Tables

**14 per sample** in `tables/<sample_id>/`, **7 cross-sample** in
`tables/cross-sample/`. Every table carries a leading `sample` column.

`tables/cross-sample/06-compare-with-mgatk.xlsx` collects all of them into
**23 sheets** and is the one to send out. `00_Samples`, `01_Criteria` and
`02_Definitions` make it self-contained.

| Cross-sample file | Contents |
| --- | --- |
| `07-sample-overview.tsv` | the registry plus every headline count per sample |
| `07-arm-yield.tsv` | variants retained per arm per sample |
| `07-cell-filter.tsv` | cells, detections and variants lost to the cell filter |
| `07-strand-support.tsv` | strand-correlation summary and the rho that failed to show a mechanism |
| `07-exclusion-reasons.tsv` | which mgatk cutoff rejects each S1 variant |
| `07-gate-test.tsv` | the gate test per sample, with a `testable` flag |
| `06-criteria-comparison.tsv` | the criterion table, machine-readable |

| Per-sample file | Contents |
| --- | --- |
| `02-cell-inclusion.tsv` | cell counts, detections, variants lost to the cell filter |
| `03-funnel-counts.tsv` | the three-arm funnel |
| `03-gate-crossapplied.tsv` | each gate on each S1 set |
| `03-exclusion-reasons.tsv` | the numbers behind `03d` and `03e` |
| `03-arm-combinations.tsv` | membership combinations |
| `03-variant-membership.tsv` | one row per variant: arm flags, exclusion reason, carrier AF, mgatk vmr/strand |
| `04-af-bins.tsv` | AF composition per arm |
| `04-af-bins-own-definition.tsv` | the same, each arm on its own AF definition, with a `measured` column |
| `04-tests.tsv` | the Wilcoxon test with CI and the background rate |
| `04-measure-sensitivity.tsv` | the same test under four carrier definitions |
| `05-vmr-strand-rejected.tsv` | gate outcome by AF bin |
| `05-arm-in-mgatk-plane.tsv` | arm counts in `05d` |
| `05-arm-in-scmocha-plane.tsv` | arm counts in `05f`, split by which scMOCHA axis the variant clears |
| `05-read-support.tsv` | median alt reads per carrier cell |

`03-variant-membership.tsv` is the lookup table: to ask why one specific
variant appears in one arm and not another, find it there.

---

## 8. What this stage does not claim

- **No sensitivity, recall, or false-negative rate.** There is no orthogonal
  truth set - no matched bulk mtDNA-seq, no simulation - so every result is
  criterion-to-criterion on identical reads. Any recall number would use one
  caller as a stand-in for truth, which is circular (`D4`).
- **No mechanism for the strand floor.** rho = -0.013, P = 0.66 in the only
  sample large enough to test it, so `07c` is descriptive only (`M4`).
- **No pooled or combined test across samples.** Seven of the eight have no
  comparison group, so a stratified test collapses to GSE181279 (`M2`).
- **The seven shallow samples do not answer the editor's question.** They answer
  a simpler and different one: original mgatk returns nothing at all on them.
  Do not present the two results as the same result.
- **P = 0.0029 rests on one sample.** Stable across carrier definitions, but
  treat it as a demonstration of direction, not an established effect size.
- **Neither caller was re-run.** Both output sets already existed and came from
  the same allele-count matrices, so the two differ only in filtering logic.
  The archives come from Ting; the pipeline versions and written confirmation
  that both arms used identical allele counts in all eight samples are still
  outstanding (`DATA.md`, Unknowns).

---

## 9. Rebuilding

```bash
cd "$(git rev-parse --show-toplevel)"
bash newplots/compare-with-mgatk/run-all.sh
```

That extracts the archives if they are not already extracted, runs steps 01 to
05 for each of the eight samples, then `07-cross-sample.R` and
`06-summary-workbook.R`. To rebuild one sample:

```bash
pixi run Rscript newplots/compare-with-mgatk/01-load-harmonize.R \
  --sample=GSE181279_GSM5494116_5PPE
```

Order matters: 02-05 read step 01's cache, 04 and 05 read step 03's, 07 reads
the per-sample tables, and 06 must run last. Per-step contracts are in the
paired `NN-*.md` files.
