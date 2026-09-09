# Reading this stage

A criterion-by-criterion comparison of three mitochondrial variant-calling
arms on one sample, built to answer the Cell Metabolism editorial concern in
`EDITOR.md`:

> "...concerns regarding the reliance on mgatk, which is biased toward
> detecting higher-heteroplasmy mutations ... and missing many mutations that
> would be tolerated at lower heteroplasmy."

Start here, then open `DIAGRAM.md` for the Mermaid version of the three arms.
`AGENTS.md` says how to run the code, `PLAN.md` the design intent,
`PROGRESS.md` the current state, `DECISION.md` why each choice was made.

Sample: 7,210 cells, median per-cell MT coverage 38.9. Identity still pending;
figure subtitles read `SAMPLE_LABEL` from `config.R` and currently say
`<pending>`.

---

## 1. The one-paragraph answer

The manuscript does not use original mgatk. It uses a modified caller whose
criteria differ at four points. Two of them bear on the editor's concern, and
one runs the other way.

The filter the editor is describing is mgatk's post-hoc gate,
`vmr > 0.01 AND strand correlation > 0.65`. **scMOCHA does not apply it at
all.** Of the 171 variants that scMOCHA AF>5% reports and mgatk does not,
**149 (87%) are lost to the strand-correlation cutoff alone** - not to VMR, not
to read depth. Among variants both pipelines consider reliable, the ones
mgatk's gate discards sit at median carrier heteroplasmy 0.069 against 0.205
for the ones it keeps (Wilcoxon P = 0.0029).

Running the other way: scMOCHA requires 10 alt reads in a cell where mgatk
requires 4, which is why scMOCHA's call set is smaller at the S1 stage
(736 vs 1,247). The variants only mgatk reports carry a median of **9** alt
reads per carrier cell against **48** for the shared ones, so what that
criterion removes is weak read support, not low heteroplasmy.

---

## 2. The three arms

| Arm | What it is | Variants |
| --- | --- | --- |
| **Original mgatk** | cell filter, then the VMR/strand gate | **230** |
| **scMOCHA variant call** | what the caller emits; no AF filter, no VMR, no strand filter | **736** |
| **scMOCHA AF>5%** | the call plus the downstream gate in `src/06.1-collect-variants-new.R` | **216** |

Arm 3 is a strict subset of arm 2. The 5% AF rule belongs to the downstream
analysis, not to variant calling, which is why it is a separate arm.

### Where they differ

| Step | Original mgatk | scMOCHA variant call | scMOCHA AF>5% |
| --- | --- | --- | --- |
| Cell inclusion | mean MT cov > 10 (6,724 / 7,210) | none (7,210) | none |
| Confident cell | fwd >= 2 AND rev >= 2 | plus fwd+rev >= 10 | same |
| Variant retention | n_cells_conf >= 3 | same | same |
| Reliability gate | vmr > 0.01 AND strand r > 0.65 | **none** | blacklist, then >= 10 cells at AF >= 0.05 with depth >= 10 |

---

## 3. How heteroplasmy is measured, and why it matters

The single thing most likely to be misread. Every measure below is a "median
AF across the cells carrying the variant". They differ only in **which cells
count as carriers**, and that choice changes the answer.

| Column | Carrier definition | Status |
| --- | --- | --- |
| `mean` (either caller) | all 7,210 cells | measures **prevalence**, not heteroplasmy |
| `af_gate_*` | cell AF >= 0.05 and depth >= 10 | **censored at 0.05** |
| `af_loose_*` | depth >= 10 and >= 2 alt reads | **background-dominated** |
| `af_carrier_*` | alt count inconsistent with the background error rate | **adopted** |

`mean` is total alt reads over total coverage across every cell, so a variant
at 60% heteroplasmy in 12 cells scores below 0.002.

`af_gate_*` cannot see below 0.05, because 0.05 is the floor inside its own
carrier definition. Using it for an AF distribution deletes the region under
study.

`af_loose_*` removes the AF floor but replaces it with one so weak that
background reads take over. The measured background alt rate is **0.026%**, so
two alt reads in a cell of depth 1,000 is expected by chance. Variant
`1801A>G` has 19 genuine carriers above 5% and 876 background cells; its
`af_loose_median` is 0.0017, which is the background rate, not the variant.

`af_carrier_*` keeps a cell when its alt count is inconsistent with that
background rate under a binomial test, Bonferroni-corrected across all
9,005,290 cell-by-variant observations. The threshold scales with depth, so it
excludes background without imposing an AF floor.

**`04f-measure-sensitivity` reports the headline test under all four
definitions**, so the choice is auditable rather than asserted:

| Carrier definition | median rejected | median passed | P |
| --- | --- | --- | --- |
| AF >= 0.05 (censored) | 0.103 | 0.107 | 0.59 |
| >= 2 alt reads (background-dominated) | 0.019 | 0.053 | 0.043 |
| **error model (adopted)** | **0.069** | **0.205** | **0.0029** |
| max across error-model carriers | 0.82 | 1.00 | 0.0036 |

The direction is identical in all four. Only the censored measure fails to
resolve it, and it fails for a structural reason rather than a biological one.
Both wrong measures were adopted first and both returned a clean null; see
`D12` and `D17` in `DECISION.md`.

### Why an AF>5% variant can still have a carrier median below 0.05

The two rules count different cells. The gate asks for **at least 10** cells
above 5%; the carrier median is taken over **every** cell that genuinely
carries the variant, including the ones below 5%. A variant with 19 cells at
7% and 40 more real carriers at 1-3% passes the gate and has a median under
0.05. That is a real pattern, not an artifact - and it is the pattern the
editor's question is about. Six of the 216 AF>5% variants sit below 0.01 on
this measure; all six have a carrier maximum between 0.09 and 0.53.

---

## 4. Figure guide

20 PDFs in `figures/`. Panel letters follow the step that produced them.

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

`05d` is the figure to put next to `03d`: one shows the count, the other shows
the geometry behind it.

`05d` and `05e` drop variants with no mgatk `vmr`, so "scMOCHA AF>5% only" is
168 there against 171 in `03b`. The three missing variants are ones mgatk never
proposed as candidates.

---

## 5. Tables

13 files in `tables/`. `06-compare-with-mgatk.xlsx` collects all of them into
13 sheets and is the one to send out; `00_Criteria` and `01_Definitions` make
it self-contained.

| File | Contents |
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
| `05-read-support.tsv` | median alt reads per carrier cell |
| `06-criteria-comparison.tsv` | the criterion table, machine-readable |

`03-variant-membership.tsv` is the lookup table: to ask why one specific
variant appears in one arm and not another, find it there.

---

## 6. What this stage does not claim

- **No sensitivity, recall, or false-negative rate.** There is no orthogonal
  truth set - no matched bulk mtDNA-seq, no simulation - so every result is
  criterion-to-criterion on identical reads. Any recall number would use one
  caller as a stand-in for truth, which is circular.
- **One sample.** The editor's criticism is about the whole resource.
  P = 0.0029 is stable across carrier definitions but rests on a single sample;
  treat it as a demonstration of direction, not as an established effect size.
- **Neither caller was re-run.** Both output sets already existed and came from
  the same allele-count matrices, so the two differ only in filtering logic.

---

## 7. Rebuilding

```bash
cd "$(git rev-parse --show-toplevel)"
for s in 01-load-harmonize 02-cell-inclusion 03-variant-funnel-overlap \
         04-heteroplasmy-spectrum 05-vmr-strand 06-summary-workbook; do
  pixi run Rscript "newplots/compare-with-mgatk/${s}.R" || break
done
```

Order matters: 02-05 read step 01's cache, 04 and 05 read step 03's, and 06
reads the `.tsv` files the others wrote. The whole stage takes about a minute.
Per-step contracts are in the paired `NN-*.md` files.
