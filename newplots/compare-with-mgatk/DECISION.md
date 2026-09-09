# DECISION: scMOCHA-updated mgatk vs original mgatk variant calling

Stage: `newplots/compare-with-mgatk/`
Append-only. Each entry: what was chosen, what it beat, why, evidence, date.

---

## D1 — Compare at the criterion level, not by re-running the callers

**Date:** 2026-09-09
**Chosen:** Use the two existing output sets in
`/research/rgs01/home/clusterHome/cliu68/project/scmocha/compare-with-mgatk`
and compare the filtering criteria applied to the same underlying allele
counts.
**Beat:** Re-running original mgatk from BAMs to regenerate the comparison.
**Why:** Both output sets already exist and were produced from the same
`cell.{A,T,C,G}.txt.gz` allele-count matrices, so the two callers differ only
in their filtering logic. Re-running adds cluster cost and introduces version
drift without adding information.
**Evidence:** Both `cell.variant_stats.tsv.gz` (25,746 rows) and
`cell.variant_stats_mgatk_original.tsv.gz` (25,580 rows) carry the identical
14-column schema and share the same candidate-variant construction step; the
counts differ only where the criteria differ.

---

## D2 — Report the criterion that makes scMOCHA more stringent, not only the two that make it more permissive

**Date:** 2026-09-09
**Chosen:** Include criterion C2 (`fwd + rev alt reads >= 10`) in the
comparison, with a dedicated panel (02c) showing its cost as a raised minimum
detectable AF per cell, and a test (05c) of whether what it removes is
low-read-support or genuinely low-heteroplasmy.
**Beat:** Presenting only C1 (no cell filter) and C4 (relaxed strand
correlation), the two changes that favour the rebuttal.
**Why:** Post-C3 variant counts are 736 for scMOCHA versus 1,247 for original
mgatk. A reviewer will compute that ratio immediately. A comparison that does
not explain it is not credible, and the explanation is itself a strength if the
removed calls are read-support artefacts.
**Evidence:** Verified counts from
`cell.cell_heteroplasmic_df.tsv.gz` (7,210 x 736) and
`cell.cell_heteroplasmic_df_mgatk_original.tsv.gz` (6,724 x 1,247).

---

## D3 — Treat the VMR + strand-correlation filter as the locus of the editor's criticism

**Date:** 2026-09-09
**Chosen:** Make step 04e / 05 the load-bearing analysis: the AF distribution
of the variants that mgatk's `vmr > 0.01 & strand_correlation > 0.65` gate
rejects but scMOCHA's reliability gate keeps.
**Beat:** Framing the rebuttal around the cell-inclusion filter (C1), or around
the downstream homoplasmic/heteroplasmic/somatic classification (C7).
**Why:** The editor's specific claim is a bias toward higher-heteroplasmy
variants. Among the criteria that differ, VMR combined with a high
strand-correlation floor is the one with a mechanistic link to that bias: a
low-heteroplasmy variant contributes few alt reads per cell, so its
forward/reverse read correlation is noise-dominated and it fails a 0.65 floor
regardless of whether it is real, while its across-cell variance is small so it
also fails `vmr > 0.01`. C1 removes cells, not preferentially low-AF variants,
and C7 is a different question entirely.
**Evidence:** The 0.65 and 0.01 lines are hard-coded in
`newplots/thecode/original-mgatk-variant-calling.py`.
**Status:** confirmed by D11.

---

## D11 — scMOCHA applies no VMR and no strand-correlation filter (corrects the plan)

**Date:** 2026-09-09
**Chosen:** State criterion C4 as "mgatk applies `vmr > 0.01 &
strand_correlation > 0.65`; scMOCHA applies neither", and represent scMOCHA's
reliability gate as C6: position blacklist, then >= 10 cells at cell AF >= 0.05
with cell depth >= 10.
**Beat:** The earlier draft of this plan, which claimed scMOCHA uses
`vmr > 0.01 & strand_correlation > 0.3`.
**Why:** The user stated that `02-update-high-quality-variant.py` is not part of
the variant-calling path, and the code confirms it. Presenting a 0.3 cutoff as
a scMOCHA calling criterion would have been factually wrong in a response to an
editor, and would have weakened the rebuttal: "we lowered the threshold" is a
much weaker answer than "we do not apply that filter at all".
**Evidence:** `grep 'strand|vmr' src/06.1-collect-variants-new.R` returns no
filter — the only hits are `CUTOFF_HETEROPLASMIC` and the
`cluster.cell_heteroplasmic_df.tsv.gz` read. The production classification
applies `af >= 0.05 & depth >= 10` then `n >= 10` cells (`fn_n_cell`),
`total_reads >= 10` at cluster level (`fn_homo_hete`), and the
`POS_RNA_EDITING` / `POS_MISSALIGNMENT_ERROR` blacklists. Repository-wide grep
shows `high_quality_variant.tsv` is read only by
`preprocessing/high_quality_variant/03-collect-hqv.R` and `04-hqv-flextable.R`,
which build a summary flextable under `analysis/zzz/hqv`.
**Consequence:** Q3 in `PLAN.md` is superseded, not answered. Step 03's
cross-applied funnel and step 05's scatter now compare a filter that one caller
applies and the other does not, rather than two settings of the same filter.

---

## D4 — No orthogonal truth set, so no sensitivity claim

**Date:** 2026-09-09
**Chosen:** Frame every result as criterion-to-criterion on identical reads.
Do not report sensitivity, recall, or false-negative rate.
**Beat:** Claiming scMOCHA "recovers N% more true low-heteroplasmy variants".
**Why:** The intermediate directory contains no matched bulk mtDNA-seq, no
simulated data, and no external validated variant list. Any recall number would
be computed against one of the two callers used as a stand-in for truth, which
is circular.
**Evidence:** Directory listing of the intermediate root: only mgatk allele
counts, coverage, the two stats tables, the three AF matrices, and the two RDS
objects.
**Status:** provisional — reverses if matched bulk or simulated data exists
(Q6).

---

## D5 — Sample label carried as a single constant

**Date:** 2026-09-09 (Q1 approved as recommended)
**Chosen:** `SAMPLE_LABEL` in `config.R`, set to `"<pending>"`.
**Beat:** Blocking implementation until the GSE/GSM/SRR identifier is known.
**Why:** The identifier affects only figure captions, not any computation. One
constant means one edit relabels every output.

---

## D6 — Build on one sample, extend later

**Date:** 2026-09-09 (Q2 approved as recommended)
**Chosen:** Run the stage on the single sample that has paired output, and make
extension to more samples a separate dated campaign.
**Beat:** Making the scripts multi-sample from the start.
**Why:** Only one sample currently has a paired `*_mgatk_original.tsv.gz`.
Generalising ahead of that is speculative work against an unknown layout.

---

## D7 — Intermediates cached outside the git tree

**Date:** 2026-09-09 (Q4 approved as recommended)
**Chosen:** `CACHEDIR = <intermediate root>/derived/`.
**Beat:** `~/tmp/cmp-mgatk/`, which is cleared between sessions, and the repo
tree, which must not carry large binaries.
**Why:** The joined variant table and the position-filtered coverage cache are
too large for git and too expensive to rebuild every session.

---

## D8 — New colours appended to the existing track colour file

**Date:** 2026-09-09 (Q5 approved as recommended)
**Chosen:** Append `color_caller` and `color_af_bin` to
`high-res/00-colors.R`, additive only.
**Beat:** A new `newplots/color.R`.
**Why:** `newplots/` scripts already source `high-res/00-colors.R` through
`HIGHRESDIR`. A second colour file would split the source of truth. No existing
colour is modified, so manuscript figures are unaffected.

---

## D9 — Coverage extracted by position-filtered streaming

**Date:** 2026-09-09
**Chosen:** Filter `cell.coverage.txt.gz` to candidate variant positions with
`awk` during decompression, before loading into R.
**Beat:** `fread` on the whole file (about 1.2e8 rows) and subsetting in R.
**Why:** Only order 1e3 of the 16,569 positions are ever needed, so streaming
cuts the row count by roughly two orders of magnitude and keeps the step inside
a normal interactive memory budget.
**Note:** This became necessary only after D11 made scMOCHA's gate depend on
per-cell depth, which the AF matrices do not carry.

---

## D10 — Criterion constants copied verbatim from the production script

**Date:** 2026-09-09
**Chosen:** `config.R` copies `CUTOFF_HETEROPLASMIC`, `CUTOFF_MIN_READS`,
`CUTOFF_NOTRELIABLE`, `POS_RNA_EDITING`, and `POS_MISSALIGNMENT_ERROR`
verbatim from `src/06.1-collect-variants-new.R`, with a comment naming the
source.
**Beat:** Sourcing `src/06.1-collect-variants-new.R` directly, or re-deriving
the values from `src/06.4-somatic.md`.
**Why:** That script is not a config file; sourcing it would execute the whole
pipeline. The `.md` is documentation and can drift from the code. Copying from
the code with the source named makes any future drift visible in review.

---

## D12 -- Heteroplasmy measured on carrier cells with no AF floor

**Date:** 2026-09-09
**Chosen:** Measure heteroplasmy as `af_obs_median`, the median AF over cells
with depth >= 10 carrying >= 2 alt reads and **no AF floor**.
**Beat:** Two alternatives, both of which were implemented first and both of
which returned a null result.

1. The mgatk `mean` column (population AF). Rejected because it is total alt
   reads over total coverage across every cell, so a variant at 60%
   heteroplasmy in 12 of 7,210 cells scores below 0.002. It ranks variants by
   prevalence, not heteroplasmy.
2. `af_carrier_median`, the median AF over cells passing the scMOCHA detection
   rule (AF >= 0.05, depth >= 10). Rejected because that definition carries an
   AF floor of 0.05, so the measure is censored at exactly the boundary of the
   region the editor is asking about.

**Why it matters:** the choice inverted the headline result.

| Measure | Rejected by mgatk gate | Passes mgatk gate | P |
| --- | --- | --- | --- |
| population `mean` | 0.00138 (n=215) | 0.00140 (n=54) | 0.82 |
| `af_carrier_median`, floored at 0.05 | 0.103 | 0.107 | 0.59 |
| `af_obs_median`, no floor | **0.0191** | **0.0545** | **0.044** |

Both censored measures said mgatk's gate is unbiased with respect to
heteroplasmy. The uncensored measure shows the discarded variants sit at
roughly a third the heteroplasmy of the retained ones. Neither wrong version
errored, warned, or looked odd.

**Evidence:** `newplots/compare-with-mgatk/tables/04-tests.tsv`, and the
per-bin rejection counts in `05-vmr-strand-rejected.tsv`, where the lowest AF
bin is rejected at 76/84 = 90.5% against 11/19 = 58% in the 0.05-0.10 bin.
**Consequence:** `config.R` now carries both measures with the rule for each,
and `config.R.md` states it. Any new AF panel uses `af_obs_*`.

---

## D13 -- Report the counter-direction result rather than only the supporting one

**Date:** 2026-09-09
**Chosen:** Keep panel 05c and the `05-read-support.tsv` table, which show that
the variants only mgatk reports carry a median of 9 alt reads per carrier cell
against 48 for shared variants, and keep panel 02c, which shows scMOCHA's
per-cell detection floor is 10/depth against mgatk's 4/depth.
**Beat:** Presenting only the funnel and the gate test, which favour scMOCHA.
**Why:** scMOCHA's final set is smaller than mgatk's on this sample (216 vs
230) and its S1 set is much smaller (736 vs 1,247). A reviewer computes that
ratio immediately. The read-support panel is what turns it from a weakness into
an argument, and omitting the cost panel would make the comparison advocacy
rather than analysis.
**Evidence:** `05-read-support.tsv`; `03-funnel-counts.tsv`.

---

## D14 -- Cell-level panels use density, not counts, on a log axis

**Date:** 2026-09-09
**Chosen:** Panel 04d plots kernel density of cell-level AF per caller.
**Beat:** A histogram or frequency polygon with a log-scaled count axis.
**Why:** The two sets carry very different numbers of cell observations, so raw
counts compare size rather than shape, and empty bins are `-Inf` on a log axis,
which raised a warning on every run.
**Evidence:** `log-10 transformation introduced infinite values` from
`scale_y_log10` in the run log; gone after the change.

---

## D15 -- Three arms, not two

**Date:** 2026-09-09
**Chosen:** Compare three arms: original mgatk, the scMOCHA variant call, and
scMOCHA after its downstream AF>5% gate.
**Beat:** Two arms, with the AF>5% gate folded into "scMOCHA".
**Why:** The user pointed out that AF>5% is applied in the downstream analysis,
not during variant calling. Folding it in misrepresents what the caller emits
and makes scMOCHA look more restrictive than it is: the call returns 736
variants, the AF>5% analysis set is 216. Both numbers matter, and a reviewer
asking "what does your caller output" must get 736.
**Evidence:** `src/06.1-collect-variants-new.R` applies the AF and cell-count
rules; `newplots/thecode/scmocha-mgatk-variant-calling.py` does not.

---

## D16 -- Exclusion attributed to the first failed criterion

**Date:** 2026-09-09
**Chosen:** Assign each excluded variant to the first criterion it fails in
that arm's own order, with "fails both VMR and strand r" kept as its own
category.
**Beat:** Reporting every criterion each variant fails, which double-counts and
makes the bars unreadable; and picking one arbitrarily when several fail.
**Why:** The bars must sum to the set size for the figure to be interpretable
as "where the difference comes from". Keeping the both-fail case explicit
stops the strand-correlation bar from absorbing variants that VMR would also
have removed, which is what makes the 87% figure trustworthy.
**Evidence:** 171 scMOCHA AF>5%-only variants split 149 / 17 / 5 across
strand-only, both, and VMR-only; 185 mgatk-only variants split 117 / 68.
**Result:** 87% of what mgatk misses is lost to the strand-correlation cutoff
alone. This is the cleanest single number the stage produces.

---

## D17 -- Carriers defined by a background error model, not by a fixed floor

**Date:** 2026-09-09
**Chosen:** A cell carries the variant when its alt count is inconsistent with
the measured background error rate under a binomial test, Bonferroni-corrected
over all 9,005,290 cell-by-variant observations. Per-variant heteroplasmy is
the median AF over those cells (`af_carrier_*`).
**Beat:** The `>= 2 alt reads` floor adopted in D12, and before that the
AF >= 0.05 floor and the population `mean` column.
**Why:** D12 fixed the censoring but introduced the opposite error. A fixed
read floor does not scale with depth, so at the measured background rate of
0.026% two alt reads in a cell of depth 1,000 is expected by chance and
qualifies as a carrier. The median then describes the background rather than
the variant. This surfaced as an apparent contradiction the user spotted in
panel 04c: 29% of the scMOCHA AF>5% arm appeared to have carrier AF below 0.01,
for a set gated on having 10 cells above 0.05.
**Evidence:** `1801A>G` -- 19 cells at AF >= 0.05 with median 0.067, but 895
cells with >= 2 alt reads whose median is 0.0017. Across the AF>5% arm the
median carrier count was 26 under the gate rule and 144 under the 2-read rule.
After the change, only 6 of 216 variants sit below 0.01, and all six have a
carrier maximum between 0.09 and 0.53, which is a real pattern rather than an
artifact.
**Threshold choice:** tested at 1e-2, 1e-3, 1e-5, 1e-8 and Bonferroni
(5.55e-9). The direction is identical throughout and the separation strengthens
as the definition tightens, which is what dilution by background predicts.
Bonferroni was taken because it needs no tuning parameter.
**Background estimate:** computed over cells nowhere near carrying the variant
and **including the zero-alt cells**. Excluding them, as the first attempt did,
inflates the rate from 0.026% to 0.059%.
**Consequence:** the headline test moved from P = 0.043 to P = 0.0029, medians
0.069 rejected against 0.205 passed. Panel 04f now reports all four definitions
so the choice is auditable rather than asserted.

---

## D18 -- Report the sensitivity, not just the adopted measure

**Date:** 2026-09-09
**Chosen:** Ship `04f-measure-sensitivity` and
`04-measure-sensitivity.tsv` with every carrier definition tried, including the
two that were rejected.
**Beat:** Reporting only the adopted measure.
**Why:** Three different definitions were adopted in sequence and the P-value
moved from 0.59 to 0.043 to 0.0029. A result that moves that much with an
analyst choice has to carry the choice with it, or it is not evidence. The
direction is identical under all four, which is the claim that actually holds.
**Evidence:** `04-measure-sensitivity.tsv`.

---

## D19 -- The 5% restriction applies to the scMOCHA AF>5% arm only

**Date:** 2026-09-09
**Chosen:** Panel `04g` measures each arm on the rule that defines it: the
scMOCHA AF>5% arm on cells at AF >= 0.05, and original mgatk and the scMOCHA
variant call on all carrier cells, because neither applies an AF rule. `04c`
stays on the common measure.
**Beat:** Two alternatives. Replacing `04c` outright, and applying the 5%
restriction symmetrically to all three arms.
**Why:** **User override.** The symmetric version was implemented first and
recommended on comparability grounds. The user rejected it: the 5% rule is a
property of the AF>5% arm alone and imposing it on the other two describes
something neither caller does. That is correct as a description of what each
arm delivers, which is what the panel is for.
**Cost, carried on the figure:** the three bars no longer share a measure. The
AF>5% bar cannot fall below 0.05 while the other two can, so part of its
apparent shift toward higher heteroplasmy is definitional. Under the common
measure the arm medians are 0.058 mgatk, 0.018 call, 0.100 AF>5%; under the
symmetric 5% restriction they were 0.122, 0.098, 0.125. The subtitle of `04g`
states that it is not a like-for-like comparison and points to `04c`, and the
output table carries a `measured` column naming the measure used per arm.
**Status:** settled by the user; the caveat is documentation, not a reopening.
