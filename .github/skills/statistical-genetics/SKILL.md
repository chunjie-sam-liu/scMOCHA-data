---
name: statistical-genetics
description: "Apply statistical-genetics conventions when writing or reviewing association analysis code. Use when fitting or reviewing a QTL, GWAS, EWAS, meQTL, eQTL, or interaction model; choosing a threshold or multiple-testing correction; selecting covariates; handling dosages, effect alleles, or strand; deciding a cis window; running ancestry-stratified or meta-analysis passes; fine-mapping or colocalizing; or judging whether a result is real. Covers the genome-build rule, the effect-allele contract, variant and sample QC, collinearity and latent-factor screening, selection-bias traps (screen-then-refit, winner's curse, double dipping), calibration (inflation, QQ, negative controls, permutation), molecular-QTL specifics (sample identity from SNP probes, probe-SNP artifacts, compositional cell fractions, technical batch, sex chromosomes), heterogeneity across ancestries, and the validity review that code review never catches."
---

# Statistical genetics

Conventions for association analysis. Every other skill in this repository
guards against code that fails loudly. This one guards against code that runs
cleanly and produces a wrong conclusion, which no syntax check, output
verification, or code review will catch.

This skill is portable and names no project value. The cis window, the
significance thresholds, the covariate set, the cohort definition, and the
software choice come from the repository's `.github/instructions/` bindings and
its project charter. Those take precedence where they conflict.

Related skills: `data-verification` (schema and key checks before any join),
`analysis-pipeline` (the plan/decision files that record which of the choices
below were made), `r-figure` (QQ and Manhattan plots).

---

## 1. Non-negotiables

Each of these has produced a published erratum somewhere. None is caught by a
passing test.

1. **Never join two datasets on genomic coordinates without asserting the build
   on both sides.** Coordinates are meaningless without a build.
2. **Never assume which allele the effect is relative to.** State the effect
   allele explicitly, in a column, in every result table.
3. **Never select on a statistic and then estimate that same statistic in the
   same samples.** Selection inflates the estimate.
4. **Never report a p-value without saying how many tests it was drawn from**
   and what correction was applied.
5. **Never adjust for a covariate on the causal path between exposure and
   outcome.** It removes the effect you are measuring.
6. **Never run the full genome before the calibration checks pass on one
   chromosome.** An uncalibrated model wastes the whole run.

---

## 2. Genome build

The single cheapest source of silently wrong results.

- Carry the build in the column name or the object name, not in a comment:
  `pos_hg38`, `start_hg19`. A bare `pos` column is a defect.
- Record the build of every coordinate-bearing input in the provenance file at
  the point the input is registered, not later.
- A cross-build join must go through an explicit liftover step that reports how
  many records failed to map and why. Silent inner-join shrinkage is the
  failure mode: the join succeeds, the row count drops, and nobody looks.
- After a liftover, assert three things: the mapped fraction, that no record
  changed chromosome unexpectedly, and that strand was handled.
- Annotation packages carry their own build. Two packages for the same platform
  are frequently on different builds. Read the package's build, do not infer it
  from the platform name.

## 3. Genotype handling

Before any association model touches a genotype file:

| Property             | What to check                                                        | Why                                                                         |
| -------------------- | -------------------------------------------------------------------- | --------------------------------------------------------------------------- |
| Effect allele        | Which column the effect is relative to, on both sides of every merge | `REF/ALT` and `A1/A2` are different conventions; flipping is silent         |
| Strand               | Whether ambiguous variants (A/T, C/G) are present                    | They cannot be resolved by allele match alone                               |
| Dosage vs hard call  | Which one the model consumes                                         | Hard calls discard imputation uncertainty and bias effect sizes toward null |
| Multiallelic sites   | Whether they are split and normalized                                | An unsplit site silently mis-codes                                          |
| Indel representation | Left-alignment                                                       | The same indel written two ways fails to merge                              |
| Missingness          | Per-variant and per-sample rate                                      | Differential missingness creates false association                          |

Standard variant filters, all of which must be stated as numbers in the plan
and carried into the results provenance: minor allele frequency or minor allele
count, call rate, Hardy-Weinberg equilibrium p-value, and imputation quality.
Apply HWE filters within an ancestry group, never pooled across ancestries -
population structure violates HWE for real variants.

Sample-level: relatedness, sex-check against genotype, heterozygosity outliers,
and ancestry outliers. Decide whether relatedness is removed or modeled
(a genetic relationship matrix in a mixed model); do not do both, and do not do
neither.

## 4. Covariates

- **Justify each covariate as a confounder, not as "it improved fit".** A
  variable on the causal path is a mediator; adjusting for it attenuates the
  effect. A variable caused by both exposure and outcome is a collider;
  adjusting for it creates association where none exists.
- **Check collinearity before fitting, not after a convergence failure.**
  Report the actual statistic used and its threshold.
- **Latent factors absorb whatever they correlate with.** Surrogate variables,
  expression/methylation PCs, and similar unsupervised factors will happily
  absorb the exposure of interest. Screen every retained factor against the
  exposure and against genotype, and drop or flag the ones that are confounded
  with what is being tested. This is mandatory for an interaction model, where
  a treatment-correlated factor removes the interaction variance.
- **Interaction models: center the modifier.** Without centering, the main
  effect coefficient is the effect at modifier zero, which is usually not a
  real value. State in the plan what the coefficient means after centering.
- **An interaction needs roughly four times the sample size of a main effect**
  for the same power. Say so before proposing a genome-wide interaction scan.
- Covariates that are in the model must stay in the model when a cross-product
  term with them is estimated; residualizing them out first makes the
  cross-product inestimable.

## 5. Outcome transformation

- Bounded proportions are transformed before linear modelling; clamp before
  transforming so the boundary does not produce infinities, and state the clamp.
- Outlier handling is a decision, not a default: state the rule, the number of
  iterations, and how many observations it removed.
- When a rank-based transform is used as a sensitivity analysis, it is a
  separate run reported alongside the primary, not a replacement for it.
- Any transformation applied to the primary outcome must be identical across
  every model that will be compared.

## 6. Multiple testing

State three numbers, always together: **how many tests, which correction,
which threshold.**

- A threshold inherited from another study is only valid under that study's
  test count and correlation structure. Say where it came from.
- Correcting across a _selected_ subset is different from correcting genome
  wide. Both are legitimate; conflating them is not.
- Hierarchical designs (screen genome wide, then test a subset) need the
  correction stated at each level.
- Permutation-based thresholds are the defensible option when the tests are
  strongly correlated and the effective number of tests is unknown.

## 7. Selection bias

The traps, in order of how often they appear:

- **Screen-then-refit on the same samples.** Selecting lead pairs with a fast
  model and refitting them with a full model in the same data means the refit
  estimates are conditional on selection. The p-values are anticonservative and
  the effect sizes are inflated. This is sometimes the only affordable design -
  when it is, say so explicitly in the plan, and treat the refit as
  hypothesis-generating unless an independent replication or a
  selection-corrected estimator is applied.
- **Winner's curse.** The effect size of the top hit is biased upward. Never
  quote the discovery effect size as the effect size.
- **Double dipping.** Deriving a feature from the outcome and then testing it
  against the outcome.
- **Optional stopping on the analysis.** Changing the covariate set, threshold,
  or transformation after seeing the result is a decision that must land in the
  decision log with its date, not a quiet edit.

## 8. Calibration checks

Run these on one chromosome or one trait before the full pass. They are the
cheapest rung that can fail.

| Check                    | Passes when                    | Fails means                                                                       |
| ------------------------ | ------------------------------ | --------------------------------------------------------------------------------- |
| Genomic inflation factor | Near 1                         | Unmodeled structure or relatedness, or a real polygenic signal - distinguish them |
| QQ plot                  | On the diagonal until the tail | Systematic inflation or deflation                                                 |
| Negative control set     | Shows no signal                | A model or annotation bug, not biology                                            |
| Permuted outcome         | Produces a null distribution   | The pipeline itself creates association                                           |
| Effect-size distribution | Centered at zero               | Allele coding is flipped somewhere                                                |
| Known positive control   | Recovers a published signal    | The pipeline is not measuring what you think                                      |

A negative control is the highest-value check per line of code and is almost
always skipped. Choose a set where the effect being tested cannot exist by
construction, and assert it is null.

**Genomic inflation is not the pass criterion for a cis-QTL scan.** A cis
molecular QTL scan (meQTL, eQTL, pQTL) has a high fraction of true positives,
so lambda well above 1 is expected and a lambda near 1 can mean the scan is
underpowered. For a cis scan the calibration checks that decide are the
negative-control set, the permuted outcome, and the effect-size distribution;
lambda and the QQ plot are decisive for the trans scan and for any
interaction scan, where true signal is sparse.

## 8.1 Molecular QTL specifics

Rules that apply when the outcome is a molecular phenotype measured on an
array or by sequencing (methylation, expression, protein). Each is a failure
mode the generic checks above do not catch.

- **Sample identity before the first model.** A molecular dataset and a
  genotype dataset are linked by a sample identifier that was typed by a
  person. Before any QTL fit, verify identity with data: the array's own
  genotyping probes (the SNP-probe betas an EPIC pipeline exports) or a
  sequencing-derived genotype call against the WGS/array genotypes of the
  supposedly same individual. Report concordance per sample, name the
  threshold, and remove or reassign every discordant sample. A swap does not
  inflate lambda; it silently deletes true cis signal and creates trans
  artifacts. Sex-vs-genotype is a subset of this check, not a substitute.
- **Probe-SNP artifacts.** A variant inside the probe body, at the extension
  base, or at the CpG itself changes hybridization or the measured site, so a
  cis association there is technical, not regulatory. Every cis result table
  carries a flag column derived from the operative probe mask plus a computed
  distance from the variant to the probe interval, and the calibration report
  states what fraction of leads is flagged. Do not drop flagged pairs
  silently; flag them, and report with and without.
- **Cell-fraction covariates are compositional.** Estimated cell fractions
  sum to one, so including all K with an intercept is exactly collinear and
  the fit either fails or silently drops one. The plan states which coding is
  used: K-1 fractions with the reference cell named, or a CLR/ILR transform.
  Check the fraction estimates against the exposure and against genotype like
  any other latent factor (section 4).
- **Technical batch is a named covariate, not only a latent factor.** Slide,
  plate, processing date, and array version are known and should be listed
  explicitly with their confounding against the exposure reported. When two
  array versions or platforms are pooled, the samples measured on both are the
  only direct measurement of the cross-platform shift; a pooled model must
  show that this shift is absorbed before the full scan.
- **Sex chromosomes are a separate decision.** X and Y probes are not pooled
  with autosomes under one model and one threshold. The plan says one of:
  excluded; sex-stratified; or sex-adjusted with X dosage coded explicitly.
  Y is male-only by construction. A filter that deliberately keeps sex-chromosome
  probes has not made this decision yet; it still needs a `Q`.
- **Outcome transform and outlier rule are per-feature.** Clamping, logit /
  M-value transform, and outlier removal run per CpG or per gene, and the
  per-feature N after outlier removal goes into the result table (section
  11), not a single cohort N.

## 9. Ancestry and heterogeneity

- Pooled analysis with principal components is weighted toward the largest
  group. Reporting it alone hides everything about the smaller groups.
- A stratified pass plus a heterogeneity test uses all samples and is better
  powered than significance-testing a small stratum on its own.
- Report the heterogeneity statistic, not just the per-stratum p-values.
- A small stratum is hypothesis-generating. Say that in the results text, not
  only in the plan.
- Linkage disequilibrium references must match the ancestry of the samples.
  A mismatched reference is a common source of spurious fine-mapping and
  colocalization results.

## 10. Fine-mapping and colocalization

- Linkage disequilibrium must come from the analysis samples themselves
  wherever possible. Reference-panel LD mismatched to the sample is the main
  cause of implausible credible sets.
- Sample overlap between the two traits breaks the standard colocalization
  assumption. Either use a method that models the overlap, or require a shared
  high-posterior variant from independent fine-mapping of both sides.
- Report credible set size alongside the posterior. A posterior of 0.9 spread
  over 300 variants is a different claim than one over 3.

## 11. Reporting

Every association result table carries, at minimum: the effect allele, the
other allele, the effect size and its standard error, the sample size actually
used in that test, the allele frequency in the samples used, the test count and
correction, and the software and version. A result table without the effect
allele is not interpretable and not shareable.

## 12. Statistical-validity review

Use this when reviewing an analysis, in addition to the code review. Each
question has a right answer that can be pointed at in the plan or the code; "it
looks fine" is not one.

1. What is the model, written out, including every term?
2. What is the unit of observation, and are the observations independent? If
   not, how is the dependence modeled?
3. Which allele is the effect relative to, and where is that recorded?
4. What build is every coordinate on?
5. How many tests, what correction, what threshold, and where did the threshold
   come from?
6. Was anything selected on a statistic and then re-estimated in the same data?
7. Which covariates are confounders, which are mediators, and is any of them a
   collider?
8. Do the latent factors correlate with the exposure being tested?
9. What is the genomic inflation factor, and what does the negative control
   show?
10. What is the sample size of the smallest cell actually contributing to the
    reported estimate?
11. If this is an interaction, is the modifier centered, and what does the main
    effect coefficient mean?
12. What would make this result wrong, and was that checked?
13. For a molecular QTL: was sample identity verified with data before the
    fit, and how many samples failed? Are probe-SNP artifacts flagged in the
    result table? How are cell fractions coded, and how are sex chromosomes
    handled?

Findings from this review belong in the decision log as `D` entries, because
each one settles a fork that a later session would otherwise reopen.
