# Reference: domain checklists

The questions a brainstorm must raise, by analysis type, and the parameters
that always need a justification row.

This file lists **what to ask**. It does not say what the answer is -- that is
the analysis design, and it belongs in the brainstorm's `Q` blocks. The rules
for judging an answer correct live in the `statistical-genetics` skill; the
schema and key checks live in `data-verification`.

Use the universal list every time, then the one or two blocks that match.

---

## A. Universal — any cohort-scale omics analysis

Ask all of these before the first stage is planned.

**Cohort**

- What defines the analysis cohort, and which single file will own that
  definition? Two competing definitions is the most common source of
  irreproducible counts.
- What is the analysis `n` after intersecting every required input? This is
  usually one join of tables that already exist, costs minutes, and every power
  statement and resource request depends on it. Compute it before planning
  anything else.
- Are there repeated measures per subject? If so, which model absorbs them, and
  does the chosen tool support that at all?
- Are there technical replicates, lab controls, or cross-platform duplicates?
  They are QC assets, not samples, and must be excluded explicitly.
- Which exclusions are QC and which are scientific? They need separate counts.

**Identity and joins**

- What is the unit of each input: subject, sample, draw, array, library, cell?
  A subject identifier used as a sample key silently collapses rows.
- What is the join key on each side, in its exact written form?
- Is that key actually unique in the table in hand, not just in the filtered
  table the producer saw?
- What does a row-count assertion on both sides of every merge look like?

**Coordinates and builds**

- What genome build is each coordinate-bearing input on, and where is that
  recorded? A bare position column is a defect.
- Does any step cross builds? If so, the liftover is a named step that reports
  its unmapped count.
- Chromosome naming: `chr1` or `1`? A mismatch usually yields zero results with
  no error.
- Coordinate base: BED is 0-based half-open, most everything else is 1-based.

**Confounding**

- Which covariates are mandatory, which are optional, and what is the source
  for each claim?
- Is any candidate covariate on the causal path between exposure and outcome?
- Which latent factors will be estimated, from what, and how many?
- Is any covariate compositional (summing to a constant)? One category must be
  dropped or the design is singular.
- What technical batch structure exists, and is it confounded with the exposure?

**Scale and compute**

- What is the test count, and therefore the correction?
- What is the largest single file, and does any step have to hold it in memory?
- What is the wall time of one unit, measured on a smoke run?

**Deliverables**

- What is the smallest result that would answer the question? Plan to produce
  that first.
- What would falsify the result? A negative control, a permutation, a held-out
  set.

---

## B. GWAS

- Imputation panel and its build; imputation quality threshold and the metric
  it is on.
- MAF or MAC cutoff, and which one -- MAC is the honest choice at small `n`.
- Hardy-Weinberg threshold, applied **within** ancestry group, never pooled.
- Relatedness: threshold, method, and whether relateds are removed or modelled.
- Population structure: principal components, a mixed model, or both -- and if
  both, why that is not double-counting.
- Which allele the effect is relative to, stated as a column in the output.
- Ambiguous strand variants (A/T, C/G): kept, dropped, or resolved by frequency.
- Sex chromosomes: separate pass, and which dosage coding.
- Trait transformation, and whether a case-control imbalance needs a
  saddlepoint or Firth correction.
- Genomic inflation: what lambda would be acceptable, and what the QQ plot must
  look like before the full genome runs.
- Replication or meta-analysis: which cohort, exact-lead or proxy matching, and
  the heterogeneity statistic reported.

---

## C. Molecular QTL (eQTL, meQTL, pQTL, sQTL, caQTL)

- Cis window, and the argument for that width in this modality. The 1 Mb
  convention comes from eQTL and is wide for methylation.
- Cis and trans thresholds, each backed by the actual test count rather than a
  citation.
- Permutation scheme: per-feature, beta-approximation, number of permutations.
- Which engine is primary, which is the cross-check, and what a disagreement
  between them would mean. Two independent engines is evidence; three is file
  formats.
- Feature-level covariates: latent factors from the molecular matrix, and
  whether any of them is itself genetically driven. A genotype-associated
  latent factor absorbs the signal being looked for.
- Whether externally computed factors may be fed to the engine at all, or
  whether its built-in factor mode must be used -- and what its built-in mode
  does **not** do.
- Feature normalization and transformation, and its direction of effect on
  heteroscedasticity.
- Platform or assay version effects, and whether they are a covariate, a
  stratification, or a separate pass meta-analyzed.
- Probe or feature artifacts: a polymorphism under a probe manufactures a QTL
  out of nothing and is the most common artifact in this field. Which mask
  file, which build, which allele frequency cutoff.
- Sample identity: the modality's own genotype-like probes as an independent
  check that the molecular sample and the genotype are the same person.
- Cell composition, which is compositional and usually the dominant confounder
  in a bulk tissue.
- What the engine does **not** support -- fine-mapping, colocalization,
  interaction with more than one factor, region-based enrichment -- and what
  external tool fills each gap.
- Interaction models: the term must be mean-centered, or the main effect is the
  effect at modifier = 0, not the average effect.
- Screening then refitting the same samples is selection bias. If the design
  does it, the brainstorm says so and the results say so.

---

## D. EWAS and methylation arrays

- Beta or M-value for the model, and the clamp applied before the logit.
- Probe universe: which mask union, which manifest, which build. Two manifests
  for the same platform are frequently on different builds.
- Cross-reactive probes, SNP-overlapping probes, and sex-chromosome probes:
  removed, flagged, or kept.
- Replicate probe groups on newer platforms: collapsed how, and what happens to
  a group whose members disagree about the mask.
- Normalization recipe, and whether a second recipe exists as a sensitivity
  track. If two finished matrices disagree in probe count, which is primary and
  what the intersection count is.
- Outlier handling: the rule, the number of iterations, and whether it is per
  probe or per sample.
- Cell composition: reference-based or reference-free, which reference, per
  sample and never per subject.
- Detection p-value and bead-count thresholds, and the sample-level pass rate.
- Batch: slide, array position, plate, scan date, and which of them survives
  into the model.
- Predicted sex versus recorded sex, and which one the analysis trusts.

---

## E. Bulk RNA-seq differential expression

- Quantifier and reference annotation version, both recorded.
- Gene-level or transcript-level, and the summarization method.
- Filtering rule for low-expression features, stated as counts and sample
  fraction.
- Normalization, and whether composition bias is plausible.
- Dispersion model, and whether the design has enough replicates to estimate it.
- Covariates versus latent factors, and the risk of removing the effect of
  interest.
- Shrinkage of effect sizes, and which estimator.
- Independent filtering and the multiple-testing method.
- Effect-size threshold as well as a p-value threshold, and its justification.
- Downstream enrichment: which background set, which ontology version.

---

## F. Single-cell

- Which cells are cells: droplet calling method, ambient RNA, doublets.
- QC thresholds per modality, chosen per sample or globally, and why.
- Integration method, what it assumes, and what a failure looks like.
- Clustering resolution, and what makes a resolution the right one here.
- Annotation: reference-based, marker-based, or manual -- and who adjudicates.
- Pseudobulk versus cell-level testing. Cell-level tests treat cells as
  independent replicates and inflate everything.
- Number of subjects, not number of cells, is the replication unit for any
  between-group claim.
- Compositional analysis: cell-type proportions sum to one.

---

## G. Integration, colocalization, and Mendelian randomization

- Which summary statistics, which build, which effect allele, which `n`.
- Sample overlap between the two datasets, and what it does to the method.
- Colocalization prior probabilities, and the sensitivity of the conclusion to
  them.
- Fine-mapping: which LD reference, and whether it matches the study ancestry.
  A mismatched LD panel produces confident, wrong credible sets.
- Instrument selection for MR, and whether it uses the same samples as the
  outcome -- winner's curse.
- Pleiotropy tests planned, not chosen after seeing the result.

---

## H. Parameters that always need a justification row

Any of these appearing in a plan without a row in the brainstorm's parameter
table is an unexamined default.

| Parameter class       | Examples                                                     |
| --------------------- | ------------------------------------------------------------ |
| Inclusion thresholds  | MAF, MAC, call rate, HWE, imputation quality, detection p    |
| Window and distance   | cis window, LD window, promoter definition, TSS distance     |
| Significance          | genome-wide, cis, trans, FDR method, permutation count       |
| Transformation        | log, logit, rank-inverse-normal, clamp bounds, scaling       |
| Latent structure      | number of factors, variance explained cutoff, selection rule |
| Filtering of features | expression floor, variance ranking, mask union               |
| Model form            | fixed versus random, interaction centering, link function    |
| Compute               | cores, memory per slot, wall limit, array width, chunk size  |

Each row answers four things: the value, why this value, what changes if it
moves, and where the value came from. "The default" is a valid answer only when
it names the tool and version whose default it is.

---

## I. The silent failure modes

These produce zero errors, a clean exit, and a wrong or empty result. Every
brainstorm that touches the relevant data type lists the ones in scope under
`## Traps already paid for`, with the assertion that would catch each.

| Failure                                         | Symptom                                    |
| ----------------------------------------------- | ------------------------------------------ |
| Chromosome naming mismatch between two inputs   | zero pairs tested, no error                |
| Sample identifiers that do not intersect        | silently drops to zero samples             |
| Genome build mismatch                           | wrong features tested, no error            |
| 0-based versus 1-based coordinate confusion     | every window shifted one base              |
| Transposed matrix orientation                   | tool reads samples as features             |
| Effect allele flipped between tool versions     | every sign reversed                        |
| Uncentered interaction term                     | main effect means something else           |
| A key assumed unique that is not                | a left join inflates row count             |
| Exploded lookup table summed without collapsing | counts exceed the total                    |
| Pre-allocated output file                       | file size tells you nothing about progress |
