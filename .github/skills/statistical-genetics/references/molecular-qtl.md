# Reference: molecular QTL specifics

The full form of `SKILL.md` section 8.1. These rules apply when the outcome is a
molecular phenotype measured on an array or by sequencing — methylation,
expression, protein. Read this before fitting any meQTL, eQTL, pQTL, sQTL, or
caQTL model, and when reviewing one.

Each is a failure mode the generic calibration checks do not catch. None of them
inflates the genomic inflation factor, so every one of them produces a clean QQ
plot and a wrong answer.

---

## A. Sample identity before the first model

A molecular dataset and a genotype dataset are linked by a sample identifier
that was typed by a person. Before any QTL fit, verify identity with data: the
array's own genotyping probes (the SNP-probe betas an EPIC pipeline exports) or
a sequencing-derived genotype call against the WGS/array genotypes of the
supposedly same individual.

Report concordance per sample, name the threshold, and remove or reassign every
discordant sample. A swap does not inflate lambda; it silently deletes true cis
signal and creates trans artifacts. Sex-vs-genotype is a subset of this check,
not a substitute.

---

## B. Probe-SNP artifacts

A variant inside the probe body, at the extension base, or at the CpG itself
changes hybridization or the measured site, so a cis association there is
technical, not regulatory.

Every cis result table carries a flag column derived from the operative probe
mask plus a computed distance from the variant to the probe interval, and the
calibration report states what fraction of leads is flagged. Do not drop flagged
pairs silently; flag them, and report with and without.

---

## C. Cell-fraction covariates are compositional

Estimated cell fractions sum to one, so including all K with an intercept is
exactly collinear and the fit either fails or silently drops one. The plan
states which coding is used: K-1 fractions with the reference cell named, or a
CLR/ILR transform.

Check the fraction estimates against the exposure and against genotype like any
other latent factor (`SKILL.md` section 4).

---

## D. Technical batch is a named covariate, not only a latent factor

Slide, plate, processing date, and array version are known and should be listed
explicitly with their confounding against the exposure reported.

When two array versions or platforms are pooled, the samples measured on both
are the only direct measurement of the cross-platform shift; a pooled model must
show that this shift is absorbed before the full scan.

---

## E. Sex chromosomes are a separate decision

X and Y probes are not pooled with autosomes under one model and one threshold.
The plan says one of: excluded; sex-stratified; or sex-adjusted with X dosage
coded explicitly. Y is male-only by construction.

A filter that deliberately keeps sex-chromosome probes has not made this
decision yet; it still needs a `Q`.

---

## F. Outcome transform and outlier rule are per-feature

Clamping, logit / M-value transform, and outlier removal run per CpG or per
gene, and the per-feature N after outlier removal goes into the result table
(`SKILL.md` section 11), not a single cohort N.
