# PLAN: scMOCHA-updated mgatk vs original mgatk variant calling

Stage: `newplots/compare-with-mgatk/` (flat stage)
Written: 2026-09-09
Status: **awaiting approval — do not implement**

---

## 1. Goal

Answer the Cell Metabolism editorial concern recorded in `EDITOR.md`:

> "...concerns regarding the reliance on mgatk, which is biased toward
> detecting higher-heteroplasmy mutations that are more likely to be subject
> to functional selection, and missing many mutations that would be tolerated
> at lower heteroplasmy."

The rebuttal is that the manuscript does **not** use original mgatk. It uses a
modified caller whose criteria differ at three points, two of which directly
remove the high-heteroplasmy bias the editor names. This stage produces the
figures and tables that quantify that difference on real data.

The stage must be honest in both directions: one of the three criterion changes
makes scMOCHA **more** stringent, and the analysis reports that too.

---

## 2. Criterion audit (already verified from source)

Source files audited: `newplots/thecode/original-mgatk-variant-calling.py`,
`newplots/thecode/scmocha-mgatk-variant-calling.py`,
`src/06.1-collect-variants-new.R` (the production classification script), and
`src/06.4-somatic.md`.

| ID | Stage of the caller | Original mgatk | scMOCHA-updated | Direction |
| --- | --- | --- | --- | --- |
| C1 | Cell inclusion | drop cells with mean MT coverage <= 10 | no cell filter; all barcodes retained | scMOCHA **more permissive** |
| C2 | Confident detection in one cell | `fwd_alt >= 2 AND rev_alt >= 2` | `fwd_alt >= 2 AND rev_alt >= 2 AND (fwd_alt + rev_alt) >= 10` | scMOCHA **more stringent** |
| C3 | Variant retention | `n_cells_conf_detected >= 3` | `n_cells_conf_detected >= 3` (identical) | same |
| C4 | Reliability gate | `vmr > 0.01 AND strand_correlation > 0.65` | **no VMR filter, no strand-correlation filter** | different basis |
| C5 | Cell-by-variant AF matrix | only variants passing C3 are written | `*_raw.tsv.gz` with all candidate variants is also written | scMOCHA **more permissive** |
| C6 | Reliability gate, scMOCHA side | none beyond C4 | position blacklist (mis-alignment + RNA editing), then >= 10 cells at cell AF >= 0.05 with cell depth >= 10 | different basis |
| C7 | Downstream classification | none (mgatk stops at the variant list) | homoplasmic / heteroplasmic / somatic by cross-cell-type AF at 0.05 / 0.95, cluster reads >= 10 | different question |

**C4 is the centre of the rebuttal.** VMR (variance-to-mean ratio) combined
with a high strand-correlation floor is exactly the combination that selects
clonally expanded, higher-heteroplasmy variants: a low-heteroplasmy variant
carries few alt reads per cell, so its forward/reverse read correlation is
noise-dominated and it fails a 0.65 floor whether or not it is real, and its
across-cell variance is small so it fails `vmr > 0.01`.

scMOCHA does not apply that filter at all. It replaces it with C6, a
reliability gate built on read depth, cell count, and a position blacklist.
None of those conditions on how variable a variant is across cells, so none of
them preferentially removes low-heteroplasmy variants.

C2 is the counter-direction and is reported as such: it trades away
few-read singleton calls for per-cell AF that is actually measurable. The
analysis must show that what C2 removes is low-**read-support** calls, not
low-**heteroplasmy** calls. If the data show otherwise, that is the finding and
it goes in `DECISION.md`.

### Not part of the variant-calling path

`preprocessing/high_quality_variant/02-update-high-quality-variant.py` does
apply `vmr > 0.01 & strand_correlation > 0.3`, but it is **not** used for
variant calling. Its output `high_quality_variant.tsv` is consumed only by
`03-collect-hqv.R` and `04-hqv-flextable.R`, which build a summary flextable
under `analysis/zzz/hqv`. `src/06.1-collect-variants-new.R` never reads it, and
neither `vmr` nor `strand_correlation` appears anywhere in that script. This
stage must not present the 0.3 cutoff as a scMOCHA calling criterion.

---

## 3. Data inventory (verified 2026-09-09)

Root: `/research/rgs01/home/clusterHome/cliu68/project/scmocha/compare-with-mgatk`

| File | Content | Verified shape |
| --- | --- | --- |
| `cell.variant_stats.tsv.gz` | scMOCHA per-variant stats | 25,746 variants x 14 cols |
| `cell.variant_stats_mgatk_original.tsv.gz` | original mgatk per-variant stats | 25,580 variants x 14 cols |
| `cell.cell_heteroplasmic_df.tsv.gz` | scMOCHA cell x variant AF, post-C3 | 7,210 cells x 736 variants |
| `cell.cell_heteroplasmic_df_mgatk_original.tsv.gz` | original cell x variant AF, post-C3 | 6,724 cells x 1,247 variants |
| `cell.cell_heteroplasmic_df_raw.tsv.gz` | scMOCHA cell x variant AF, all candidates | 7,210 cells x 25,746 variants |
| `cell.depthTable.txt` | per-cell mean MT coverage | 7,210 cells |
| `cell.{A,T,C,G}.txt.gz` | per-cell strand-resolved allele counts | `pos,barcode,fwd,rev` |
| `cell.coverage.txt.gz` | per-cell per-position total coverage | `pos,barcode,cov` |
| `MT_refAllele.txt` | rCRS reference allele per position | 16,569 rows |
| `cell.rds`, `cell.signac.rds` | mgatk SummarizedExperiment / Signac object | not needed for this stage |

Columns in both `variant_stats` files (identical schema):
`position, nucleotide, variant, vmr, mean, variance, n_cells_conf_detected,
n_cells_over_5, n_cells_over_10, n_cells_over_20, n_cells_over_95,
max_heteroplasmy, strand_correlation, mean_coverage`

Applying scMOCHA's C6 gate needs per-cell depth at each variant position, which
the AF matrices do not carry. It comes from `cell.coverage.txt.gz`
(`pos,barcode,cov`; about 1.2e8 rows uncompressed). Step 01 streams that file
through `awk`, keeping only the positions of candidate variants in the union of
the two rule sets (order 1e3 positions, so order 1e7 rows), before it ever
reaches R.

Numbers already measured that anchor the figures:

- 7,210 cells total; **486 (6.74%)** have mean MT coverage <= 10 and are
  discarded by C1 in original mgatk.
- Per-cell mean MT coverage: min 0.25, q05 8.62, median 38.92, q95 222.82.
- Post-C3 variant counts: scMOCHA 736, original mgatk 1,247.

The 736 < 1,247 gap is entirely attributable to C2 and is the number the
analysis has to explain rather than hide.

---

## 4. Steps

Flat stage. Every script ships a sibling `.md` (Purpose / Inputs / Outputs /
Run / Verify).

### `config.R` + `config.R.md`

Stage constants and paths only. Holds: input root, cache root, figure root,
table root, the criterion constants (`CUTOFF_STRAND_MGATK = 0.65`,
`CUTOFF_VMR_MGATK = 0.01`, `CUTOFF_NCELLS_CONF = 3`, `CUTOFF_CELL_MEANCOV = 10`,
`CUTOFF_ALT_READS = 10`, `CUTOFF_HETEROPLASMIC = 0.05`, `CUTOFF_MIN_READS = 10`,
`CUTOFF_NOTRELIABLE = 10`, `POS_RNA_EDITING`, `POS_MISSALIGNMENT_ERROR`, AF bin
breaks), and the shared `fn_*` helpers. The blacklists and the four scMOCHA
cutoffs are copied verbatim from `src/06.1-collect-variants-new.R` so the two
never drift silently. No colours, no hex.

The two rule sets compared throughout:

| | Rule set A, original mgatk | Rule set B, scMOCHA |
| --- | --- | --- |
| cells | mean MT cov > 10 | all cells |
| per-cell confident call | fwd >= 2 & rev >= 2 | fwd >= 2 & rev >= 2 & fwd+rev >= 10 |
| variant retention | n_cells_conf >= 3 | n_cells_conf >= 3 |
| reliability gate | vmr > 0.01 & strand_corr > 0.65 | position blacklist, then >= 10 cells with cell AF >= 0.05 & cell depth >= 10 |

### `01-load-harmonize.R`

Read both `variant_stats` tables and `depthTable`, full-outer join on
`variant`, tag each variant with which caller proposed it and which stage it
survives under each rule set. Extract per-cell depth at candidate variant
positions from `cell.coverage.txt.gz` (streamed, position-filtered) so rule set
B's C6 gate can be evaluated. Cache to the stage cache dir as qs2.

Output: `<cache>/01-variant-joined.qs2`, `<cache>/01-cell-depth.qs2`,
`<cache>/01-cell-pos-coverage.qs2`.

### `02-cell-inclusion.R` (criterion C1)

- (a) Histogram of per-cell mean MT coverage on log10 x, cutoff line at 10,
  cells split kept / dropped by original mgatk.
- (b) Bar: cells retained, scMOCHA 7,210 vs original mgatk 6,724.
- (c) Minimum detectable AF as a function of per-cell coverage under each
  criterion (mgatk floor `4 / depth`, scMOCHA floor `10 / depth`), overlaid on
  the observed coverage distribution. This panel is where the C2 cost is stated
  explicitly.
- (d) For the 486 discarded cells: how many variant detections at cell
  AF >= 0.05 with cell depth >= 10 exist in them, i.e. what the cell filter
  throws away.

Figures `02a`..`02d`, table `02-cell-inclusion.tsv`.

### `03-variant-funnel-overlap.R` (criteria C3, C4, C6)

Variant counts at each stage under each rule set:

| Stage | Rule set A, mgatk | Rule set B, scMOCHA |
| --- | --- | --- |
| S0 candidates | 25,580 | 25,746 |
| S1 `n_cells_conf_detected >= 3` | 1,247 | 736 |
| S2 reliability gate | `+ vmr > 0.01 & strand_corr > 0.65` | `+ blacklist, >= 10 cells at AF >= 0.05 & depth >= 10` |

- (a) Funnel bar chart of the two rule sets side by side.
- (b) Venn of final S2 sets (both / mgatk-only / scMOCHA-only).
- (c) Cross-applied funnel: rule set A's reliability gate applied to scMOCHA's
  S1 set and vice versa, to separate how much of the final difference comes
  from C2 and how much from C4 versus C6.

Figures `03a`..`03c`, table `03-funnel-counts.tsv`.

### `04-heteroplasmy-spectrum.R` (the editor's question)

The load-bearing figure. For the three variant classes (both / mgatk-only /
scMOCHA-only) and for the S1 vs S2 transition:

- (a) ECDF of population-level AF (`mean`) on log10 x, one curve per class.
- (b) Density of `max_heteroplasmy`.
- (c) Stacked bar: variants per AF bin (`<0.01`, `0.01-0.05`, `0.05-0.10`,
  `0.10-0.20`, `0.20-0.50`, `>=0.50`) per rule set.
- (d) Cell-level view from `cell_heteroplasmic_df_raw`: distribution of nonzero
  per-cell AF restricted to each rule set's variant list.
- (e) The direct test of the editor's claim: AF distribution of the variants
  that mgatk's `vmr > 0.01 & strand_corr > 0.65` gate **rejects** but scMOCHA's
  C6 gate **keeps**, versus the AF distribution of variants both keep. Wilcoxon
  test, with the shift and an effect size reported.

Figures `04a`..`04e`, table `04-af-bins.tsv`, `04-tests.tsv`.

### `05-vmr-strand.R` (criterion C4 detail)

VMR and strand correlation are **computed** by both callers but **acted on**
only by mgatk. Plotting scMOCHA's variants in the same VMR / strand-correlation
plane shows exactly which of them mgatk's gate would have discarded, and what
their heteroplasmy is.

- (a) VMR vs strand-correlation scatter, one panel per rule set, mgatk's
  `strand_corr = 0.65` and `vmr = 0.01` lines drawn on both, points coloured by
  AF bin. scMOCHA's panel carries the same lines annotated as "not applied".
- (b) AF-bin composition of the region mgatk's gate rejects, tested against the
  region it accepts. This is the quantitative form of "biased toward higher
  heteroplasmy".
- (c) Read-support view: alt-read count per confidently detected cell for
  mgatk-only variants versus shared variants, testing whether C2 removes
  low-read-support calls rather than low-heteroplasmy calls.

Figures `05a`..`05c`, table `05-vmr-strand-rejected.tsv`.

### `06-summary-workbook.R`

One `.xlsx` via `openxlsx2` assembling: criteria table (section 2 of this
plan, machine-readable), cell inclusion, funnel counts, AF bins, statistical
tests, recovered-variant table. Plus a rendered criteria table figure for
direct use in the response letter.

Output: `tables/06-compare-with-mgatk.xlsx`, figure `06-criteria-table`.

---

## 5. Conventions this stage follows

- R only, run through pixi: `pixi run Rscript newplots/compare-with-mgatk/NN-*.R`
  from the repository root. Conda is not installed on this host; `pixi` is the
  only environment.
- `library(jutils)` first, `dotenv()` for `.env`, `fs::path()` for paths,
  `saveplot()` for figures, `GetoptLong` for CLI args (verified available).
- Script headers use the `@AUTHOR / @CONTACT / @DATE / @DESCRIPTION / @VERSION`
  block. No AI attribution.
- The intermediate data root is read as input only; nothing is written back
  into it except under a new `derived/` subdirectory (see Q4).
- Scratch goes to `~/tmp/cmp-mgatk/`, never `/tmp`.

---

## 6. Questions, resolved 2026-09-09

All six were approved as recommended; recorded as D5-D10 in `DECISION.md`.
Q3 was then superseded by the correction in section 2.

| # | Question | Approved answer |
| --- | --- | --- |
| Q1 | Sample identity for figure labels | `SAMPLE_LABEL` constant in `config.R`, set to `"<pending>"` until supplied; one edit relabels every figure |
| Q2 | One sample or more | this sample now; extension to more samples becomes a dated follow-up campaign |
| Q3 | scMOCHA strand-correlation cutoff | **superseded** — scMOCHA applies no VMR and no strand-correlation filter; the 0.3 in `02-update-high-quality-variant.py` is a side branch, not the calling path |
| Q4 | Intermediate cache location | `<intermediate root>/derived/`, bound to `CACHEDIR` in `config.R` |
| Q5 | Colour palette location | append `color_caller` and `color_af_bin` to `high-res/00-colors.R`, additive only |
| Q6 | Low-heteroplasmy truth set | none exists; criterion-to-criterion only, no sensitivity claim, limitation stated in every legend |

---

## 7. What this stage will not do

- Not re-run either caller. Both output sets already exist.
- Not touch `preprocessing/`, `analysis/`, or `src/`.
- Not modify anything under the intermediate data root except a new `derived/`
  subdirectory.
- Not claim a sensitivity improvement that the data cannot support (see Q6).
