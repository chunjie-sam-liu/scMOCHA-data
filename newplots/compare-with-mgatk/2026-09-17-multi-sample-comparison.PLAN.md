# PLAN: extend the mgatk comparison from one sample to five

Stage: `newplots/compare-with-mgatk/` (flat stage)
Campaign: `2026-09-17-multi-sample-comparison`
Written: 2026-09-17
Status: **approved 2026-09-17, all recommendations accepted; implemented**

Parent: `PLAN.md` (stage charter). This campaign does not change the criterion
audit or the rebuttal argument; it changes the evidence base from one sample to
five and reworks the stage layout to carry that.

---

## 1. Goal

`PROGRESS.md` records the standing weakness of this stage:

> **One sample only.** The editor's criticism is about the whole resource.
> Extending to more samples needs paired `*_mgatk_original.tsv.gz` outputs,
> which currently exist for this sample alone.

Four more samples with paired outputs now exist. Close that gap: run the same
criterion-by-criterion comparison on all five, and add a cross-sample layer
that states what holds across samples rather than what held in one.

---

## 2. Data inventory (verified 2026-09-17, read from the archives)

Root: `${ISILON_BASE}/compare-with-mgatk`. Five archives, each with the same
17-file flat layout and the same 14-column `variant_stats` schema.

| Archive | Chemistry | Cells | Cells mean cov > 10 | S0 candidates | S1 scMOCHA | S1 mgatk | **S2 mgatk** |
| --- | --- | ---: | ---: | ---: | ---: | ---: | ---: |
| `GSE149689_GSM4509019_3PV3.zip` | SC3Pv3 | 721 | 332 | 1,897 | 20 | 19 | **0** |
| `GSE163314_GSM4976997_3PV2.zip` | SC3Pv2 | 7,949 | 274 | 7,077 | 9 | 4 | **0** |
| `GSE163668-GSM4995445_5PR2.zip` | SC5P-R2 | 191 | 170 | 2,484 | 20 | 21 | **0** |
| `GSE181279-GSM5494116_5PPE.zip` | SC5P-PE | 7,210 | 6,724 | 25,746 | 736 | 1,247 | 230 |
| `GSE271107_GSM8369876_3PV3.zip` | SC3Pv3 | 8,645 | 1,490 | 7,540 | 2 | 6 | **0** |

S1 is `n_cells_conf_detected >= 3`. S2 mgatk is S1 plus `vmr > 0.01 AND
strand_correlation > 0.65`.

### 2.1 The finding that forces the design

**Original mgatk retains zero variants in four of the five samples.** Not
because of missing data: every mgatk S1 variant in those four samples fails the
`strand_correlation > 0.65` floor, and the values are not near it. The highest
strand correlation observed in any mgatk S1 variant across the four samples is
**0.59** (`2803A>T`, GSE163668); most sit between `-0.3` and `+0.4`, and many
are negative. Both `vmr` and `strand_correlation` are populated, so this is a
real rejection, not an `NA` artefact.

These four samples are shallow 10x libraries: median per-cell MT coverage is far
below GSE181279, and only 4% (GSE163314) to 89% (GSE163668) of cells even clear
mgatk's own coverage filter. That is exactly the regime in which a
forward/reverse read-correlation floor degenerates, because each cell
contributes too few alt reads for the correlation to be estimable.

Two consequences, both of which the implementation has to handle:

1. **The within-sample contrast that answers the editor does not exist in four
   of five samples.** Step 04's Wilcoxon splits variants passing the scMOCHA
   gate by mgatk's gate call. In four samples that split has one level, so the
   test is undefined and several panels have an empty group. The scripts
   currently assume both groups are non-empty and will error.
2. **The cross-sample layer, not the per-sample layer, is where the result
   lives.** Per-sample variant counts of 2, 9, 20, 20 cannot carry a
   distribution comparison. What they can carry is a yield statement and a
   pooled mechanism analysis.

### 2.2 Second verified defect: the mgatk cell filter is not reproducible from `depthTable`

`cell.depthTable.txt` stores mean coverage **rounded to two decimals**.
GSE163314 has four cells recorded as exactly `10.0`; two of them are present in
`cell.cell_heteroplasmic_df_mgatk_original.tsv.gz` and two are not. So
`mean_coverage > 10` recomputed from that file gives 274 kept cells while mgatk
itself kept 276, and step 01's

```r
stopifnot(setequal(cell_depth[kept_by_mgatk == TRUE, barcode], barcodes_mgatk))
```

fails on that sample. The other four samples agree exactly.

---

## 3. Design

### 3.1 Sample registry in `config.R`

One table, the single source of truth for which samples exist:

| sample_id | archive | gse | gsm | chemistry |
| --- | --- | --- | --- | --- |
| `GSE149689_GSM4509019_3PV3` | `GSE149689_GSM4509019_3PV3.zip` | GSE149689 | GSM4509019 | SC3Pv3 |
| `GSE163314_GSM4976997_3PV2` | `GSE163314_GSM4976997_3PV2.zip` | GSE163314 | GSM4976997 | SC3Pv2 |
| `GSE163668_GSM4995445_5PR2` | `GSE163668-GSM4995445_5PR2.zip` | GSE163668 | GSM4995445 | SC5P-R2 |
| `GSE181279_GSM5494116_5PPE` | `GSE181279-GSM5494116_5PPE.zip` | GSE181279 | GSM5494116 | SC5P-PE |
| `GSE271107_GSM8369876_3PV3` | `GSE271107_GSM8369876_3PV3.zip` | GSE271107 | GSM8369876 | SC3Pv3 |

`sample_id` normalises the `-` in two archive names to `_` so one token is safe
in a path, a factor level, and a file name. `archive` keeps the real filename.
Chemistry levels match the existing `color_chemistry` keys.

This retires stage `Q1` / `D5`: `SAMPLE_LABEL <- "<pending>"` is replaced by a
per-sample label built from the registry, so no figure subtitle says
`<pending>` again.

### 3.2 New step `00-extract-archives.sh`

Unzip each archive into `${ISILON_BASE}/compare-with-mgatk/samples/<sample_id>/`,
extracting **only the seven files this stage reads**:

```
cell.variant_stats.tsv.gz
cell.variant_stats_mgatk_original.tsv.gz
cell.cell_heteroplasmic_df.tsv.gz
cell.cell_heteroplasmic_df_mgatk_original.tsv.gz
cell.cell_heteroplasmic_df_raw.tsv.gz
cell.depthTable.txt
cell.coverage.txt.gz
```

`cell.rds`, `cell.signac.rds`, the four allele-count matrices and the PNGs are
skipped. That is about 1.3 GB extracted instead of 3.6 GB (measured after the
run). Idempotent: a sample whose seven files already exist is skipped, so
re-running costs nothing.

Sources `config.sh`-style environment activation per the repository contract
(`.env` then `pixi shell-hook`), with the bash-completion strip already recorded
in the environment pitfalls.

The already-extracted flat `cell.*` files at the input root are GSE181279. They
are left untouched; nothing will read them after this change.

### 3.3 Per-sample invocation

`01` to `05` gain one argument and otherwise keep their bodies:

```bash
pixi run Rscript newplots/compare-with-mgatk/01-load-harmonize.R --sample=GSE149689_GSM4509019_3PV3
```

`GetoptLong`, direct call, validated against `SAMPLES$sample_id` with the valid
list printed on error. A driver `run-all.sh` loops samples x steps and then runs
the cross-sample steps, so the whole stage is still one command.

Rejected alternative: looping inside each script. It would force every 300-line
body into a function or an extra indent level, make a single-sample re-run
impossible, and hold five samples' caches in one R session.

### 3.4 Output layout

| What | Now | After |
| --- | --- | --- |
| Cache | `<indir>/derived/` | `<indir>/derived/<sample_id>/` |
| Figures | `figures/` | `figures/<sample_id>/` and `figures/cross-sample/` |
| Tables | `tables/` | `tables/<sample_id>/` and `tables/cross-sample/` |

The existing flat `figures/` and `tables/` are GSE181279 outputs that the new
run reproduces under `figures/GSE181279_GSM5494116_5PPE/`. They are retired by
rename to `figures_stale-2026-09-17/` and `tables_stale-2026-09-17/`, never
deleted, and the new run recreates the original paths (see `Q1`).

### 3.5 Robustness to degenerate arms

Required by section 2.1, and checked on every panel that groups or tests:

- `04e` Wilcoxon and `04f` sensitivity table: run only when both gate-call
  groups have at least 3 variants. Otherwise emit the group medians with
  `p_value = NA`, log a warning, and skip the panel.
- `03b` Venn, `03d` / `03e` exclusion bars, `05c` read support: an empty set or
  an empty group produces a stated "no variants" placeholder panel rather than
  an error or a silently blank PDF.
- Every per-sample table gains a leading `sample` column so the cross-sample
  step can bind them without re-deriving provenance.

A degenerate arm is a **result**, not a failure, and the tables must record the
zero rather than omit the row.

### 3.6 Cell inclusion taken from the caller, not recomputed

`kept_by_mgatk` becomes membership in the barcode set of
`cell.cell_heteroplasmic_df_mgatk_original.tsv.gz`, which is mgatk's own answer.
The `mean_coverage > CUTOFF_CELL_MEANCOV` rule is kept as a reported criterion
and cross-checked: the symmetric difference is logged, and the assertion becomes
"differs only for cells whose rounded coverage equals the cutoff" instead of
`setequal`. See section 2.2 and `Q4`.

### 3.7 New step `07-cross-sample.R`

The load-bearing step after this change. Reads the per-sample tables and
produces:

- **07a Yield.** Variants retained by each arm in each sample, grouped bars.
  The headline: original mgatk retains 0 in four of five samples, scMOCHA
  retains 2 to 216.
- **07b Cell filter.** Fraction of cells discarded by mgatk's coverage filter
  per sample, with the detections and variants lost with them.
- **07c Mechanism, pooled.** All mgatk S1 variants from all five samples
  (about 1,297), strand correlation against per-cell alt-read support, with the
  0.65 floor drawn. This is the cross-sample demonstration that the floor tracks
  read depth rather than variant truth, and it is the analysis that replaces the
  within-sample test that four samples cannot support.
- **07d Exclusion reasons.** Which mgatk cutoff rejects each S1 variant, stacked
  per sample. Expected to be dominated by strand correlation everywhere.
- **07e Heteroplasmy of what the gate rejects**, per sample, on the error-model
  carrier measure, with per-sample n. Where both gate-call groups exist
  (GSE181279 only, on current data) the Wilcoxon is shown; elsewhere only the
  rejected group is drawn and labelled as having no comparison group.

No pooled p-value across samples. With one sample contributing both groups, a
stratified test reduces to that sample, and Fisher's method over four undefined
tests is not a combination of anything. The cross-sample claim is the yield
table plus the pooled mechanism panel; it is not a meta-analysis (`Q6`).

### 3.8 `06-summary-workbook.R` becomes cross-sample

One workbook, one sheet per topic, every sheet carrying a `sample` column, built
by binding the five per-sample tables. Adds a `00_Samples` sheet from the
registry and a `01_Yield` sheet from `07a`. Replaces five separate workbooks,
which nobody would open five of.

### 3.9 Colours

`color_sample`, keyed by `sample_id`, added to `high-res/00-colors.R` - the one
colour file this track sources. Derived from the existing `color_chemistry`
anchors so a sample's colour stays tied to its chemistry, with the two SC3Pv3
samples separated by lightness. No hex is written in stage code.

---

## 4. Work items

| # | File | Change |
| --- | --- | --- |
| 1 | `high-res/00-colors.R` | add `color_sample` |
| 2 | `config.R` + `.md` | `SAMPLES` registry, sample-aware `stage_paths()` / `stage_inputs()`, `fn_sample_label()`, cross-sample paths, drop `SAMPLE_LABEL <- "<pending>"` |
| 3 | `00-extract-archives.sh` + `.md` | new, idempotent, seven files per sample |
| 4 | `01-load-harmonize.R` + `.md` | `--sample`; cell filter from the mgatk AF matrix; per-sample cache |
| 5 | `02-cell-inclusion.R` + `.md` | `--sample`; per-sample dirs; `sample` column in the table |
| 6 | `03-variant-funnel-overlap.R` + `.md` | `--sample`; empty-set guards on 03b/03d/03e |
| 7 | `04-heteroplasmy-spectrum.R` + `.md` | `--sample`; guard the Wilcoxon and the sensitivity table |
| 8 | `05-vmr-strand.R` + `.md` | `--sample`; empty-group guards |
| 9 | `06-summary-workbook.R` + `.md` | cross-sample workbook |
| 10 | `07-cross-sample.R` + `.md` | new |
| 11 | `run-all.sh` | driver |
| 12 | `AGENTS.md`, `README.md`, `DIAGRAM.md` | run commands, layout, new pitfalls, five-sample reading guide |
| 13 | `PROGRESS.md` | campaign index row; supersede the "one sample only" open item |

---

## 5. Verification

1. `Rscript -e 'parse(...)'` on every changed `.R`; `bash -n` on every `.sh`.
2. `grep -nP '[^\x00-\x7F]'` over stage sources returns nothing.
3. Extraction: seven files present and non-empty for all five samples; re-run is
   a no-op.
4. Smoke rung: full chain on `GSE271107_GSM8369876_3PV3` (2 variants, the most
   degenerate case) before any other sample. If the guards hold there they hold
   everywhere.
5. Then `GSE181279_GSM5494116_5PPE` and confirm its numbers reproduce the values
   in `PROGRESS.md` (7,210 cells, 486 dropped, S1 736 / 1,247, S2 230, gate test
   0.069 vs 0.205, P = 0.0029). A mismatch means the refactor changed behaviour.
6. Then the remaining three, then `06` and `07`.
7. Figure count and table row counts checked against the registry: every sample
   present in every cross-sample table, zeros included.

---

## 6. Risks

- **Runtime and memory.** Step 01 streams two coverage files over 200 MB
  gzipped. Handled one sample per process, as now. Estimated a few minutes per
  large sample.
- **Reproducing GSE181279 exactly.** Item 5 above is the guard. The cell-filter
  change in 3.6 is a no-op on that sample (verified: the two sets agree).
- **Scope.** This campaign does not re-audit criteria, does not touch
  `src/06.1-collect-variants-new.R`, and adds no new statistical claim beyond
  those in section 3.7.

---

## 7. Questions

Each has a recommended answer. Approving the plan without comment approves the
recommendations.

### Q1 - Retire the existing flat `figures/` and `tables/`?

The 22 PDFs and 15 tables at `figures/` and `tables/` are GSE181279 outputs,
which the new run reproduces under a per-sample subdirectory.

**Recommended:** rename to `figures_stale-2026-09-17/` and
`tables_stale-2026-09-17/`. Nothing is deleted, the old response-letter figures
stay reachable, and the new run recreates the original paths.
**Alternative:** leave them in place, so the directory mixes one sample's flat
output with five samples' subdirectories.

### Q2 - One sample per invocation, or loop inside each script?

**Recommended:** `--sample=<id>` on steps 01 to 05, plus `run-all.sh`. Keeps
bodies unchanged, allows a single-sample re-run, and bounds memory.
**Alternative:** loop inside each script, one command, larger diff, no
single-sample re-run.

### Q3 - Is the cross-sample step in scope for this campaign?

Per-sample variant counts of 2 to 20 cannot carry the comparison alone.

**Recommended:** yes, add `07-cross-sample.R` now; it is where the five-sample
claim is actually made.
**Alternative:** per-sample outputs only, cross-sample deferred - which leaves
five sets of figures and no statement about the resource.

### Q4 - Define `kept_by_mgatk` from the mgatk AF matrix instead of recomputing it?

Section 2.2: `depthTable` is rounded, so the recomputed filter disagrees with
mgatk on GSE163314 by two cells.

**Recommended:** take the kept set from mgatk's own output and keep the
coverage rule as the reported criterion, logging the symmetric difference.
**Alternative:** relax the assertion to a tolerance and keep recomputing - still
reports a cell count that mgatk did not use.

### Q5 - How should a sample with zero variants in an arm be presented?

**Recommended:** the panel is drawn with a stated "no variants in this arm"
annotation, and the table carries the explicit `0` row. A zero here is the
finding.
**Alternative:** skip the figure entirely, which makes a missing file
indistinguishable from a failed run.

### Q6 - Any pooled statistical test across samples?

**Recommended:** none. Report per-sample results plus the pooled mechanism panel
`07c`. Four samples have no comparison group, so a stratified test collapses to
GSE181279 and a combined p-value would be one p-value wearing a meta-analysis
label.
**Alternative:** report a stratified van Elteren test anyway, with the caveat.

### Q7 - Create `DATA.md` at the repository root for these five archives?

The repository has no `DATA.md` today and no `.github/instructions/` bindings,
so this would establish the convention rather than follow it.

**Recommended:** yes, a minimal `DATA.md` covering the five archives, their
origin, date, the env variable that exposes the root, and the stage that
consumes them - five rows, since the provenance of the four new archives is
currently recorded nowhere.
**Alternative:** record the inventory inside this campaign's files only, and
leave the repository-level convention for a later decision.

### Q8 - Where did the four new archives come from, and when?

Needed for `Q7` and for the figure captions. The files landed on 2026-09-11 but
the producer is not recorded anywhere in the repository.

**Recommended answer needed from the user:** who produced them, with what
pipeline version, and whether the scMOCHA and mgatk arms were run on identical
allele counts as they were for GSE181279 (stage `D1` depends on that).
