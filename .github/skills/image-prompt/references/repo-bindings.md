# Reference: bindings for the AA GWAS repository

Project-specific values the general skill deliberately leaves open. Everything
here is read off files in this repository; re-verify before reuse, because the
color files do change.

---

## A. Color file resolution order

Resolve every visual role to a hex from a file that the real plots already use,
in this order. Stop at the first hit and record the source in the brief's
palette table.

1. The stage's own config, e.g.
   [src_newwgs/04.03-gwas-replication/00-config.R](../../../../src_newwgs/04.03-gwas-replication/00-config.R).
   A stage that defines its own palette wins, because its PDFs already use it.
2. [src_newwgs/color.R](../../../../src_newwgs/color.R) for new-WGS track work.
3. [src/color.R](../../../../src/color.R) for original-track work.

Steps 2 and 3 are byte-identical today, so the choice only changes the source
path recorded in the brief — and matters the day they diverge. Cite the file
belonging to the track the figure's results came from.

Never invent a hex. If a role has no color yet, add it to the appropriate
`color.R` first (see the `palette` skill), then bind the figure to it.

## B. Neutral roles

The same in every figure in this project, so briefs stay comparable.
Each neutral palette row may cite this `references/repo-bindings.md` file as its
source in the brief.

| Role | Hex |
| --- | --- |
| Text | `#2B3238` |
| Rules, card outlines, hairlines | `#C9D1D6` |
| Panel tint fill | `#F4F6F7` |
| Background | `#FFFFFF` |
| Missing / not applicable | `#CCCCCC` |

## C. NPG palette in use

Both `color.R` files draw from the same NPG family, so figures mix cleanly.

| Hex | Named in `color.R` as | Typical figure role |
| --- | --- | --- |
| `#E64B35` | `color_highlight` | Accent, alarm, the thing that failed |
| `#4DBBD5` | `color_pheno_hist` | Echo / continuous outcomes |
| `#00A087` | `color_age_hist` | Second entity, adequate / passing |
| `#3C5488` | `color_scatter_point` | First entity, discovery |
| `#F39B7F` | `color_pheno_int` | Third entity, intermediate band |
| `#8491B4` | `color_sex_male` | Deliberately neutral bar |
| `#B09C85` | `color_age_dx_hist` | Low-emphasis category |
| `#374E55` | `color_group$nhw` | NHW cohort |

## D. Cohort colors: two conventions coexist

This is the single most common source of a wrong figure. Pick the one the
underlying results used and say which in the brief.

**Global, `src_newwgs/color.R` (`color_newgroup`)** — used by metadata,
phenotype, and cohort-construction figures.

| Cohort | Hex on the figure | As written in `color.R` |
| --- | --- | --- |
| `AA-Primary` | `#E64B35` | `color_group$aa`, the literal string `"red"` |
| `AA-Independent` | `#F39B7F` | `#F39B7F` |
| `NHW` | `#374E55` | `color_group$nhw`, `#374E55FF` |

**Stage-local, `src_newwgs/04.03-gwas-replication/00-config.R`
(`COHORT_COLORS`)** — used by every replication figure and inherited by stage
`04.03.02`.

| Cohort | Hex |
| --- | --- |
| `AA-Primary` | `#3C5488` |
| `AA-Independent` | `#00A087` |
| `NHW` | `#B09C85` |

Two documented substitutions to apply without asking. Cite this
`references/repo-bindings.md` file as the source for the substituted binding:

- `color_group$aa` is the literal R string `"red"`. Use `#E64B35` in figures; it
  is the NPG red already used as `color_highlight` and reproduces far better in
  print than `#FF0000`.
- `color_group$nhw` is `"#374E55FF"` with an alpha suffix. Drop the trailing
  `FF` in a prompt; image models treat an eight-digit hex as a typo.

## E. Known collisions

Check for these every time. Two roles sharing a hex may never appear in the same
panel; say so in the brief and let it fix the layout.

| Hex | Colliding roles | Resolution used before |
| --- | --- | --- |
| `#00A087` | `AA-Independent` (stage-local) and "power >= 0.80, adequately powered" | Cohort colors in panel a only; the power ramp in panels b and d only |
| `#F39B7F` | `AA-Independent` (global) and "power 0.20-0.50, low" | Do not mix the global cohort trio with a power ramp |
| `#E64B35` | `AA-Primary` (global), CM, and the generic accent | Inside a cohort figure, choose a different existing accent role from the stage or project palette and record that role's source |
| `#B09C85` | `NHW` (stage-local `COHORT_COLORS`) and `color_age_dx_hist`, the low-emphasis category | Never put a stage-local NHW mark and an age-at-diagnosis category in the same panel |

## F. The 8 echo traits

Use the short names on figures, in this order, matching `AGENTS.md`.

`EF`, `GLS`, `Ee`, `LAVI`, `LVESV`, `LVEDV`, `LVMI`, `RWT`

Never put the raw column names (`Priority_EF_Value_measure`) on a figure.

## G. Where the rendered image goes

| Figure | Path |
| --- | --- |
| Cross-cutting, project-level | `imgs/<name>.png` |
| Stage-level or question-level | `results*/NN-stage/figures/<brief-stem>.png` |

Name the raster after its brief: `2026-08-24-my-question.DIAGRAM.md` renders to
`2026-08-24-my-question.DIAGRAM.png`. Never overwrite an accepted image; rename
the old one `<name>_stale-YYYY-MM-DD.png`.

## H. Case studies

These briefs are evidence and case studies, not contract templates. `SKILL.md`
wins wherever they differ from the current rules. Read the relevant case study
for its evidence, but do not copy obsolete structure or counts.

| File | Archetype | What remains useful |
| --- | --- | --- |
| [src_newwgs/04.03.02-gwas-replication/2026-08-24-aa-independent-replication-power.DIAGRAM.md](../../../../src_newwgs/04.03.02-gwas-replication/2026-08-24-aa-independent-replication-power.DIAGRAM.md) | Method schematic with supporting result panels | This brief predates the deterministic counting rule. Do not copy its element count. It is useful only for its verified-number inventory, recomputation, collision note, and render log. |
| [src_newwgs/03.01-prepare-model/DIAGRAM.md](../../../../src_newwgs/03.01-prepare-model/DIAGRAM.md) | Study design | This brief predates `## Audience`, `## Panel inventory`, and `## Honest scope`. It is useful only for label provenance and its documented source correction. |

Its predecessor, `power_stale-2026-08-25.DIAGRAM.md`, is kept in the same
directory as the retirement example. The superseding brief records each change
and its evidence under `## What supersedes the first draft`. In a new compliant
brief, keep that record after the required `## Audience` section rather than
copying the historical section order.

## I. Historical figure prompts

34 plans still carry an inline flowchart prompt block, spread over 16 different
heading spellings (`## Flowchart prompt`, `## 9. Flowchart image-generation
prompt (for gemini / gpt-image-2)`, `## 12. Pipeline Flowchart —
Image-Generation Prompt`, ...), so search for a heading containing both
"Flowchart" and "prompt" rather than one exact string.
`src/05-gwas-local-ancestry/12-figure-spec.md` uses a third naming convention.
These are read-only history: read them when resuming a stage, do not
rename them, do not migrate them in bulk. Extract one into a `DIAGRAM.md` only
when that specific figure is being rewritten.
