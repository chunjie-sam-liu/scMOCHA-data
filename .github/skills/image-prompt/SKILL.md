---
name: image-prompt
description: Write image-generation prompts for scientific figures in the Nature Reviews Genetics visual style, rendered by gpt-image-2, nano banana, or a similar model. Use when producing a pipeline flowchart, study-design figure, method schematic, result schematic, or conceptual overview for a paper, talk, slide, or stage plan; when creating or editing a DIAGRAM.md; when a PLAN.md needs its figure brief; when a prompt renders as a flowchart, dashboard, slide, or data table instead of a journal figure; when a render has garbled formulas, dropped labels, wrong colors, or invented numbers; or when deciding whether a generated image may stand in for a real plot. Covers the DIAGRAM.md brief format, the Mermaid workflow view, the NRG visual language and glyph table, the five archetypes, the word budget, and palette binding.
---

# Scientific figure prompts

Produce a **brief**, not a prompt string. The brief is a repository file that
records what the figure must say, which numbers are real, which palette it
inherits, what the reader sees as a Mermaid workflow, and only then the literal
text to paste into the image model. The prompt is the last section, not the
whole document.

The brief serves two readers. A person opens it and reads the logic from the
`## Structure preview` Mermaid block without rendering anything. The image
model receives only the fenced block under `## Prompt` and renders a figure in
the _Nature Reviews Genetics_ style described in
[references/nrg-style.md](./references/nrg-style.md). Both views are built from
the same `## Canvas inventory`, so they cannot drift apart.

An image model cannot be trusted to invent factual content. Every factual or
project-specific label, identifier, scientific or data number, cohort name,
threshold, and hex code is copied from a file that already exists in the
repository. Authored headings, explanatory synthesis, and layout choices such
as aspect ratio, panel proportions, and axis framing are design choices, not
sourced facts; they still must not imply unsupported findings.

Related skills: `analysis-pipeline` (where the brief lives inside a stage),
`palette` (choosing colors), `data-verification` (proving the numbers are real).

**Project bindings.** This skill is portable and names no project value. The
concrete color files, cohort palettes, neutral roles, trait short names, and
rendered-image paths come from the repository's own bindings file, normally
`.github/instructions/diagram-bindings.instructions.md` or the repository
`AGENTS.md`. Where the repository and this skill disagree, the repository wins.
If the repository has no bindings file, ask for the values instead of inventing
them.

---

## 1. File convention

The brief is its own file, next to the plan it belongs to.

```
src/NN-stage-name/
  PLAN.md                                  links to DIAGRAM.md
  PROGRESS.md
  DECISION.md
  DIAGRAM.md                               stage-level figure: the whole stage
  {date}-{short-title}.PLAN.md
  {date}-{short-title}.DIAGRAM.md          figure for that one question
  {date}-{short-title}.PROGRESS.md
  {date}-{short-title}.DECISION.md
```

- Stage-level figure -> `DIAGRAM.md`.
- Figure for one specific question -> `{date}-{short-title}.DIAGRAM.md`, sharing
  the exact `{date}-{short-title}` stem of its `.PLAN.md` / `.PROGRESS.md` /
  `.DECISION.md` set.
- More than one figure for the same question -> `{stem}.{figure-key}.DIAGRAM.md`,
  where `{figure-key}` is lowercase kebab-case naming the figure
  (`power`, `cohort-flow`).
- Never inline a prompt inside `PLAN.md`. `PLAN.md` carries a one-line link
  under a `## Figure` heading and nothing else.
- Never name it `figure-spec.md`, `NN-figure-spec.md`, or `flowchart.md`. Those
  are historical names; leave existing ones alone, do not create new ones.

Existing plans that still hold an inline `## Flowchart prompt` block are
read-only history. Do not migrate them in bulk. Extract one only when you are
already rewriting that figure.

## 2. Required brief structure

Write these sections in this order. Skipping a section is how a figure ends up
with a wrong number in it, or with a function name on the canvas.

| Section                                | Content                                                                                                                                                                                                                                               |
| -------------------------------------- | ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| Title + one paragraph                  | What the figure renders, for which model, and which script or plan it belongs to. Link the companions.                                                                                                                                                |
| `## Audience`                          | The first H2 after the title and introductory paragraph. State `paper` or `discussion`, the one primary archetype, and the word and number count per panel; point to `## Canvas inventory`. See section 3 rules 1 and 2.                              |
| `## What the figure must communicate`  | Numbered list of the claims a reader must leave with. Name which one is the punchline and where on the canvas it sits.                                                                                                                                |
| `## The formula` _(if any)_            | ASCII-only, in a fenced block, exactly as it must appear on the page, plus a short annotation block defining every symbol.                                                                                                                            |
| `## Worked example` _(method figures)_ | One real unit carried through every stage of the calculation, as a table.                                                                                                                                                                             |
| `## Verified number inventory`         | Every factual or project-specific label, identifier, scientific or data number, cohort name, threshold, and hex in the prompt or the caption, with the file it was read from, the date, and the rounding rule where a value is rounded on the canvas. |
| `## Palette`                           | Role -> hex -> source file table: the four to six hues, what each one means, and the neutrals. Flag any two roles sharing a hex and say how the layout keeps them apart. Name the one emphasis point and its device.                                  |
| `## Canvas inventory`                  | Per panel: each shape or glyph with its fill role and exact quoted label, each arrow with its quoted process label, each annotation line; then that panel's word count and number count against rule 2.                                               |
| `## Structure preview`                 | A Mermaid diagram of the same content, so a reader follows the workflow and the logic in the editor without rendering an image. Required. See rule 12.                                                                                                |
| `## Caption`                           | The legend a reader gets beneath the image: exact counts, thresholds, software and versions, job IDs, and every value the canvas shows rounded. Code-level detail lives here, never on the canvas. See rule 13.                                       |
| `## Honest scope`                      | What the generated raster may and may not be used for.                                                                                                                                                                                                |
| `## Prompt`                            | The primary literal block to paste, fenced as ```text`: Full for `paper`, Discussion for `discussion`.                                                                                                                                                |
| Required lower-density prompt          | `## Leaner variant` for `paper`; `## One-message variant` for `discussion`. Additional tiers are optional. See the density ladder in the style contract.                                                                                              |

## 3. Hard rules

1. **Declare the audience first.** After the brief title and introductory
   paragraph, its first H2 is `## Audience`, naming one of two modes. This is
   not a preference: it sets the budget in rule 2 and decides which prompt tier
   is the primary deliverable. It reports the per-panel word and number counts
   and points to `## Canvas inventory`.
   - `discussion` — a figure for sharing, a talk, a message thread, or a
     conversation with a collaborator. One idea, read in ten seconds, still
     legible projected or at half width. **This is the default.**
   - `paper` — a manuscript figure. May carry a full derivation, a worked
     example, and several result panels, because the reader has the caption and
     unlimited time.
2. **Respect the word budget.** An NRG figure is sparse; the caption does the
   explaining. Before writing any prompt prose, list every string that will
   appear in each panel in `## Canvas inventory`, count its words and its
   numbers, and record the totals in `## Audience`.
   - A word is any whitespace-separated token in a quoted label, including
     numbers and units. `"n = 24"` is three words.
   - A number is any numeric value shown on the canvas. `"n = 24"` is one
     number; `"1.05 M"` is one number.
   - Panel titles count; panel letters and the caption do not.

   Count formulas separately. Each independently displayed equation, model
   relation, or numeric calculation expression is one formula. Wrapping one
   expression across lines does not add another formula; a symbol glossary with
   no operators does not count. Only formulas intended to appear on the canvas
   count.

   | Mode         | Words per panel | Numbers per panel | Formulas per panel | Panels                          |
   | ------------ | --------------- | ----------------- | ------------------ | ------------------------------- |
   | `discussion` | about 25        | 1                 | 1                  | 1 to 2                          |
   | `paper`      | about 40        | 3                 | 2                  | up to 4, plus at most one inset |

   Over budget means moving detail to `## Caption`, cutting content, or
   splitting the figure into two. It never means shrinking the type, adding
   "make the text clearer" to the prompt, or packing values into a table on the
   canvas. A value shown rounded on the canvas (`~860 k`, `1.05 M`) has its
   exact form and rounding rule in the inventory and the caption.

3. **Verified factual content only.** Every factual or project-specific label,
   identifier, scientific or data number, cohort name, threshold, and hex in
   the prompt is recorded in the inventory with its source path and read date.
   Never round a number carried over from a narrative, an earlier draft, or
   another document. If a source document disagrees with the data, put the
   correction in the brief with the evidence, and use the corrected value on
   the page. Authored headings, explanatory synthesis, aspect ratio, panel
   proportions, and axis framing are documented design choices rather than
   sourced facts.
4. **Palette bound to a source.** Take project-role hexes from the color file
   of the track the figure belongs to, the same file its real plots source.
   A stage never owns a palette of its own; when the repository bindings name
   a frozen legacy exception, they say so. Every palette row requires a
   source; neutral rows and documented repository-specific substitutions may
   cite the repository bindings file. If two roles
   resolve to the same hex, they may never appear in the same panel; state that
   constraint in the brief and let it drive the layout. A collision override
   must use another existing palette role and record its exact source; never
   invent a replacement hex.
5. **Formulas stay legible and always carry an ASCII fallback.** Numbers,
   variable names, cohort names, trait names and variant IDs are ASCII, always.
   Typeset subscripts are allowed only where the source ASCII has an underscore;
   the brief must hold that same formula in plain ASCII so the fallback is one
   edit away. Greek letters, summation or product symbols, a square-root glyph
   and LaTeX markup are still banned: those garble. Use the function name from
   the code (`pnorm`, not `Phi`) and spell out multiplication with `*`. Define
   every retained symbol in a symbolic formula. A fully numeric discussion
   example may omit a symbol glossary only when no undefined symbols remain.
6. **Opening style block, complete.** Every prompt starts with the full
   opening block from the style contract: flat filled shapes with no outline,
   four to six muted hues, the glyph sentence filled in for this figure,
   dark-grey arrows with process labels beside them, sentence-case sans-serif
   in two sizes, and the forbidden forms (gradients, shading, shadows, 3D,
   glow, texture, photographic elements, cartoon characters, twisted DNA
   helices, laboratory scenes, small caps, letter-spacing, monospace). See
   [references/style-contract.md](./references/style-contract.md) section A.
7. **Honest scope.** A generated image is a schematic. It may illustrate a
   method, a design, or the shape of a result. It may not replace a real plot.
   Name the real output file the figure must not be substituted for.
8. **Message before layout.** Write `## What the figure must communicate` before
   writing any prompt text. The layout exists to serve the numbered claims; if
   a panel serves none of them, delete the panel. A `discussion` figure carries
   at most three claims.
9. **One archetype per figure.** Pick from section 5 and follow its skeleton. Do
   not fuse a method schematic and a pipeline flowchart into one canvas. A
   figure has one primary archetype. Supporting panels are allowed only where
   that archetype's skeleton permits them and only when they serve a listed
   claim. A method panel with a permitted supporting result panel is still one
   method-schematic archetype, not two archetypes.
10. **Complete prompt tiers.** Every brief contains primary `## Prompt` plus the
    next lower density rung. For `paper`, `## Prompt` is Full and
    `## Leaner variant` is required. For `discussion`, `## Prompt` is
    Discussion and `## One-message variant` is required. Additional tiers are
    optional. Every standalone prompt, including every variant, repeats the
    complete opening negative-constraint block and closing impression block;
    abbreviated variants are forbidden.
11. **Glyphs are semantic, flat, and from the table.** Simplified flat
    biological glyphs are what make a genetics figure recognisable, and they
    are allowed: a CpG lollipop, an array grid, a chromosome bar, an IDAT page,
    a sample circle, a read stack, a matrix. Each is a single flat fill in the
    hue of its class, placed where a reader would otherwise need a word to know
    what kind of thing sits there. Photographic, 3D, cartoon, or laboratory
    imagery is forbidden, as is any glyph that only illustrates a word already
    present. The table and the forbidden list are in
    [references/nrg-style.md](./references/nrg-style.md) section B4. Name the
    two or three glyphs the figure uses in the opening block.
12. **Ship a Mermaid structure preview.** Every brief carries a
    `## Structure preview` section holding one Mermaid block built from
    `## Canvas inventory`: one `subgraph` per panel, one node per entity or
    glyph, one labelled edge per process arrow, the exact quoted label
    strings, and `classDef` fills bound to the same hexes as `## Palette`. It
    renders in the editor and on GitHub, so a reader follows the workflow and
    the logic of the figure without generating an image, and the structure,
    wording, and color plan are checked before a render is paid for. It is the
    reader's view of the logic, never the deliverable and never a substitute
    for the rendered figure. Template and rules:
    [references/structure-preview.md](./references/structure-preview.md).
13. **Nothing from the code on the canvas.** Function calls, package names,
    versions, object classes, file paths, job IDs, run times, and parameter
    syntax never appear in a panel. The canvas names the thing in plain words
    ("Raw intensities", "Detection P-value filter", "Mapped to hg38"); the
    `## Caption` section carries the exact calls, thresholds, versions, and
    counts. The translation table is in
    [references/nrg-style.md](./references/nrg-style.md) section D.

## 4. Workflow

1. **Set the audience.** `discussion` or `paper`. After the title and
   introductory paragraph, write `## Audience` as the first H2; it fixes the
   budget for every later step.
2. **Name the claims.** Write the numbered `must communicate` list. Mark the
   punchline and the canvas position where the eye lands last. Past three
   claims in a `discussion` figure, start a second figure.
3. **Harvest the numbers.** Read the real outputs; build the inventory table
   with paths and the read date. Recompute rather than copy where cheap.
4. **Translate code into canvas words.** For each step, file, and count the
   script produces, decide what the canvas says (a plain-word label, a glyph,
   one rounded number) and what the caption says (the call, the version, the
   exact value). Use the table in
   [references/nrg-style.md](./references/nrg-style.md) section D. Everything
   that lands in the caption column goes into `## Caption`, not the prompt.
5. **Bind the palette.** Choose four to six hues from the project color file,
   one per kind of thing. Check for collisions. Name the emphasis point and
   whether it is carried by a warm hue, a heavier outline, or size.
6. **Pick the archetype** and its canvas from
   [references/archetypes.md](./references/archetypes.md).
7. **Write the canvas inventory and count it.** Per panel: every shape or
   glyph with its quoted label, every arrow with its quoted process label,
   every annotation line. Count words and numbers against rule 2. Cut or move
   to the caption here, before any prompt prose exists.
8. **Draw the Mermaid structure preview** from the inventory and read it as a
   collaborator would. It is the cheapest check that can fail: if the workflow
   is not clear from the Mermaid, it will not be clear from the render either.
   Fix the structure here, then continue.
9. **Write the caption.** Exact counts, thresholds, software and versions, the
   rounding rule for every rounded canvas value.
10. **Write the prompt** using the style contract: opening block with the
    glyph sentence filled in, layout sentence, panel blocks in reading order,
    closing impression block.
11. **Write the required lower-density prompt** in the same pass, not after
    the first render fails. Repeat the complete opening and closing blocks in
    every standalone tier.
12. **Render and inspect.** Check against section 6.
13. **Link it** from `PLAN.md` under `## Figure`.

## 5. Archetypes

| Archetype           | Renders                                          | Default canvas                            |
| ------------------- | ------------------------------------------------ | ----------------------------------------- |
| Pipeline flowchart  | A stage's data flow, inputs to outputs           | landscape 3:2, single row or two branches |
| Study design        | Cohorts, roles, outcome families, model matrix   | landscape 3:2, columns left to right      |
| Method schematic    | How one calculation works, formula on the page   | landscape 16:9 or wide panel a            |
| Result schematic    | The shape of a result, expected against observed | landscape 3:2, 2-4 small panels           |
| Conceptual overview | The idea, for a talk opener; no paths, no counts | landscape 16:9                            |

The archetype names the shape; the audience names the density. Supporting
panels do not create a second archetype when the primary archetype's skeleton
permits them and they serve a listed claim. A method schematic for `discussion`
keeps the pipeline and one worked example and drops the inset, the annotation
layer, and every result panel.

Layout skeletons, per-archetype failure modes, what each must never contain, and
how to trim each one to the `discussion` budget:
[references/archetypes.md](./references/archetypes.md).

## 6. Inspecting the render

Do not accept the first image because it looks good. Check, in this order:

1. **Numbers.** Read every number in the image against the inventory table.
   Models silently alter digits. A wrong number is a retraction risk.
2. **Formulas.** Confirm no Greek letter or invented or altered symbol appeared,
   that typeset subscripts occur only where the source ASCII has underscores,
   and that operator order survived.
3. **Spelling of labels.** Cohort names, trait names, and variant IDs are
   frequently corrupted (`AA-Independant`, `echo.RWT` becoming `echo.RVVT`).
4. **Palette.** Sample the fills; confirm each hue means what the palette table
   says and emphasis appears only at the named point.
5. **Collisions.** Confirm no two roles sharing a hex ended up in one panel.
6. **Dropped content.** Walk `## Canvas inventory` and confirm every shape,
   glyph, arrow label, and annotation is present. Silent omission is the most
   common failure and it always hits the densest region, so check that region
   first.
7. **Style.** Run the seven-point NRG check in
   [references/nrg-style.md](./references/nrg-style.md) section E: no outlines
   on filled shapes, four to six hues, sentence case, word count, glyphs
   present and flat, thin grey arrows, no depth. A render that looks like a
   flowchart or a dashboard fails here even when every number is right.
8. **Density.** If any text is unreadable at 50% zoom, drop a rung on the ladder
   rather than adding a re-render instruction.

If two rounds of editing do not fix a garbled formula or label, stop editing and
drop a rung on the density ladder. Text density is almost always the cause.

## 7. Forbidden patterns

- Inventing a number, a sample size, or a p-value to fill a box.
- Using a generated raster as a manuscript figure that reports real data.
- Writing prompt prose before `## Audience` and the word count exist.
- Function calls, package names, versions, object classes, paths, job IDs, or
  run times on the canvas. They go in `## Caption`.
- Rows of numbers in a box. One or two numbers per panel, rounded for display;
  the exact values go in the caption.
- Dashboard vocabulary in a prompt: card, role bar, small-caps, letter-spaced,
  monospaced, strip, chip, pill, hairline. These words are what pull the render
  toward a docs site.
- Outlined boxes on white as the main construct. Shapes are filled, no stroke.
- Packing a worked example onto one long line. Stack it as short rows, one
  stage per row; a `discussion` figure may not use the single-line form at all.
- Greek letters, summation or product symbols, or LaTeX anywhere in a prompt.
- Photographic, 3D, or cartoon imagery; twisted DNA helices; laboratory
  scenes; any glyph placed where no biological thing sits.
- Gradients, drop shadows, 3D extrusion, glow, or paper texture.
- Hard-coded pixel dimensions in the prompt body; state aspect ratio in words so
  the same brief works across models.
- A new inline flowchart prompt inside `PLAN.md`.
- Shipping the Mermaid structure preview as the figure, or pasting Mermaid
  source into an image prompt.
- Writing the prompt before the Mermaid preview has been drawn and read.
- Deleting a superseded brief. Rename it `<name>_stale-YYYY-MM-DD.DIAGRAM.md`.

---

## References

- [references/nrg-style.md](./references/nrg-style.md) — what a Nature Reviews
  Genetics figure actually looks like: canvas, shapes, colour, the glyph table,
  arrows, typography, the word budget, the code-to-canvas translation table,
  and the style checklist for a render.
- [references/style-contract.md](./references/style-contract.md) — the reusable
  opening and closing blocks, the layout sentence, panel-block rules, the
  entity / process / group vocabulary, formula rules, the density ladder, and
  how gpt-image-2 and nano banana differ in practice.
- [references/structure-preview.md](./references/structure-preview.md) — the
  Mermaid workflow view: template, binding rules, and how to view it.
- [references/archetypes.md](./references/archetypes.md) — the five archetypes
  with layout skeletons and failure modes.

Project values (color files, hue sets, neutral roles, case studies) are not in
this skill. Read the repository bindings file named under "Project bindings"
above.
