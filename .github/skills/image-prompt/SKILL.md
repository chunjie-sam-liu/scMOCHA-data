---
name: image-prompt
description: Write image-generation prompts for scientific figures in the Nature Reviews Genetics visual style, rendered by gpt-image-2, nano banana, or a similar model. Use when producing a pipeline flowchart, study-design figure, method schematic, or conceptual overview for a paper, talk, or slide; when creating or editing a DIAGRAM.md; when a PLAN.md needs its figure brief; when a prompt renders as a flowchart, dashboard, or data table instead of a journal figure; when a render has garbled formulas, dropped labels, or wrong colors; or when deciding whether a generated image may stand in for a real plot. Covers the DIAGRAM.md brief format, the NRG visual language, the five archetypes, and palette binding.
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

Write these sections in this order. Skipping one is how a figure ends up with a
wrong number in it, or with a function name on the canvas.

1. **Title + one paragraph** -- what it renders, for which model, and which
   script or plan it belongs to. Link the companions.
2. `## Audience` -- `paper` or `discussion`, the one primary archetype, and the
   per-panel word and number counts. Rules 1-2.
3. `## What the figure must communicate` -- numbered claims; name the punchline
   and where on the canvas it sits.
4. `## The formula` _(if any)_ -- ASCII only, fenced, plus a glossary defining
   every symbol.
5. `## Worked example` _(method figures)_ -- one real unit carried through every
   stage of the calculation, as a table.
6. `## Verified number inventory` -- every fact with the file it was read from,
   the date, and the rounding rule where the canvas rounds. Rule 3.
7. `## Palette` -- role -> hex -> source file; flag any two roles sharing a hex
   and say how the layout keeps them apart. Rule 4.
8. `## Canvas inventory` -- per panel: every shape, glyph, arrow, and annotation
   with its exact quoted label, then that panel's counts against rule 2.
9. `## Structure preview` -- the Mermaid view of the same content. Rule 12.
10. `## Caption` -- exact counts, thresholds, software versions, job IDs, and
    every value the canvas shows rounded. Rule 13.
11. `## Honest scope` -- what the generated raster may and may not be used for.
12. `## Prompt` -- the primary literal block to paste: Full for `paper`,
    Discussion for `discussion`.
13. The required lower-density variant: `## Leaner variant` for `paper`,
    `## One-message variant` for `discussion`. Rule 10.

What each section must contain, in full:
[references/brief-rules.md](./references/brief-rules.md) section A.

## 3. Hard rules

Thirteen rules govern every brief. The rule is stated here; the procedure for
satisfying it -- the counting definitions, the palette collision protocol, the
formula fallback contract, the prompt-tier rules -- is
[references/brief-rules.md](./references/brief-rules.md). Read that before
writing the first brief, and again whenever a render comes back over budget,
mis-coloured, or with garbled type.

1. **Declare the audience first.** The brief's first H2 is `## Audience`, naming
   `discussion` (**the default**: one idea, read in ten seconds, legible at half
   width) or `paper` (a manuscript figure, whose reader has the caption and
   unlimited time). It sets the rule 2 budget and decides which prompt tier is
   the primary deliverable.
2. **Respect the word budget.** List every canvas string in
   `## Canvas inventory`, count its words, numbers, and formulas, and record the
   totals in `## Audience`.

   | Mode         | Words per panel | Numbers per panel | Formulas per panel | Panels                          |
   | ------------ | --------------- | ----------------- | ------------------ | ------------------------------- |
   | `discussion` | about 25        | 1                 | 1                  | 1 to 2                          |
   | `paper`      | about 40        | 3                 | 2                  | up to 4, plus at most one inset |

   Over budget means moving detail to `## Caption`, cutting content, or
   splitting the figure in two. It never means shrinking the type or packing
   values into a table on the canvas.

3. **Verified factual content only.** Every factual or project-specific label,
   identifier, number, cohort name, threshold, and hex is recorded in the
   inventory with its source path and read date. Never round a number carried
   over from a narrative or an earlier draft.
4. **Palette bound to a source.** Project-role hexes come from the color file of
   the track the figure belongs to, the same file its real plots source. Every
   palette row cites a source. Two roles resolving to the same hex may never
   appear in one panel. Never invent a replacement hex.
5. **Formulas stay legible and always carry an ASCII fallback.** No Greek, no
   summation or product symbols, no square-root glyph, no LaTeX -- those garble.
   Use the function name from the code (`pnorm`, not `Phi`), spell out
   multiplication with `*`, and define every retained symbol.
6. **Opening style block, complete.** Every prompt starts with the full opening
   block from [references/style-contract.md](./references/style-contract.md)
   section A.
7. **Honest scope.** A generated image is a schematic and may not replace a real
   plot. Name the real output file it must not be substituted for.
8. **Message before layout.** Write `## What the figure must communicate` before
   any prompt text, and delete any panel that serves none of its claims. A
   `discussion` figure carries at most three claims.
9. **One archetype per figure.** Pick from section 5 and follow its skeleton.
   Never fuse a method schematic and a pipeline flowchart onto one canvas.
10. **Complete prompt tiers.** Primary `## Prompt` plus the next lower rung:
    `## Leaner variant` for `paper`, `## One-message variant` for `discussion`.
    Every variant repeats the complete opening and closing blocks; abbreviated
    variants are forbidden.
11. **Glyphs are semantic, flat, and from the table.** Name the two or three the
    figure uses in the opening block. Table and forbidden list:
    [references/nrg-style.md](./references/nrg-style.md) section B4.
12. **Ship a Mermaid structure preview.** `## Structure preview` is built from
    `## Canvas inventory` and checked before a render is paid for. It is the
    reader's view of the logic, never the deliverable. Template:
    [references/structure-preview.md](./references/structure-preview.md).
13. **Nothing from the code on the canvas.** No function calls, package names,
    versions, file paths, job IDs, or parameter syntax in a panel. The canvas
    names the thing in plain words; `## Caption` carries the exact calls,
    thresholds, and counts. Translation table:
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

- [references/brief-rules.md](./references/brief-rules.md) — the thirteen hard
  rules in full: the word/number/formula counting definitions, the palette
  collision protocol, the ASCII formula fallback contract, and the prompt-tier
  rules.
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
