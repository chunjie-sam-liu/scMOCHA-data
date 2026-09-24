# Reference: the thirteen brief rules, in full

`SKILL.md` section 3 states each rule in one or two lines. This file holds the
procedure for satisfying it: the counting definitions, the palette collision
protocol, the formula fallback contract, and the prompt-tier rules. Read it
before writing the first brief in a repository, and again whenever a render
comes back over budget, mis-coloured, or with garbled type.

---

## A. The required brief structure

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

---

## 1. Declare the audience first

After the brief title and introductory paragraph, its first H2 is `## Audience`,
naming one of two modes. This is not a preference: it sets the budget in rule 2
and decides which prompt tier is the primary deliverable. It reports the
per-panel word and number counts and points to `## Canvas inventory`.

- `discussion` — a figure for sharing, a talk, a message thread, or a
  conversation with a collaborator. One idea, read in ten seconds, still
  legible projected or at half width. **This is the default.**
- `paper` — a manuscript figure. May carry a full derivation, a worked
  example, and several result panels, because the reader has the caption and
  unlimited time.

---

## 2. Respect the word budget

An NRG figure is sparse; the caption does the explaining. Before writing any
prompt prose, list every string that will appear in each panel in
`## Canvas inventory`, count its words and its numbers, and record the totals in
`## Audience`.

- A word is any whitespace-separated token in a quoted label, including numbers
  and units. `"n = 24"` is three words.
- A number is any numeric value shown on the canvas. `"n = 24"` is one number;
  `"1.05 M"` is one number.
- Panel titles count; panel letters and the caption do not.

Count formulas separately. Each independently displayed equation, model
relation, or numeric calculation expression is one formula. Wrapping one
expression across lines does not add another formula; a symbol glossary with no
operators does not count. Only formulas intended to appear on the canvas count.

| Mode         | Words per panel | Numbers per panel | Formulas per panel | Panels                          |
| ------------ | --------------- | ----------------- | ------------------ | ------------------------------- |
| `discussion` | about 25        | 1                 | 1                  | 1 to 2                          |
| `paper`      | about 40        | 3                 | 2                  | up to 4, plus at most one inset |

Over budget means moving detail to `## Caption`, cutting content, or splitting
the figure into two. It never means shrinking the type, adding "make the text
clearer" to the prompt, or packing values into a table on the canvas. A value
shown rounded on the canvas (`~860 k`, `1.05 M`) has its exact form and rounding
rule in the inventory and the caption.

---

## 3. Verified factual content only

Every factual or project-specific label, identifier, scientific or data number,
cohort name, threshold, and hex in the prompt is recorded in the inventory with
its source path and read date. Never round a number carried over from a
narrative, an earlier draft, or another document. If a source document disagrees
with the data, put the correction in the brief with the evidence, and use the
corrected value on the page. Authored headings, explanatory synthesis, aspect
ratio, panel proportions, and axis framing are documented design choices rather
than sourced facts.

---

## 4. Palette bound to a source

Take project-role hexes from the color file of the track the figure belongs to,
the same file its real plots source. A stage never owns a palette of its own;
when the repository bindings name a frozen legacy exception, they say so. Every
palette row requires a source; neutral rows and documented repository-specific
substitutions may cite the repository bindings file.

If two roles resolve to the same hex, they may never appear in the same panel;
state that constraint in the brief and let it drive the layout. A collision
override must use another existing palette role and record its exact source;
never invent a replacement hex.

---

## 5. Formulas stay legible and always carry an ASCII fallback

Numbers, variable names, cohort names, trait names and variant IDs are ASCII,
always. Typeset subscripts are allowed only where the source ASCII has an
underscore; the brief must hold that same formula in plain ASCII so the fallback
is one edit away. Greek letters, summation or product symbols, a square-root
glyph and LaTeX markup are still banned: those garble. Use the function name
from the code (`pnorm`, not `Phi`) and spell out multiplication with `*`. Define
every retained symbol in a symbolic formula. A fully numeric discussion example
may omit a symbol glossary only when no undefined symbols remain.

---

## 6. Opening style block, complete

Every prompt starts with the full opening block from the style contract: flat
filled shapes with no outline, four to six muted hues, the glyph sentence filled
in for this figure, dark-grey arrows with process labels beside them,
sentence-case sans-serif in two sizes, and the forbidden forms (gradients,
shading, shadows, 3D, glow, texture, photographic elements, cartoon characters,
twisted DNA helices, laboratory scenes, small caps, letter-spacing, monospace).
See [style-contract.md](./style-contract.md) section A.

---

## 7. Honest scope

A generated image is a schematic. It may illustrate a method, a design, or the
shape of a result. It may not replace a real plot. Name the real output file the
figure must not be substituted for.

---

## 8. Message before layout

Write `## What the figure must communicate` before writing any prompt text. The
layout exists to serve the numbered claims; if a panel serves none of them,
delete the panel. A `discussion` figure carries at most three claims.

---

## 9. One archetype per figure

Pick from `SKILL.md` section 5 and follow its skeleton. Do not fuse a method
schematic and a pipeline flowchart into one canvas. A figure has one primary
archetype. Supporting panels are allowed only where that archetype's skeleton
permits them and only when they serve a listed claim. A method panel with a
permitted supporting result panel is still one method-schematic archetype, not
two archetypes.

---

## 10. Complete prompt tiers

Every brief contains primary `## Prompt` plus the next lower density rung. For
`paper`, `## Prompt` is Full and `## Leaner variant` is required. For
`discussion`, `## Prompt` is Discussion and `## One-message variant` is
required. Additional tiers are optional. Every standalone prompt, including
every variant, repeats the complete opening negative-constraint block and
closing impression block; abbreviated variants are forbidden.

---

## 11. Glyphs are semantic, flat, and from the table

Simplified flat biological glyphs are what make a genetics figure recognisable,
and they are allowed: a CpG lollipop, an array grid, a chromosome bar, an IDAT
page, a sample circle, a read stack, a matrix. Each is a single flat fill in the
hue of its class, placed where a reader would otherwise need a word to know what
kind of thing sits there. Photographic, 3D, cartoon, or laboratory imagery is
forbidden, as is any glyph that only illustrates a word already present. The
table and the forbidden list are in [nrg-style.md](./nrg-style.md) section B4.
Name the two or three glyphs the figure uses in the opening block.

---

## 12. Ship a Mermaid structure preview

Every brief carries a `## Structure preview` section holding one Mermaid block
built from `## Canvas inventory`: one `subgraph` per panel, one node per entity
or glyph, one labelled edge per process arrow, the exact quoted label strings,
and `classDef` fills bound to the same hexes as `## Palette`. It renders in the
editor and on GitHub, so a reader follows the workflow and the logic of the
figure without generating an image, and the structure, wording, and color plan
are checked before a render is paid for. It is the reader's view of the logic,
never the deliverable and never a substitute for the rendered figure. Template
and rules: [structure-preview.md](./structure-preview.md).

---

## 13. Nothing from the code on the canvas

Function calls, package names, versions, object classes, file paths, job IDs,
run times, and parameter syntax never appear in a panel. The canvas names the
thing in plain words ("Raw intensities", "Detection P-value filter", "Mapped to
hg38"); the `## Caption` section carries the exact calls, thresholds, versions,
and counts. The translation table is in [nrg-style.md](./nrg-style.md) section D.
