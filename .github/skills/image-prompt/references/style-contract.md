# Reference: the Nature-style contract

The reusable text that makes a generated figure look like it came out of
_Nature Reviews Genetics_ rather than out of a diagramming tool. Copy the
blocks verbatim and fill the angle-bracket slots.

The visual language these blocks encode is described in
[nrg-style.md](./nrg-style.md). Read it once before writing a first prompt,
and read section E of it whenever a render comes back looking like a flowchart,
a slide, or a dashboard.

---

## A. Opening style block

Every prompt starts with this. It is the single highest-leverage paragraph.
The glyph sentence is filled per figure from the glyph table in
`nrg-style.md` section B4; never leave it generic and never leave it out.

```text
A flat two-dimensional scientific schematic in the visual style of a Nature
Reviews Genetics figure: an illustrated, colourful, precisely drawn editorial
diagram, not a slide, not a dashboard, not a flowchart of boxes. <orientation>,
<aspect ratio> aspect ratio, white background, no outer frame, generous even
whitespace.

Shapes are rounded rectangles filled with flat muted colour and no outline.
Related shapes may sit on a very pale tinted rounded field of their own hue.
Four to six coordinated muted hues, each standing for one kind of thing and
used consistently across every panel. No gradients, no shading, no shadows, no
3D, no glow, no texture.

Simplified flat biological glyphs are used where they stand for a thing:
<two or three glyphs from the glyph table, each described in one clause>. Each
is a single flat fill in the hue of its class. No twisted DNA helices, no
photographic or 3D elements, no cartoon characters, no laboratory scenes.

Connectors are thin dark-grey #2B3238 lines with small solid triangular heads,
all the same size. Process labels sit beside the arrow in small dark-grey text,
not inside a box.

Typography: one neutral sans-serif such as Helvetica, sentence case throughout,
two sizes only, dark grey #2B3238 on light fills and white on dark fills. No
small caps, no letter-spacing, no monospace, no all-caps. Very few words: short
labels only.

Panel letters are lowercase bold "a", "b", "c" at the top-left of each panel.
```

Drop the last sentence for a single-panel figure. Never drop the shapes
paragraph or the typography paragraph: those two are what separate an NRG
figure from a flowchart.

## B. Layout sentence

One sentence immediately after the style block, before any panel description.
State proportions as percentages of the canvas, never in pixels.

```text
Layout: one wide panel a occupying the full width of the top 55 percent of the
canvas, and three equal-width panels b, c, d in a single row across the bottom
45 percent.
```

## C. Panel blocks

One block per panel, in reading order, each opening with
`PANEL <letter>, titled "<panel title>".` Then describe, in this order:

1. The overall shape (a left-to-right chain, two bars, a matrix, one curve).
2. Axes and their end labels, if it is a plot.
3. Each shape or glyph with its fill hex and its exact label text in double
   quotes.
4. Each arrow, with the process label that sits beside it in double quotes.
5. Any single line of small annotation text.

Rules that matter more than they look:

- **Quote every string that must appear.** Unquoted label text gets paraphrased.
- **Give a hex for every fill.** "Blue" becomes a different blue each render.
- **Say the corner radius and arrowhead size once**, in the opening block, not
  per element.
- **Describe position relatively** ("to the right of the array", "in the lower
  right of the plot area"), never in coordinates.
- **Name the scale relationship** when it carries meaning: "a very short bar,
  roughly one tenth the height of the left bar" survives; "26 versus 242.8"
  does not.
- **Labels name the thing, never the code.** "Raw intensities", not
  `RGChannelSet`; "Detection P-value filter", not `detectionP()`. Function
  names, package names, versions, job IDs, run times, paths and object classes
  belong in `## Caption`, never in a panel block.
- **Stack, never run on.** A worked example or a chain of intermediate values
  goes one stage per row. A single long line of tokens is where the model
  drops a segment and misspells a word, and both are easy to miss because the
  line still looks full.
- **Count words before you write.** Use the budget in `SKILL.md` section 3
  rule 2: about 25 words per panel for `discussion`, about 40 for `paper`. A
  panel that needs more is a panel whose explanation belongs in the caption.

## D. Shape, glyph and arrow vocabulary

Three constructs carry an NRG schematic. Reuse this wording so figures across
a project match.

**Entity** - a thing that exists: a cohort, a file, a matrix, a locus.

```text
a rounded rectangle filled flat <HEX> with no outline, containing the label
"<NAME>" in white
```

or, where the glyph table has a form for that kind of thing,

```text
a flat <glyph description> in <HEX>, labelled beneath in small dark-grey text
"<NAME>"
```

**Process** - something done to an entity. It is never a box.

```text
a thin #2B3238 arrow from <A> to <B> with the label "<VERB PHRASE>" in small
dark-grey text sitting just above the arrow
```

**Group** - related entities that belong together.

```text
sitting together on a very pale rounded field, a 10 percent tint of <HEX>, with
no outline
```

The glyph table, the forbidden glyph list, and the rule for when a glyph is
semantic rather than decorative are in `nrg-style.md` section B4. A glyph is
allowed where a reader would otherwise need a word to know what kind of thing
sits there; it is not allowed as an illustration of a word already present.

## E. Formula rules

Formulas are the most fragile part of any scientific figure prompt.

**Always ASCII, no exceptions:** numbers, variable names, cohort names, trait
names, variant IDs, function names, and every quoted label. These are what a
reader checks against the data, so a corrupted character here is a factual
error, not a cosmetic one.

**Still banned outright:** Greek letters, summation and product symbols, the
square-root glyph, LaTeX markup. Current models garble all of them.

**Allowed, with a condition:** typeset subscripts only where the source ASCII
has an underscore. `SE_expected` set as SE with a subscript reads better than
the flat form, and current models render it cleanly. The brief must hold the
same formula in plain ASCII, so falling back costs one edit turn rather than a
rewrite. Inspection rejects invented or altered subscripts, not valid
subscripts corresponding to source underscores.

- Use `sqrt( )`, `abs( )`, `*`, `/`, and explicit parentheses.
- Use the function name that appears in the code (`pnorm`, not `Phi`), so a
  reader can grep the script.
- Write variables with an underscore suffix in the brief: `SE_disc`, `N_rep`,
  `f_disc`. The model may set the suffix as a subscript; the meaning is
  unchanged and the fallback is already written down.
- Put every formula inside a fenced block in the brief first, then paste that
  exact text into the prompt.
- Add a symbol annotation block of short lines under every symbolic formula,
  one retained symbol per line, so the figure is self-contained. A fully
  numeric discussion example may omit the glossary only when no undefined
  symbols remain.
- Explain any term that looks missing. A factor that cancels is absent by
  design, and a reader will otherwise assume it was dropped by mistake.
- A `discussion` figure gets one formula, in one panel. If the argument needs
  two, it is a `paper` figure.
- In every prompt that contains a formula, place this safeguard immediately
  before the closing impression block:

```text
Every number, variable name, cohort name and quoted label is plain ASCII exactly
as written. No Greek letters, no summation symbols, no LaTeX styling. Subscripts
are acceptable only where an underscore appears in the text above.
```

If a formula does come back garbled, re-issue that panel with this line
appended rather than rewriting the brief:

```text
Set all formulas in plain ASCII exactly as written, in the same sans-serif as
the labels, with no subscripts or superscripts of any kind.
```

## F. Closing impression block

Last paragraph of every prompt. It fixes the details that per-element
instructions cannot.

```text
Overall impression: a warm, colourful, calm and precise illustration from a
Nature Reviews journal. Flat colour carries the shapes, whitespace carries the
grouping, glyphs carry the biology, and the caption carries the explanation.
The figure is sparse enough to read in ten seconds. Emphasis falls on <the one
place> through <a warmer hue among cool neighbours / a heavier outline / size>,
and nowhere else. No title text outside the panel letters and their panel
headings.
```

Naming exactly where emphasis is allowed, and by which device, is what keeps
the figure calm. Without it the model spreads saturated colour across every
panel and the punchline stops being visible.

## G. Density ladder

The rung is chosen by the audience, not by whether the last render failed.
Write the rung the audience calls for plus the one below it in the same pass.
For `paper`, `## Prompt` is Full and a complete `## Leaner variant` is required.
For `discussion`, `## Prompt` is Discussion and a complete
`## One-message variant` is required. Additional tiers are optional. Every
standalone tier, including every variant, repeats the complete opening style
block from section A and the complete closing impression block from section F;
never abbreviate a variant by referring back to another prompt.

| Rung        | Scope                                                                                    | Canvas      | Budget per panel    | Primary for                                     |
| ----------- | ---------------------------------------------------------------------------------------- | ----------- | ------------------- | ----------------------------------------------- |
| Full        | All panels, worked example, annotations                                                  | 3:2         | 40 words, 3 numbers | `paper`                                         |
| Leaner      | The primary panel or panels plus at most one supporting panel that carries the punchline | 3:2         | 40 words, 3 numbers | `paper`, when the full rung garbles             |
| Discussion  | One or two panels, one formula, stacked worked example, no inset                         | 3:2 or 16:9 | 25 words, 1 number  | `discussion`                                    |
| One-message | Single panel, one fact or one minimal method relation, no supporting panels              | 16:9        | 25 words, 1 number  | `discussion`, for a thread or a projected slide |

Going down a rung means removing panels and annotations, never shrinking type.
Never respond to a garbled or truncated render by adding "make the text clearer"
to the same prompt: the cause is density, and the prompt has to lose content.

A figure that started as `paper` and is now needed for a meeting does not get
resized. Render the discussion rung from the same brief; the numbers are already
verified, so the second render costs nothing but the prompt.

## H. Practical differences between models

Both models take the same brief. These are phrasing differences, not capability
claims; confirm current size and aspect-ratio support in the vendor
documentation before hard-coding anything.

| Concern        | Guidance                                                                                                                                                                                       |
| -------------- | ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| Canvas size    | State the aspect ratio in words ("landscape orientation, 3:2 aspect ratio"). Do not put pixel dimensions in the prompt body; set them in the API call or the UI.                               |
| Long prompts   | Both accept long structured prompts. Keep one blank line between panel blocks; unbroken walls of text lose their tail.                                                                         |
| Text rendering | Short quoted strings survive best. Long sentences inside a box get paraphrased. Break box content into several short quoted rows.                                                              |
| Iterating      | Prefer an edit turn over a re-render: state what must stay identical, then the one change ("keep the layout, colors and all other labels identical; change only the label in panel b to ..."). |
| Re-rolling     | If an edit turn introduces a new error, go back to the original prompt and change the brief instead of stacking edit turns.                                                                    |
| Determinism    | Neither is deterministic. Two renders of the same prompt differ. Keep the brief as the source of truth, and store the accepted image alongside it.                                             |

## I. Where to store the rendered image

- The accepted raster goes under `imgs/` or the stage's results figure
  directory, named after the brief (`2026-08-24-my-question.DIAGRAM.png` for
  `2026-08-24-my-question.DIAGRAM.md`).
- Never overwrite an accepted image. Version it or rename the old one
  `<name>_stale-YYYY-MM-DD`.
- Never let a generated raster replace a real plot produced by code. Name the
  real file in the brief's `## Honest scope` section so the boundary is on the
  record.

## J. The Nature visual language

The full description lives in [nrg-style.md](./nrg-style.md): canvas, shapes,
colour, the glyph table, arrows, typography, word budget, plots inside
schematics, depth, and emphasis. This section keeps only the two things a
prompt author reaches for mid-render.

### J1. Tints are not new colours

A bound hue is used at two strengths: full for the shape, and a pale wash of
about 10 percent for the field a related group sits on. State the level in
words: `a very pale 10 percent tint of #4E79A7`. Do not compute the tint and do
not write a new hex for it. A tint of a bound hex needs no extra row in
`## Palette` and does not violate hard rule 4. A genuinely different hue always
does.

### J2. What to fix when a render misses the style

| Render looks                       | Cause                                              | Fix                                                                                |
| ---------------------------------- | -------------------------------------------------- | ---------------------------------------------------------------------------------- |
| Flowchart, wireframe, "PowerPoint" | Outlined boxes on white                            | Filled rounded shapes with no outline; group by pale field or whitespace           |
| Dashboard, docs site, infographic  | Cards with side bars, small caps, monospace, chips | Remove those words from the prompt; sentence case; plain sans                      |
| Grey and lifeless                  | One accent on neutrals                             | Four to six muted hues, each meaning one thing                                     |
| Data table with boxes around it    | Too many words and numbers on the canvas           | Cut to the word budget; move exact values to `## Caption`                          |
| Not recognisably genetics          | No glyphs                                          | Add the two or three glyphs from the table that name what kind of thing sits where |
| Cheap, clip-art, 3D                | Wrong glyph forms                                  | Repeat the forbidden list; name the specific glyph to remove                       |
| Busy, no focal point               | Saturated colour everywhere                        | Name the one emphasis point and device in the closing block                        |

If a render fails the first two rows, grep the prompt for `card`, `role bar`,
`small-caps`, `monospaced`, `strip`, `chip`, `pill`, `hairline` and remove
every occurrence before re-rendering. If two rounds do not help, the cause is
density rather than styling: drop a rung on the ladder in section G.
