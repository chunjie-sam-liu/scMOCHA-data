# Reference: what a Nature Reviews Genetics figure actually looks like

A description of the _Nature Reviews Genetics_ (NRG) house style for
schematics, written so it can be pasted into a brief and checked against a
render. It describes character, not project values: every hex still comes from
the project color file named in the repository bindings.

Read section A first. It explains why a prompt could satisfy every rule in the
contract in force before 2026-09-05 and still render as a software dashboard
rather than a journal figure. Briefs written under that contract are history;
read them, do not copy their prompt vocabulary.

---

## A. Why earlier renders missed the style

The earlier contract encoded a **product-documentation / dashboard** look:
white cards with a coloured bar on the left edge, letter-spaced small-caps role
labels, monospaced code rows, count pills, one accent colour on grey. That is
the visual language of a SaaS docs site. It is clean, but it is not NRG.

| Earlier contract said                                      | NRG actually does                                                         |
| ---------------------------------------------------------- | ------------------------------------------------------------------------- |
| White card, `#C9D1D6` outline, thick vertical role bar     | Rounded rectangle filled with flat colour, **no outline**, no role bar    |
| Letter-spaced small-caps `"1 READ"`, `"INPUT"`             | Sentence case, plain sans, no letter-spacing, no small caps               |
| Monospaced rows: `minfi::detectionP(rg_set)`               | No code on the canvas, ever; the function name goes in the caption        |
| Three rows of numbers per step, `1,051,539 addresses x 24` | One short label, at most one number, rounded (`~1.05 M probes`, `n = 24`) |
| Neutral grey everything + one accent                       | Four to six coordinated muted hues; every class of thing owns a hue       |
| No DNA, no cells, no people, no glyphs                     | Stylised flat biological glyphs **are** the identity                      |
| Strips, chips, pills, hairline dividers                    | Not present; grouping is done by a pale tinted field or by whitespace     |
| 25 elements and 2 formulas per panel                       | Roughly 40 words per panel; the caption carries the explanation           |
| `"Overall impression: restrained, editorial"`              | Warm, colourful, illustrative, and still flat and precise                 |

Three root causes, in order of impact:

1. **Text density.** A figure that reads as a data table cannot look like NRG
   regardless of styling. NRG figures are sparse; the caption does the work.
2. **Glyph ban.** Forbidding every biological glyph removes the one thing that
   makes a genetics figure recognisably a genetics figure. NRG uses flat,
   simplified chromosomes, CpG lollipops, array grids, read stacks, and
   sample silhouettes constantly; what it never uses is 3D, photographic, or
   cartoon versions of them.
3. **UI vocabulary.** Cards with role bars, small-caps labels, and monospace
   rows are design-system idioms. Replacing them with filled rounded shapes and
   plain sentence-case labels changes the render more than any other edit.

---

## B. The visual identity

### B1. Canvas and layout

- White background, no outer frame, no panel boxes. Panels are compositions
  separated by whitespace, not by rules.
- Lowercase bold panel letters `a`, `b`, `c` at the top-left of each panel;
  nothing else is bold except an occasional short heading inside a panel.
- Landscape by default. Reading order left to right, then top to bottom.
- A group of related shapes may sit on a **pale tinted rounded field** (about
  10 percent tint of that group's hue) with no outline. That field is the only
  grouping device; never an outlined box.
- Whitespace is generous and even. Elements do not touch panel edges.

### B2. Shapes

- Entities are **rounded rectangles filled with flat colour and no stroke**.
  Text sits inside in white (on mid or dark fills) or dark grey (on pale fills).
- One corner radius across the figure, small relative to the shape.
- Shapes are sized to their label with even padding; widths match within a
  column.
- Processes are not boxes. A process is a labelled arrow, or a short label
  sitting on the arrow, or a glyph of the thing being done.

### B3. Colour

- **Colourful but muted.** Medium saturation, medium lightness, never neon and
  never pure primaries. Think teal, soft blue, coral, mustard, lavender, sage.
- **Four to six hues per figure**, each bound to one class of thing (a cohort,
  an array version, a data type, a stage). The same hue means the same thing
  across every panel.
- Each hue is used at two strengths: full for the shape, a pale tint for its
  grouping field or for a "secondary" state. No third strength, no gradient.
- Text is dark grey, not black. Arrows are mid or dark grey, not coloured,
  unless the arrow colour itself carries meaning.
- There is no single "accent". Emphasis is carried by a **warm hue against
  cool neighbours** (a coral shape among teals and blues), or by a heavier
  outline on one shape, or by size.
- Hex values come from the project color file. This section binds character;
  the bindings file binds numbers. If the project file has only cool hues,
  say so in the brief and choose emphasis by outline or size instead.

### B4. Glyph vocabulary

Flat, single-fill, simplified glyphs that stand for a thing at the place it
occurs. They are the recognisable NRG element and they are allowed. Each is
drawn in the hue of the class it belongs to.

| Thing                | Glyph                                                                                          |
| -------------------- | ---------------------------------------------------------------------------------------------- |
| DNA / genomic region | Two thin parallel lines or a flat single ribbon; **never a twisted 3D helix**                  |
| Chromosome           | Rounded vertical bar with a centromere pinch; a band shows a locus                             |
| CpG methylation      | Lollipop on a DNA line: stem plus filled circle = methylated, open circle = unmethylated       |
| Methylation array    | Small grid of squares, some filled, on a rounded plate                                         |
| Raw data file / IDAT | Rounded page shape with one folded corner, flat                                                |
| Sample / individual  | Simple filled circle, or a minimal rounded head-and-shoulders silhouette; a row means a cohort |
| Sequencing reads     | Short horizontal bars stacked under a genomic line                                             |
| Matrix / table       | Grid of cells with a darker header row and column                                              |
| Cluster / heatmap    | Small grid with cells in two or three tints                                                    |
| Distribution         | One smooth filled curve, no ticks                                                              |
| Genotype / variant   | Filled circle on the DNA line, or a short vertical tick                                        |
| Tissue / cell        | Single rounded blob with one nucleus circle, flat                                              |
| Filter / QC step     | Funnel outline, or a set of circles with some greyed out                                       |
| Model / equation     | Never a glyph; one short plain-text relation, ASCII                                            |

Forbidden, still: photographic textures, 3D rendering, twisted helices,
cartoon characters, laboratory scenes, microscopes, pipettes, lightbulbs,
magnifying glasses, badge-style checkmarks, and any glyph that decorates a
place where nothing biological is happening.

Rule of thumb: a glyph is allowed where a reader would otherwise need a word to
know **what kind of thing** sits there. It is not allowed as an illustration
of the word already present.

### B5. Arrows and connectors

- Thin lines, about one point, mid or dark grey, with small solid triangular
  heads. All heads the same size.
- Straight or gently curved; curves used to route around content, never for
  flourish.
- Dashed line for an optional or excluded path.
- T-bar for inhibition or exclusion.
- A process label sits **on or beside** the arrow in small dark-grey text, not
  in a box.

### B6. Typography

- One neutral sans-serif (Helvetica, Arial, or similar). Do not name a
  display or coding font.
- **Sentence case everywhere.** No small caps, no all-caps, no letter-spacing.
- Two sizes: label and small annotation. Panel letters are bold; almost nothing
  else is.
- No monospace anywhere on the canvas.
- Italic only for gene names.
- Dark grey text on white or pale fills; white text on mid and dark fills.

### B7. Numbers and words

- A panel carries **roughly 40 words**. If it needs more, the extra goes in the
  caption or the paired `.md`.
- One or two numbers per panel, the ones that carry the message. Rounded for
  display (`~860 k CpGs`, `n = 24`); the exact value and the rounding rule are
  recorded in the brief's inventory.
- Labels name the **thing**, not the code that made it: "Raw intensities", not
  `RGChannelSet`; "Detection P-value filter", not `detectionP()`.
- Software names, function calls, job IDs, run times, version strings, file
  paths and object classes never appear on the canvas.

### B8. Plots inside schematics

A plot in a schematic is a **shape**, not data: a smooth curve, three or four
dots, one dashed reference line. Axis lines are thin grey with at most the two
end labels; no tick marks, no gridlines, no legend box. If the reader needs the
numbers, the figure is a real plot and belongs to the r-figure skill.

### B9. Depth

Strictly flat. Overlap implies order. The only permitted "many" device is a
flat offset duplicate behind a shape. No shadows, no bevels, no glow.

### B10. Emphasis

The eye is led by three devices only: a warm hue among cool ones, a slightly
heavier outline on one shape, and size. Never by a saturated red, a glow, a
starburst, or a bold sentence.

---

## C. Prompt text that encodes this

These are the opening and closing blocks now used by `style-contract.md`
sections A and F. They are repeated here so the reasoning in section B sits
next to the text it produced.

### C1. Opening block

```text
A flat two-dimensional scientific schematic in the visual style of a Nature
Reviews Genetics figure: an illustrated, colourful, precisely drawn editorial
diagram, not a slide, not a dashboard, not a flowchart of boxes. <orientation>,
<aspect ratio>, white background, no outer frame, generous even whitespace.

Shapes are rounded rectangles filled with flat muted colour and no outline.
Related shapes may sit on a very pale tinted rounded field of their own hue.
Four to six coordinated muted hues, each standing for one kind of thing and
used consistently across every panel. No gradients, no shading, no shadows, no
3D, no glow, no texture.

Simplified flat biological glyphs are used where they stand for a thing:
<list the two or three glyphs this figure uses, from the glyph table>. Each is
a single flat fill in the hue of its class. No twisted DNA helices, no
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

### C2. Closing block

```text
Overall impression: a warm, colourful, calm and precise illustration from a
Nature Reviews journal. Flat colour carries the shapes, whitespace carries the
grouping, glyphs carry the biology, and the caption carries the explanation. The
figure is sparse enough to read in ten seconds. Emphasis falls on <the one
place> through <a warmer hue among cool neighbours / a heavier outline / size>,
and nowhere else. No title text outside the panel letters and their panel
headings.
```

---

## D. Translating pipeline content into figure content

Briefs in this repository start from scripts, so the raw material is code and
counts. This table is how that material becomes NRG content.

| Raw material from the script                | On the canvas                                      | In `## Caption`                     |
| ------------------------------------------- | -------------------------------------------------- | ----------------------------------- |
| `minfi::read.metharray.exp(...)`            | IDAT page glyph -> array glyph, label "Read IDATs" | the function call and its arguments |
| `RGChannelSet 1,051,539 addresses x 24`     | array glyph, label "Raw intensities", `n = 24`     | the exact probe count               |
| `detectionP > 0.01, 5 percent rule`         | funnel glyph, label "Detection P-value filter"     | the thresholds                      |
| `0 samples, 3,052 probes flagged`           | small grey annotation "~3 k probes removed"        | the exact numbers per array version |
| `preprocessQuantile` -> `GenomicRatioSet`   | matrix glyph, label "Normalised beta values"       | the method name                     |
| hg38 manifest join, 99.995 percent          | chromosome glyph, label "Mapped to hg38"           | the exact fraction                  |
| `getSnpBeta`, `getQC`, `getSex`             | three small glyphs on a side branch, one word each | the calls and the counts            |
| LSF job id, runtime, R and package versions | never                                              | always                              |

A useful test before writing a prompt: read each label aloud to a collaborator
who has not seen the script. If the label needs the script to make sense, it is
caption material.

---

## E. Inspecting a render for NRG style

Check, after the factual checks in `SKILL.md` section 6:

1. **Outlines.** Filled shapes have no stroke. If every box has a grey border,
   the render is a flowchart.
2. **Hue count.** Four to six, each meaning one thing, consistent across
   panels. A grey figure with one red is the old contract, not NRG.
3. **Case.** No small caps, no all-caps labels, no letter-spacing, no
   monospace.
4. **Words.** Roughly 40 per panel. Count them.
5. **Glyphs.** Present where a kind of thing needs naming, flat, in-palette,
   and none of the forbidden forms.
6. **Arrows.** Thin, grey, small identical heads.
7. **Depth.** Nothing cast a shadow.

If a render fails 1 to 3, the prompt still contains dashboard vocabulary. Find
and remove the words "card", "role bar", "small-caps", "monospaced", "strip",
"chip", "pill" from the prompt before re-rendering.
