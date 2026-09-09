# Reference: figure archetypes

Five shapes cover essentially every figure a computational-biology project
needs. Pick one primary archetype. Do not fuse two on one canvas — a method
schematic bolted onto a pipeline flowchart reads as neither. Supporting panels
are allowed only where the selected skeleton permits them and only for claims
listed in the brief; they do not create a second archetype.

Each entry gives the canvas, the layout skeleton, what carries the message, what
the figure must never contain, and the failure mode that shows up in practice.
The skeletons below are the `paper` form. The last section trims each one to the
`discussion` budget of about 25 words and one number per panel, counted with
the rules in `SKILL.md`.

In every skeleton, an entity is a filled shape or a glyph with a short label,
and a process is a labelled arrow. Nothing is an outlined box, and no label
names a function, package, version, or file.

---

## 1. Pipeline flowchart

**Renders** a stage's data flow: what goes in, what it becomes, what comes out.
This is the default figure for a stage-level `DIAGRAM.md`.

**Canvas** landscape, 3:2. Single left-to-right chain when the stage is linear;
two parallel chains converging when it is not.

**Skeleton**

```
input glyphs (left)                      e.g. IDAT pages, one per array version
  -> labelled arrow "Read intensities"
    -> entity glyph "Raw intensities"    e.g. an array grid, n = 24
      -> labelled arrow "Detection P-value filter"
        -> labelled arrow "Quantile normalisation"
          -> entity glyph "Beta matrix"  e.g. a matrix with a header row
            -> labelled arrow "Map to hg38"
              -> entity glyph "hg38 positions"  e.g. a chromosome bar
side branch: one dashed arrow from an entity down to a small glyph, one label
```

Each array version, cohort, or branch owns one hue and keeps it across the
chain. Steps are the arrows, not boxes; each arrow carries one short verb
phrase.

**Carries the message**

- The chain of things the data becomes, in order, each named in plain words.
- The one or two counts that matter: what went in, what came out.
- Hue separates roles handled differently (array versions, discovery versus
  lookup, QC'd versus raw).

**Never contains** function names, package versions, script numbers, file
paths, parameter syntax, run times, or any step that produces no artifact.
Those go in `## Caption`, where the reader who needs them will look.

**Failure mode** the model merges adjacent steps when there are more than about
six in a row, and it turns steps into boxes when the prompt describes them as
nouns. Keep the chain to four or five entities, phrase every step as a verb
phrase on an arrow, and group a long chain into two or three labelled phases
on pale fields.

---

## 2. Study design

**Renders** cohorts and their roles, outcome families, and the model matrix.
Answers "what was analysed, in whom, with which covariates".

**Canvas** landscape, 3:2. Columns left to right: cohorts, outcomes, models.

**Skeleton**

```
column 1  cohort glyphs, one per cohort, each a row of sample circles or one
          filled shape in the cohort hue, with role and n beneath; the
          asymmetry between roles is the point
column 2  one common model grammar in ASCII plus compact outcome-specific term
          rows; if three distinct equations are essential, split them across
          panels so no panel contains more than two formulas
column 3  the model matrix: rows M0, S1..Sn; columns = outcome families;
          cells = filled / empty, with the primary row visually distinguished
```

**Carries the message**

- The primary model must be unmistakably separated from the sensitivity models,
  or a reader counts them as independent scans.
- If the equations differ in exactly one term, say so on the page and mark the
  term. That single difference is usually the reason the figure exists.
- A cell that is deliberately absent needs its reason on the page, not in a
  caption.

**Never contains** a cohort that is excluded from the design being shown, or a
covariate list abbreviated differently across the equations.

**Failure mode** the matrix loses rows or transposes. Keep it under about eight
rows and four columns; state the row and column counts in the prompt.

---

## 3. Method schematic

**Renders** how one calculation works, with the formula on the page and one
real unit carried through it.

**Canvas** landscape 16:9 as a standalone figure, or a full-width panel a above
smaller supporting result panels that each serve a listed claim. Those panels
remain part of the method-schematic skeleton rather than a second archetype.

**Skeleton**

```
one short question line across the top: the question the method answers
  stage 1 glyph   the input entity, with its one real value
  stage 2 arrow   the transform, with the formula beside it
  stage 3 glyph   the target entity, with its one real value
  stage 4 arrow   the output quantity, with the formula beside it
  optional inset  a small smooth curve showing the transform's behaviour, with
                  the worked example plotted as a single dot
one pale tinted band beneath: the worked example as stacked short rows, one
stage per row, ending in the actual outcome
```

**Carries the message**

- The question line at the top states the conditioning clause. A conditional
  quantity whose condition lives only in the caption will be misread as
  absolute.
- The worked example is what makes the figure reproducible. One real unit, real
  values, all the way to the answer, including when the answer is a failure.
- Mark which inputs are transported and which are carried over unchanged.

**Never contains** a second method, or a summary statistic that the worked
example does not touch.

**Failure mode** garbled formulas. Fix by shortening the formula, not by asking
for clearer text. If two rounds fail, drop to the Discussion rung and retain
only the one essential formula.

---

## 4. Result schematic

**Renders** the shape of a result: a distribution, an expected-against-observed
comparison, a calibration check.

**Canvas** landscape 3:2, two to four small panels in one row. Never more than
four.

**Skeleton**

```
panel  distribution   one smooth filled curve, a dashed threshold line with a
                      short label, the one number that matters beneath
panel  comparison     two bars on a shared baseline, the ratio between them
                      described in words, the departure carried by the warmer
                      hue
panel  calibration    a few dots against a dashed 1:1 diagonal, a solid
                      reference line at the null, one short annotation
```

**Carries the message**

- Bar height ratios and point coordinates must be described in words as well as
  numbers; the model reproduces the described geometry, not the numbers.
- Say where the punchline sits. Bottom-right is where the eye lands last.
- A caption line under the axis is the place to kill the most likely
  misreading ("expected = sum of individual probabilities, not a count").

**Never contains** real per-observation data claiming to be a plot. This
archetype is the one that most often gets mistaken for a real figure, so its
`## Honest scope` section is mandatory and must name the real output file.

**Failure mode** the model "fixes" the result by making the points climb the
diagonal. State explicitly that the points lie on the null line and nowhere near
the diagonal.

---

## 5. Conceptual overview

**Renders** the idea, for a talk opener or a grant figure. No file paths, no
step numbers, no counts.

**Canvas** landscape 16:9.

**Skeleton**

```
three to five large labelled regions, left to right, joined by broad arrows,
each region holding one short phrase and at most one flat glyph from the
glyph table that names what kind of thing the region is about
```

**Carries the message** one sentence. If it needs two, it is the wrong
archetype.

**Never contains** numbers, software names, file formats, or anything that will
be out of date next month.

**Failure mode** decoration creeps back in - twisted helices, cells with
organelles, microscopes, scientists. Repeat the forbidden list and name the
specific forms to exclude.

---

## Choosing between them

| The reader should leave knowing        | Archetype           |
| -------------------------------------- | ------------------- |
| What this stage does to the data       | Pipeline flowchart  |
| Who was analysed and under which model | Study design        |
| How this number was computed           | Method schematic    |
| What the answer looks like             | Result schematic    |
| Why this project exists                | Conceptual overview |

If two rows describe separate primary reader goals, the figure is doing two
jobs. Split it into two briefs. A supporting panel explicitly permitted by the
selected skeleton does not create a second primary goal.

---

## Trimming to the `discussion` budget

About 25 words and one number per panel, one formula, one or two panels. Cut in
this order until the count fits. Stop as soon as it does; do not keep cutting.

| Cut                                                           | Why it goes first                             |
| ------------------------------------------------------------- | --------------------------------------------- |
| Every result panel except the one that is the punchline       | The method or the result, not both            |
| The inset plot                                                | It restates what the main shape already shows |
| The annotation layer: leader lines, callouts, QC statistics   | These answer a reviewer, not a colleague      |
| Secondary formulas, keeping the one that defines the quantity | Two formulas is a `paper` figure              |
| Legend rows beyond four                                       | Merge or label in place                       |
| Every number except the one the figure is about               | The rest live in the caption                  |

What never gets cut from the Discussion rung: the question the figure answers,
the one formula or the one comparison that is the point, and the real numbers
attached to them. The lower One-message rung may compress a method figure to
one fully defined symbolic formula or one fully numeric relation plus its
outcome; it never keeps supporting panels.

Per archetype:

- **Pipeline flowchart.** Collapse the chain to three or four entities with
  one labelled arrow between each. Keep the input count and the output count;
  drop every intermediate number.
- **Study design.** Keep the cohort row and the model matrix. Use one common
  model grammar plus compact outcome-specific term rows, or keep exactly one
  symbolic equation and say in words how the others differ. Define every
  retained symbol.
- **Method schematic.** Keep the question line and a compressed chain with
  one worked example integrated into its stages. Merge stages or remove
  duplicated labels until the word count fits; do not repeat the same values
  in a separate example band. Drop the inset and every supporting result
  panel. Drop the symbol annotation block only when the example is fully
  numeric and leaves no undefined symbols; otherwise retain definitions for
  every symbolic variable.
- **Result schematic.** Keep one panel. Almost always the expected-against-
  observed comparison, because a ratio between two bars survives any shrink.
- **Conceptual overview.** Already at budget. If it is not, it was never a
  conceptual overview.

A `paper` brief that will also be discussed may add the Discussion tier as an
optional standalone prompt, rendered from the same verified inventory. It does
not replace the required Full `## Prompt` and complete Leaner variant.
