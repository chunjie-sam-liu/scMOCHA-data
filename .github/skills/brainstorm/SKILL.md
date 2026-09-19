---
name: brainstorm
description: "Run and record a pre-plan discussion in a BRAINSTORM.md or {date}-{short-title}.BRAINSTORM.md. Use when the shape of the work is not yet clear and a PLAN.md would be premature: scoping a new analysis, choosing between tools or methods, evaluating software nobody on the project has used before, mapping what blocks what and who owns each blocker, sizing an unfamiliar dataset, comparing parameter settings, or turning a broad charter into an ordered work queue. Aimed at bioinformatics and statistical genetics - GWAS, QTL, eQTL, meQTL, EWAS, RNA-seq, single-cell, methylation, multi-omics - where the fork taken changes the analysis meaning rather than just the code. Covers the verified-state rule, the multi-round question lifecycle from OPEN to RESOLVED, the Mermaid diagram library, parameter tables that record why a value was chosen, the blocker and trap inventories, and the handoff into PLAN.md and DECISION.md."
---

# Brainstorm

A brainstorm is a **thinking document**: the shared, durable record of a
discussion held before anyone knows enough to write a plan. It exists because
the person asking often does not have the whole picture, and because the
technique or the software is new to the project. Its job is to build that
picture out of verified facts, name every fork that changes the analysis
meaning, and settle those forks over as many rounds as it takes.

It survives the chat session. That is the entire point. A decision argued for
an hour in chat and never written down gets re-argued next week.

A brainstorm **authorizes nothing**. It creates no scripts, submits no jobs,
and invalidates no outputs, so it does not trigger the Plan-First Gate. It is
what runs _before_ the gate, so that the plan written afterwards has no open
guesses in it.

This skill is portable and names no project value. The track layout, the stage
numbering, the env variables, the output roots, and the note directory come
from the repository's `.github/instructions/` bindings or its `AGENTS.md`.
Those win where they conflict.

Related skills: `analysis-pipeline` (the `PLAN.md` / `PROGRESS.md` /
`DECISION.md` set this hands off to), `markdown-doc` (choosing between Mermaid,
table, checklist, code block, and prose for each section), `statistical-genetics`
(the validity questions a brainstorm must raise), `data-verification` (how every
state claim is verified), `image-prompt` (when a diagram graduates into a real
figure), `lsf-resource-planner` and `long-running-jobs` (when a fork is about
runtime).

---

## 1. When to brainstorm

Brainstorm when **any** of these is true:

- The work has more than one plausible shape and nobody can yet say which.
- A tool, method, or file format is new to the project and its capabilities
  have to be established before a plan can assume them.
- The request is a charter, a paper section, or a paragraph of intent, and it
  has to become an ordered queue of stages.
- Progress is blocked on something outside the repository -- a collaborator's
  file, an unconfirmed genome build, an unapproved cohort definition -- and it
  is not yet clear what can proceed anyway.
- Two finished pipelines disagree and which one becomes primary is undecided.
- A number that every downstream decision depends on has never been computed.

Do **not** brainstorm when the shape is already clear. Writing a brainstorm for
a settled question is the slow way to write a plan.

| Document        | Answers                                     | Written                     |
| --------------- | ------------------------------------------- | --------------------------- |
| `BRAINSTORM.md` | What is the shape of this, and which forks? | before a plan exists        |
| `PLAN.md`       | What will be built, step by step?           | after the forks are settled |
| `DECISION.md`   | What was chosen, over what, and why?        | alongside the plan          |
| `PROGRESS.md`   | What has actually run?                      | before the first run        |
| `AGENTS.md`     | How do I work here without breaking it?     | continuously                |
| `DIAGRAM.md`    | What does the published figure show?        | when a figure is needed     |

A brainstorm is not a `DECISION.md`. It holds the argument, the options, the
evidence, and the reversals. `DECISION.md` holds one line per settled choice.

---

## 2. File convention

Scope decides the location.

| Scope                                          | Path                                              |
| ---------------------------------------------- | ------------------------------------------------- |
| One stage that already exists                  | `src/NN-stage/{date}-{short-title}.BRAINSTORM.md` |
| Scoping a stage or track that has no plan yet  | `src/NN-stage/BRAINSTORM.md`                      |
| Spans stages, spans tracks, or is project-wide | `{note-dir}/{date}-{short-title}.BRAINSTORM.md`   |

- `{date}` is `YYYY-MM-DD`, the day the brainstorm is opened, and never changes
  when later rounds are added.
- `{short-title}` is lowercase kebab-case, 2-5 words, naming the **question**
  (`meqtl-roadmap`, `which-normalization`), not the stage.
- The uppercase `.BRAINSTORM.md` suffix is what makes it greppable and what
  separates it from a scratch note.
- `{note-dir}` is whatever the repository bindings name as the home for
  cross-cutting documents. Ask if the bindings name none; never invent a
  directory.
- A cross-cutting brainstorm is **always** dated. Only a stage may hold a bare
  `BRAINSTORM.md`, and only one.
- Never park a brainstorm in a directory that a test runner, a build, or a
  pipeline globs. A thinking document must not become an input.

When a stage-scoped brainstorm turns out to be cross-cutting, move it and leave
nothing behind. Two copies diverge within a day.

---

## 3. Required structure

Write these in this order. Each one exists because leaving it out is how a
brainstorm turns into a confident, wrong plan.

| Section                                       | Content                                                                                                         |
| --------------------------------------------- | --------------------------------------------------------------------------------------------------------------- |
| Title + status header                         | Date opened, `Status`, last updated, and one line naming the question.                                          |
| `## What this file is and is not`             | That it authorizes nothing, triggers no gate, and which charter or plan it feeds.                               |
| `## 0. TL;DR`                                 | The blockers and the one thing to start now. Readable in thirty seconds. Written last, placed first.            |
| `## 1. Where we are, verified {date}`         | Three tables: done, exists but empty, not started. **Every row carries a source path and a read date.**         |
| `## 2. The whole picture`                     | One Mermaid diagram of the landscape: done, blocked, to build, and what depends on what. See section 6.         |
| `## 3. Blockers, in detail`                   | One subsection per blocker: what it blocks, who resolves it, what can proceed meanwhile.                        |
| `## 4. Questions to settle before <x>`        | One `Q` block per fork that changes the analysis meaning. See section 5.                                        |
| `## 5. What each stage should do`             | A table of stage -> question it answers -> key output -> blocked by. Then prose for the hard ones only.         |
| `## 6. Technical detail`                      | Parameters, tool capabilities, formats. Only where the value is contested or the tool is unfamiliar. Section 7. |
| `## 7. Scale reality check`                   | The sizes, counts, and runtimes the plan will be sized against, each marked measured or estimated.              |
| `## 8. Traps already paid for, do not re-buy` | Numbered pitfalls carried forward from `AGENTS.md`, finished stages, and failed runs.                           |
| `## 9. Proposed order of work`                | Ordered by what unblocks the most, not by stage number. Grouped by what each group waits on.                    |
| `## 10. What I need from you`                 | A decision table: number, decision, and **what happens by default if the user says nothing**.                   |
| `## 11. Discussion log`                       | Append-only, one row per round: date, what changed, which `Q` moved. See section 5.3.                           |

The headings above are the exact strings the template writes, numbering
included; keep them verbatim so a section can be found by grep. Sections may be
dropped only when they are genuinely empty, and then the heading stays with the
word `none` under it. A missing `## 3. Blockers, in detail` reads as "nothing
blocks this"; an absent one reads as "nobody checked".

The full copyable skeleton is in
[references/brainstorm-template.md](./references/brainstorm-template.md).

---

## 4. Hard rules

1. **Every state claim carries a source and a date.** A row saying
   "718,863 probes x 7,680 arrays" names the file it was read from and the day
   it was read. A claim with no source is a guess wearing a table cell. Use the
   cheapest check that answers the question -> `data-verification`.
2. **Never state a tool's capability from memory.** Run its `--help`, read its
   version string, or cite its documentation with the version. Builds differ:
   a mode described in the upstream docs is frequently absent from the
   installed binary, and a plan built on an absent mode fails late and
   expensively. Record what the tool **cannot** do as carefully as what it can.
3. **Never quote a parameter value without its justification and its source.**
   See section 7.
4. **Distinguish measured from estimated, every time.** A measured number names
   its source. An estimated number says so in the same cell and says what will
   replace it.
5. **Name the fork, do not resolve it silently.** Anything that changes cohort
   inclusion, covariate set, transformation, threshold, model family, reference
   panel, or which tool is primary is a `Q` block with options, not a sentence
   in prose.
6. **Never delete a rejected option or a superseded recommendation.** The
   alternative a choice beat is exactly what `DECISION.md` needs later. Amend
   by appending; strike nothing out.
7. **The brainstorm never writes code and never runs anything that writes.**
   Read-only inspection is expected and encouraged; edits are not.
8. **Mermaid blocks take their color from the documentation palette**, never
   from an invented hex. Section 6.
9. **ASCII only.** No emoji, no Greek, no LaTeX. `>=`, `<=`, `+/-` in prose and
   in labels.
10. **Update it in the round it changes.** A brainstorm that lags the
    discussion by two rounds is misinformation, not history.

---

## 5. Questions: the multi-round lifecycle

A `Q` block is the unit of a brainstorm. It opens undecided and closes only
after discussion, usually over several rounds, sometimes after new evidence is
fetched between rounds.

### 5.1 The block

```markdown
### Q3 — One array per person, or all of them? [OPEN]

Genotype is per person; methylation is per array, and 812 persons have two or
more surviving draws. Neither engine takes repeated measures on one genotype.

| Option                                    | Consequence                                                  |
| ----------------------------------------- | ------------------------------------------------------------ |
| (a) One array per person, by a fixed rule | Simple, matches both tools, loses ~1,200 arrays              |
| (b) All arrays, person as a random effect | Correct, but neither tool supports it without a custom refit |
| (c) All arrays, ignore the duplication    | **Wrong.** Inflates every test statistic. Not an option      |

**Recommended: (a)** for the main atlas, with the selection rule written down
once and never re-picked. Keep the extra draws as a longitudinal sensitivity
set.

**Evidence:** `src/02-demo/PROGRESS.md`, read 2026-09-15.
```

Rules for the block:

- One sentence of context, then **why it changes the analysis meaning**. A fork
  that only changes code style is not a `Q`.
- Options as a table: the option, and its consequence or cost. Two is the
  minimum; an option listed only to be dismissed is still worth a row, because
  it stops the same idea returning next round.
- Exactly one **Recommended**, with the reason. Never hedge across two.
- **Evidence** naming the files the recommendation rests on, with read dates.
- A `Q` that blocks nothing says so, and says which stage plan will inherit it.

### 5.2 Status, and how a `Q` closes

The status tag sits in the heading and is the only place status is recorded.

| Tag                       | Meaning                                                                    |
| ------------------------- | -------------------------------------------------------------------------- |
| `[OPEN]`                  | Posed, not yet discussed, or still contested                               |
| `[RESOLVED YYYY-MM-DD]`   | Settled. The block gains a `**Resolution:**` paragraph                     |
| `[DEFERRED -> <file> Qn]` | Cannot be settled here; carried into a plan with the recommendation intact |
| `[DROPPED YYYY-MM-DD]`    | No longer a question, with the one-line reason it stopped being one        |

A resolution paragraph states the option chosen, who chose it, the date, and --
where the recommendation was overridden -- the reason the user gave. It is
appended below the recommendation, never written over it.

```markdown
**Resolution 2026-09-16:** (a), selection rule = earliest passing draw.
Overrode the "highest median intensity" variant because draw order is
reproducible from `demo_spine.parquet` and intensity depends on the QC pass
that produced it.
```

A `Q` resolved without the user is marked `provisional` in its resolution line,
exactly as in `DECISION.md`, and is re-raised in `## What I need from you`.

### 5.3 Rounds

Discussion is iterative and the file must show it. After each round, append one
row:

```markdown
## Discussion log

| Round | Date       | What changed                                                               |
| ----- | ---------- | -------------------------------------------------------------------------- |
| 1     | 2026-09-15 | Opened. Q1-Q8 posed, three blockers named                                  |
| 2     | 2026-09-16 | Q1 and Q3 resolved; Q5 recommendation reversed after reading the tool docs |
| 3     | 2026-09-18 | Q4 deferred to `src/03.01-prep-geno-pca/PLAN.md` Q1                        |
```

Between rounds, go and get the evidence rather than arguing from memory: read
the file, run the tool's `--help`, compute the one number everything depends
on. A round that changes no evidence usually changes no minds.

The brainstorm reaches `SETTLED` when every `Q` is `RESOLVED`, `DEFERRED`, or
`DROPPED`. Then, and only then, write the plan.

---

## 6. Mermaid

Diagrams are the reason a brainstorm is readable. Use several, each answering a
different question, rather than one diagram carrying everything.

The direction rule, the documentation palette, and one copyable template per
question are owned by `markdown-doc`
[references/mermaid-patterns.md](../markdown-doc/references/mermaid-patterns.md).
Read it before drawing. Three rules a brainstorm adds on top:

1. **A number on a node came from a file**, never from memory -- rule 4.1
   applies inside a diagram exactly as it does in the prose.
2. **Status is carried by grouping and edge style as well as by fill**: solid
   for required, dotted for optional or sensitivity, `-.blocks.->` with a label
   for a blocked dependency. A brainstorm that needs a real scientific palette
   has outgrown itself and belongs in a `DIAGRAM.md` -> `image-prompt`.
3. **Several small diagrams, not one big one.** A brainstorm carrying every
   question in a single graph is unreadable; one diagram per question in the
   library table is the shape.

---

## 7. Technical depth

The brainstorm is where an unfamiliar tool or parameter is made legible. Go
deep only where the value is contested or the behaviour is surprising; a
parameter with one obvious setting needs no table.

### 7.1 Parameters

Never write a bare value. Every contested parameter gets a row:

| Parameter  | Proposed | Why this value                                          | What changes if it moves                          | Source                        |
| ---------- | -------- | ------------------------------------------------------- | ------------------------------------------------- | ----------------------------- |
| cis window | +/- 1 Mb | Matches GoDMC and the published atlases we compare to   | Narrower cuts the test burden; wider finds little | `Analysis.PLAN.md` s5.5       |
| MAF        | >= 0.01  | Below this the per-variant power is not worth the tests | Lower inflates the variant count and the runtime  | cohort n, computed 2026-09-15 |

- **Why this value** is the column that justifies the setting. "The default" is
  an acceptable answer only when it is stated as such, with the tool and
  version whose default it is.
- **What changes if it moves** is what makes the row useful next year.
- **Source** is a file path, a tool version, or a citation. Never "standard
  practice".
- A parameter that cannot be settled before seeing data says so, and names the
  cheap check that will settle it -- one chromosome, one trait, one array
  index.

### 7.2 Unfamiliar software

When the project has not used a tool before, the brainstorm establishes the
ground truth, in this order:

1. **What is installed** -- the exact version, from the binary, not the docs.
2. **What modes it has** -- from `--help` on that binary.
3. **What it does not do** -- the modes the plan assumes and the build lacks.
   This is the highest-value paragraph in the whole document.
4. **What its input contract is** -- orientation, coordinate base, index base,
   chromosome naming, which allele the effect refers to. These are the silent
   failure modes; every one of them can produce zero results with no error.
5. **What it costs to run** -- measured on a smoke run, not guessed
   -> `lsf-resource-planner`.
6. **What would validate it** -- simulated data with known truth, scored,
   before any production run depends on it.

Where two engines will be run, say plainly which is primary and which is the
cross-check, and what a disagreement between them would mean. Three engines
that produce the same evidence as two is three sets of file formats to keep in
sync.

Per-analysis question sets -- GWAS, QTL, EWAS, RNA-seq, single-cell,
multi-omics -- and the parameters that always need justification:
[references/domain-checklists.md](./references/domain-checklists.md).

---

## 8. Handoff

A settled brainstorm becomes a plan. It does not become the plan.

1. **The file stays where it is**, read-only, with
   `Status: SETTLED YYYY-MM-DD -> <plan file>` in its header. Never rename it
   `_stale-`, never fold it into the plan, never delete it. It is the only
   record of what the choices beat.
2. **Each `RESOLVED` `Q` becomes a `D` entry** in the target `DECISION.md`, in
   the form "chose X over Y because Z", with `Evidence:` pointing at the
   brainstorm file and the `Q` number. It does **not** reappear as a `Q` in
   `PLAN.md`; re-opening a settled fork is the failure this whole convention
   exists to prevent.
3. **Each `DEFERRED` `Q` becomes a `Q`** in the plan that first needs it, with
   the brainstorm's recommendation carried over verbatim as the recommended
   answer.
4. **`PLAN.md` links back** under a `## Brainstorm` heading, one line.
5. **`## Traps already paid for` is checked against `AGENTS.md`.** Any trap not
   already there is added, one line each, in the same change set.
6. **Blockers owned by a person do not migrate.** They stay in the brainstorm
   and are repeated in the plan's scope section as a dependency, so the plan
   does not silently assume a file that has not arrived.

One brainstorm often feeds several plans. That is normal and is why it is not
renamed into one of them.

---

## 9. Forbidden patterns

- A state claim, count, or size with no source path and read date.
- A tool capability, flag, or default quoted from memory instead of from the
  installed binary or its versioned docs.
- A number carried over from an earlier draft, a chat message, or another
  document without being re-read from the source.
- A fork resolved in prose instead of a `Q` block.
- Deleting a rejected option, or editing a recommendation in place after it was
  discussed.
- A hex in a brainstorm Mermaid diagram that is not one of the documentation
  palette classes.
- Committing a rendered PNG or SVG of a brainstorm diagram.
- Implementation, file creation, or job submission as a side effect of a
  brainstorm.
- A brainstorm with no Mermaid at all, or with one diagram carrying every
  question.
- A brainstorm in a directory that a test runner or pipeline globs.
- Emoji, anywhere.

---

## References

- [references/brainstorm-template.md](./references/brainstorm-template.md) --
  the full copyable skeleton, every required section with its shape.
- [markdown-doc/references/mermaid-patterns.md](../markdown-doc/references/mermaid-patterns.md)
  -- one copyable template per diagram role, the `TD` direction rule, and the
  documentation palette. Owned by `markdown-doc`, not by this skill.
- [references/domain-checklists.md](./references/domain-checklists.md) --
  the questions a brainstorm must raise per analysis type, and the parameters
  that always need a justification row.
