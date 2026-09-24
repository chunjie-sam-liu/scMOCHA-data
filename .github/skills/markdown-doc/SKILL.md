---
name: markdown-doc
description: "Choose the representation for any Markdown file, after deciding its class: an operational document another session resumes (PLAN.md, PROGRESS.md, DECISION.md, BRAINSTORM.md, AGENTS.md, DATA.md, a paired script .md, a debugging note) or a reference document that travels without the repository (a README shipped with a data release, a data dictionary, a USAGE file). Use when creating or editing either; when a document has drifted into walls of prose or repeats the same attributes per item; when a Mermaid diagram renders too small to read or needs meaningful color; when choosing between a diagram, a table, a task list, a code block, and prose; or when a Markdown doc must be rendered to one self-contained HTML file with Mermaid diagrams rendered. Covers the class split, the release-README skeleton, the representation table, the TD direction rule, the documentation palette, the resume-first spine, the measured-vs-projected rule, which skill owns each template, and the bundled pandoc+Mermaid HTML-export script."
---

# Markdown documents

A Markdown file in this repository is read under time pressure, usually by somebody who did not write it. Two very different readers exist, and the first decision is which one the file serves (section 0):

- An **operational** reader -- a colleague picking up a stalled run, or an agent session with no memory of yesterday -- needs the current state fast. Optimize for **fast comprehension, operational usefulness, and resumption**, in that order; anything that makes this reader reconstruct state from narrative is a defect.
- A **reference** reader -- somebody who received a data release, a data dictionary, or a usage guide, often without the repository -- cannot ask and cannot follow a link. Optimize for **completeness, self-containment, and correctness against the shipped files**; anything that sends this reader to a file they do not have is a defect (section 11).

This skill owns **how a document is represented**: which of Mermaid, table, task list, code block, or prose carries each piece of information, and what shape an operational or a reference document takes. It does not own the required sections of any specific document type; section 10 says who does.

---

## 0. Decide the document class first

| Class | Examples | Reader | Shape |
| --- | --- | --- | --- |
| Operational | `PLAN.md`, `PROGRESS.md`, `DECISION.md`, `BRAINSTORM.md`, `AGENTS.md`, a paired script `.md`, a debugging or migration note | the next session, inside the repository | current state first; sections 5, 6 and 8 |
| Reference | a `README.md` shipped inside or beside a data release, a data dictionary, a `*.USAGE.md` handed to another team, a methods appendix | somebody with the files but not the repository | complete and self-contained; section 11 |

The deciding question is **where the reader will be**. A README that describes a release is a reference document even though producing the release was an operational event: the progress file records the run, the README describes the product. A file that tries to be both fails the reference reader first -- it opens with a status block, describes the release as a diff against the previous one, stamps numbers `measured`, and links to plans the recipient cannot open. Keep the run in the progress file and the product in the README.

`DATA.md` sits between the two: it is read inside the repository, so links are fine, but every row must stand on its own. `data-catalog` owns it.

---

## 1. Choose the representation before writing

**Do not default to prose.** Prose is one of five choices, and it is the right one only for the things the other four cannot carry. Decide per section, from the structure of the information rather than from habit.

| The information is | Use |
| --- | --- |
| Relational or directional: workflows, pipelines, job dependencies, DAGs, execution order, architecture, data flow, state machines, sequence interactions, branching or parallel work | Mermaid |
| Repeated attributes a reader will scan or compare: job status, inputs and outputs, parameters, dimensions and counts, expected vs observed, dependencies and their reasons, validation results, runtimes, config differences, one row per debugging attempt | Table |
| Actionable progress: what is done and what is not | Task list |
| Meant to be copied or executed: commands, code, configuration, logs, queries, file excerpts | Code block |
| Sequential or categorical, with nothing to compare across items | Simple list |
| Interpretation, rationale, caveats, conclusions, operational warnings | Concise prose |

The test for a table is **repetition**: the moment the same three or four attributes are stated for a third item, the paragraphs are a table that has not been written yet. The test for a diagram is **tracing**: if the reader has to hold "A feeds B, which also needs C" in their head, draw it.

---

## 2. Mermaid

Reach for a diagram when it communicates dependency or flow better than a sentence can -- not to decorate a page. Every block obeys the same rules: `flowchart TD` by default, under about fifteen nodes, ASCII labels broken with `<br/>` and, in an operational document, carrying live state when it is known (stage name, job id, array range, state, runtime, output dimensions) -- a reference document's labels carry counts and shapes only -- one sentence above saying what to take away, prose below that says what the diagram cannot rather than repeating it, and every number on a node read from a file rather than recalled.

**Color means something and is never invented.** One fixed documentation palette supplies it: `done` / `run` / `pend` / `todo` / `fail` / `dead` for a submitted pipeline, `input` / `step` / `output` / `fork` / `blocked` / `ext` for a structural diagram. One set per diagram, `classDef` copied verbatim, only the classes used declared. That palette is **documentation chrome, not a figure palette**: it encodes no cohort, trait, ancestry, or outcome, so it is never re-derived and never added to a track color file. A fill never carries information the label does not also carry, and any other hex is a defect. -> `palette`

One documented exception to all of the above: the `## Structure preview` block in a `DIAGRAM.md` stands in for a printed landscape figure, so it follows its archetype's direction and its own theme block. -> `image-prompt`

### Where the rules and templates live

The paragraphs above are the summary. The direction rule with its escape hatches, the edge vocabulary, the `classDef` blocks, and one copyable template per question live in the libraries -- read one before drawing anything non-trivial.

| Need | Library |
| --- | --- |
| Direction, the full palette, and one template per question | [references/mermaid-patterns.md](./references/mermaid-patterns.md) |
| Submitted pipeline, its live state, and the classes it uses | `analysis-pipeline` [plan-and-progress.md](../analysis-pipeline/references/plan-and-progress.md) section D.1 |

---

## 3. Tables

Use tables aggressively where they cut reading time, and never merely to decorate. A good column answers a concrete question the reader actually has:

- What is it?
- What state is it in?
- What does it depend on, and why does that dependency exist?
- What was measured, and what was expected?
- What action does this row require?

```markdown
| Step | Job ID | State | Depends on | Reason |
| --- | --- | --- | --- | --- |
| matrix | 323076536 | RUN | `numdone(323076535,*)` | reads the exported GDS |
| release | 323076538 | PEND | `numdone(323076536,*)` | packages every matrix stem |
```

- **Keep cells short.** A cell carrying two clauses is interpretation; move it to a sentence below the table.
- **One grain per table.** Mixing per-job rows with per-file rows forces the reader to re-read the header for every row.
- **A table of one row is a sentence.** A table of two columns where the second is always the same is a list.
- Put the long interpretation below the table, not inside it.

---

## 4. Task lists

Actionable progress is a checklist, never a paragraph. Each item carries the evidence that closed it, so that "done" is checkable.

```markdown
- [x] 01-export -- done 2026-09-16, 2 of 2 tasks exit 0, 7,680 columns out
- [ ] 02-matrix -- submitted, job 323076536 [1-10]
- [ ] 03-release -- waiting on 02
```

An unchecked box with no note is indistinguishable from work nobody has looked at. Say whether it is queued, blocked, or not started.

---

## 5. Measured state over narrative history

An operational document's job is to say what is true **now**. History stays only where it explains the present. A reference document handles numbers differently -- it names the shipped file each count can be re-derived from, and never carries a projected number (section 11).

- **Distinguish measured from projected, every time.** A measured number names where it was read from; a projected one says so in the same cell and says what will replace it.
- **Never silently replace a measured value with an earlier estimate**, and never carry a number forward from a chat message or an older draft without re-reading the source.
- Keep these pairs visibly apart: completed vs running vs queued; expected vs observed; current failure vs historical failure; required dependency vs incidental execution order.
- Prefer rewriting a section to appending a correction below a stale one. The exception is a rejected option or a superseded decision, which is appended to and never overwritten -- that is what the decision log needs later.

---

## 6. Structure an operational document for resumption

This applies to progress, implementation, migration, debugging, and execution documents, including the progress file of a release campaign -- anything a second person or session has to pick up. The reader must reach the current state without reading the history. It does **not** apply to the README that ships with a release, which is a reference document (section 11): a resume block, job ids, and links to the plan are exactly what that reader cannot use.

Near the top it answers, in order: what is complete, what is running now, what is waiting and on what, what is still required from the operator, what happens next, how the state is checked, and what to do when something fails. The spine that answers them:

| Element | Include when |
| --- | --- |
| Title and current scope | always |
| Links to the plan, decision record, issue, or spec | always |
| Short resume / current-state summary | always, near the top |
| Mermaid graph of the active workflow or dependencies | relationships matter, or any step has been submitted |
| Compact table of jobs, stages, artifacts, or status | more than two of anything |
| Task list of completed and remaining work | always |
| Verified measured results, marked as measured | any number was produced by a run |
| Commands to inspect, verify, resume, or operate the work | always |
| Caveats, failure behavior, operational warnings | an operator could do the wrong thing without them |
| Error history | it explains the current state or prevents re-debugging |
| Key paths and artifacts needed to continue | always |

Sections stay in a fixed order so a reader can scan for one; organize by section and never by chronology.

---

## 7. Emphasis

Bold is for facts the reader must not miss, and it stops working the moment it is used for emphasis in general. Reserve it for:

- what is currently running, and what is blocked;
- what must not be re-run, and what has already been submitted;
- destructive or irreversible actions;
- a surprising measured result;
- failure behavior that changes what the operator should do.

Two or three bold spans per screen is the working budget. A paragraph with four is a paragraph with none.

---

## 8. Keep failure history, compactly

Do not delete a failure that explains an environmental quirk, a workaround, or why a design is the way it is -- it gets paid for twice otherwise. But summarize it; a document is not a log file.

When an investigation ran several experiments, the experiments are a table and the conclusion is prose:

```markdown
| Attempt | Change                       | Result                            |
| ------- | ---------------------------- | --------------------------------- |
| 1       | baseline                     | exit 2 in 12 s, no program output |
| 2       | forced the job shell to bash | same failure                      |
| 3       | stripped the completion loop | exit 0, program ran               |
```

Then one short paragraph: what the cause turned out to be, what the fix was, and what it means for the next run.

**Never dump raw logs** unless the raw text itself is the evidence. When it is, quote the few lines that matter inside a code block, not the file.

---

## 9. Avoid mechanical formatting

The objective is not to maximize visual elements. A document where every section is a table is as hard to read as one that is all prose. Take the simplest representation in the section 1 table that preserves the important information.

Forbidden patterns:

- A wall of prose carrying repeated attributes that a table would show at a glance.
- A table with one row, or with a column that never varies.
- A diagram large enough that the reader traces crossing edges, or `LR` and wide enough that the renderer shrinks its labels to unreadable.
- A hex outside the documentation palette, or a fill carrying a meaning the node label does not also carry.
- A diagram or table with no sentence saying what to take away.
- Bold used for general emphasis rather than for a must-not-miss fact.
- A progress document whose current state can only be found by reading its history in order.
- A number with no indication of whether it was measured or projected.
- A raw log pasted in full where three lines were the evidence.
- A reference document that links to, or only makes sense with, a file the recipient will not have.
- A release README written as a status report or as a diff against the previous release, or one that silently drops a section its predecessor carried.
- A column dictionary copied from a plan or an earlier README instead of read from the shipped header.
- A number in a reference document that no measurement script printed, or a whole-matrix claim (no missing values, a value range) resting on a spot check or a log line.
- Emoji, anywhere.

### Source line breaks belong to the formatter

Markdown collapses a single newline inside a paragraph into a space, so where the source breaks changes nothing for the reader and everything for the next person editing the file. The repository's prettier config owns it through `proseWrap`: under `"preserve"` whatever breaks you type are frozen and never normalized, under `"never"` they are collapsed the first time the file is formatted -- and that collapse silently invalidates every multi-line anchor anyone had taken from it.

**Never hand-wrap.** Read `.prettierrc` before writing a `.md` and match it: one line per paragraph and per list item under `"never"`, the wrapping already in the file under `"preserve"`. Soft-wrap in the editor is the display fix; a hard line break is not. A file that predates the config is not reflowed as a side effect of an unrelated edit -- that is its own change, run and reviewed on its own. -> `pixi-env` section 5

---

## 10. Who owns which document

This skill governs representation. The required sections, naming, and tier rules for each document type belong to the skill or rule that owns it -- follow that one for structure, and this one for how each part is rendered.

| Document | Owner |
| --- | --- |
| `PLAN.md`, `PROGRESS.md`, `DECISION.md` and their dated tier-2 sets | `analysis-pipeline`, section 1.1 and `references/plan-and-progress.md` |
| `BRAINSTORM.md` | `brainstorm`, section 3 |
| `DIAGRAM.md` | `image-prompt` |
| `AGENTS.md` | `copilot-instructions.md`, "Track Guide" |
| `DATA.md` | `data-catalog` |
| Paired `config.sh.md` / `config.R.md` / `NN-step.md` | `analysis-pipeline`, section 3 |
| Release `README.md`, data dictionary, `*.USAGE.md` handed on | this skill, section 11 and `references/reference-docs.md` |

Where an owner prescribes a section, that section still obeys the rules above: its job table is a table, its pipeline graph is Mermaid, its verify block is a code block.

---

## 11. Reference documents that travel with data

A reference document is read by somebody who has the files and nothing else: no repository, no plan, no chat history, nobody to ask. Sections 1-4 and 7 still apply -- tables for repeated attributes, Mermaid for flow, code blocks for recipes -- but the document is judged by four different questions. Where each kind of content comes from, the measurement script, the section skeleton, the table and diagram sets, the recipe contract, the surprises to hunt for, and the pre-ship checklist are in [references/reference-docs.md](./references/reference-docs.md); read it before writing one.

**Is it self-contained?**

- No link to a file outside the shipped bundle, and no sentence whose meaning depends on one. Inline what the reader needs instead: the definition, the rule, the count, the reason, the caveat.
- Every path is a code span, never a link, even inside the bundle: the same text has to work in the repository copy and in the bundle copy. A companion file outside the bundle is named by a path relative to a stated root and by whoever owns it.
- The repository copy and the bundle copy are one document, kept identical.

**Is it complete?** The reader cannot ask, so length is not a defect here. A release README carries, per version or unit and with totals: the composition of what went in, every gate and funnel with its counts, the removal reasons exactly as they appear in the files, the feature filter with every count, the column dictionary of every shipped table, the defects the data itself carries, and recipes that were run against the shipped files. When a new release replaces an old one, take the section structure and the history from the previous release's README -- never its numbers -- update every section, and give a reason in the history section for any section that goes.

**Is it organized by the product, not by the run?** Start from the previous release's README when there is one, keeping its section numbers, and otherwise from the fixed skeleton of `reference-docs.md` section C: identity, a quick start carrying the surprises that silently corrupt an analysis, contents, how it was made, inputs, QC, what was removed and why, tiers, the feature filter, parameters, the files, recipes, the column dictionary, caveats, related files, release history, provenance. History goes near the end and never becomes the frame: the document describes what the release **is**, not how it differs from the last one.

**Is it grounded in the shipped files, after the build?**

- Finalize it only once the release exists. Every count, shape, size, range, and column name comes from a read-only measurement script run against the shipped files, and its labelled output is kept until review -- never from a plan, a proposal, a progress file, an earlier README, or memory, and never from a log line when the shipped files can answer. A plan names columns the build may not emit.
- Rules, thresholds, and definitions come from the producing code at the recorded commit; a reason comes from the decision record and is restated inline in one sentence. Reading the plan and the decisions is right; copying their numbers is not.
- A property of a whole matrix -- no missing values, a value range -- comes from a full scan. A spot check or a log line proves only what it covered.
- Say once, near the top, when and from what the numbers were read, then name the source file beside each table ("Counted from `keep/<version>_sample_keep.tsv`"). Do not stamp individual numbers `measured`; a number from outside the bundle names its source in the same sentence. Never carry a projected number: if something is still projected, the document is not ready to ship.
- Hunt for what the data itself will do to a user -- missing values, stale columns, keys that are not unique, rows that are aggregates, features that do not measure what the rest do -- and put each in the quick start and the caveats with its count. A caveat nobody measured is a guess.

Three presentation differences from an operational document:

- No resume block, job ids, task list, or `Open` section. Build verification belongs in the provenance section at the end.
- The one end-to-end overview diagram may exceed the fifteen-node guide when it is grouped into subgraphs and still reads top to bottom. Every other diagram keeps the limit, reads the role classes as `mermaid-patterns.md` section C.3 does for a finished product, and draws every removal as a funnel (section D.10 there).
- The bold budget of section 7 is per section rather than per screen: a quick-start list of must-not-miss facts may bold each one.

---

## 12. Rendering to a single-file HTML

When a Markdown document must travel outside the repository as one file -- attached to an email, handed to a collaborator without the repository, or archived beside a data release -- render it once to a self-contained `.html`: no external stylesheet, script, or CDN reference, so it opens correctly years later on a machine with no network access. This applies most often to reference documents (section 11), but works for any Markdown file.

Use the bundled script rather than improvising a pandoc invocation each time. A hand-rolled command silently drops Mermaid diagrams -- pandoc emits a fenced mermaid block as an inert code block, not a rendered figure -- unless the block-format fix and a Mermaid renderer are both wired in.

```bash
pixi run bash .github/skills/markdown-doc/scripts/render_html.sh path/to/README.md
# writes path/to/README.html beside the source; pass a second argument to redirect it elsewhere.
```

What it does, in order:

- Converts with `pandoc --standalone --embed-resources`, folding in the shared stylesheet (`assets/doc.css`) so every document this skill renders looks the same.
- Runs a Lua filter (`scripts/mermaid-codeblock.lua`) that turns every fenced mermaid block into `<pre class="mermaid">` -- the markup Mermaid's renderer expects -- instead of the plain code block pandoc emits by default.
- If the document contains at least one Mermaid block, embeds a vendored copy of mermaid.js (`assets/mermaid.min.js`) directly in the file and calls `mermaid.run()` on load. The library is vendored, not linked from a CDN, so the diagrams render with the browser's own JavaScript engine and no network call -- open the file offline and every diagram still draws itself.
- Reports the byte size and the Mermaid block count, and warns if the rendered count does not match the source.

**Verify by opening the file in an actual browser**, not an editor preview or a mail client's HTML preview -- both can disable script execution, which leaves the Mermaid source sitting as inert text. A rendered diagram is an inline `<svg>` inside `<pre class="mermaid">`; `grep -c '<svg'` is a cheap proxy when a browser is not at hand, though it only proves the script ran, not that the diagram is legible.

Read [references/html-export.md](./references/html-export.md) for the asset contract, how to refresh the vendored Mermaid version, and troubleshooting.

---

## References

- [references/mermaid-patterns.md](./references/mermaid-patterns.md) -- the direction rule, the documentation palette, and one copyable template per question a diagram can answer, including the funnel and aggregation diagrams a release needs.
- [references/reference-docs.md](./references/reference-docs.md) -- for a reference document: where each kind of content comes from, the measurement script, the release-README skeleton, the table and diagram sets, the recipe contract, the surprises to hunt for, bundle mechanics, and the pre-ship checklist.
- [references/html-export.md](./references/html-export.md) -- the single-file HTML export contract behind `scripts/render_html.sh`: the asset list, how to refresh the vendored Mermaid version, and troubleshooting.
