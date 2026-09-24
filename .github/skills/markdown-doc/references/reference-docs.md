# Reference: documents that travel with data

The procedure behind `SKILL.md` section 11: where each statement may come from, how every number is measured, what a data-release README contains and in what order, which tables and diagrams carry it, and what to check before it ships. Everything here is portable; the concrete paths, versions, units, and tiers come from the project.

Three terms are used throughout:

| Term | Means |
| --- | --- |
| version | a variant of the product whose items are processed and tabulated separately: a platform, an assay chemistry, a panel |
| unit | one processing run inside a version whose output ships as its own product; a version with one run has one unit, and some items of a version may be in no unit |
| primary unit | the unit most analyses use; the project names it |

Read it before writing or rewriting a release README, a data dictionary, or a usage guide that leaves the repository.

---

## A. Where each kind of content comes from

A reference document states four kinds of content, and each has one legitimate source. Reading a plan or a decision record for a rule or a reason is correct; copying a number out of one is not.

| Content | Read it from | Never from |
| --- | --- | --- |
| Counts, shapes, sizes, value ranges, column names, key formats | the shipped files, through a measurement script (section B) | a plan, a progress file, a proposal, an earlier README, a chat message, or a log line when the shipped files can answer |
| Rules, thresholds, parameters, definitions | the producing code at the recorded commit, and the shipped provenance file | memory, a generated stub, or a plan the code may not implement |
| Why a choice was made | the decision record or plan, restated inline as the rule and its one reason | a link to that record |
| What ran for this release, and how the build was verified | the run's own records -- provenance, runtime files, the producer's verification output -- each named in the sentence | a progress file's summary of them |
| Genome build, coordinate base, chromosome naming | the coordinate column names and the producing code; if the bundle records none, say where it was read | the platform's usual default |
| A publication behind a rule | a plain-text citation: first author, year, journal, and the DOI as plain text | a link |
| What earlier releases contained | their READMEs and decision records, and their files if still on disk; each number labelled with its release | an old number presented as current |

When the provenance does not say whether the working tree was clean at the recorded commit, compare the producing code at that commit with the working copy, and read the rule from the commit.

A number the bundle cannot reproduce -- a value from a run log, or one measured during an earlier build -- appears only with its source in the same sentence: "from the normalize run log of 2026-01-31", "measured during the v1 build". Arithmetic on measured numbers is fine when the sentence says what it adds.

"The previous README" means the README of the previous release of the same product, whether its repository copy or the one in the previous bundle. Take its **structure** and its **history** from it; never take a current number from it. When one exists, start from it: keep its sections and their numbers, so a reader of both finds the same thing in the same place, and add any part of the section C skeleton it lacks. Without a predecessor, use section C as it stands. The rule text of a generated stub is a lead to check in the code, like a plan, never a source.

---

## B. Measure before writing

Numbers come out of a script, not out of memory. Write one read-only measurement script per document in the repository's scratch root (never the system `/tmp`), have it print one labelled value per line, and keep its output log until the document has been reviewed. Every number in the document is then either in that log, arithmetic on it, or a number that names another source. -> `data-verification`

```text
version_A items | 412
version_A person_id unique | 390
version_A drop_reason | detection_failure | 3
```

The steps, in order:

1. **Inventory the bundle.** List every file with its logical size and count the files per directory. Say whether a quoted footprint is logical bytes or allocated blocks; a network filesystem can report twice the logical size.
2. **Read every header** before counting anything: text-table headers, parquet schemas, matrix dimension names. Save the column lists; the column dictionary is generated from them (section E).
3. **Treat blank and missing alike** in text files. Many writers write a missing value as an empty cell and many readers return it as an empty string, so an `is missing` test silently matches nothing.
4. **Count** everything in B.1.
5. **Scan values in full.** Missing cells, rows or columns that are entirely missing, and the value range come from a scan of every cell, read in blocks along the storage order (whole columns for a column-major format). A range or a "no missing values" claim taken from a first block, a sample, or a log line is not a property of the matrix: scan in full, or say exactly what was checked. Scan every matrix in one format in full; for the other formats compare shape and ids, and say what value check was made. A full scan is long work; launch it the way the project launches long jobs.
6. **Run every recipe** and keep its output (section G).
7. **Reproduce every published selection** the document describes -- a cohort list, a one-per-person rule -- from the shipped columns alone. State an exact match as such; investigate a mismatch before writing, because it usually exposes a stale column.
8. **Compare metadata with its current source of truth** -- the current sample manifest, the current genotype freeze. A shipped column that disagrees is stale: count the rows, say when the shipped copy was taken, and make it a caveat.
9. **Check the identities** in B.2. Only then write.

### B.1 What to count

| Area | Count |
| --- | --- |
| Items and keys | rows per version and per unit; for every key: rows, unique non-missing values, missing values |
| Composition | the full value table of every categorical column with a small fixed vocabulary -- group, tissue or label, sex, study arm -- per version, with totals; for a column with many levels (a slide, a plate, a batch), the number of levels and the range of level sizes |
| Removal | every gate's count per version; the full value table of every reason column, with the exact strings; items carrying more than one reason |
| Flags | each flag over all items, over shipped items, and over the primary unit |
| Feature filter | each source list's size and its matches per version; the union; the overlaps that change a decision; the sequential funnel; every universe the files use |
| Aggregated features | groups, members, the size distribution, groups whose members disagree |
| Products | the shape of every matrix in every format and of every sidecar |
| Values | missing cells, all-missing rows and columns, minimum and maximum, per matrix, from the full scan |
| Recipes | the output of each |
| Companions | the row count of every file outside the bundle that the document describes |

### B.2 Identities to check before writing

- Every funnel adds up: items in = kept + removed at every step, and the units of a version plus its items in no unit sum to the version total.
- Every reason table sums to its row count, or the document says how many rows carry more than one reason.
- A tier or product declared a subset of another is one: the anti-join is empty.
- One quantity read from two files agrees: a keep-list count with the matrix rows, a sidecar with the matrix dimensions, each storage format with the others.
- Counts on different universes are never compared without naming both.

---

## C. The section skeleton of a data-release README

Use this order, or the predecessor's (section A). Keep the numbering stable, so a caveat can say "section 7.4": leave out a part the product does not have, but when the previous release's README had a section this release no longer needs, keep its number with one line saying why it is empty. The contents list is plain text; an anchor link fails the link check of section J. Every section carries its tables and diagrams; a sentence that summarizes a table is not the table.

| # | Section | Carries |
| --- | --- | --- |
| -- | Title and identity | a two-column identity table (field, value): release version; generated date, from the provenance file; producer and software versions; commit; items processed (entering the first step) per version; items shipped (in a product) per unit; deliverables (product stems, and the files per stem); footprint, logical and allocated; genome build (section A); what it supersedes. Then one sentence that the file is self-contained and when its numbers were read, and a plain contents list |
| 0 | Quick start | a runnable minimal load that asserts the sidecars align; any exception the user must apply at once; a numbered list of the surprises of section H that are present, each bold, one sentence, with its count and a section pointer |
| 1 | What is in this directory | the tree as a code block, with counts and sizes; the directory map diagram; a table of companion files that are not in the bundle and where they sit |
| 2 | How it was made, end to end | the overview diagram; the manual gates; a table of which steps ran for this release and why |
| 3 | Inputs and their validation | the validation chain diagram; the composition table per version with totals; the value tables of the categorical calls; special subsets as a short list with counts |
| 4 | Per-item QC and its flags | the flag-flow diagram; how each flag that is not self-explanatory is computed |
| 5 | What was removed, and why | the gate table with the exact reason strings; the per-unit and whole-bundle funnel tables; the criteria table; the flag counts; one funnel diagram per version; the removed ids; every selection rule with its counts and how to reproduce it; the measured QC signals |
| 6 | The tiers | the tier diagram; the tier table; what each tier is and is not |
| 7 | The feature filter | the criteria table naming the column that shows each; the funnel diagram per version; the reason table with totals; the universe note; the source-list table with matches per version and overlaps; feature classes; aggregated features |
| 8 | Processing parameters | the unit diagram; the parameter table, values from the code; which values are conventions and which were fitted; scale warnings |
| 9 | The files | the inventory (stem, tier, shape, size per format); formats and how to read each; sidecars; the full value scan; small slices or samples |
| 10 | Recipes | section G |
| 11 | Column dictionary | section E |
| 12 | Caveats | numbered; each opens with a bold sentence, gives its count, and says what to do |
| 13 | Related, not in the bundle | derived products with their own funnel; companion lists; parallel pipelines that are not interchangeable |
| 14 | Release history | one table across releases, then one paragraph per transition; every older number labelled with its release |
| 15 | Provenance and integrity | checksum commands with entry counts and what they cover; the provenance keys; the build verification and its result; the re-checks made while writing; build notes: defects found and fixed, known cosmetic defects |

Removed ids: a table when there are about thirty or fewer, a collapsible `<details>` block holding the table up to a few hundred, and beyond that a recipe that lists them.

---

## D. Table patterns

Per-version composition, with a total only where the quantity adds. A dash marks a cell that does not add, never a missing measurement. Typographic dashes and minus signs are fine in tables and prose; diagram labels and code stay ASCII (`mermaid-patterns.md` A.4).

```markdown
|                          | Version A | Version B |     Total |
| ------------------------ | --------: | --------: | --------: |
| Items                    |       412 |       188 |       600 |
| Unique persons           |       390 |       170 |       541 |
| Persons in both versions |         — |         — |        19 |
| Group X / Group Y        | 300 / 112 |   90 / 98 | 390 / 210 |
```

Gates, one row per gate, with the string the gate writes:

```markdown
| Gate | Removes | Version A | Version B | `drop_reason` it writes |
| --- | --- | --: | --: | --- |
| QC | a failed measurement | 3 | 2 | `detection_failure` |
| Biospecimen | an unusable sample type | 5 | 1 | `unusable biospecimen` |
```

A funnel, removals as negative rows and the shipped count in bold:

```markdown
| Level                  |   Items |
| ---------------------- | ------: |
| In the source          |     600 |
| minus biospecimen gate |      −6 |
| minus QC gate          |      −5 |
| **Shipped**            | **589** |
| In no product          |      11 |
```

A reason vocabulary, with the exact strings, an `empty` row, and rows that sum to the table:

```markdown
| `drop_reason`          | Version A | Version B | Means            |
| ---------------------- | --------: | --------: | ---------------- |
| empty                  |       404 |       185 | shipped          |
| `detection_failure`    |         3 |         2 | QC gate          |
| `unusable biospecimen` |         5 |         1 | biospecimen gate |
```

A feature filter has three different counts, and a column says which one it holds: **marginal** (flagged by this source), **union** (flagged by any), and **sequential** (removed at this step, after the earlier steps). Put the marginal counts and the union in the source table, and the sequential removals in the funnel.

The inventory is `| Stem | Tier | Features x items | parquet | tsv.gz | gds |`, with the unit in the sentence above it ("GB of 10^9 bytes, logical"). The value scan is `| Matrix | Cells scanned | Missing cells | Non-finite | Rows entirely missing | Columns entirely missing | Min | Max |`, plus a column for the missing cells outside the all-missing rows when there are any; a column that is zero for every matrix may go when the sentence above says so.

---

## E. The column dictionary

- One subsection per table family. State its row count, its column count, and how the families relate: a sidecar is the annotation plus one column, a keep list is the annotation plus the decision columns.
- List **every** column, generated from the header lists saved in section B, and group them by purpose when there are many.
- Then a meaning table for each column whose meaning is not obvious from its name: its unit, its coding, which of two alternatives is preferred, any known staleness.
- For each categorical or reason column, its full vocabulary with counts.
- Say which expected columns do **not** exist, and where the value lives instead.
- A trivial format -- a checksum file, key/value provenance, a runtime record -- gets one line.

---

## F. Diagrams

A funnel is both a diagram and a table: the diagram shows where items split, the table carries the exact counts and reason strings. Both take their numbers from the same measurement log, and the table is authoritative.

Draw each diagram below when the product has the part it describes. Node labels carry counts and shapes, never job ids, states, or runtimes. The role classes are in `mermaid-patterns.md` section C.3. A funnel merges criteria that act at the same point into one removal node, so it stays under fifteen nodes; the funnel table carries the split.

| Diagram | Answers | Pattern |
| --- | --- | --- |
| Directory map | what in the bundle is derived from what | D.4 |
| End-to-end overview | inputs to products, grouped by phase | D.1 subgraphs; may exceed fifteen |
| Validation chain | how inputs were admitted | D.10 as a chain |
| Flag flow | which flags remove items and which only report | D.10 |
| Item funnel, one per version | where items were removed and which unit the rest went into | D.10 |
| Tiers | which tiers exist, which ship, how each is derived | D.10, unshipped tiers as `ext` |
| Feature funnel | where features were removed; both versions as subgraphs while the pair stays under fifteen nodes, else one per version | D.10 |
| Filter sources | which sources are operative, report-only, or retired | D.10 |
| Aggregation | which features were merged into one row | D.11 |
| Processing units | which runs produced which products | D.10 |
| Derived products | what was built from the bundle, outside it | D.10 |

---

## G. Recipes

- Run every recipe against the shipped bundle before the document ships; a count in a comment is that run's output.
- Write paths relative to the bundle root, with one placeholder for the root.
- Name the packages once, at the top of the section. Write the recipes in the audience's language, or the producer's when the audience is unknown; add a second language only where the audience uses it.
- Warn in a comment where a line loads a large matrix into memory, and show the column-pruned read beside it.
- Cover at least: reading a product with its sidecars and asserting alignment; reading the text tables with blank treated as missing; applying every quick-start exception; reading a subset of items without loading the whole product; rebuilding any tier that is not shipped; reproducing each published selection; listing the items each gate removed; any transform the product omits for some tier; prototyping on the small slice; and a covariate table naming the column behind each suggested covariate.

---

## H. Surprises to hunt for

Each surprise found becomes a quick-start item and a caveat, with its count. Look for every one; one that is absent needs no sentence.

| Surprise | How it shows up |
| --- | --- |
| Keys that are not unique | a person or sample id repeated across items, or missing on some rows |
| Features with no value | rows or columns entirely missing in the full scan |
| Features that are not the analyte | control, variant, or off-type features, usually by an id prefix |
| Rows that are aggregates | a member-count column above 1 |
| Stale metadata | a column that disagrees with its current source of truth |
| Encoding traps | blank written for missing; types lost in a text format |
| Different scales | products processed in separate runs |
| Different universes | one filter counted on different feature sets |
| Coordinates | genome build, coordinate base, chromosome naming |
| A spot check presented as a full property | a range or a missing-value claim that came from a sample or a log line |
| A shipped file that points outside | a provenance value naming an unshipped plan; explain the reference inline |
| A design that looks like a defect | an unshipped tier, a slice without a unit in its name; confirm it from the code and say so |

---

## I. Bundle mechanics

- **One document, two copies.** Keep the canonical copy in the repository beside the producer and copy it into the bundle root at release; the two stay identical. Say whether the checksum manifest covers it.
- **Copying it in changes nothing else in the bundle.** The stubs and checksum manifests stay as the build wrote them. An older README at the same path is renamed with the project's stale suffix, never overwritten. -> `data-result-layout`
- **Paths are code spans, never links**, even to files inside the bundle: a link that resolves in the bundle breaks in the repository copy, and the reverse.
- **Companion files** are named by a path relative to a stated root and by the track or person that owns them, never by an absolute cluster path.
- **Generated stubs** in the bundle, such as per-version READMEs written by the producer, are superseded by this document. Check their numbers against the measurement log; a disagreement is a finding, not a choice.
- **A defect found while measuring** goes into the caveats with its count, and to the producer's owner as a separate report, through the channel the project uses for defects. The README never promises a fix.

---

## J. Before it ships

Four of the checks are commands. Run them on the finished file; a command kept inside a table cell needs its pipes escaped, and the escaped copy silently matches nothing.

```bash
grep -nE '\]\(' README.md                               # links: must print nothing
grep -niE 'projected|TBD|pending|job [0-9]+' README.md  # run-state words: read every hit
grep -nE '[`( ]/[A-Za-z]+/' README.md                    # absolute paths: only placeholders may remain
cmp README.md "$bundle_root/README.md"                    # once the bundle copy exists: must print nothing
```

| Check | Passes when |
| --- | --- |
| No links | the first command prints nothing |
| No run-state words | every hit of the second command is removed, or is history that names its release |
| No machine paths | every hit of the third command is a placeholder such as `/path/to/<bundle>` |
| Copies identical | the fourth command prints nothing |
| Every number sourced | each number is in the measurement log, is arithmetic on it, or names its source in the sentence |
| Repeated numbers agree | a count stated in several sections -- identity table, quick start, a funnel, a caveat -- is the same in each |
| Tree and inventory match the bundle | every path in the tree exists, and every file in the bundle is covered |
| Dictionary matches the headers | each dictionary subsection lists exactly the saved header of its family |
| Tables add up | the identities of B.2, re-checked on the finished tables |
| Recipes ran | every count in a recipe comment equals that recipe's output |
| Diagrams | palette hexes only, ASCII labels, one sentence above each, and a preview without an error box; an agent that cannot preview says so in its report |
| Predecessor coverage | every section of the previous release's README exists here, or the history section says why it went |
| On disk | run the project's formatter on this file only, then re-read the disk copy |
