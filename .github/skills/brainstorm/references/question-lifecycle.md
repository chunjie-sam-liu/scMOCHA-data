# Reference: the question lifecycle and technical depth

The worked detail behind `SKILL.md` sections 5 and 7. Read this when opening the
first `Q` block of a brainstorm, when closing one, when a round adds evidence,
or when an unfamiliar tool or a contested parameter has to be made legible.

---

## A. The `Q` block, worked

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

---

## B. Status, and how a `Q` closes

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

---

## C. Rounds

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

---

## D. Parameters

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

---

## E. Unfamiliar software

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
[domain-checklists.md](./domain-checklists.md).
