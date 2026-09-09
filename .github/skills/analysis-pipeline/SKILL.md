---
name: analysis-pipeline
description: "Author, run, and monitor multi-step analysis pipelines organized as numbered stage directories (src/NN-stage/) with paired script + .md files. Use when adding or modifying a pipeline stage, porting a lab pipeline, writing a PLAN.md, PROGRESS.md, or DECISION.md, scaffolding config.sh/.env/env-manager activation, running the plan -> approve -> implement -> syntax-check -> progress -> run -> verify cycle, or resuming a pipeline across sessions. Covers script headers, stage layout, paired documentation, plan-first gating, progress tracking, and the decision log that carries each choice, the alternative it beat, and why between sessions. Applies to HPC/cluster bioinformatics (GWAS, ancestry, RNA-seq, methylation, single-cell) and any cohort-scale or batch data pipeline."
---

# Analysis pipeline stages

Authoring and execution guide for multi-step pipelines that follow the
numbered-stage convention. The rules below are the minimum contract every
stage must satisfy.

This skill is portable and names no project value. The env files, output roots,
color files, scheduler, and legacy naming exceptions come from the repository's
`.github/instructions/` bindings, its `copilot-instructions.md`, or its
`AGENTS.md`. Those take precedence where they conflict.

Related skills: `long-running-jobs` (tmux / cluster arrays),
`data-verification` (schema + output checks), `data-result-layout`
(where outputs go), `image-prompt` (the stage's `DIAGRAM.md` figure brief).

---

## 1. Stage layout

A stage lives inside a **track**: any top-level analysis directory the
repository bindings list -- `src/`, `src_<variant>/`, `pipeline/`,
`workflow/`, or a directory named after its pipeline. A track is either
**staged**, holding numbered `NN-stage-name/` directories (`NN` for a
top-level stage, `NN.MM` for a sub-stage, e.g. `04.03-gwas-echo-cm-top-hits`),
or **flat**, where the track directory itself is the single stage and
everything below sits at its root. `src/NN-stage/` throughout this skill means
"a stage in either layout". Each stage contains:

- `PLAN.md` + `PROGRESS.md` + `DECISION.md` — stage-level overview; see
  section 1.1
- `{date}-{short-title}.PLAN.md` + `{date}-{short-title}.PROGRESS.md` +
  `{date}-{short-title}.DECISION.md` — one set per specific question; see
  section 1.1
- `DIAGRAM.md` — optional figure brief for the stage; a question-scoped figure
  is `{date}-{short-title}.DIAGRAM.md`. Owned by the `image-prompt` skill.
- `AGENTS.md` — track guide: run commands, entry points, pitfalls, known drift,
  key output paths, and a pointer to the active plan/progress/decision set.
  Required at the track root and optional per stage; the rules are in
  `copilot-instructions.md` under "Track Guide". Update it in the same change
  set when a run command or output path changes. When it disagrees with the
  active progress file, the progress file wins.
- `config.sh` + `config.sh.md` — required once the stage has shell or cluster
  scripts; sourced by every one of them
- `config.R` + `config.R.md` — the R-side twin, required once the stage has R
  steps; holds `STAGE`, the stage constants, `stage_paths()`, and the shared
  `fn_*` helpers, and is sourced by every R step. Never a color literal: source
  the track `color.R` instead. A stage with both kinds of step ships both
  configs. Stages predating this use a single `config.md`, or `00-config.R` +
  `00-config.md`; leave those names alone.
- `NN-step-name.{sh,lsf,sbatch,R,py}` + sibling `NN-step-name.md` — **one task per file**
- Cluster wrappers (`.lsf` / `.sbatch`) submit array jobs; the matching
  `.sh` / `.R` / `.py` runs the single unit of work. Never inline the unit
  work inside the wrapper.

Numbering, `PLAN.md`, the directory listing, and the documented execution
order must all agree.

## 1.1 Plan, progress, and decision files (two tiers)

Three files travel together at each tier. The plan holds the intent, the
progress file holds the state, and the decision file holds the choices: at each
fork, what was chosen, what it beat, and the one reason why. The first two let
a new session resume the work; the third stops it from re-opening a settled
choice or proposing an alternative that was already rejected.

**Tier 1 — stage level, one of each, always:**

- `PLAN.md` — stage charter: goal, data inventory, step index, cross-cutting
  decisions that hold for every campaign in the stage.
- `PROGRESS.md` — stage dashboard: current state in one screen, plus an index
  table linking every tier-2 campaign with its status.
- `DECISION.md` — stage decision log: the cross-cutting choices, numbered
  `D1..Dn`, each with the alternative it beat and the evidence.

**Tier 2 — one specific question, one set per question:**

- `{date}-{short-title}.PLAN.md` — the design for that question only, with its
  own `Q1..Qn` blocks; approved before implementation.
- `{date}-{short-title}.PROGRESS.md` — run state for that question only: job
  IDs, output counts, error log, resume block.
- `{date}-{short-title}.DECISION.md` — decisions for that question only,
  including the outcome of each `Q` once approved or overridden.

Naming: `{date}` is `YYYY-MM-DD` of the day the plan is written;
`{short-title}` is lowercase kebab-case, 2-5 words, naming the _question_
(`gwas-replication-exact-match`), not the stage. The `.PLAN.md` /
`.PROGRESS.md` / `.DECISION.md` suffix is uppercase — that is what separates a
live tier-2 file from a historical lowercase `*.progress.md`.

Which tier to write:

- New stage, or a stage that exists to answer exactly one question with no
  follow-up expected → the tier-1 set carries it; no tier-2 files.
- Any later question in that stage — a revision, a rerun with changed
  decisions, a new sub-analysis → a new tier-2 set. Do not retro-split the
  first question out of the tier-1 files; leave them as the stage overview and
  add the campaign to the `PROGRESS.md` index table.
- Never open a dated section inside `PROGRESS.md` or `DECISION.md` for a new
  question. A new question is a new set of files.

Stages created before this convention use lowercase `plan.md`,
`plan.<aimodel>.md`, and dated `{date}-{short-name}.progress.md`, and have no
decision file at all. Those are historical, read-only records: read them when
resuming a stage, but never rename them and never write new run state into
them. New work in such a stage adds the tier-1 set (if missing) plus a tier-2
set for the new question, alongside the historical files. A decision file
added to an existing stage starts at `D1` with the first decision made that
day; it does not reconstruct earlier decisions from chat that no longer exists.

`src_<variant>/NN-stage/` (e.g. a second cohort or a re-analysis) and a flat
pipeline track (`pipeline/`, `<name>_pipeline/`) follow the same contract as
`src/NN-stage/`. Each track has its own `color.R` and its own `config.*`; a
script never sources those of another track.

## 2. Script headers

Every script carries the same metainfo block, comment-adapted per language,
placed after the shebang and after any scheduler directives. Never write
"Generated by <AI tool>" attribution.

```text
# @AUTHOR: <Name>
# @CONTACT: <email>
# @DATE: YYYY-MM-DD
# @DESCRIPTION: One-sentence purpose. Usage line if it takes args.
# @VERSION: v0.1.0
```

Full shell / LSF / SLURM / R / Python templates and the `config.sh` +
`config.R` contracts:
[references/script-templates.md](./references/script-templates.md).

## 3. Paired `.md` for every script

Every `.sh / .lsf / .sbatch / .R / .py` ships a sibling `.md` with these
sections:

- **Purpose** — one paragraph
- **Inputs** — bullet list of exact paths or file patterns
- **Outputs** — bullet list of exact paths or file patterns
- **Run** — one or more copy-pasteable commands
- **Verify** — 1-3 exact shell commands that confirm the step succeeded

The `.md` documents the script's _contract_. Do NOT write `.md` files that
describe what changed in a given session.

One exception: `config.sh` and `config.R` are not executable, so their docs
drop **Run** and **Verify**. They also append `.md` to the full filename
(`config.sh.md`, `config.R.md`) rather than the stem, since both would
otherwise claim `config.md`.

## 4. Plan-first gate (mandatory)

Applies when BOTH hold:

- The work creates or modifies scripts inside any track the bindings list:
  `src/NN-stage/`, `src_*/NN-stage/`, or a step of a `pipeline/`, `workflow/`,
  or flat pipeline directory, AND
- It creates a new stage, adds a new step to an existing stage, or would
  invalidate existing outputs under any output root the bindings list for that
  track: `results/`, `data/intermediate/`, or a flat track's own data root. A
  file under a track's env-bound root counts exactly like one under
  `results/`.

A new stage or a new step always triggers the gate, even though it invalidates
nothing. When the change only edits an existing step, decide by naming which
already-produced files it invalidates. "None" means act directly -- that covers
a bug fix inside one existing script, re-running an existing step, and
doc/config edits.

1. **Plan** — write the plan file for this question (section 1.1: tier-1
   `PLAN.md` for a new stage, tier-2 `{date}-{short-title}.PLAN.md` for a later
   question). Lock contested design decisions in numbered Q&A blocks
   (Q1, Q2, ...) with a recommended answer for each. Do NOT guess; do NOT
   answer the questions in chat instead of in the plan. Create the matching
   decision file in the same step: the alternatives rejected while writing the
   plan are its first `D` entries (section 6).
2. **Review gate** — present the plan and wait for explicit approval.
   **Do NOT implement until the user approves.** On approval, record each
   `Q` outcome as a `D` entry — especially where the user overrode the
   recommendation, with the reason they gave.
3. **Implement** — write scripts in numbered order, each with its `.md`.
4. **Syntax gate** — run before any execution (section 5).
5. **Progress file** — create the matching progress file before running
   anything, and link it from the stage `PROGRESS.md` index (section 6).
6. **Run** — execute sequentially; diagnose and fix errors inline, then re-run.
   When a fix could take more than one form, the form chosen and what it beat
   is a `D` entry; a choice taken without asking the user is marked
   `provisional`.
7. **Verify** — run the Verify block from each step's own `.md` (section 3),
   confirm the outputs exist, are fresh, and have sane row counts, then print
   the key results to the user. A Verify block that is never executed is not a
   check.

Plan structure and the two-tier file convention:
[references/plan-and-progress.md](./references/plan-and-progress.md). The
stage's figure brief is a separate `DIAGRAM.md`; see the `image-prompt` skill.

## 5. Syntax gate

Run before execution, every time. **When the repository has a test suite, it is
the first rung** -- run it before the syntax parse, because it catches whole
classes of drift a parse cannot, such as a wrapper's hard-coded array falling
out of step with a config constant. Check for one before assuming there is none.

**The formatter and linter are the second rung.** When the repository ships a
format task or a lint config, run them on the changed files before the parse:
they are cheaper than a run and a formatting diff is a code change that belongs
in the same change set. The repository bindings name the commands; never run a
whole-repository format from a skill when a changed-files form exists.

```bash
find src/NN-stage -maxdepth 1 \( -name '*.sh' -o -name '*.lsf' -o -name '*.sbatch' \) -print0 | xargs -0 -r -n1 bash -n
pixi run Rscript -e 'for (f in list.files("src/NN-stage", "\\.R$", full.names=TRUE)) tryCatch(parse(f), error=function(e) stop(f, ": ", conditionMessage(e)))'
```

Python equivalent: `pixi run python -m compileall -q src/NN-stage` in a pixi
repository; `uv run python -m compileall -q src/NN-stage` only where there is
no `pixi.toml`.

A parse check cannot catch an argument validated against a constant that is
defined later, or a wrapper array that disagrees with the config. Only invoking
the driver, or a test that compares the two, will. Add the test rather than
relying on the smoke run to find it.

Prefer the fastest meaningful test. Never launch an expensive pipeline when
a syntax check or a single-unit smoke run would catch the bug first.

## 6. Progress and decision files (mandatory)

**Progress file — before the first run.** Create the progress file that
matches the plan file for this question (section 1.1) before the first run. It
is the single source of truth for the current run state, for both the user and
any future agent session.

Required sections: Resume block, Todo list, Job table, Output file counts,
Error log, Key paths. A tier-2 progress file carries all of them for its own
question; the stage `PROGRESS.md` keeps the one-screen current state plus the
campaign index table and never duplicates a tier-2 job table.

**Update immediately** — never batch or defer — after: a step completes, a job
is submitted, an error is diagnosed, a script/config is modified, a job is
re-submitted, or the session is about to end. When a tier-2 campaign starts or
finishes, update its row in the stage `PROGRESS.md` index in the same edit.

**Decision file — with the plan.** Create the decision file that shares the
plan's stem in the same step as the plan itself (section 4, step 1). It is a
log of decisions, not of reasoning: one numbered entry `Dn` per choice, each
with four fields — `Instead of`, `Why`, `Evidence`, `Status` — and nothing
else. The test for an entry is that it can be written as "chose X over Y
because Z"; narrative, hypotheses considered, and tool calls fail that test
and stay out.

Entries are dated, numbered in the order made, and append-only. A reversal is
a new entry that marks the old one `superseded by Dn`; a choice taken without
asking the user is `provisional` until confirmed.

**Update immediately** the moment: a choice is made between alternatives while
planning, implementing, or debugging; the user approves or overrides a `Q`; a
choice is taken without the user; an earlier choice is reversed. When the
progress file changes for the same reason, update both in one edit. A tier-2
decision that holds for the whole stage is copied into the stage `DECISION.md`
in that edit.

**Resume order:** `PROGRESS.md` for where things stand, `DECISION.md` for what
is settled and what is still provisional, then `PLAN.md` for the design as
needed.

Templates, the what-goes-where table, and the multi-session check-fix-submit
cycle: [references/plan-and-progress.md](./references/plan-and-progress.md).

## 7. Execution rules

- Test one unit before the full array: submit index `[1-1]`, verify its output
  is non-empty and correct, then submit the full range.
- The inline / tmux / cluster thresholds are owned by the `long-running-jobs`
  skill. Never poll with `sleep`.
- Verify real outputs after an array finishes (`ls <pattern> | wc -l`, row
  counts) — a DONE status does not prove the work happened.
- When a script or config is fixed mid-pipeline, re-run a downstream step if it
  consumes the changed file AND its outputs are missing, or older than the
  change (compare `stat -c '%Y %n'` on both).

## 8. Porting an external pipeline

- Reuse the source repo's helpers by absolute path
  (`Rscript ${LAB_DIR}/R/NN_*.R ...`) rather than re-implementing them.
- Mirror the source step numbering inside one submitter `.sh` per multi-step
  chain; do not split each external step into its own file unless it buys
  clarity.
- Keep environment-specific paths (reference panels, tool installs) in
  `config.sh`, never hard-coded in step scripts.

## 9. Forbidden patterns

- Deleting a superseded step script, stage directory, or its outputs. Rename
  it `<name>_stale-YYYY-MM-DD` instead; see `data-result-layout`.
- `bash -c "..."` sub-shells inside agent tool calls.
- Globally installed conda/mamba envs for R or Python; always use the
  project's pixi / uv / conda env.
- Comments that restate what the code does or which task added them. Comments
  only for non-obvious _why_: hidden constraints, workarounds, invariants.
- `.md` files that explain a diff instead of a script contract.
- A decision file that narrates reasoning, lists hypotheses, or restates the
  plan. It holds only entries of the form "chose X over Y because Z", with
  evidence.
- Emoji anywhere in code, commits, or docs.
- AI-tool attribution lines in script headers.

---

## References

- [references/script-templates.md](./references/script-templates.md) — shell,
  LSF, SLURM, R, Python headers and the `config.sh` contract.
- [references/plan-and-progress.md](./references/plan-and-progress.md) —
  plan structure, progress template, decision-log template and what-goes-where
  table, the two-tier file convention, multi-session run/monitor cycle.
- [references/pipeline-pattern.md](./references/pipeline-pattern.md) — worked
  GWAS / local-ancestry stage as the canonical 10-file shape.
