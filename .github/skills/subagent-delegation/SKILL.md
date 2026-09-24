---
name: subagent-delegation
description: "Coordinate work across explorer, worker, runner, monitor, and reviewer subagents instead of doing everything inline. Use when a task spans more than two files, when the right files are not yet known, when tracing how code works or what calls what, when a script/test/build/cluster submission must run, when a submitted job must be checked or its outputs verified, when non-trivial edits need review, or when several independent areas must be understood in parallel. Covers the roles including the user-invocable coordinator, automatic delegation triggers, the subagent prompt contract, parallel vs sequential rules, and model selection."
---

# Subagent delegation

Delegate scoped subtasks instead of doing heavy reading, wide searching, and
command execution inline. Keep the coordinating context light.

Subagents live in `.github/agents/`. Delegate with the `runSubagent` tool; it
must be enabled in the tools picker for delegation to work. If delegation is
unavailable, continue inline and keep the same role boundaries.

---

## 1. Roles

| Agent         | Scope                                                                                            |
| ------------- | ------------------------------------------------------------------------------------------------ |
| `coordinator` | Plans and delegates; user-invocable only, the model cannot start it                              |
| `explorer`    | Read-only mapping, dependency tracing, entrypoint discovery, execution-path summaries            |
| `worker`      | Minimal scoped edits after the relevant path is understood                                       |
| `runner`      | Targeted command execution with concrete output reporting                                        |
| `monitor`     | What happened after a job was submitted: state, outputs, error logs, progress and decision files |
| `reviewer`    | Correctness, regression, statistical-validity, and scope audit                                   |

`runner` and `monitor` split one job in two. `runner` submits and reports the
submission; a `bsub` returns in a second and tells you nothing about the work.
`monitor` is what comes back hours later to read `bjobs`, the per-index `.err`
files, and the real outputs, and to update `PROGRESS.md`. A repository whose
work is mostly cluster arrays uses `monitor` more than any other agent.

`coordinator` is the entry point for a whole multi-step task and is the only
one the model cannot start on its own (`disable-model-invocation: true`). The
other five are `user-invocable: true` as well, so the user can call any of them
directly. `coordinator.agent.md` points at sections 2-6 of this skill rather
than copying them, so this file is the only place they are written.

The editor also ships a built-in `Explore` agent. Prefer `explorer` for work in
this repository, since it carries the repo's conventions; use `Explore` only for
a generic throwaway lookup.

### What a subagent inherits

Measured 2026-09-21 with a smoke-test delegation, so no session has to guess at
it or re-derive it. Delegation works; what does not carry over is listed on the
right.

| Inherited                                                                        | Not inherited                                                            |
| -------------------------------------------------------------------------------- | ------------------------------------------------------------------------ |
| `.github/copilot-instructions.md`, in full                                       | the `<agents>` list, so a subagent cannot delegate further               |
| The `<skills>` list: name, description, path                                     | every `SKILL.md` body, which stays lazy-loaded                           |
| The `.agent.md` body, injected as `<modeInstructions>`                           | the calling conversation, in any form                                    |
| The `.github/instructions/` file list, with `applyTo` and description            | every `.instructions.md` body                                            |
| The `tools:` frontmatter, enforced: `explorer` truly has no edit or execute tool | `tool_search`, so any deferred tool offered to a subagent is unreachable |

A subagent therefore knows the rules and the name of every skill, but has read
nothing. Name the skill it must load and the files it must open; never assume
it has already seen them, and never send it after a deferred tool.

## 2. Automatic triggers

Start a subagent without being asked when any of these hold:

| Signal                                                                             | Action                          |
| ---------------------------------------------------------------------------------- | ------------------------------- |
| "where is X", "how does X work", "what calls X", "trace"                           | `explorer`                      |
| Task touches more than 2 files, or the right files are unknown                     | `explorer` first                |
| A code edit is needed and the target path is known                                 | `worker`                        |
| A script, test, build, install, or `bsub`/`sbatch` must run, or its output checked | `runner`                        |
| "check on the job", "is it done", a job finished, DONE but no output               | `monitor`                       |
| Non-trivial edits are complete                                                     | `reviewer`                      |
| An association model was written or changed                                        | `reviewer`, statistical pass    |
| Several independent areas must be understood                                       | multiple `explorer` in parallel |

Do NOT delegate a single-file typo fix, a one-line config change, or a question
already answerable from loaded context.

## 3. Sequence by task type

- **Analysis only**: `explorer`
- **Execute only**: `explorer` -> `runner`
- **Implement only**: `explorer` -> `worker` -> `runner` (syntax gate only:
  test suite if one exists, formatter, `bash -n`, the lint task) -> `reviewer`.
  `worker` cannot run commands, so without this step nobody runs the gate
- **Implement and run**: `explorer` -> `worker` -> `runner` -> `reviewer`
- **Submit to the cluster**: `explorer` -> `worker` -> `runner` (submits, then
  the turn ends) -> `monitor` (next turn, once the job has had time to run)

A cluster task is not finished when `runner` returns. A `bsub` exit code says
the job was accepted, nothing more. Report the job ID and stop; never close the
task out on a submission, and never poll. The next step is `monitor`, in a
later turn.

Summarize each subagent result before assigning the next step. Report the
synthesized outcome to the user, not raw subagent transcripts.

## 4. Prompt contract

A subagent starts with a clean context and sees only the prompt it is given.
Every delegation prompt must state:

1. The concrete goal, in one sentence.
2. Known file paths, commands, or symbols to start from.
3. Hard constraints — what not to touch, what not to run.
4. The exact expected output — file list, execution path, diff summary, exit
   code, or findings by severity.

Never forward the whole conversation. Never ask a subagent to infer intent.

The reply comes back only to the caller; the user never sees a subagent's
output. Anything worth showing has to be restated in your own message.

**Ask for the report in English.** A subagent inherits the language rule from
`copilot-instructions.md` and answers in Chinese by default, but its report is
read by you and rewritten before anyone sees it, so that translation is pure
waste and one more place a number can shift.

Pass along any hard environment rule the subagent could violate — for example,
that scratch files go to the repository's own `tmp/` root as
`"${tmpdir}/<task-name>/"`, never `/tmp` and never inside the repository tree,
falling back to `$HOME/tmp/<task-name>/` only after reporting that the
repository has no `tmp/` root.

## 5. Parallel vs sequential

- Run read-only subagents in parallel when their subtasks are independent.
- Run sequentially when one step's output feeds the next.
- Never run two `worker` subagents on the same file concurrently.

## 6. Failure and retry loop

- `runner` reports the failure once and stops.
- Summarize the failure before reassigning.
- `worker` uses the failure evidence to make the next minimal fix.
- `runner` reruns only after code has changed.
- Do not let `runner` enter unbounded retry loops. After two failed fix
  attempts on the same error, stop and report to the user.

## 7. Model selection

Each subagent's model is pinned in its own `.agent.md` frontmatter:
`coordinator` on `claude-opus-5`, and `explorer`, `worker`, `runner`,
`monitor`, `reviewer` on `claude-sonnet-5`. The `runSubagent` `model` parameter
overrides that for one call; use it only when a specific model is genuinely
required.

**The two fields take different spellings.** Frontmatter takes the model id
(`claude-sonnet-5`). The `runSubagent` `model` parameter takes the display name
`"<Model Name> (<vendor>)"` — `"Claude Sonnet 5 (copilot)"`. Passing the id
there fails outright (verified 2026-09-21), so omit the parameter and let the
frontmatter decide unless you have the exact display name.

The frontmatter is the contract, not a verified fact about routing. Which model
actually served a call cannot be confirmed from the local session store, so
never claim it afterwards and never build logic on an assumed fallback.
