---
name: long-running-jobs
description: "Launch and track commands that run longer than ~60 seconds without blocking or polling. Use when starting a build, install, training run, large R/Python job, whole-genome pass, or any LSF/SLURM submission; when parallelizing independent work across chromosomes, traits, samples, files, or parameter sweeps; when wrapping a command in tmux or nohup; when submitting a multi-stage pipeline as one dependency chain so it runs unattended overnight; when resuming a chained pipeline from the stage that broke; when asked to check progress on a running job; or when a job array reports DONE but produced no output. Covers the no-sleep-polling rule, tmux session patterns, LSF/SLURM array wrappers, the DAG submitter and its dependency expressions, bounded local parallelism, how to check status once via progress file, queue, logs and real outputs, and the two silent array-wrapper bugs."
---

# Long-running jobs

Rules for launching work that outlives a single tool call, and for keeping the
agent from blocking on it.

This skill is portable and names no project value. Which scheduler the site
runs, which queue to use, and the default resource sizes come from the
repository's `.github/instructions/` bindings or `AGENTS.md`.

---

## 1. Never poll

**Never use `sleep` to wait for completion.** After launching a long-running
command, report the session/job ID and an estimated runtime, then stop. No
loops, no polling, no "checking back in 30 seconds". The user or a later
session resumes from the active progress file or the queue.

If a command finishes within a few seconds (e.g. a quick `bsub` submission),
just proceed.

## 2. Decide: inline, tmux, or cluster

This table is the source of truth for the thresholds. Other files point here
rather than redefining them.

| Work                                                                    | Where             |
| ----------------------------------------------------------------------- | ----------------- |
| ≤ ~60 s: quick CLI, syntax check, single-unit smoke run, `ls` / `head`  | inline            |
| > ~60 s, single node: builds, installs, model fits, large R/Python jobs | tmux              |
| Many independent units, or needs cluster resources                      | LSF / SLURM array |

Every launch redirects its output to `logs/<stage>/`: a tmux run to
`logs/<stage>/<session>.log`, a cluster task to
`logs/<stage>/<NN-step>_%J_%I.out` and `.err`. Create the directory before
launching. Never write a run log next to the script or under `data/` or
`results/`. -> `data-result-layout`

If tmux is unavailable, fall back to
`nohup <full-command> > logs/<stage>/<name>.log 2>&1 & echo $!` and report that
PID. Record the PID in the progress file: it is the only handle a later session
has on that process, because `tmux ls` will never show it.

## 3. Single tmux session

One session per launch, named after the job. The trailing `; exit` closes the
session on completion, so a session that still exists means its job is still
running.

```bash
SESSION=<descriptive-kebab-case-name>
STAGE=<NN-stage-name>
mkdir -p "logs/$STAGE"
tmux has-session -t "$SESSION" 2>/dev/null \
  && { echo "session $SESSION is still running; wait or pick another name"; exit 1; }
tmux new-session -d -s "$SESSION"
tmux send-keys -t "$SESSION" "<full-command> > logs/$STAGE/$SESSION.log 2>&1; exit" C-m
```

**Never `has-session ... || new-session ...` and then `send-keys` anyway.** If
the session already exists and its pane has a foreground process, `send-keys`
types the command into that process's stdin instead of running it. The launch
looks successful and nothing starts. Refuse to reuse a live session instead.

After launching, always print:

```bash
tmux attach -t $SESSION
tail -f logs/$STAGE/$SESSION.log
```

The trailing `; exit` also means `tmux ls` shows only jobs still running, so
note that the session disappears when the job finishes.

## 4. Parallelize the right dimension

Run independent units concurrently instead of in a serial loop. Valid parallel
dimensions: per-chromosome shards, per-trait runs, per-sample or per-batch
preprocessing, per-file jobs, per-parameter sweeps.

Do NOT parallelize dependent steps — QC -> PCA -> GWAS must stay ordered.
Parallelize only _within_ one independent step.

Bounded local parallelism. Pass each unit as an **argument**; never interpolate
it into the `sh -c` string:

```bash
printf '%s\n' unit1 unit2 unit3 \
  | xargs -P 4 -I{} sh -c '<command> --input "$1" > "logs/<stage>/$1.log" 2>&1' _ {}
```

`sh -c '<command> --input {} ...'` is command injection: a unit named
`u3;touch X` runs `touch X`. The trailing `_ {}` passes the unit as `$1`, where
it stays a literal string.

`-P N`, or `parallel -j N` (GNU parallel, provided by pixi), bounds concurrency
so the node's CPU and RAM are not exhausted. Put the pipeline in the stage's
`.sh` and launch that script with the section 3 tmux pattern. Do not nest this
one-liner inside `send-keys`, where the quoting becomes unmaintainable.

## 5. Cluster array wrapper

One array over the independent dimension, one log per task:

```bash
#BSUB -J myjob[1-N]
#BSUB -o logs/<stage>/myjob_%J_%I.out
#BSUB -e logs/<stage>/myjob_%J_%I.err
#BSUB -n 4
#BSUB -R "rusage[mem=16000]"
#BSUB -W 8:00
#BSUB -q <queue>

TASKS=(unit1 unit2 ... unitN)
unit=${TASKS[$((LSB_JOBINDEX-1))]}

<command> --input "$unit" --out data/intermediate/NN-stage/$unit.result ; status=$?
echo "[done] $unit (exit $status)"
exit $status
```

The resource directives are not optional: `-n` cores, `-R "rusage[mem=]"` MB per
core, `-W` wall time, `-q` queue. Size `-W` above the worst expected task, since
LSF kills the task at the limit. `logs/<stage>/` must already exist at submit
time; `config.sh` creates it.

This skill does not choose those numbers. Where the site has a
`lsf-resource-planner` skill, load it to pick `-n`, `mem`, `-q`, and the `%K`
throttle before writing the header; otherwise take them from the repository
bindings.

SLURM equivalent: `#SBATCH --array=1-N` and `${SLURM_ARRAY_TASK_ID}`.

Anchor any repo path inside the wrapper on `${LS_SUBCWD:-$PWD}` (LSF) or
`${SLURM_SUBMIT_DIR:-$PWD}` (SLURM), never on `${BASH_SOURCE[0]}`, and submit
from the repository root. Full wrapper contract:
`analysis-pipeline/references/script-templates.md` sections B and C.

**Test before full submission**: submit `[1-1]` first, verify its output file
exists and is non-empty, then submit the full range. This prevents burning
cluster hours on N tasks that all fail on the same config or column-name bug.

## 6. Two silent array bugs

Either one makes an entire array report DONE while producing nothing:

- **Bash-reserved array name.** `GROUPS` is a readonly special variable holding
  the caller's Unix group IDs, so `GROUPS=(a b c)` is silently ignored and
  `${GROUPS[i]}` returns a numeric GID. Avoid `PIPESTATUS`, `BASH_SOURCE`,
  `FUNCNAME`, `UID`, `EUID`, `PPID` too. Use `TASKS`, `UNITS`, or `NEWGROUPS`.
- **Swallowed exit code.** Any command after the worker (even a trailing
  `echo`) makes the task exit 0, so the scheduler records DONE on failure.
  Capture `status=$?` immediately and `exit $status`.

Because a masked failure passes the DONE check, always verify real outputs
after an array finishes (`ls <pattern> | wc -l`, row counts), never the DONE
count.

## 7. Chain the whole pipeline in one submission

Submitting stage N, waiting for it, then submitting stage N+1 wastes every hour
between "stage N finished" and "somebody noticed it finished". Submit the
**whole** graph at once and let the scheduler hold each stage PEND until its
prerequisites succeed. The pipeline then runs unattended overnight.

Do this only after the one-unit smoke rung in section 5 has passed for every
stage. A chain multiplies the cost of a config bug by the number of stages.

### The dependency expression

Hold the pipeline as a **DAG**, not a list. Stages fan out and fan in, and a
stage chained behind something it does not read sits on the critical path for
nothing.

| Parent               | Expression                | Means                                       |
| -------------------- | ------------------------- | ------------------------------------------- |
| a job array          | `numdone(<jid>,*)`        | every element DONE, i.e. every one exited 0 |
| a single job         | `done(<jid>)`             | that job DONE                               |
| two parents (fan-in) | `numdone(A,*) && done(B)` | both                                        |

**Never `done(<jid>[*])`.** In LSF that is the element-wise form: element `i` of
the child waits for element `i` of the parent, and it requires the two arrays to
be the same size. It does not mean "all elements succeeded", and it silently
does the wrong thing when the sizes differ. It is the most common mistake in
LSF chaining advice.

SLURM equivalent: `--dependency=afterok:<jobid>`, which for an array already
means every task; `aftercorr` is the element-wise form.

### The helper, and where the DAG lives

**The Mermaid block in the campaign's `PROGRESS.md` is the DAG.** Not a picture
of it — the actual submission input.
[scripts/submit_lsf_pipeline.sh](./scripts/submit_lsf_pipeline.sh) parses it,
submits every stage in topological order, and prints the same block back with
job IDs and live classes for you to paste in.

That is deliberate. A separate spec file is a fourth thing that drifts from the
code, the plan, and the progress file. Here the graph the user reads at 7 a.m.
is provably the graph that ran, because it would not have submitted otherwise.
It is also written and approved during planning, when the DAG is decided.

```bash
# from the repository root
.github/skills/long-running-jobs/scripts/submit_lsf_pipeline.sh \
  --dry-run --graph <track>/<date>-<title>.PROGRESS.md
```

Graph contract — the first ` ```mermaid ` fence in the file:

| Element                                                  | Meaning                                           |
| -------------------------------------------------------- | ------------------------------------------------- |
| node id                                                  | the stage id, `[A-Za-z0-9_]` only                 |
| a `<br/>` field holding a path ending `.lsf` / `.sbatch` | makes the node a **stage**                        |
| any other node                                           | display-only artifact box; dropped with its edges |
| `ext_<jobid>`                                            | a job already in the queue                        |
| `parent --> child`                                       | a dependency, with or without an edge label       |
| `:::class`                                               | display state only; never read as input           |

```
  matrix["matrix<br/><track>/run_matrix.lsf [1-10]<br/>323076536 run 0/10"]:::run
```

So an artifact node like `OUT["release/v2<br/>10 matrix stems"]` stays in the
picture and is simply not submitted. The block round-trips: feeding the printed
block back in reproduces the same plan.

Two fallbacks, neither of which should be maintained by hand:

- positional `.lsf` files — a linear chain, for a quick two- or three-step run.
- `--spec <file.tsv>` — three columns `id script deps`, for a pipeline that has
  no `PROGRESS.md` yet or for testing.

What the helper guarantees, and why each one matters:

- Every `.lsf` is read, its `#BSUB -J` parsed, and its `-o` / `-e` log directory
  checked **before the first `bsub`**. A missing log directory fails the job at
  dispatch, by which time the whole chain is already queued behind it.
- Array-ness comes from the parsed `-J`, not from the array spec in the node
  label, so a stale label cannot swap `numdone` for `done`.
- The topological order is computed and a cycle is refused.
- Submission stops at the first `bsub` failure and prints what was already
  submitted, so nothing is ever queued behind an unparsed job ID.
- `--dry-run` prints the plan and submits nothing.
- It writes nothing back into the markdown and never runs `bkill`. Updating
  `PROGRESS.md` is a normal edit; killing queued work needs a human.

**Always `--dry-run` first** and read the `waits for` lines. That is the only
cheap check that the graph is the graph you meant.

### `done()` fails closed, not loudly

If any element of a parent array EXITs, `numdone(parent,*)` can never become
true, so the child sits PEND **forever** instead of failing. The scheduler does
not report this as an error and nothing cleans it up. A chain left overnight
after a mid-pipeline failure therefore shows PEND, not EXIT, and the queue looks
healthy. This is why the state file and the progress graph exist: the queue
alone does not tell you the chain is dead.

### Resume from the breakpoint, never from the start

When a stage fails mid-chain, fix it and restart **at the break**, leaving the
already-verified stages untouched.

1. **Find the break.** `--status <state.tsv>` re-reads the graph, queries the
   scheduler, and marks the failed stage plus every downstream job now stranded
   PEND behind it.
2. **Read the failures** — the first two or three `.err` files, not all of them.
3. **Kill the stranded jobs.** They can never run. `--status` prints the exact
   `bkill` line; review it and run it yourself. The helper never kills anything.
4. **Fix** the script with the minimal change.
5. **Rerun only the failed indices**, not the whole array:
   `bsub -J "meth-matrix[3,7]" < <track>/run_matrix.lsf`. A command-line `-J`
   overrides the one in the file.
6. **Requeue the tail** behind that rerun:

   ```bash
   submit_lsf_pipeline.sh --graph <PROGRESS.md> \
     --after <rerun-jobid> --from <next-stage>
   ```

   `--from` keeps the named stage and everything downstream and drops the rest;
   `--skip` drops named stages and rewires their children onto their
   dependencies.

Depending on the rerun array alone is correct: the elements that already
finished are DONE and verified, and the new array covers exactly the ones that
were not.

**`--after` only works inside the scheduler's record-retention window.** Once a
finished job's record is purged, its ID cannot be depended on. The helper
refuses rather than queueing a dependency that can never resolve; verify the
outputs on disk and start a fresh root instead.

### Record it in the progress file

The graph is already in `PROGRESS.md`; that is what was submitted. After each
submission and each `--status` check, paste the printed block back over the old
one so the classes are current. The six state classes and the update rule are
owned by `analysis-pipeline/references/plan-and-progress.md` section D.1.

The auto-written state file under `logs/lsf-pipeline/` is the machine record
that `--status` reads. Nobody maintains it; record its path in the progress
file's job table so the next session can find it.

## 8. Check status when asked

"Check progress" is a single read, not a poll. Do the passes below once, report,
update the progress file, and stop.

1. **Active progress file** — the campaign's `PROGRESS.md` names the submitted
   job IDs and the expected output counts. Read it first; everything else is
   checked against it.
2. **Queue** — use the scheduler the repository bindings name. On LSF:
   `bjobs -A` for array summaries, `bjobs -l <jobid>` for one job including
   elapsed time, `bqueues` for queue state. On SLURM: `squeue -u $USER` and
   `sacct -j <jobid>`. For a chained pipeline, `submit_lsf_pipeline.sh --status
<state.tsv>` does this pass for every stage at once and flags the ones
   stranded behind a failure (section 7).
3. **tmux** — `tmux ls` lists only sessions whose job is still running
   (section 3); `tail -n 50 logs/<stage>/<session>.log` shows where it is. For
   a `nohup` fallback launch there is no session, so check the PID the progress
   file recorded with `ps -p <PID>` instead; an empty `tmux ls` says nothing
   about it.
4. **Logs** — for a failing array read the first few failures, not all of them:
   `logs/<stage>/<step>_<jobid>_<idx>.err`.
5. **Real outputs** — `ls <pattern> | wc -l` against the expected count. This is
   the check that decides, not the DONE count (section 6).

## 9. Estimate the runtime

Give a rough range when reporting a launch (e.g. "~30-40 min/task",
"~2-3 h total"), based on input size, prior runs of a similar step
(`bjobs -l <jobid>` elapsed time), and the wall-time header as an upper bound.
Precision is not needed; it just tells the user whether to wait.
