# Reference: chaining a whole pipeline in one submission

The full procedure behind `SKILL.md` section 7. Read this before writing a
`--graph` submission, when a chained pipeline breaks mid-way, or when deciding
what the Mermaid graph contract allows.

Site-specific dependency facts (record retention, name-vs-ID resolution,
re-evaluation lag) belong to the scheduler binding and to
`lsf-resource-planner`; this file owns the method.

---

## A. The dependency expression

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

---

## B. The helper, and where the DAG lives

**The Mermaid block in the campaign's `PROGRESS.md` is the DAG.** Not a picture
of it — the actual submission input.
[../scripts/submit_lsf_pipeline.sh](../scripts/submit_lsf_pipeline.sh) parses
it, submits every stage in topological order, and prints the same block back
with job IDs and live classes for you to paste in.

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

---

## C. `done()` fails closed, not loudly

If any element of a parent array EXITs, `numdone(parent,*)` can never become
true, so the child sits PEND **forever** instead of failing. The scheduler does
not report this as an error and nothing cleans it up. A chain left overnight
after a mid-pipeline failure therefore shows PEND, not EXIT, and the queue looks
healthy. This is why the state file and the progress graph exist: the queue
alone does not tell you the chain is dead.

---

## D. Resume from the breakpoint, never from the start

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

---

## E. Record it in the progress file

The graph is already in `PROGRESS.md`; that is what was submitted. After each
submission and each `--status` check, paste the printed block back over the old
one so the classes are current. The six state classes and the update rule are
owned by `analysis-pipeline/references/plan-and-progress.md` section D.1.

The auto-written state file under `logs/lsf-pipeline/` is the machine record
that `--status` reads. Nobody maintains it; record its path in the progress
file's job table so the next session can find it.
