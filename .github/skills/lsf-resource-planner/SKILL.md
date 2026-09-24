---
name: lsf-resource-planner
description: 'Size and place an LSF job on this site before writing the .lsf file. Use when creating or editing any .lsf wrapper or bsub command; when choosing -n, -R "rusage[mem=]", -M, -W, or -q; when a job needs a lot of memory or many cores on one host; when a job was killed for exceeding memory; when an array pends too long or would be too wide; when deciding between many narrow tasks and one wide task; when partitioning work into array tasks or chunks; when writing a -w dependency expression or chaining stages into one submission; or when asked which queue or node dispatches fastest. Covers the wide-and-thin vs narrow-and-fat shape decision, the per-slot memory trap, the queue decision table, array throttling, and chunking. Site-bound: carries a dated inventory of this cluster''s nodes, queues, and limits.'
---

# LSF resource planner

How to choose `-n`, `-R "rusage[mem=]"`, `-W`, and `-q` on this site, and how to
place an array so it dispatches instead of sitting in the queue.

This skill is **site-bound**, not portable: it carries a measured inventory of
this cluster. `long-running-jobs` owns the launch-and-never-poll rules and the
array wrapper contract; `analysis-pipeline` owns the script template. This skill
only decides the numbers that go into those headers.

The decision path is below. The evidence, recipes, templates, and symptom tables
live in the references, so load only the one the task needs:

| Reference                                                           | Holds                                                                          |
| ------------------------------------------------------------------- | ------------------------------------------------------------------------------ |
| [cluster-inventory.md](./references/cluster-inventory.md)           | full node and queue tables, raw config evidence, refresh commands              |
| [shape-and-queue.md](./references/shape-and-queue.md)               | the concurrency measurements, both shape recipes, the 500 GB and `short` cases |
| [array-and-headers.md](./references/array-and-headers.md)           | the dispatch probe, throttle and split recipes, chunking, filled headers       |
| [chaining-and-diagnosis.md](./references/chaining-and-diagnosis.md) | `-w` dependency semantics and the full diagnosis table                         |

---

## 1. Eight site facts that change how bsub is written

Verified on `hpcf_research_cluster`, IBM Spectrum LSF 10.1.0.14, 2026-09-04.

| Fact                                        | Setting                                                                       | Consequence                                                                                                                                |
| ------------------------------------------- | ----------------------------------------------------------------------------- | ------------------------------------------------------------------------------------------------------------------------------------------ |
| Memory unit is MB                           | `LSF_UNIT_FOR_LIMITS=MB`                                                      | `rusage[mem=16000]` means 16000 MB                                                                                                         |
| **Reservation** is per slot                 | `RESOURCE_RESERVE_PER_TASK=Y`                                                 | scheduling reserves `n` x `mem`. Governs _placement only_                                                                                  |
| **The hard limit comes from a site esub**   | `LSB_ESUB_METHOD="hpcf"`                                                      | `esub.hpcf` rewrites every submission: `RLIMIT_RSS = max(rusage_mem, -M, app_limit) * n`. This, not `rusage` itself, is what kills the job |
| The limit kills, job-wide                   | `LSB_MEMLIMIT_ENFORCE=y`, `LSB_JOB_MEMLIMIT=y`                                | over the limit is a kill, not a suspend, and usage is summed across **all** the job's processes                                            |
| Cores are pinned                            | `LSB_RESOURCE_ENFORCE="cpu gpu"`, `CPUSET_CPU_EXCLUSIVE`, `affinity[core(1)]` | 1 slot = 1 exclusive physical core; threads beyond `-n` contend on the same cores and make the job slower                                  |
| `-W` is the only wall-clock limit           | `ABS_RUNLIMIT=Y`; no CPU queue sets `RUNLIMIT`                                | `-W 8:00` is 8 wall-clock hours, not CPU-factor scaled. Omit `-W` and the job has no wall-clock limit at all                               |
| **`short` limits CPU time, not wall clock** | `CPULIMIT = 30` on `short`, no `RUNLIMIT`                                     | see section 5. A fully loaded `-n 8` job dies after ~6 wall-clock minutes                                                                  |
| Arrays cap at 4000                          | `MAX_JOB_ARRAY_SIZE=4000`                                                     | an index range wider than 4000 is rejected at submit time                                                                                  |

**Reservation and the hard limit are different things.**
`RESOURCE_RESERVE_PER_TASK=Y` only decides how much memory the scheduler holds
on a host while placing the job. What kills the job is `LSB_SUB_RLIMIT_RSS`,
which the site esub computes as `max(rusage_mem, -M, app_profile_limit) *
ncores`. That conversion is this site's policy script, not portable LSF
behaviour; code path and evidence in
[references/cluster-inventory.md](./references/cluster-inventory.md) section 3.
Three consequences:

- Requesting no memory at all gets `default_mem() = 2499` MB per slot plus a
  stderr warning. Never rely on it.
- `-M` is not ignored, it is a **floor**: an `-M` larger than `rusage[mem]`
  raises the kill threshold while leaving the reservation small. Use it only
  when the peak is genuinely uncertain, and know the node is then
  oversubscribed relative to what was reserved.
- The `* ncores` multiplication happens in the esub, so changing `-n` moves the
  kill threshold even when `rusage[mem]` is untouched.

## 2. Estimate the footprint

Two numbers, in this order.

**Cores `n`** is the number of threads the code will actually run, not a wish.
Because cores are pinned, over-requesting wastes slots and makes the job pend
longer, and under-requesting throttles the process to fewer cores than it spawns.

- Single-threaded R or Python: `-n 1`, or `-n 2` when a large `data.table` or GC
  pass benefits from one spare core.
- Threaded BLAS, `data.table`, `mclapply`, `future`: `-n` = the thread count, and
  set it from LSF rather than hard-coding it:

  ```bash
  export OMP_NUM_THREADS=${LSB_DJOB_NUMPROC:-1}
  export MKL_NUM_THREADS=${LSB_DJOB_NUMPROC:-1}
  export OPENBLAS_NUM_THREADS=${LSB_DJOB_NUMPROC:-1}
  ```

**Total memory `M`** comes from evidence, in descending order of trust:

1. `bjobs -o "jobid jobindex max_mem memlimit" -d <jobid>` on a previous run of
   the same step. This is the only real answer.
2. A one-unit smoke run at a deliberately generous limit, then read `max_mem`.
3. An estimate from input size: for R, peak is usually 3-6x the size of the
   largest object held in memory, more if the code copies rather than mutates.

Add **30 percent headroom** over the observed peak, never less. The limit is
enforced by kill, so an under-estimate costs a whole run.

## 3. Convert to per-slot memory

```
mem_per_slot_MB = ceil( M_GB * 1024 / n )
```

Write that into `-R "rusage[mem=<mem_per_slot_MB>]"`. Always state the total in
a comment on the same line so the next reader does not have to redo the math:

```bash
#BSUB -n 8
#BSUB -R "rusage[mem=64000]"     # 8 x 64000 MB = 512000 MB = 500 GB, reserved and enforced
```

That single value sets both numbers: the scheduler reserves `n x mem`, and the
site esub turns the same figure into the `n x mem` kill threshold (section 1).
Add `-M` only to raise the kill threshold above the reservation when the peak is
genuinely uncertain; an `-M` smaller than `rusage[mem]` has no effect, because
the esub takes the max of the two.

## 4. Classify the shape

The unit count `N`, the cores per task `n`, and the memory per task decide
everything that follows. Two shapes need opposite treatment; get this wrong and
the rest of the plan is wrong.

| Shape              | `N`   | `n` per task      | Memory per task | Goal                       |
| ------------------ | ----- | ----------------- | --------------- | -------------------------- |
| **Wide and thin**  | > 50  | 1-4               | < 32 GB         | **maximize concurrency**   |
| **Mixed**          | 10-50 | 4-8               | 32-100 GB       | concurrency, but throttled |
| **Narrow and fat** | 1-10  | >= 8, or > 100 GB | large           | **place precisely**        |

**Width is the throughput dial, not memory.** A task needs `n` free cores on
one host, and every doubling of `-n` roughly halves how many tasks can start:
measured on `rhel8_cpu`, `-n 1` placed 1049 tasks where `-n 64` placed 7. On
`standard`, `bjobs -p` shows cores blocking ~30 hosts against memory blocking
4, so cutting `-n` buys far more dispatch speed than cutting memory. Measured
table and evidence:
[references/shape-and-queue.md](./references/shape-and-queue.md) section A.

### Split or bundle

One question decides between the two shapes: **is the memory shared?**

| Situation                                           | Do this                                                    |
| --------------------------------------------------- | ---------------------------------------------------------- |
| Each unit loads its own data, nothing shared in RAM | **Split.** N independent tasks at `-n 1` or `-n 2`         |
| All threads read one large in-memory object         | **Bundle.** One task, `-n` = thread count, `span[hosts=1]` |
| Units share a read-only file on disk                | **Split.** The filesystem is the shared layer, not RAM     |
| `N` is huge and a unit runs under ~2 min            | **Split, then chunk.** Narrow tasks, several units each    |

The common mistake is bundling because the _total_ is large. 8 units x 64 GB is
not a 512 GB job unless all 512 GB must be resident at once.

### Recipe per shape

**Wide and thin** -> `-n 1` or the true thread count, never padded; no
`span[hosts=1]`; no `%K` on `standard`; chunk if a unit runs under ~2 min.

**Narrow and fat** -> match the pool's RAM-per-core ratio, `span[hosts=1]` is
mandatory, move to `large_mem` or `large_core_count`, target a node class with
`select[maxmem>...]` rather than `-m <host>`, and let reservation gather slots.

Both recipes in full, with the ratios and the `select[]` caveats:
[references/shape-and-queue.md](./references/shape-and-queue.md) section C.

## 5. Choose the queue

Pick by **total memory on one host** and **cores on one host**. Do not default
to `standard`.

| Condition                                      | Queue              | Node pool                     | Typical pending                              |
| ---------------------------------------------- | ------------------ | ----------------------------- | -------------------------------------------- |
| `M` < 500 GB and `n` <= 64                     | `standard`         | 200 noderome + 2 nodecn       | heavy, 13 pending per running                |
| same, but urgent and under 400 concurrent jobs | `priority`         | same                          | light, 0.5 pending per running               |
| runtime <= 30 min                              | `short`            | same                          | light, but `CPULIMIT`, see the warning below |
| I/O bound, not CPU bound                       | `heavy_io`         | same                          | none                                         |
| **`M` >= 500 GB**, `n` <= 80                   | **`large_mem`**    | 19 nodelmr, 2.8-3.9 TB each   | **none**                                     |
| `n` > 64 on one host, `M` <= 2.8 TB            | `large_core_count` | 4 nodelcr, 128 cores / 2.9 TB | **none**                                     |
| `M` > 3.9 TB or `n` > 128 on one host          | `superdome`        | 3 nodesd, 224 cores / 11.8 TB | worst on the cluster, 17 pending per running |
| interactive shell                              | `interactive`      | 6 noderome                    | none, but only 8 jobs per user               |

The 500 GB boundary is the queue's own stated purpose ("For jobs that needs over
500GB of memory to run"). `large_mem` has `USERS = all` and no minimum-memory
admission check, so nothing stops a small job from landing there. Do not use
that -- it is empty precisely because people respect the boundary, and the pool
is only 19 nodes. `large_core_count` (4 nodes) and `superdome` (3 nodes) take
single large jobs, never arrays.

**`short` limits CPU time, not wall clock**, and that time sums across threads.
The real allowance is ~45 CPU-minutes, so a fully loaded `-n 8` job dies after
~6 wall-clock minutes; `short` is only useful for `-n 1` or `-n 2` work. **No
CPU queue sets `RUNLIMIT`**, so a job submitted without `-W` has no wall-clock
limit at all. `-W` is never optional, including on `short`.

The worked 500 GB comparison, the `short` CPU-time derivation, and the
`RUNLIMIT` verification:
[references/shape-and-queue.md](./references/shape-and-queue.md) sections D-E.

## 6. Probe dispatch speed before submitting

Run the probe once before writing the `-q` line. It is a read, not a poll, so it
does not violate the no-polling rule in `long-running-jobs`.

Read it as **`P/R` is the wait, `FREESLOT` is the headroom.** A queue with `P/R`
near 0 and non-zero `FREESLOT` dispatches now; `standard` above 10 means hours
to days for a wide array. Queue pressure moves hourly, so probe again before
re-submitting.

The probe script:
[references/array-and-headers.md](./references/array-and-headers.md) section A.

## 7. Check the caps, then size the array

**Caps that exist on this site** (verified, see the inventory for evidence):

| Cap                          | Value                                                                                                                                  | Scope                           |
| ---------------------------- | -------------------------------------------------------------------------------------------------------------------------------------- | ------------------------------- |
| Running job slots per user   | **4000**                                                                                                                               | cluster-wide, `busers`          |
| Array index range            | **4000**                                                                                                                               | per array, `MAX_JOB_ARRAY_SIZE` |
| Pending jobs / pending slots | effectively unlimited                                                                                                                  | `MAX_PEND_JOBS` = 2^31-1        |
| Jobs per user per queue      | `priority` 400, `short` 300, `large_mem` 300, `large_core_count` 300, `heavy_io` 500, `interactive` 8, `compbio` 3000; `standard` none | `JL/U` in `bqueues`             |
| CPU or memory quota per user | **none**                                                                                                                               | only GPU has `PER_USER` limits  |

There is no CPU-hour or memory quota. The real limits are the 4000 running
slots, the queue `JL/U`, and fairshare. LSF dynamic priority falls with
**running** slots and accumulated run time, not with pending jobs, so a large
pending backlog costs nothing and LSF self-limits as your tasks start.

**Wide and thin: do not throttle `standard`.** No `JL/U`, no fairshare penalty
for pending. Submit the full range and let LSF place whatever fits; a `%K` here
caps only your own throughput. The ceiling is `min(4000 / n, 4000)`.

```bash
#BSUB -J stage-step[1-3000]      # -n 1, standard, no %K
```

**Narrow and fat, or a small pool: throttle** with the `%K` suffix, where
`K = min( JL_U(queue), floor(4000 / n), floor(0.25 * FREESLOT(queue) / n) )`.
The third term keeps the job from taking more than a quarter of a pool's free
slots.

**When the unit count `N` exceeds 4000**, do not submit several arrays back to
back. Chunk instead: keep the array narrow, give each task a contiguous slice,
and size `-W` for the **chunk** rather than for one unit. Chunking also
amortizes per-task startup, which matters below about two minutes per unit.

The split-across-queues recipe, the per-queue throttle table, and the chunking
template:
[references/array-and-headers.md](./references/array-and-headers.md)
sections B-D.

## 8. Write the script

Use the template in `analysis-pipeline/references/script-templates.md` section B
unchanged, and fill the header from steps 2 to 7. Only the resource lines are
this skill's business.

Filled headers for both shapes, with the `span[hosts=1]` and `-W` caveats:
[references/array-and-headers.md](./references/array-and-headers.md) section E.
The array traps themselves are `long-running-jobs` sections 5 and 6.

## 9. Smoke test, measure, retune

The cheapest rung first, always.

1. Submit `[1-1]` with a generous `mem`. Verify the output file is fresh and
   non-empty. -> `data-verification`
2. Read what it actually used:

   ```bash
   bjobs -o "jobid:10 jobindex:6 stat:6 slots:5 memlimit:9 max_mem:10 run_time:12" -d <jobid>
   ```

3. Set the real value: `mem_per_slot = ceil(max_mem_GB * 1.3 * 1024 / n)`.
4. Submit the full range with the tuned header.

Over-requesting is not free. A real case from this repository: an array ran
`-n 8 -R "rusage[mem=8000]"`, so a 63 GB limit, and peaked at 13.1 GB. It
reserved five times what it used and pended five times longer than it needed to.
If its units were independent, `-n 1 -R "rusage[mem=18000]"` would have run the
same work with roughly eight times the concurrency.

## 10. Chain stages into one submission

`long-running-jobs` section 7 owns the method and the DAG submitter. Four site
facts decide whether a chain works here:

- **Finished job records are purged after 24 h** (`CLEAN_PERIOD = 86400`), so
  `-w` on a job that finished more than a day ago can never resolve. Verify its
  outputs on disk and start a fresh root instead.
- **A name dependency resolves to the newest job of that name**
  (`JOB_DEP_LAST_SUB = 1`) -- a trap after a smoke `[1-1]` under the same `-J`.
  **Depend on job IDs, not names.**
- **The whole-array success condition is `numdone(<jid>,*)`.** `done(<jid>[*])`
  is element-wise pairing between two equal-sized arrays, not a whole-array
  test; `numexit(<jid>, >0)` is the matching "did anything fail".
- **Pending jobs are free.** Dynamic priority falls with running slots, not
  pending ones, so there is no reason to hold stages back and submit them one
  at a time.

The `man bsub` wording, the dependency re-evaluation lag, and the evidence:
[references/chaining-and-diagnosis.md](./references/chaining-and-diagnosis.md)
sections A-B.

## 11. Diagnose

| Symptom                 | Command                          | What it means                                                      |
| ----------------------- | -------------------------------- | ------------------------------------------------------------------ |
| Pending for a long time | `bjobs -p <jobid>`               | hosts blocked on cores vs on memory; act on whichever dominates    |
| Killed unexpectedly     | `bhist -l <jobid> \| grep TERM_` | `TERM_MEMLIMIT` / `TERM_RUNLIMIT` / `TERM_CPULIMIT` name the limit |
| Child PEND, parent done | `bjobs -l <child>`               | one EXIT element makes `numdone(parent,*)` unsatisfiable, forever  |

A wide array reporting DONE is not proof. Count real outputs, per
`long-running-jobs` section 6. Full symptom table:
[references/chaining-and-diagnosis.md](./references/chaining-and-diagnosis.md)
section C.
