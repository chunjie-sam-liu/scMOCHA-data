---
name: lsf-resource-planner
description: 'Size and place an LSF job on this site before writing the .lsf file. Use when creating or editing any .lsf wrapper or bsub command; when choosing -n, -R "rusage[mem=]", -M, -W, or -q; when a job needs a lot of memory or many cores on one host; when a job was killed for exceeding memory; when an array pends too long or would be too wide; when deciding whether to split work into many narrow tasks or bundle it into one wide task; when partitioning large-scale work into array tasks or chunks; or when asked which queue or node dispatches fastest. Covers the wide-and-thin vs narrow-and-fat shape decision, the per-slot memory trap, the queue decision table, the dispatch-speed probe, per-user and per-queue caps, array throttling, and chunking. Site-bound: carries a dated inventory of this cluster''s nodes, queues, and limits.'
---

# LSF resource planner

How to choose `-n`, `-R "rusage[mem=]"`, `-W`, and `-q` on this site, and how to
place an array so it dispatches instead of sitting in the queue.

This skill is **site-bound**, not portable: it carries a measured inventory of
this cluster. `long-running-jobs` owns the launch-and-never-poll rules and the
array wrapper contract; `analysis-pipeline` owns the script template. This skill
only decides the numbers that go into those headers.

Full node table, queue table, raw config evidence, and refresh commands:
[references/cluster-inventory.md](./references/cluster-inventory.md).

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

Two of these are easy to get wrong.

**Reservation and hard limit are different things, set by different mechanisms.**
`RESOURCE_RESERVE_PER_TASK=Y` is stock LSF and only decides how much memory the
scheduler holds on a host while placing the job. The value that kills the job is
`LSB_SUB_RLIMIT_RSS`, which the site esub computes as
`max(rusage_mem, -M, app_profile_limit) * ncores`. Stock `LSB_JOB_MEMLIMIT=y`
would **not** by itself turn `rusage[mem]` into a limit: IBM ties that parameter
to `-M` / queue `MEMLIMIT`, and no queue here sets `MEMLIMIT`. The
`rusage` -> limit conversion is this site's policy script, a local convention
rather than portable LSF behaviour. Code path and evidence:
[references/cluster-inventory.md](./references/cluster-inventory.md) section 3.

Practical consequences of the esub formula:

- Requesting no memory at all gets `default_mem() = 2499` MB per slot plus a
  stderr warning. Never rely on it.
- `-M` is not ignored, it is a **floor**: an `-M` larger than `rusage[mem]`
  raises the kill threshold while leaving the reservation small. Use it only
  when the peak is genuinely uncertain, and know the node is then
  oversubscribed relative to what was reserved.
- The `* ncores` multiplication happens in the esub. Changing `-n` moves the
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
#BSUB -R "rusage[mem=64000]"     # 8 x 64000 MB = 512 GB reserved and enforced
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

### Width costs concurrency, superlinearly

A task needs `n` free cores **on one host**. Measured on `rhel8_cpu`,
2026-09-03:

| `-n` | hosts with that many cores free | tasks that could start now |
| ---- | ------------------------------- | -------------------------- |
| 1    | 47                              | **1049**                   |
| 2    | 41                              | 518                        |
| 4    | 34                              | 250                        |
| 8    | 26                              | 118                        |
| 16   | 18                              | 52                         |
| 32   | 14                              | 21                         |
| 64   | 7                               | **7**                      |

Every doubling of `-n` roughly halves how many tasks can start. `-n 1` places
**150x** more tasks than `-n 64`. The absolute numbers move hourly; the shape of
the curve does not. Refresh with the command in the inventory.

**So `-n` is the throughput dial for a wide array, not `mem`.** Drop `-n` to
what the code truly uses. Two `-n 1` tasks that start now beat one `-n 2` task
that never does.

### Split or bundle

One question decides between the two shapes: **is the memory shared?**

| Situation                                           | Do this                                                    |
| --------------------------------------------------- | ---------------------------------------------------------- |
| Each unit loads its own data, nothing shared in RAM | **Split.** N independent tasks at `-n 1` or `-n 2`         |
| All threads read one large in-memory object         | **Bundle.** One task, `-n` = thread count, `span[hosts=1]` |
| Units share a read-only file on disk                | **Split.** The filesystem is the shared layer, not RAM     |
| `N` is huge and a unit runs under ~2 min            | **Split, then chunk.** Narrow tasks, several units each    |

The common mistake is bundling because the _total_ is large. 8 units x 64 GB is
not a 512 GB job unless all 512 GB must be resident at once. Split gives 8 tasks
that each need one host with 64 GB free, which is nearly every host. Bundled, it
is one task needing 8 free cores **and** 512 GB on a single host, which right
now is 26 hosts instead of 47.

### Recipe per shape

**Wide and thin** -> get more tasks running:

1. `-n 1`, or the true thread count if higher. Never pad.
2. Omit `span[hosts=1]`; it is meaningless at `-n 1` and only adds a constraint.
3. Keep `rusage[mem]` tight but do not agonize: memory almost never blocks on
   rome (section 6).
4. Do **not** throttle on `standard`. It has no `JL/U`, and pending jobs cost no
   fairshare -- only _running_ slots lower dynamic priority. Throttling here only
   caps your own throughput.
5. If a unit runs under ~2 min, chunk (section 7). Scheduling overhead otherwise
   exceeds the work.
6. Consider splitting the range across queues (section 7).

**Narrow and fat** -> get the one task placed:

1. Match the pool's RAM-per-core ratio so the node is not left half-idle:
   rome ~15.7 GB/core, nodelmr up to ~62, nodelcr ~23, nodesd ~54.
2. `span[hosts=1]` is mandatory for anything shared-memory.
3. Move to `large_mem` or `large_core_count` (section 5). Both are near-idle,
   and a wide `-n` there competes with far fewer jobs.
4. Target a node class with `select[]` when the memory needs it, scoped by the
   queue's own host group:

   ```bash
   #BSUB -q standard
   #BSUB -R "select[maxmem>1900000]"   # the 36 fat rome nodes, 1.9 TB each
   ```

   `maxmem` is the host's total RAM; `mem` is what is free right now. Prefer
   `maxmem` -- a `mem` predicate re-evaluates and can leave the job pending on a
   transient. Never pin with `-m <host>`.

5. Let reservation work. `standard` accumulates slots for up to 7200 s
   (`MAX_RESERVE_TIME`), `large_mem` and `large_core_count` for 1800 s, so a wide
   task does not starve. `bjobs -l` shows `Reserved <N> job slots` while it
   gathers.

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
that. It is empty precisely because people respect the boundary, and the pool is
only 19 nodes.

**Worked example, the 500 GB case.**

| Option             | Directive                                  | Verdict                                                                        |
| ------------------ | ------------------------------------------ | ------------------------------------------------------------------------------ |
| Naive              | `-q standard -n 32 -R "rusage[mem=16000]"` | 512 GB, but 32 slots behind ~114000 pending jobs                               |
| Better on standard | `-q standard -n 4 -R "rusage[mem=128000]"` | 512 GB on 4 slots; dispatches far sooner, memory is not the constraint on rome |
| Correct            | `-q large_mem -n 8 -R "rusage[mem=64000]"` | 512 GB, matches the queue's purpose, pending is 0                              |

Take the third when the work genuinely needs 500 GB. Take the second when the
real need is 300-450 GB and `large_mem` would be an abuse.

**Two queues that are not for arrays.** `large_core_count` has 4 nodes and
`superdome` has 3. A wide array on either will starve everyone including itself.
Send single large jobs there; send arrays to `standard`, `priority`, or
`large_mem`.

### `short` is a CPU-time queue, not a 30-minute wall-clock queue

`short` sets `CPULIMIT = 30` and **no `RUNLIMIT`**. Its own `DESCRIPTION` string
says "hard run time limit of 30 minutes" and that description is wrong. Three
consequences:

1. **It limits CPU time, not wall clock.** An I/O-bound job that sleeps on the
   filesystem for hours accrues almost no CPU time and is never killed.
2. **The 30 is normalised to a reference host.** `bqueues -l short` reports
   `30.0 min of svlprhpc02`, and `svlprhpc02` has `cpuf 60.0` while every compute
   node has `cpuf 40.0`. The real allowance on a compute node is
   `30 x 60/40 = 45` CPU-minutes.
3. **CPU time sums across threads.** `LSF_HPC_EXTENSIONS="CUMULATIVE_RUSAGE"`
   accumulates usage over all the job's processes, so a fully loaded multi-core
   job burns the allowance `n` times faster.

Safe wall-clock budget on `short`, assuming full CPU utilisation:

| `-n` | CPU-minutes allowed | wall-clock before kill |
| ---- | ------------------- | ---------------------- |
| 1    | ~45                 | ~45 min                |
| 2    | ~45                 | ~22 min                |
| 4    | ~45                 | ~11 min                |
| 8    | ~45                 | ~6 min                 |

So `short` is only useful for `-n 1` or `-n 2` work. Still set `-W` on it: `-W`
is the only wall-clock guard anywhere in the CPU queues, and a job that stalls
on I/O will otherwise sit forever without ever tripping `CPULIMIT`.

**No CPU queue sets `RUNLIMIT`.** Verified for `standard`, `priority`, `short`,
`large_mem`, `large_core_count`, `superdome`, `heavy_io`, `interactive`, and
`compbio`. The six queues that do set one are all cryoem or GPU queues. A job
submitted without `-W` therefore has no wall-clock limit at all, which is why
`-W` is not optional.

## 6. Probe dispatch speed before submitting

Run this once before writing the `-q` line. It is a read, not a poll, so it does
not violate the no-polling rule in `long-running-jobs`.

```bash
printf '%-18s %7s %7s %6s %9s %8s\n' QUEUE PEND RUN P/R FREESLOT FREEHOST
for q in standard priority short large_mem large_core_count heavy_io interactive; do
  read -r pend run <<<"$(bqueues -w "$q" 2>/dev/null | awk '
    NR==1{for(i=1;i<=NF;i++){if($i=="PEND")p=i; if($i=="RUN")r=i}}
    NR==2{print $p, $r}')"
  hg=$(bqueues -l "$q" 2>/dev/null | awk '/^HOSTS:/{print $2}' | tr -d '/')
  read -r fs fh <<<"$(bhosts -w "$hg" 2>/dev/null | awk 'NR>1 && $4!="-"{f=$4-$6; s+=f; if(f>0 && $2=="ok") h++} END{print s+0, h+0}')"
  printf '%-18s %7s %7s %6s %9s %8s\n' "$q" "$pend" "$run" \
    "$(awk -v p="$pend" -v r="$run" 'BEGIN{printf "%.1f",(r>0?p/r:p)}')" "$fs" "$fh"
done
```

Read it as: **`P/R` is the wait, `FREESLOT` is the headroom.** A queue with
`P/R` near 0 and non-zero `FREESLOT` dispatches now. `standard` at `P/R` above
10 means hours to days for a wide array.

**On `standard`, slots are the bottleneck and memory is not.** Of 220 rome/cn
nodes, 210 have at least 800 GB free while only about 1900 of 14197 slots are
free. LSF states this directly in `bjobs -p`:

```
Affinity resource requirement cannot be met because there are not enough
processor units to satisfy the job affinity request: 30 hosts;
Job's requirements for resource reservation not satisfied (Resource: mem): 4 hosts;
```

30 hosts blocked on cores, 4 on memory. So on `standard`, **cutting `-n` buys
far more dispatch speed than cutting memory.** Reach for a smaller `n` with a
larger `mem` per slot before reaching for a different queue.

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
for pending. Submit the full range and let LSF place whatever fits. A `%K` here
caps only your own throughput.

```bash
#BSUB -J stage-step[1-3000]      # -n 1, standard, no %K
```

The ceiling is `min(4000 / n, 4000)`: the per-user running-slot cap and the
array index cap. At `-n 1` both are 4000.

**Split the range across queues when the backlog is long.** `priority`
dispatches roughly 26x sooner than `standard` per job but caps at 400 jobs per
user; `short` takes 300 more if every task finishes inside 30 minutes. Same
script, three submissions, three job IDs to record in `PROGRESS.md`:

```bash
bsub -q priority < step.lsf      # edit -J to [1-400],    fills within minutes
bsub -q short    < step.lsf      # edit -J to [401-700],  tasks under 30 min only
bsub -q standard < step.lsf      # edit -J to [701-3000], absorbs the remainder
```

Only worth the extra bookkeeping when `standard` shows `P/R` above ~10 and the
array is large. Never split a _dependent_ range this way.

**Narrow and fat, or a small pool: throttle.** `%K` is the suffix on the array
spec:

```
K = min( JL_U(queue), floor(4000 / n), floor(0.25 * FREESLOT(queue) / n) )
```

The third term keeps the job from taking more than a quarter of a pool's free
slots. Round down to something readable.

| Queue              | `n` | `JL/U` | `4000/n` | `0.25*free/n` | Use                             |
| ------------------ | --- | ------ | -------- | ------------- | ------------------------------- |
| `standard`         | 1   | none   | 4000     | ~260          | none, submit the full range     |
| `standard`         | 2   | none   | 2000     | ~240          | none, or `%200` if being polite |
| `standard`         | 8   | none   | 500      | ~60           | `%50`                           |
| `large_mem`        | 8   | 300    | 500      | ~22           | `%20`                           |
| `large_core_count` | 64  | 300    | 62       | ~1            | not an array queue              |

**When the unit count `N` exceeds 4000**, do not submit several arrays back to
back. Chunk instead: keep the array narrow and give each task a contiguous slice.

```bash
#BSUB -J stage-step[1-400]%200

CHUNK=25                       # 400 tasks x 25 units = 10000 units
TASKS=( ... )                  # never name this GROUPS; see long-running-jobs
start=$(( (LSB_JOBINDEX - 1) * CHUNK ))
fail=0
for (( i = start; i < start + CHUNK && i < ${#TASKS[@]}; i++ )); do
  unit="${TASKS[$i]}"
  <command> --input "$unit" ; status=$?
  if (( status != 0 )); then echo "[fail] $unit (exit $status)"; fail=1; fi
done
exit $fail
```

Chunking also amortizes per-task startup, which matters when a task is shorter
than about two minutes. Size `-W` for the **chunk**, not for one unit.

## 8. Write the script

Use the template in `analysis-pipeline/references/script-templates.md` section B
unchanged, and fill the header from steps 2 to 7. Only the resource lines are
this skill's business.

**Wide and thin** -- 3000 units, each single-threaded and about 6 GB:

```bash
#BSUB -J <stage>-<step>[1-3000]
#BSUB -o logs/<stage>/<step>_%J_%I.out
#BSUB -e logs/<stage>/<step>_%J_%I.err
#BSUB -n 1
#BSUB -R "rusage[mem=8000]"      # 1 x 8000 MB = 8 GB per task
#BSUB -W 2:00
#BSUB -q standard
```

No `%K`, no `span[hosts=1]`. At `-n 1` about 1000 tasks can start immediately
and LSF backfills the rest as slots free up.

**Narrow and fat** -- 2 units, each needing 8 threads over one 512 GB object:

```bash
#BSUB -J <stage>-<step>[1-2]%2
#BSUB -o logs/<stage>/<step>_%J_%I.out
#BSUB -e logs/<stage>/<step>_%J_%I.err
#BSUB -n 8
#BSUB -R "rusage[mem=64000]"     # 8 x 64000 MB = 512 GB total
#BSUB -R "span[hosts=1]"
#BSUB -W 48:00
#BSUB -q large_mem
```

- `span[hosts=1]` is required for anything shared-memory. Without it LSF may
  spread the slots across hosts and the threads cannot see each other.
- `-W` above the worst expected task, since the task is killed at the limit, but
  not absurdly above it: a tight `-W` helps backfill scheduling place the job.
- The three array traps (`BASH_SOURCE`, reserved array names, swallowed exit
  code) are in `long-running-jobs` section 6. They still apply.

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

## 10. Diagnose

| Symptom                                     | Command                                   | What it means                                                                                                                          |
| ------------------------------------------- | ----------------------------------------- | -------------------------------------------------------------------------------------------------------------------------------------- |
| Pending for a long time                     | `bjobs -p <jobid>`                        | reports hosts blocked on cores vs on memory; act on whichever dominates                                                                |
| Array running far fewer tasks than expected | recheck `-n` against the section 4 table  | width, not memory, is usually the cap                                                                                                  |
| Array status                                | `bjobs -A <jobid>`                        | PEND / RUN / DONE / EXIT per array                                                                                                     |
| Killed unexpectedly                         | `bhist -l <jobid> \| grep TERM_`          | `TERM_MEMLIMIT` = over the esub's `n` x `mem` threshold; `TERM_RUNLIMIT` = over `-W`; `TERM_CPULIMIT` = over `short`'s CPU-time budget |
| Reserved vs used                            | `bjobs -o "... memlimit max_mem" <jobid>` | the retune input from step 9                                                                                                           |
| Queue changed                               | rerun the step 6 probe                    | queue pressure moves hourly                                                                                                            |

A wide array reporting DONE is not proof. Count real outputs, per
`long-running-jobs` section 6.
