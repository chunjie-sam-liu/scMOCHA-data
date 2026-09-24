# Reference: shape evidence and queue selection

The measured evidence behind the shape decision in `SKILL.md` section 4 and the
queue table in section 5. Read this when the shape is not obvious, when a queue
choice is being defended, or before using `short`.

Verified on `hpcf_research_cluster`, IBM Spectrum LSF 10.1.0.14. The absolute
numbers move hourly; the shape of every curve below does not. Refresh commands:
[cluster-inventory.md](./cluster-inventory.md).

---

## A. Width costs concurrency, superlinearly

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
**150x** more tasks than `-n 64`.

**So `-n` is the throughput dial for a wide array, not `mem`.** Drop `-n` to
what the code truly uses. Two `-n 1` tasks that start now beat one `-n 2` task
that never does.

**On `standard`, slots are the bottleneck and memory is not.** Of 220 rome/cn
nodes, 210 have at least 800 GB free while only about 1900 of 14197 slots are
free. LSF states this directly in `bjobs -p`:

```
Affinity resource requirement cannot be met because there are not enough
processor units to satisfy the job affinity request: 30 hosts;
Job's requirements for resource reservation not satisfied (Resource: mem): 4 hosts;
```

30 hosts blocked on cores, 4 on memory. Cutting `-n` buys far more dispatch
speed than cutting memory. Reach for a smaller `n` with a larger `mem` per slot
before reaching for a different queue.

---

## B. Why bundling a large total is the common mistake

8 units x 64 GB is not a 512 GB job unless all 512 GB must be resident at once.
Split gives 8 tasks that each need one host with 64 GB free, which is nearly
every host. Bundled, it is one task needing 8 free cores **and** 512 GB on a
single host, which right now is 26 hosts instead of 47.

---

## C. Recipe per shape

### Wide and thin -- get more tasks running

1. `-n 1`, or the true thread count if higher. Never pad.
2. Omit `span[hosts=1]`; it is meaningless at `-n 1` and only adds a constraint.
3. Keep `rusage[mem]` tight but do not agonize: memory almost never blocks on
   rome (section A).
4. Do **not** throttle on `standard`. It has no `JL/U`, and pending jobs cost no
   fairshare -- only _running_ slots lower dynamic priority. Throttling here only
   caps your own throughput.
5. If a unit runs under ~2 min, chunk. -> [array-and-headers.md](./array-and-headers.md)
6. Consider splitting the range across queues. -> [array-and-headers.md](./array-and-headers.md)

### Narrow and fat -- get the one task placed

1. Match the pool's RAM-per-core ratio so the node is not left half-idle:
   rome ~15.7 GB/core, nodelmr up to ~62, nodelcr ~23, nodesd ~54.
2. `span[hosts=1]` is mandatory for anything shared-memory.
3. Move to `large_mem` or `large_core_count`. Both are near-idle, and a wide
   `-n` there competes with far fewer jobs.
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

---

## D. The 500 GB boundary, worked

`large_mem` has `USERS = all` and no minimum-memory admission check, so nothing
stops a small job from landing there. Do not use that. It is empty precisely
because people respect the boundary, and the pool is only 19 nodes.

| Option             | Directive                                  | Verdict                                                                        |
| ------------------ | ------------------------------------------ | ------------------------------------------------------------------------------ |
| Naive              | `-q standard -n 32 -R "rusage[mem=16000]"` | 500 GB, but 32 slots behind ~114000 pending jobs                               |
| Better on standard | `-q standard -n 4 -R "rusage[mem=128000]"` | 500 GB on 4 slots; dispatches far sooner, memory is not the constraint on rome |
| Correct            | `-q large_mem -n 8 -R "rusage[mem=64000]"` | 500 GB, matches the queue's purpose, pending is 0                              |

All three reserve the same 512000 MB; only the slot count and the queue differ.
Remember that `mem` is MB and 1 GB is 1024 MB, so 64000 is not 64 GB of a
512 GB total -- it is 62.5 GB of a 500 GB one. Take the third when the work
genuinely sits at or above 500 GB. Take the second when the real need is
300-450 GB and `large_mem` would be an abuse.

**Two queues that are not for arrays.** `large_core_count` has 4 nodes and
`superdome` has 3. A wide array on either will starve everyone including itself.
Send single large jobs there; send arrays to `standard`, `priority`, or
`large_mem`.

---

## E. `short` is a CPU-time queue, not a 30-minute wall-clock queue

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
