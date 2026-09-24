# Reference: dispatch probe, array sizing, and header templates

The procedures behind `SKILL.md` sections 6, 7, and 8. Read this when sizing an
array, choosing whether to throttle, chunking work that exceeds 4000 units, or
writing the resource lines of a `.lsf`.

---

## A. Probe dispatch speed before submitting

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

---

## B. Wide and thin: do not throttle `standard`

No `JL/U`, no fairshare penalty for pending. Submit the full range and let LSF
place whatever fits. A `%K` here caps only your own throughput.

```bash
#BSUB -J stage-step[1-3000]      # -n 1, standard, no %K
```

The ceiling is `min(4000 / n, 4000)`: the per-user running-slot cap and the
array index cap. At `-n 1` both are 4000.

### Split the range across queues when the backlog is long

`priority` dispatches roughly 26x sooner than `standard` per job but caps at 400
jobs per user; `short` takes 300 more if every task finishes inside 30 minutes.
Same script, three submissions, three job IDs to record in `PROGRESS.md`:

```bash
bsub -q priority < step.lsf      # edit -J to [1-400],    fills within minutes
bsub -q short    < step.lsf      # edit -J to [401-700],  tasks under 30 min only
bsub -q standard < step.lsf      # edit -J to [701-3000], absorbs the remainder
```

Only worth the extra bookkeeping when `standard` shows `P/R` above ~10 and the
array is large. Never split a _dependent_ range this way.

---

## C. Narrow and fat, or a small pool: throttle

`%K` is the suffix on the array spec:

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

---

## D. Chunking, when `N` exceeds 4000

Do not submit several arrays back to back. Keep the array narrow and give each
task a contiguous slice.

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

---

## E. Header templates

Use the script template in `analysis-pipeline/references/script-templates.md`
section B unchanged, and fill the header from `SKILL.md` steps 2 to 7. Only the
resource lines are this skill's business.

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

**Narrow and fat** -- 2 units, each needing 8 threads over one 500 GB object:

```bash
#BSUB -J <stage>-<step>[1-2]%2
#BSUB -o logs/<stage>/<step>_%J_%I.out
#BSUB -e logs/<stage>/<step>_%J_%I.err
#BSUB -n 8
#BSUB -R "rusage[mem=64000]"     # 8 x 64000 MB = 512000 MB = 500 GB total
#BSUB -R "span[hosts=1]"
#BSUB -W 48:00
#BSUB -q large_mem
```

- `span[hosts=1]` is required for anything shared-memory. Without it LSF may
  spread the slots across hosts and the threads cannot see each other.
- `-W` above the worst expected task, since the task is killed at the limit, but
  not absurdly above it: a tight `-W` helps backfill scheduling place the job.
- The array traps still apply: the `BASH_SOURCE` anchoring rule is in
  `long-running-jobs` section 5, and the two silent array bugs (reserved array
  name, swallowed exit code) are in its section 6.
