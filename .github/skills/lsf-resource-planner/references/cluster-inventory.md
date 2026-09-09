# Cluster inventory: hpcf_research_cluster

Measured snapshot. **Hardware and site config change rarely; queue pressure
changes hourly.** Trust the tables in sections 1-4 without re-checking. Re-probe
section 5 numbers with the command in `SKILL.md` step 6 before every submission.

| Field         | Value                               |
| ------------- | ----------------------------------- |
| Cluster       | `hpcf_research_cluster`             |
| Scheduler     | IBM Spectrum LSF Standard 10.1.0.14 |
| Master        | `splprhpc07`                        |
| Default queue | `standard`                          |
| Measured      | 2026-09-03                          |

---

## 1. CPU node families

From `lshosts`, joined against `bmgroup -r -w <hostgroup>`.

| Family                | Count | Cores/node | RAM/node   | RAM/core | Queues                                             |
| --------------------- | ----- | ---------- | ---------- | -------- | -------------------------------------------------- |
| `noderome` (standard) | 162   | 64         | 1003.3 G   | ~15.7 G  | `standard` `priority` `short` `heavy_io` `compbio` |
| `noderome` (fat)      | 36    | 64         | 1.9 T      | ~30 G    | same                                               |
| `noderome` (one-off)  | 1     | 64         | 940.3 G    | ~14.7 G  | same                                               |
| `nodecn`              | 2     | 64         | 1.4 T      | ~22 G    | same                                               |
| `nodelmr`             | 11    | 64         | **3.9 T**  | ~62 G    | `large_mem`                                        |
| `nodelmr`             | 6     | 80         | 2.9 T      | ~37 G    | `large_mem`                                        |
| `nodelmr`             | 1     | 80         | 2.8 T      | ~35 G    | `large_mem`                                        |
| `nodelmr`             | 1     | 40         | 1.4 T      | ~35 G    | `large_mem`                                        |
| `nodelcr`             | 4     | **128**    | 2.9 T      | ~23 G    | `large_core_count`                                 |
| `nodesd`              | 3     | **224**    | **11.8 T** | ~54 G    | `superdome`                                        |
| `nodecem`             | 9     | 40         | 750.3 G    | ~19 G    | `cryoem_cpu`                                       |
| `dragen`              | 5     | 32         | ~503 G     | ~15.7 G  | `dragen`                                           |

Total in `rhel8_cpu` (the `standard` pool): **201 hosts, 12864 cores**.

GPU pools, for completeness only: `nodegpu201-251` A100 (64 or 128 cores,
1-1.4 TB, `a100_40g` or `a100_80g`), `nodegpu252` H100 (64 cores, 2.9 TB),
`dgxac*` (128 cores, 1-1.9 TB). GPU work is out of scope for this skill.

## 2. Queues

| Queue              | PRIO | `JL/U`               | Host group          | Members         | Hard limits       |
| ------------------ | ---- | -------------------- | ------------------- | --------------- | ----------------- |
| `standard`         | 50   | none                 | `rhel8_cpu`         | 200 rome + 2 cn | none              |
| `priority`         | 100  | 400                  | `rhel8_cpu`         | same            | none              |
| `short`            | 100  | 300                  | `rhel8_cpu`         | same            | `CPULIMIT` 30 min |
| `heavy_io`         | 50   | 500 (queue MAX 2000) | `rhel8_cpu`         | same            | none              |
| `interactive`      | 200  | 8                    | `rhel8_interactive` | 6 noderome      | none              |
| `large_mem`        | 60   | 300                  | `rhel8_lm`          | 19 nodelmr      | none              |
| `large_core_count` | 60   | 300                  | `large_core`        | 4 nodelcr       | none              |
| `superdome`        | 100  | none                 | `superdome`         | 3 nodesd        | none              |
| `compbio`          | 50   | 3000                 | `rhel8_cpu`         | 200 rome + 2 cn | none              |

**No CPU queue sets `RUNLIMIT`, `MEMLIMIT`, or `PROCLIMIT`.** Verified across all
nine queues above. `bqueues -l | grep -c PROCLIMIT` returns 0 and
`grep -cE "^ *MEMLIMIT" lsb.queues` returns 0. Six queues in `lsb.queues` do set
`RUNLIMIT` (`cryoem` 72:0, `cryoem_rtx` 72:0, `cryoem_cpu` 240:0,
`kellogg_gpu` 240:0, `gpu_short` 30, `gpu_interactive` 48:00) but all are cryoem
or GPU queues, outside the scope of this skill.

The only limit any CPU queue sets is `CPULIMIT = 30` on `short`, and it is a
**CPU-time** limit, not a wall-clock one. See section 3.1.

Both `large_mem` and `large_core_count` are `USERS = all` with no minimum-memory
admission check. The "over 500GB" in their descriptions is a convention, not an
enforced rule. Respect it anyway; see `SKILL.md` step 5.

## 3. Site configuration

From `${LSF_ENVDIR}/lsf.conf` (`LSF_ENVDIR=/hpcf/lsf/lsf_prod/conf`) and
`bparams -a`.

```
LSF_UNIT_FOR_LIMITS=MB
LSB_MEMLIMIT_ENFORCE=y
LSB_JOB_MEMLIMIT=y
LSB_RESOURCE_ENFORCE="cpu gpu"
LSB_ESUB_METHOD="hpcf"
LSF_HPC_EXTENSIONS="CUMULATIVE_RUSAGE HOST_RUSAGE"

RESOURCE_RESERVE_PER_TASK = Y
ABS_RUNLIMIT              = Y
MAX_JOB_ARRAY_SIZE        = 4000
MAX_PEND_JOBS             = 2147483647
MAX_PEND_SLOTS            = 2147483647
DEFAULT_QUEUE             = standard
```

Every CPU queue carries:

```
DEFAULT_EXTSCHED = CPUSET[CPUSET_OPTIONS=CPUSET_CPU_EXCLUSIVE]
RES_REQ          = select[rhel8] affinity[core(1)]
```

so one slot is one exclusive physical core.

### Per-slot memory, and where the hard limit really comes from

Two separate mechanisms, easy to conflate.

**Reservation** is stock LSF: `RESOURCE_RESERVE_PER_TASK = Y` makes the
scheduler hold `n x rusage_mem` on the chosen host. It affects placement only.

**The hard limit** is set by a site esub. `lsf.conf` has
`LSB_ESUB_METHOD="hpcf"`, so `esub.hpcf` runs on every submission and calls
`update_mem()` from `${LSF_SERVERDIR}/hpcf_esub_funs.py` (that is the module
`esub.hpcf` imports; the dated `hpcf_esub_funs_*.py` files beside it are older
copies and none currently match it):

```python
hpcf_mem_lim = max( hpcf_mem_req, user_mem_lim, app_mem_lim )
...
if( 'LSB_SUB_NUM_PROCESSORS' in lsb_params ):
    ncores = int(lsb_params['LSB_SUB_NUM_PROCESSORS'])

lsb_params_new['LSB_SUB_RLIMIT_RSS']  = hpcf_mem_lim * ncores
lsb_params_new['LSB_SUB_RLIMIT_SWAP'] = hpcf_swp_lim
```

`hpcf_mem_req` is parsed from the user's `rusage[mem=...]`, `user_mem_lim` is
`-M`, `app_mem_lim` comes from an application profile. `LSB_SUB_RLIMIT_RSS` is
the `-M` value, so the esub **synthesises `-M` from `rusage[mem]` and multiplies
it by `-n`**.

Why this matters:

- `LSB_JOB_MEMLIMIT=y` alone does not do this. IBM ties that parameter to `-M`
  and queue `MEMLIMIT`; it only changes the limit from per-process to job-wide.
  No queue here sets `MEMLIMIT`, so without the esub a plain `rusage[mem]`
  submission would carry **no** memory limit.
- `LSB_MEMLIMIT_ENFORCE=y` decides the reaction: exceed the limit and the job is
  killed rather than suspended.
- This is a local policy script and can change without notice. Re-read
  `${LSF_SERVERDIR}/hpcf_esub_funs.py` if memory kills start behaving
  unexpectedly.

Default when nothing is requested:

```python
def default_mem():
    ### obscure number close to 3800 to make sure we know it is being set by THIS algorithm
    return 2499
```

accompanied by `HPCF: WARNING! No Memory was requested!` on stderr.

Arithmetic verified on two real jobs with `bjobs -o memlimit`:

| Submitted                              | `n` | esub computes           | Resulting `MEMLIMIT` |
| -------------------------------------- | --- | ----------------------- | -------------------- |
| `-n 10 -M 50000 -R "rusage[mem=5000]"` | 10  | `max(5000, 50000) * 10` | 488 G                |
| `-n 2 -R "rusage[mem=24000]"`, no `-M` | 2   | `max(24000, 0) * 2`     | 46.8 G               |

The second row is the proof that `rusage` alone produces a kill threshold here.

### `short` limits CPU time, not wall clock

Raw config from `lsb.queues`:

```
QUEUE_NAME   = short
CPULIMIT      = 30
DESCRIPTION  = High priority queue with a hard run time limit of 30 minutes
```

There is no `RUNLIMIT`; the `DESCRIPTION` is inaccurate.

`bqueues -l short` renders the limit as `30.0 min of svlprhpc02`, so it is
normalised by CPU factor:

| Host                              | `cpuf` |
| --------------------------------- | ------ |
| `svlprhpc02` (reference)          | 60.0   |
| `noderome*`, `nodelmr*` (compute) | 40.0   |

The allowance on a compute node is `30 x 60/40 = 45` CPU-minutes. With
`LSF_HPC_EXTENSIONS="CUMULATIVE_RUSAGE"` that CPU time is summed across every
process the job starts, so a fully loaded `-n k` job exhausts it in roughly
`45/k` wall-clock minutes. An I/O-bound job may never trip it at all.

## 4. User limits

| Limit                      | Value                 | Source                                        |
| -------------------------- | --------------------- | --------------------------------------------- |
| Running job slots per user | **4000**              | `busers <user>` MAX column                    |
| Pending jobs / slots       | unlimited in practice | `MAX_PEND_JOBS`, `MAX_PEND_SLOTS` = 2^31-1    |
| CPU quota                  | **none**              | no `PER_USER` slot limit in `lsb.resources`   |
| Memory quota               | **none**              | no `PER_USER` memory limit in `lsb.resources` |
| Per-queue job limits       | see section 2 `JL/U`  | `bqueues`                                     |

`blimits -u <user>` returns "No resource usage found" for CPU work. Every
`PER_USER` block in `lsb.resources` targets `ngpus_physical`; the non-GPU blocks
are all `SLOTS = 0` reservations that fence off other groups' hardware
(`cab_gpu_nodes`, `yugrp_reserve`, `dragen*`, `lmicro`, `cryoem_reserve`,
`compbio_tracker`, `mgmtHost`) and do not cap an ordinary user.

Both `large_mem` and `large_core_count` use `FAIRSHARE=USER_SHARES[[default,1]]`,
equal shares for everyone, so a low-usage user dispatches quickly there.
`standard` shares are weighted (`appdpcbauto` 5, `compbio@` 2, default 1).

## 5. Queue pressure snapshot, 2026-09-03

**Re-probe before every submission.** These numbers move hourly.

| Queue              | PEND   | RUN  | P/R  | Free slots | Verdict                             |
| ------------------ | ------ | ---- | ---- | ---------- | ----------------------------------- |
| `standard`         | 113629 | 8553 | 13.3 | 1938       | congested, hours to days            |
| `priority`         | 579    | 1262 | 0.5  | 1938       | good, capped at 400 jobs/user       |
| `short`            | 78     | 5    | 15.6 | 1938       | 30-min jobs only                    |
| `large_mem`        | **0**  | 594  | 0.0  | 710        | **dispatches immediately**          |
| `large_core_count` | **0**  | 272  | 0.0  | 240        | **dispatches immediately**, 4 nodes |
| `heavy_io`         | **0**  | 47   | 0.0  | 1938       | dispatches immediately              |
| `interactive`      | **0**  | 204  | 0.0  | 180        | 8 jobs/user                         |
| `superdome`        | 7859   | 472  | 16.6 | 328        | worst on the cluster                |
| `compbio`          | 3604   | 160  | 22.5 | 1938       | worse than `standard`               |

### Slots bind, memory does not

Free-memory distribution across the 220 `rhel8_cpu` hosts:

| Free memory | Hosts |
| ----------- | ----- |
| >= 1500 G   | 36    |
| 800-1500 G  | 174   |
| 400-800 G   | 9     |
| < 32 G      | 1     |

210 of 220 hosts had at least 800 GB free while only 1938 of 14197 slots were
free. `bjobs -p` on a pending array confirms the same thing:

```
Affinity resource requirement cannot be met because there are not enough
processor units to satisfy the job affinity request: 30 hosts;
Job's requirements for resource reservation not satisfied (Resource: mem): 4 hosts;
```

`large_mem`, `large_core_count`, and `superdome` nodes were nearly idle on
memory: `nodelmr21-31` at 3.8 T free, `nodelcr01-04` at 2.8 T free,
`nodesd01-03` at 11.5 T free.

### Placement capacity by job width

A task needs `n` free cores on **one** host, so width, not total free slots,
decides how much of an array can start.

| `-n` | `rhel8_cpu` hosts qualifying | tasks placeable | `rhel8_lm` tasks placeable |
| ---- | ---------------------------- | --------------- | -------------------------- |
| 1    | 47                           | 1049            | 700                        |
| 2    | 41                           | 518             | -                          |
| 4    | 34                           | 250             | 169                        |
| 8    | 26                           | 118             | 80                         |
| 16   | 18                           | 52              | 35                         |
| 32   | 14                           | 21              | 14                         |
| 64   | 7                            | 7               | 3                          |

Each doubling of `-n` roughly halves placeable tasks. The absolute numbers move
hourly; the halving relationship is structural.

Refresh:

```bash
bhosts -w rhel8_cpu | awk 'NR>1 && $4!="-" && $2=="ok"{f=$4-$6; if(f>0) print f}' \
 | awk '{a[NR]=$0} END{
     split("1 2 4 8 16 32 64", K, " ");
     for(j=1;j<=7;j++){k=K[j]+0; c=0; s=0;
       for(i=1;i<=NR;i++) if(a[i]>=k){c++; s+=int(a[i]/k)}
       printf "  -n %-3s : %3d hosts, %5d tasks placeable\n", k, c, s}
   }'
```

### Host targeting with `select[]`

`maxmem` is the host's total RAM in MB, `mem` is what is free right now. Both
work in `-R "select[...]"`. Verified counts (cluster-wide, **not** queue-scoped;
the queue's own host group narrows it further at submit time):

| Predicate                         | Hosts matching               |
| --------------------------------- | ---------------------------- |
| `select[rhel8 && maxmem>1900000]` | 62 (36 fat rome + GPU nodes) |
| `select[rhel8 && maxmem>900000]`  | 288                          |
| `select[rhel8 && mem>1500000]`    | 23                           |

Prefer `maxmem` over `mem`: a `mem` predicate re-evaluates each scheduling cycle
and can leave a job pending on a transient. Check a predicate before using it
with `bhosts -R "select[...]" | tail -n +2 | wc -l`.

## 6. Refresh commands

Run these only when the hardware or site config is suspected to have changed.

```bash
lsid                                  # cluster and LSF version
lshosts -w                            # hardware per host
bqueues -w                            # queues, PEND/RUN, JL/U
bqueues -l <queue>                    # host group, RES_REQ, limits
bmgroup -r -w <hostgroup>             # expand a host group to hostnames
bhosts -w <hostgroup>                 # per-host slot usage and state
lsload -w                             # live free memory per host
bparams -a                            # RESOURCE_RESERVE_PER_TASK, ABS_RUNLIMIT, array cap
grep -iE "LSF_UNIT|MEMLIMIT|RESOURCE_ENFORCE|ESUB" ${LSF_ENVDIR}/lsf.conf
busers <user>                         # per-user MAX running slots
blimits -u <user>                     # limits that currently apply
```

Re-derive the memory-limit rule if `LSB_ESUB_METHOD` changes:

```bash
grep -n "RLIMIT_RSS'\] = \|hpcf_mem_lim = max\|return 2499" \
  ${LSF_SERVERDIR}/hpcf_esub_funs.py
```

Confirm which queues carry a real wall-clock limit:

```bash
for q in standard priority short large_mem large_core_count superdome \
         heavy_io interactive compbio; do
  printf '%-18s RUNLIMIT=[%s] CPULIMIT=[%s]\n' "$q" \
    "$(bqueues -l $q | grep -A1 '^ *RUNLIMIT' | tail -1 | tr -s ' ')" \
    "$(bqueues -l $q | grep -A1 '^ *CPULIMIT' | tail -1 | tr -s ' ')"
done
```

Node-family summary in one pass:

```bash
lshosts -w | awk 'NR>1 && $5!="-" {c=$1; gsub(/[0-9]+$/,"",c);
  k[c" "$5"c "$6]++} END{for(i in k) print k[i], i}' | sort -k2
```

Scratch for these probes goes to `"${tmpdir}/<task-name>/"` and is removed
afterwards. -> `data-result-layout`
