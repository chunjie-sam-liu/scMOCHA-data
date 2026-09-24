# Reference: chaining stages and diagnosing a job

The site facts behind `SKILL.md` section 10 and the symptom table behind section 11. Read this before writing a `-w` dependency expression, and when a job pends,
dies, or finishes without producing what it should.

---

## A. Site facts that decide whether a chain works here

`long-running-jobs` section 7 owns the method and the DAG submitter. This
section holds only the facts local to this cluster. Verified 2026-09-16 from
`bparams -a`, `bparams -l`, `man bsub`, and `man lsb.params`.

| Fact                                         | Setting                          | Consequence                                                                                                                                                                              |
| -------------------------------------------- | -------------------------------- | ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| Finished job records are purged after 24 h   | `CLEAN_PERIOD = 86400`           | `-w` on a job that finished more than a day ago can never resolve. Verify its outputs on disk and start a fresh root instead                                                             |
| A name dependency resolves to the newest job | `JOB_DEP_LAST_SUB = 1`           | `-w "done(meth-matrix)"` tests only the **most recently submitted** job of that name. Convenient, and a trap after a smoke `[1-1]` under the same `-J`. **Depend on job IDs, not names** |
| Dependencies are re-evaluated in batches     | `EVALUATE_JOB_DEPENDENCY = 1000` | a satisfied child does not start the instant its parent finishes; a few seconds of lag is normal, not a stuck job                                                                        |
| Arrays still cap at 4000                     | `MAX_JOB_ARRAY_SIZE = 4000`      | chaining changes nothing here; each array in the chain is capped separately                                                                                                              |

---

## B. The whole-array success condition is `numdone(<jid>,*)`

This site's `man bsub` documents `numdone(job_ID, operator number | *)` as "the
number of jobs in the DONE state satisfies the test. Use `*` (with no operator)
to specify all the jobs in the array." Its `numended`, `numexit`, `numrun` and
`numpend` siblings take the same form; `numexit(<jid>, >0)` is the matching "did
anything fail" test.

**`done(<jid>[*])` is element-wise, not whole-array.** The same man page: "Use
the `*` with dependency conditions to define one-to-one dependency among job
array elements such that each element of one array depends on the corresponding
element of another array. The job array size must be identical." Writing it to
mean "wait until the whole parent array succeeded" is wrong, and when the two
arrays differ in size LSF pairs only what it can.

`done(<jid>)` on a bare array job ID is accepted here -- the
`methylation_pipeline` release-v2 chain used it, with `323076538` held on
`done(323076536)` -- but this site's man page does not document its array
semantics. Prefer `numdone(<jid>,*)`, which does.

Pending jobs are free: `MAX_PEND_JOBS` is effectively unlimited and dynamic
priority falls with **running** slots, not pending ones. A chain that sits
queued all night costs no fairshare while it waits, so there is no reason to
hold stages back and submit them one at a time.

---

## C. Diagnose

| Symptom                                     | Command                                   | What it means                                                                                                                          |
| ------------------------------------------- | ----------------------------------------- | -------------------------------------------------------------------------------------------------------------------------------------- |
| Pending for a long time                     | `bjobs -p <jobid>`                        | reports hosts blocked on cores vs on memory; act on whichever dominates                                                                |
| Array running far fewer tasks than expected | recheck `-n` against the section 4 table  | width, not memory, is usually the cap                                                                                                  |
| Array status                                | `bjobs -A <jobid>`                        | PEND / RUN / DONE / EXIT per array                                                                                                     |
| Killed unexpectedly                         | `bhist -l <jobid> \| grep TERM_`          | `TERM_MEMLIMIT` = over the esub's `n` x `mem` threshold; `TERM_RUNLIMIT` = over `-W`; `TERM_CPULIMIT` = over `short`'s CPU-time budget |
| Reserved vs used                            | `bjobs -o "... memlimit max_mem" <jobid>` | the retune input from `SKILL.md` step 9                                                                                                |
| Queue changed                               | rerun the dispatch probe                  | queue pressure moves hourly                                                                                                            |
| Child still PEND after its parent finished  | `bjobs -l <child> \| grep -i -A2 depend`  | if the parent array had **any** EXIT element, `numdone(parent,*)` can never be satisfied and the child pends forever. Section A        |

A wide array reporting DONE is not proof. Count real outputs, per
`long-running-jobs` section 6.
