---
name: data-catalog
description: "Keep DATA.md current as the repository's data catalog: what each shared input is, what it means, where it came from, and whether it is still the one to use. Use when reading the env file for any task, when a variable is added to or changed in it, when a file arrives from a person, when a job or script produces an artifact a second stage or track will read, when a producing script changes what it writes, when an input is replaced or retired, when a row's schema facts are needed before writing a join, or when auditing the env file against DATA.md for drift. Covers the bidirectional env-vs-catalog audit, variable classification, the schema-fact block, the superseded protocol, and the no-fabricated-facts rule."
---

# Data catalog (`DATA.md`)

One `DATA.md` at the repository root is the catalog of every shared data asset.
It answers two questions that nothing else in the repository answers together:

1. **Provenance** — where did this file come from, when, from whom, and is it
   still the one to use.
2. **Description** — what _is_ this file: what one row means, which column is
   actually unique, and which columns a joining script can rely on.

It is read **before** touching data, not written afterwards as documentation. A
catalog that has fallen behind the env file is worse than no catalog, because
it is believed. Keeping it current is part of the work, not a cleanup step.

This skill is portable and names no project value. The env file's name, the
data roots, which variables are inputs versus output roots, and any existing
section layout come from the repository's `.github/instructions/` bindings or
`AGENTS.md`; where they disagree with this skill, the repository wins.

Related skills: `data-verification` (how to read the schema facts a row
records), `data-result-layout` (where a new file is written, and the
`_stale-YYYY-MM-DD` retirement rule), `analysis-pipeline` (a stage's own
outputs, which stay with the stage until a second consumer appears).

---

## 1. The env file is the watch list

The env file and `DATA.md` are two views of the same set of assets, and they
drift in **both** directions. Audit both directions, not just one.

**Forward — a variable with no row.** Match on whole words (`grep -w`). A
plain substring match silently reports a short variable as present because a
longer one contains it, and that failure is invisible:

```bash
grep -oE '^[a-z_]+=' <env-file> | tr -d '=' | while read -r v; do
  grep -qw -- "$v" DATA.md || echo "no row: $v"
done
```

This lists every variable, most of which are output roots that correctly have
no row. Read the result through the classification table below; only the ones
it marks `yes` are gaps. When the repository's bindings already list which
variables are inputs, diff against that list instead and the output is directly
actionable.

**Reverse — a file with no row.** List the actual contents of an input
directory and check each name against the catalog. This is the direction that
catches a file someone dropped in by hand, which by definition no variable
names and no script announced.

```bash
ls -1 <input-dir>/ | while read -r f; do
  grep -qF -- "$f" DATA.md || echo "no row: $f"
done
```

`-F` is not optional. Without it the filename is a regex, so `.` matches any
character, and a plain substring match reports an unlisted `foo.tsv` as present
because some row mentions `foo.tsv.gz`. Both failures are silent, and this
audit exists precisely to catch the file nobody announced.

Run this only on directories the catalog describes **per file**. A directory
that carries a single directory-level row — a large delivered tree, a rolling
package cache — is covered by that row, and running the reverse audit on it
reports every entry as missing. Catalog a directory as one row when its
contents arrive and are consumed as a unit, and per file when a script names
individual files; then audit accordingly.

**Paths resolve.** Every literal path in a row must exist on disk, and a row
that does not resolve is either a typo, a moved file, or a path written as a
glob rather than a name:

```bash
grep -oE '`<data-root>/[^`]+`' DATA.md | tr -d '`' | sort -u | while read -r p; do
  [[ -e "$p" ]] || echo "does not resolve: $p"
done
```

A `Pending` row and a row retired under section 5 are the two legitimate
non-resolving entries; everything else is a defect. `[[ -e ]]` already follows
symlinks, so no extra flag is needed here — it is `find`, `du`, and `rsync`
that need `-L` on a symlinked data root.

Run all three when the env file changes, before writing code that reads an
input, and when resuming work in a repository you have not touched recently.
Close every gap in the same change set.

### Which variables need a row

| The variable points at                                   | Row?                       |
| -------------------------------------------------------- | -------------------------- |
| A file received from a person or an outside group        | yes                        |
| A directory of external inputs                           | yes                        |
| A reference, annotation, or manifest tree                | yes                        |
| A file built here that another stage or track reads      | yes                        |
| A cross-stage data root holding files cataloged per file | no row for the root itself |
| An output root this project writes into                  | no                         |
| Scaffolding — repository root, logs, scratch, results    | no                         |

The fifth row is the one that makes the forward audit readable. A root that
gathers inputs — curated clinical tables, covariates, metadata — is neither an
external input nor an output root: the root gets no row, each file under it
does. Without that distinction most of the audit's output looks unclassifiable
and gets ignored.

Two cases the table does not cover, and both are common:

- **Pending** — the variable exists but the file does not yet. Record it with
  status `Pending`, name who owes it, and never silently omit it. An agent that
  finds the variable and no row cannot tell "not yet delivered" from
  "nobody wrote it down", and will guess.
- **Unnamed** — a data file exists that no variable points at. It still gets a
  row, with `Variable: none`. Not being in the env file is a fact about the
  file, not a reason to leave it undocumented.

An input whose origin is genuinely unknown is recorded as `UNKNOWN` and asked
about. Never infer an origin from a filename.

## 2. When a row is added or changed

Every trigger below is satisfied **in the same change set** that caused it.
Never batch catalog updates to the end of a session.

| Trigger                                                 | Action                                                                              |
| ------------------------------------------------------- | ----------------------------------------------------------------------------------- |
| A variable is added to or repointed in the env file     | Add or update the row; classify it with the table above                             |
| A file arrives from a person                            | Add a row with the sender and the delivery date                                     |
| Reference data is downloaded or re-downloaded           | Add a row with the URL and a checksum                                               |
| A second stage or track starts reading another's output | Add a row to the derived section naming the producing script                        |
| A producing script changes what it writes               | Update the row's description facts; add a `Superseded` entry if the meaning changed |
| An input is replaced or retired                         | `Superseded` entry — never an in-place edit; see section 5                          |
| A schema fact is verified while writing a join          | Record it on the row so the next session does not re-derive it                      |

**Adding a row is never a reason to edit the env file.** A missing, wrongly
named, or wrongly pointed variable is reported to the user, who changes it. The
catalog follows the env file; it does not steer it.

### Where a pipeline output belongs

A stage's own outputs belong to that stage, not to the catalog. An output earns
a `DATA.md` row **the moment a second stage or track reads it** — at that point
it is no longer a pipeline output, it is somebody's input. This is the line
between "the stage produced it" and "the project depends on it".

The row is therefore owed by the **consuming** change set, not the producing
one: the change that adds the cross-stage read adds the row, in the same change
set, before the code that reads it is written. A producer cannot know it has a
second consumer, so making it guess produces rows for outputs nobody ever reads
and misses the ones that matter.

## 3. What a row carries

Every row carries these five facts. Column layout may differ per section — an
external input records a sender, a derived input records a producing script —
but a row missing any of these is incomplete.

| Fact        | Rule                                                                                                                                                              |
| ----------- | ----------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| Path        | Repository-relative, literal, and resolvable. No ellipsis and no brace glob — `x{1,2}.parquet` is two rows, not one. The row has to be greppable and has to exist |
| Origin      | A person's name or the producing script. Never "the cluster"                                                                                                      |
| Date        | Received, downloaded, or built                                                                                                                                    |
| Variable    | The env-file variable that exposes it, or `none`                                                                                                                  |
| Consumed by | The tracks or stages that read it                                                                                                                                 |

A coordinate-bearing input carries a sixth: the **genome build**. A build that
is not known is recorded as `UNKNOWN` and confirmed with the producer before
anything joins on its positions.

## 4. The description facts

Provenance alone does not let an agent use a file. Add the facts that decide
whether code written against it is correct. Record only what was actually read,
using `data-verification` to read it.

| Fact      | Why it belongs                                                                                |
| --------- | --------------------------------------------------------------------------------------------- |
| Grain     | What one row _is_ (one sample, one array, one probe-sample pair). Decides every aggregation   |
| Key       | The column that is genuinely unique **on this file**, verified — not the one the name implies |
| Join keys | Which columns link it to other assets, and in what format (`chr1:100:A:G` vs `1_100_A_G`)     |
| Build     | For anything with coordinates                                                                 |
| Traps     | One line per gotcha that has already cost a wrong number or a failed run                      |

Keep these to one compact block or a few extra columns. **`DATA.md` is a
catalog, not a data dictionary.** When a file ships its own dictionary,
`SOURCES.tsv`, or checksum list, link it and say what it covers rather than
copying it — a duplicated provenance record forks on the first update.

## 5. Superseded, never edited

**Never silently edit a row in place when an input is replaced.** Add a
`Superseded` entry naming the replacement, the reason, and — the part that
matters — whether anything still reads the old file.

An input retired upstream but still consumed downstream is the single most
valuable thing this file records, and it is never resolved by editing
`DATA.md`. It changes cohort or analysis meaning, so it needs a plan and a
decision entry. **While anything still reads it, the file is not renamed** —
the old path stays exactly where the consumers expect it, and the `Superseded`
entry is what records that the project is knowingly running on a retired input.

Once nothing reads it, retire it the way `data-result-layout` prescribes:
rename `_stale-YYYY-MM-DD`, never delete, and record the old-path to new-path
mapping where that skill says to put it. Retiring is the one time the old row's
`Path` is updated in place — to the `_stale-` path, with a pointer to the
`Superseded` entry. That keeps the paths-resolve audit honest instead of
leaving a permanent false defect; it is bookkeeping, not a rewrite of history.

## 6. Never write a fact you did not just read

A date, size, row count, key, or build that was inferred rather than read is
worse than a blank, because the next session will trust it. Get every value
from the file itself, the env file, the producing script, or the user.

`UNKNOWN` is a legitimate and useful entry. A guess is not.

## 7. Before reporting

State what the audit actually found: which variables and files were checked,
which rows were added, and which gaps remain open with the reason each is still
open. An unresolved gap is reported, not quietly left — the value of a catalog
is that its silence means something.
