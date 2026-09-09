---
name: data-result-layout
description: "Decide where a generated file belongs: results/ for small final human-readable deliverables, data/ for everything large or intermediate, and the repository tmp/ root for scratch. Use when a script needs an output path, when creating an output directory, when choosing between a cache and a final table or plot, when versioning outputs to avoid overwriting, when retiring or superseding an output or script, when a tool needs a temp directory, when a path or directory variable should come from the env file instead of being hard-coded, or when find/du/rsync/git return a surprising result under data/ or results/. Covers the results-vs-data split, the cross-stage data roots and their env variables, the no-hard-coded-path rule, versioned output directories, the rename-to-_stale-{date} retirement rule, the never-write-to-/tmp rule, and the symlinked data/results/logs/tmp roots."
---

# Data and result layout

Keep a strict separation between `results/` and `data/`, and mirror the
analysis script structure in both. Anyone should be able to map an output file
back to the exact script that produced it, from its path alone.

**Never overwrite existing results unless explicitly requested.** Prefer
versioned output directories with a hyphen suffix (`<name>-v2`, `<name>-v3`).
The underscore form is reserved for `_stale-YYYY-MM-DD`, section 7.

This skill is portable and names no project value. The repository's own data
roots, results roots, env-file names, and pre-existing directory names come
from its `.github/instructions/` bindings or `AGENTS.md`; where they disagree
with this skill, the repository wins.

---

## 1. What goes where

- **`results/`** — ONLY final, small, human-readable outputs: plots and small
  result tables. Allowed formats: `pdf`, `png`, `svg`, `xlsx`, `tsv`, `csv`.
  These are the deliverables a person opens and reads directly. "Small" means
  under ~50 MB and openable in Excel or an image viewer; full summary
  statistics, per-variant tables, and any whole-genome-scale output belong in
  `data/` even though they are plain text.
- **`data/`** — ALL large and intermediate data: raw inputs, genotype /
  PLINK / VCF files, per-chromosome shards, model objects, and
  `.rds` / `.qs` / `.parquet` / `.fst` caches.
- **`logs/`** — every process run log, and nothing else, in one subdirectory per
  stage: `logs/<stage>/`. A cluster task writes
  `logs/<stage>/<NN-step-name>_<jobid>_<arrayindex>.out` and `.err`; a tmux
  launch writes `logs/<stage>/<session>.log`. The root is bound to `logdir` in
  the env file. Create the directory before submitting — `mkdir -p` belongs in
  `config.sh`; a scheduler that cannot write to its output path loses the log.
  Never write a run log next to the script, under `data/`, or under `results/`,
  and never promote a log to `results/`. Logs written flat under `logs/` predate
  this rule; leave them, and put new ones in `logs/<stage>/`.
- **`tmp/`** — scratch and throwaway files only, never a pipeline input or a
  deliverable. One subdirectory per task, removed when the task ends. The root
  is bound to `tmpdir` in the env file; section 6 is the full rule.

Rule of thumb: if a file is large, binary, a shard, or only consumed by a later
script rather than read by a human, it belongs in `data/`. When unsure, default
to `data/` for intermediates and promote only the final small table or plot to
`results/`.

A few directories created before this rule still hold whole-genome tables under
`results/`. Leave them where they are -- moving them would break every existing
reference -- but never add a new one.

## 2. Mirror the stage, subdivide by artifact kind

Do NOT dump every output into one flat folder. The stage directory name is what
carries through to both output roots. Inside it, subfolders are named after the
kind of artifact, not after the step number that wrote them.

```
src/NN-stage/  ->  data/intermediate/NN-stage/{cohort,geno,masks,hits,...}/
                   results/NN-stage/{figures,tables}/
```

- A new stage `src/NN-stage/` writes intermediates to
  `data/intermediate/NN-stage/` and final outputs to `results/NN-stage/`,
  reusing the numbered stage directory name verbatim.
- Name subfolders after the artifact (`cohort/`, `geno/`, `masks/`, `burden/`,
  `annotation/`, `hits/`, `thinned/`, `figures/`), not after the step. A step
  that produces a single file writes it directly under the stage directory.
- Older output directories named after the analysis rather than the stage
  number (`data/intermediate/<analysis-name>/`, `results/<analysis-name>/`)
  predate this rule. Keep writing to them for the stage that owns them; never
  rename or migrate them. Only newly created output directories use the
  numbered form. The repository bindings list which ones exist.

## 3. Cross-stage data lives at the top of the data root

An artifact that more than one stage consumes is not a stage intermediate. It
lives directly under the data root, one directory per kind, each bound to a
variable in the env file:

```
data*/geno/          genodir           genotype filesets, per-chr shards
data*/pheno/         phenodir          phenotype tables
data*/covar/         covardir          covariate tables
data*/pca/           pcadir            eigenvec / eigenval
data*/model/         modeldir          model registry, fitted null models
data*/clinical/      clinicaldir       curated clinical input
data*/annotation/    annotationdir     VEP / gnomAD annotation
data*/reference/     refdir            reference panels, GTF, external files
data*/raw/           rawdir            untouched source data
data*/gwas*/         gwasdir, ...      per-phase GWAS outputs
data*/intermediate/  intermediatedir   per-stage intermediates, section 2
```

Something used by exactly one stage stays in `data/intermediate/<stage>/`.
Promote it to a top-level directory only once a second stage needs it, and add
the matching env variable in the same change.

The table is the recommended default, not an inventory. Only the variables the
repository's env file actually defines exist; check it, and add a missing one
there rather than rebuilding the path from `repodir`.

## 4. Never hard-code an absolute path

Every path in a script is built from an env variable. The env file is the only
place an absolute path is ever written, including paths to data outside the
repository.

- Shell: the stage `config.sh` sources the env file with `set -o allexport`,
  then the script uses `"${intermediatedir}/<stage>/..."`.
- R: `jutils::dotenv()` at the top of the script, then
  `fs::path(intermediatedir, "<stage>", "file.tsv")`.
- Honor these variables instead of rebuilding the path from `repodir`, and
  append the stage subpath rather than writing to the root of the variable.
- A repo with more than one track may have one env file per track, or one
  shared env file with a distinct root variable per track. Never mix them: a
  script loads the env file the repository bindings name for its own track and
  writes under that track's own root variable, never another track's.

## 5. The output roots are symlinks

`data*/`, `results*/`, `logs/`, and `tmp/` at the repository root are symlinks
into a separate storage area, so the git work tree stays small while the bulk
data lives on another filesystem. Assume this layout in every project; confirm
with `ls -ld data results logs tmp` when it matters.

Four consequences produce a wrong answer with a zero exit code:

- `find data -maxdepth 1` returns only the link itself. `find` does not descend
  into a symlinked start point. Use `find -L data ...`, or a trailing slash
  (`find data/ ...`), which GNU `find` resolves.
- `du -sh data` reports the size of the link (`1.0K`), not the tree (`14T`).
  Use `du -shL`.
- `rsync -a results/ dest/` copies the link, not the contents. Use `rsync -aL`.
- `git status` and `git diff` never show anything under these roots. A
  generated output will not appear in a diff; that is not evidence it is
  missing.

`rm -rf results/` removes the symlink and orphans the real data, while
`rm -rf results/*` destroys the real data. Neither is allowed without an
explicit request in the current message; see section 7.

Resolve the real location with `readlink -f <path>` before quoting a path to
the user or writing it into a `PLAN.md`.

## 6. Temporary files (hard rule)

**Never write scratch, throwaway, or temporary files to `/tmp`.** Scratch goes
to the repository's own `tmp/` root, bound to `tmpdir` in the env file. It is a
symlink into the same storage area as `data/` and `results/` (section 5), so
spill files land on the filesystem sized for them instead of on `$HOME`.

- Applies to every scratch artifact: smoke-test fixtures, debug dumps, one-off
  extracts, downloaded archives, sort/join spill files, generated test data.
- Create a task-specific subfolder, `"${tmpdir}/<task-name>/"`, not loose files
  in the root.
- Set `export TMPDIR="${tmpdir}"` for tools that pick their own temp location
  (`sort`, `bcftools`, PLINK workflows, R `tempdir()`, Python `tempfile`).
- Build the path from the variable, never from a hard-coded absolute path
  (section 4).
- Clean up the task subfolder when the task finishes.
- Never write scratch files into the git work tree.

### When the repository has no `tmp/` root

No `tmpdir` in the env file, or no `tmp/` at the repository root. Do not create
one, and do not improvise a scratch path under `data/`. Instead:

1. **Tell the user** the repository has no `tmp/` root, so they can create it
   and point it at their own storage.
2. Fall back to `$HOME/tmp/<task-name>/` for that task, and
   `export TMPDIR="$HOME/tmp"`. Write `$HOME`, not `~`: a quoted `~` is not
   expanded.

The fallback is per-task, not a new default. The `tmp/` root is a symlink the
user owns, so never add `tmpdir` to the env file or create the symlink
yourself; report it and let the user do it.

## 7. Retiring outdated scripts and data (hard rule)

**Never delete an outdated script, data file, or output directory. Rename it
with a `_stale-YYYY-MM-DD` suffix**, where the date is the day it is retired,
then let the new run recreate the original path from scratch.

```
data/pheno/                 ->  data/pheno_stale-2026-08-23/
data/geno/<cohort>/         ->  data/geno/<cohort>_stale-2026-08-22/
results/NN-stage/           ->  results/NN-stage_stale-2026-07-29/
src/NN-stage/               ->  src/NN-stage_stale-2026-08-20/
```

- Underscore before `stale`, hyphen before the date: `_stale-YYYY-MM-DD`.
- Suffix the leaf name only and keep the parent path unchanged, so the archive
  sits beside the live path and is obvious in one `ls`.
- For a single file the suffix goes before the extension:
  `NN-step.R` -> `NN-step_stale-2026-08-20.R`.
- One date suffix per name. If a path goes stale a second time, the earlier
  archive already carries its own date, so never stack suffixes.
- Applies equally to `data*/`, `results*/`, and retired stage or step scripts
  under `src/` and `src_*/`.
- Write the old-path -> new-path table into the active `PLAN.md` /
  `PROGRESS.md` before renaming, so existing references stay resolvable.
- Deleting is allowed only when the user explicitly asks for it in the current
  message.
