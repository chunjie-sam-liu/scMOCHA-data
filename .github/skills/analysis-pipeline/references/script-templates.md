# Reference: script templates and `config.sh` contract

Language-specific templates for a pipeline stage. All share the same metainfo
header block; only the comment syntax and runner change.

---

## A. Shell step script

```bash
#!/usr/bin/env bash
# @AUTHOR: <Name>
# @CONTACT: <email>
# @DATE: YYYY-MM-DD
# @DESCRIPTION: One-sentence purpose. Usage: bash NN-step.sh <arg1> <arg2>
# @VERSION: v0.1.0

set -euo pipefail

# Stage config: sources the env file, activates pixi, defines every path.
source "$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/config.sh"

arg1=${1:?Usage: bash NN-step.sh <arg1> <arg2>}
arg2=${2:?Usage: bash NN-step.sh <arg1> <arg2>}
```

Resolve `config.sh` from `${BASH_SOURCE[0]}`, never from a hard-coded repo
path. Source it first, then validate positional args with `:?Usage: ...`.

This applies to a `.sh` invoked as `bash path/to/script.sh`. It does **not**
apply inside a `.lsf` or `.sbatch`; see sections B and C.

## B. LSF array wrapper

```bash
#BSUB -J stage-step[1-N]
#BSUB -o logs/NN-stage-name/NN-step_%J_%I.out
#BSUB -e logs/NN-stage-name/NN-step_%J_%I.err
#BSUB -n 4
#BSUB -R "rusage[mem=16000]"
#BSUB -W 8:00
#BSUB -q standard

# @AUTHOR: <Name>
# @CONTACT: <email>
# @DATE: YYYY-MM-DD
# @DESCRIPTION: One-sentence purpose.
# @VERSION: v0.1.0

# `bsub <` reads this file from stdin, so BASH_SOURCE is the jobspool copy.
repodir="${LS_SUBCWD:-$PWD}"

TASKS=(unit1 unit2 ... unitN)
unit=${TASKS[$((LSB_JOBINDEX-1))]}

bash "${repodir}/src/NN-stage-name/NN-step.sh" "$unit" ; status=$?
echo "[done] $unit (exit $status)"
exit $status
```

Three traps, all verified in production:

- **Never resolve a path from `${BASH_SOURCE[0]}` in a scheduler script.**
  `bsub < job.lsf` feeds the script through stdin, so the scheduler runs a copy
  under `/lsf_jobspool/...` and `BASH_SOURCE` points there. Sourcing
  `config.sh` that way fails with `No such file or directory` before any work
  starts. Use `${LS_SUBCWD:-$PWD}` (LSF) or `${SLURM_SUBMIT_DIR:-$PWD}`
  (SLURM), and submit from the repository root. A bare relative path like
  `bash src/NN-stage/NN-step.sh` happens to work from the repo root but breaks
  confusingly from anywhere else, so anchor it explicitly.
- **Never name the array with a bash-reserved variable.** `GROUPS` is readonly
  and silently ignores assignment, so `${GROUPS[i]}` returns a numeric GID.
  Avoid `PIPESTATUS`, `BASH_SOURCE`, `FUNCNAME`, `UID`, `EUID`, `PPID` too.
  Use `TASKS`, `UNITS`, or `NEWGROUPS`.
- **Propagate the worker exit code.** Any command after the worker (even
  `echo`) makes the task exit 0, so the scheduler records DONE on failure.
  Capture and re-raise as shown above.

## C. SLURM array wrapper

```bash
#SBATCH --job-name=stage-step
#SBATCH --array=1-N
#SBATCH --output=logs/NN-stage-name/NN-step_%A_%a.out
#SBATCH --error=logs/NN-stage-name/NN-step_%A_%a.err
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=08:00:00

# @AUTHOR: <Name>
# @CONTACT: <email>
# @DATE: YYYY-MM-DD
# @DESCRIPTION: One-sentence purpose.
# @VERSION: v0.1.0

repodir="${SLURM_SUBMIT_DIR:-$PWD}"

TASKS=(unit1 unit2 ... unitN)
unit=${TASKS[$((SLURM_ARRAY_TASK_ID-1))]}
```

Map the index to the unit of work with integer math. Never branch on unit-name
strings inside the wrapper. The three traps in section B apply here too, with
`SLURM_SUBMIT_DIR` in place of `LS_SUBCWD`.

## D. `config.sh` contract

Every stage with shell scripts ships a `config.sh` that:

1. Derives `repodir` from its own location, never a hard-coded path:
   `repodir=$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)`. The `../..`
   assumes `src/NN-stage/config.sh`; adjust the depth to match.
2. Sources the repo's env file with
   `set -o allexport ; source "${repodir}/.env" ; set +o allexport`. A repo with
   more than one track sources its own track's file; the repository bindings
   name them.
3. Activates pixi once for the whole process:
   `eval "$(pixi shell-hook --manifest-path "${repodir}/pixi.toml")"`. Never
   conda, mamba, or Miniforge. After this line `Rscript`, `python`, and every
   other pixi binary are on `PATH` — do not wrap them in `pixi run` again.
4. Defines every input prefix, output root, run name, and model/unit list.
5. Runs `mkdir -p` for every output directory it advertises, including
   `logs/<stage>/`, which the scheduler will not create for itself.
6. Defines `check_env()` that prints `MISSING tool: ...` / `MISSING path: ...`
   for each absent binary or file and returns non-zero on any miss.

Never enable `set -e` in `config.sh`: it is sourced, not executed, so a
non-zero test would kill the calling script.

Documented in `config.sh.md`, one row per exported variable. The R-side twin is
`config.R`; see section G.

## E. R step script (jutils style)

```r
#!/usr/bin/env Rscript
# Metainfo --------------------------------------------------------------
# @AUTHOR: <Name>
# @CONTACT: <email>
# @DATE: YYYY-MM-DD
# @DESCRIPTION: ...
# @VERSION: v0.1.0

# Reproducibility --------------------------------------------------------
# Library ----------------------------------------------------------------

suppressMessages({library(jutils)})          # ggplot2, data.table, dplyr, fs, ...

# Args -------------------------------------------------------------------

GetoptLong.options(help_style = "two-column")
verbose = FALSE
GetoptLong("verbose", "Enable verbose logging")

# Logger -----------------------------------------------------------------

log_layout(layout_glue_colors)
log_threshold(if (isTRUE(verbose)) TRACE else INFO)

# Load data --------------------------------------------------------------

library(jutils)                              # repeated on purpose, attach is idempotent
dotenv()                                     # repodir / intermediatedir / resultsdir
set.seed(9527)
suppressMessages({conflicted::conflicts_prefer(dplyr::filter, fs::path)})

# Source -----------------------------------------------------------------

source("src/NN-stage-name/config.R")          # stage constants, paths, helpers
paths <- stage_paths()
```

Section order: Metainfo / Reproducibility / Library / Args / Logger /
Load data / Source / Conn / Function / Main / Save / Session info. The
`rmeta` editor snippet emits exactly this skeleton. The `Source` section is
where the stage's `config.R` goes; see section G.

Use jutils `import()` / `export()` by default, data.table `fread` / `fwrite`
when the PLINK2 `#FID` header or exact quoting matters, and `saveplot()` for
figures; ggplot2 + patchwork for plotting, and `fs::path` for the `a / b / c`
path operator. Pull colors from the track's shared color file instead of
inlining hex codes. `jutils` is the only package ever loaded with `library()` — everything
else goes through `load_pkg(...)`. Run with `pixi run Rscript ...`. See the `jutils` skill for helper
details.

## F. Python step script

```python
#!/usr/bin/env python
# @AUTHOR: <Name>
# @CONTACT: <email>
# @DATE: YYYY-MM-DD
# @DESCRIPTION: ...
# @VERSION: v0.1.0
```

Run with `pixi run python ...` from the repository root, or `python` directly
inside a script that has sourced the stage `config.sh`; `uv run python ...`
only in a repository with no `pixi.toml`. Mirror the same argument validation
and path-from-`.env` behaviour as the R template.

## G. `config.R` contract (R stages)

A stage whose steps are R scripts ships `config.R`, the R-side twin of
`config.sh`. It is never executed on its own. Every R step in the stage sources
it, so that no step hard-codes a path, a threshold, or a column name.

```r
#!/usr/bin/env Rscript
# Metainfo ----------------------------------------------------------------
# @AUTHOR: <Name>
# @CONTACT: <email>
# @DATE: YYYY-MM-DD
# @DESCRIPTION: Shared constants + helpers for stage NN. Sourced by the NN-NN
#               scripts after library(jutils) + dotenv().
# @VERSION: v0.1.0

STAGE <- "NN-stage-name"

# Constants: thresholds, cohort lists, column contracts, trait maps.
P_GW <- 5e-8
COHORTS <- c("cohort-a", "cohort-b")
REQUIRED_COLS <- c("id", "trait", "age")

# --- Paths ----------------------------------------------------------------
# Reads Sys.getenv(), so dotenv() must already have run in the caller.
stage_paths <- function() {
  intermediatedir <- fs::path(Sys.getenv("intermediatedir"))
  list(
    interstage = intermediatedir / STAGE,
    phenodir = fs::path(Sys.getenv("phenodir")),
    outdir = fs::path(Sys.getenv("resultsdir")) / STAGE
  )
}

# --- Helpers --------------------------------------------------------------
fn_load_cohort <- function(cohort, paths) {
  f <- paths$interstage / glue::glue("{cohort}.qs2")
  if (!fs::file_exists(f)) cli::cli_abort("Not found: {.path {f}}")
  import(f)
}
```

It contains exactly four kinds of thing:

1. `STAGE` — the stage directory name, so every output path derives from one
   string.
2. Stage constants in `UPPER_SNAKE_CASE` — thresholds, cohort lists, required
   column vectors, trait maps. Anything a reviewer would want to check in one
   place.
3. `stage_paths()` — returns a named list of `fs::path` values built from
   `Sys.getenv()`. This is the only place an output path is assembled.
4. Shared helpers, `fn_` prefixed — used by two or more steps. A helper used by
   exactly one step stays in that step's `Function` section.

It must NOT:

- call `library()` — the caller attaches `jutils` first;
- call `dotenv()` — the caller does, and `stage_paths()` depends on that order;
- have side effects — no file writes, no `dir_create()`, no `set.seed()`, no
  logging at source time. Sourcing it twice must be harmless;
- contain an absolute path — those live only in the env file. -> `data-result-layout`
- define a color — every hex lives in the track `color.R`. A config that needs
  a palette aliases one (`COHORT_COLORS <- color_newgroup`). -> `palette`

How a step sources it:

```r
suppressMessages({library(jutils)})
dotenv()                                 # or dotenv("<track>.env") in a multi-track repo
source("src/NN-stage-name/config.R")
paths <- stage_paths()
```

The `source()` path is repository-relative, so every R step must be launched
from the repository root (`pixi run Rscript src/NN-stage/NN-step.R`). Never
`source("config.R")` bare: that only resolves when the working directory
happens to be the stage directory.

A stage with both shell and R steps ships **both** `config.sh` and `config.R`.
They never source each other; each is the single source of truth for its own
language, and both read the same `.env`.

### Documenting a config file

A config file's sibling doc appends `.md` to the **full** filename, not to the
stem: `config.sh` -> `config.sh.md`, `config.R` -> `config.R.md`. Step scripts
keep the stem form (`NN-step-name.R` -> `NN-step-name.md`); the configs are the
exception, because otherwise both would claim `config.md`.

`config.sh.md` documents each exported variable and the `check_env()` contract.
`config.R.md` uses **Purpose**, **How to source it**, **Constants**, **Paths**,
**Functions**.

Neither carries a `Run` or a `Verify` section, because neither config file is
executable. That is the one documented exception to the paired-`.md` structure
in `SKILL.md` section 3.

Stages written before this convention use a single `config.md`, or
`00-config.R` + `00-config.md`. Those are historical: read them, but do not
rename them. New stages use `config.sh` + `config.sh.md` and `config.R` +
`config.R.md`.
