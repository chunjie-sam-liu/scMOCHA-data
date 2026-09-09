# Reference: worked GWAS / local-ancestry pipeline pattern

A concrete example of the `analysis-pipeline` skill's rules, drawn from a
cohort-scale GWAS repository with echo, cardiotoxicity, and local-ancestry
stages. The shapes and conventions generalise to any cohort-scale HPC pipeline;
replace the `gwas-*` names with the domain-specific names for your project.

---

## A. Directory shape

```
src/NN-stage-name/
  AGENTS.md                     track guide: run commands, entry points, pitfalls
  PLAN.md                       stage charter: goal, steps, cross-cutting decisions
  PROGRESS.md                   stage dashboard + campaign index
  DECISION.md                   stage decision log, D1..Dn, append-only
  {date}-{short-title}.PLAN.md      one specific question, design + Q&A
  {date}-{short-title}.PROGRESS.md  that question's run state
  {date}-{short-title}.DECISION.md  that question's decisions
  config.sh + config.sh.md      shared paths, model, check_env()
  config.R + config.R.md        stage constants, stage_paths(), shared fn_* helpers
  01-step-name.sh + .md         one task per file
  01-step-name.lsf + .md        optional LSF array wrapper
  02-step-name.R + .md
  ...
  DIAGRAM.md                    optional figure brief, see the image-prompt skill
```

`AGENTS.md` is required at the track root and optional per stage; the other
three tier-1 files belong to every stage. A stage with only shell steps omits
the `config.R` pair, and one with only R steps omits the `config.sh` pair; a
stage with both ships both.

## B. Data layout (mirrors across stages)

```
data/intermediate/<stage-key>/
  pheno/<OUTCOME>.pheno.tsv     PLINK2 #FID IID PHENO tsv
  covar/<OUTCOME>.covar_*.tsv   PLINK2 #FID IID ... tsv
  <model>/<OUTCOME>.*.glm.{linear,logistic}   whole-genome sumstats
  <other per-stage subdirs>
results/<stage-key>/
  plots/
  comparison/
  tables/                       small filtered tables only
logs/<stage-key>/                 tmux + LSF stdout/stderr
```

Per-variant genome-wide output stays in `data/` no matter its format; only
small filtered tables and plots are promoted to `results/`. See the
`data-result-layout` skill for the threshold.

`<stage-key>` above is the legacy directory name of an existing stage
(`gwas-echo`, `gwas-cm`, `gwas-local-ancestry`). Keep using it for the stage
that already owns it; a new stage instead names its output directories after
its numbered `NN-stage` directory. The PLINK file header convention
`#FID IID ...` is non-negotiable — downstream R must
`setnames(d, sub("^#", "", names(d)))`.

## C. config.sh contract

The contract itself is in [script-templates.md](./script-templates.md)
section D. What this repo adds on top of it:

- `repodir` resolves by walking up from `${BASH_SOURCE[0]}`, so a stage works
  from any checkout. Each track sources its own env file; the repository
  bindings name them.
- Contract item 4 here means RUN_NAME, output roots, input prefixes, and the
  ancestry / model / trait lists.

A stage's `.sh` and `.lsf` scripts source `config.sh` _first_, then validate
their own positional args with `:?Usage: ...`.

## D. LSF array submitter pattern

Long chains use one submitter `.sh` that issues all `bsub` calls with
dependencies via `-w "done(${JOB_PREFIX}_prev)"`. Simpler stages use one
`.lsf` per step with `[1-N]` arrays and a sibling `.sh` that does the unit
work. Decision rule: if steps have heterogeneous resources or pure ordering
dependencies, write one chained submitter; if all units share resources,
write a single `.lsf` array.

LSF naming: `#BSUB -J stage-step[1-N]`, `-o logs/<stage>/NN-step_%J_%I.out`,
`-e logs/<stage>/NN-step_%J_%I.err`. `LSB_JOBINDEX` maps to outcome/chr via integer
math; never branch on outcome name strings inside the LSF script.

## E. R script style (jutils)

Header block first (Metainfo / Reproducibility / Library / Args / Logger /
Load data / Source / Conn / Function / Main / Save / Session info). After
`suppressMessages({library(jutils)})` these are available without any further
load: ggplot2, patchwork, data.table, dplyr, tidyr, purrr, arrow, fs, glue,
GetoptLong, logger, cli. Anything else goes through `load_pkg(...)`.

Path operator: `fs::path` lets you write `intermediatedir / "gwas-cm" / "pheno"`
instead of `file.path`. Enable it via
`conflicted::conflicts_prefer(dplyr::filter, fs::path)`.

Reading PLINK2 tsv:

```r
d <- fread(p, header = TRUE)
setnames(d, sub("^#", "", names(d)))
d[, c("FID","IID") := lapply(.SD, as.character), .SDcols = c("FID","IID")]
```

Writing PLINK2 tsv:

```r
dd <- copy(d); setnames(dd, "FID", "#FID")
fwrite(dd, p, sep = "\t", quote = FALSE, na = "NA")
```

Plot save: `jutils::saveplot(path, plot, width=, height=)`. Compose with
`patchwork`: `p1 / p2`, `wrap_plots(list_of_p, ncol=)`.

Colors: `source(repodir / "src" / "color.R")` and use named vectors like
`color_ancestry["AFR"]`, `color_sex["Female"]`. Add new palettes to the track
color file rather than inlining hex codes; see the `palette` skill.

## F. Verification gates

Always run before declaring a stage done:

```bash
find src/NN-stage -maxdepth 1 \( -name '*.sh' -o -name '*.lsf' -o -name '*.sbatch' \) -print0 | xargs -0 -r -n1 bash -n
pixi run Rscript -e 'for (f in list.files("src/NN-stage", "\\.R$", full.names=TRUE)) tryCatch(parse(f), error=function(e) stop(f, ": ", conditionMessage(e)))'
```

The per-step `.md` Verify blocks should be runnable as shell pasted by a
reviewer who hasn't read the code.

## G. Plan-first decisions

For non-trivial stages, the plan file holds Q&A blocks where contested
design choices are locked in before implementation (e.g. Q1: 2-way vs 3-way
ancestry model). When new scripts contradict the plan, update the plan in
the same change set.

## H. Tmux long-running

A whole-genome PLINK pass, or an R job that loops over 22 chromosomes, runs in
tmux. A submitter script runs in tmux only if it does real work between
`bsub` calls; a script that only issues submissions returns in seconds and
runs inline. The `long-running-jobs` skill owns the threshold.

```bash
SESSION=la-rfmix-run
STAGE=05-gwas-local-ancestry
mkdir -p "logs/$STAGE"
tmux has-session -t "$SESSION" 2>/dev/null && { echo "$SESSION still running"; exit 1; }
tmux new-session -d -s "$SESSION"
tmux send-keys -t "$SESSION" "bash src/$STAGE/02-rfmix-infer.sh > logs/$STAGE/$SESSION.log 2>&1; exit" C-m
```

The agent does NOT poll for completion; LSF + tmux are persistent.

## I. Worked grid for the LA stage

| File                                               | Triggers                                                  |
| -------------------------------------------------- | --------------------------------------------------------- |
| `01-prepare-la-target.sh`                          | dense per-chr QC, merge with plink1.9, summarize          |
| `02-rfmix-infer.sh`                                | submitter for 9-step RFMix chain, reuses lab `R/01..05`   |
| `03-summarize-global-ancestry.R`                   | `.rfmix.Q` aggregation, weighted + simple mean            |
| `04-prepare-la-inputs.R`                           | 10 outcomes × {pheno, covar_extended_la, tractor dt+keep} |
| `05-run-plink-global-la.{sh,lsf}`                  | 10-task array, linear / logistic by outcome               |
| `06-tractor-extract-tracts.{sh,lsf}`               | 220-task array (10 × 22)                                  |
| `07-run-tractor.{sh,lsf}` + `07-concat-tractor.sh` | 220-task scan + concat                                    |
| `08-post-la-gwas.R`                                | λ_GC + Manhattan/QQ with ancestry palette                 |
| `09-compare-standard-vs-la.R`                      | β scatter + hit retention vs standard extended            |
| `12-figure-spec.md`                                | historical image-prompt name; new stages use `DIAGRAM.md` |

This grid is the canonical shape for any future ancestry / GWAS extension. The
sibling `.md` required for every script (skill section 3) is elided from the
table; a real stage ships one per file.
