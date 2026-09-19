# Template inventory

What each asset is, where it lands, and what to change after the copy.

| Asset                         | Lands as                      | Tracked | Edit after init                   |
| ----------------------------- | ----------------------------- | ------- | --------------------------------- |
| `pixi.toml`                   | `pixi.toml`                   | yes     | name is substituted; prune extras |
| `pixi.minimal.toml`           | `pixi.toml` (`--minimal`)     | yes     | grow it with `pixi add`           |
| `pixi-config.toml`            | `.pixi/config.toml`           | yes     | no                                |
| `package.R`                   | `package.R`                   | yes     | uncomment the domain packages     |
| `air.toml`                    | `air.toml`                    | yes     | no                                |
| `ruff.toml`                   | `ruff.toml`                   | yes     | `target-version` per project      |
| `prettierignore`              | `.prettierignore`             | yes     | add generated trees               |
| `gitignore`                   | `.gitignore`                  | yes     | add project-specific rules        |
| `tools/format.sh`             | `tools/format.sh`             | yes     | no                                |
| `tools/hooks/pre-commit`      | `tools/hooks/pre-commit`      | yes     | no                                |
| `tools/opt/README.md`         | `tools/opt/README.md`         | yes     | one row per installed tool        |
| `tools/opt/opt.gitignore`     | `tools/opt/.gitignore`        | yes     | no                                |
| `tools/sync-lab.sh`           | `tools/sync-lab.sh`           | yes     | no; set `LAB_REPO` instead        |
| `tools/lab-exclude.gitignore` | `tools/lab-exclude.gitignore` | yes     | the Markdown allowlist            |
| `tools/hooks/pre-push`        | `tools/hooks/pre-push`        | yes     | no                                |

The last three are written only with `--lab-sync`.

Assets whose target name starts with a dot are stored without it
(`gitignore`, `prettierignore`, `opt.gitignore`) so git does not apply them
inside the skill directory. The init script renames on copy.

## Full profile: what is commented out, and why

The full `pixi.toml` keeps the whole general toolkit -- tidyverse, ggplot2 and
its extensions, arrow/duckdb/parquet I/O, the modelling stack, the core
Bioconductor infrastructure, plink/plink2/bcftools/htslib, and the air / ruff /
shfmt / prettier formatters. Commented out, each with a one-line reason in the
file:

| Block                                                    | Reason                                     |
| -------------------------------------------------------- | ------------------------------------------ |
| modelling / teaching extras (`rstanarm`, `nycflights13`) | rarely needed outside a tutorial           |
| single-cell stack (`seurat`, `signac`, `azimuth`)        | heavy solve                                |
| genome and example datasets (`BSgenome`, `EnsDb`)        | large downloads                            |
| `perl` + `perl-bioperl`                                  | only for scripts that shell out to Perl    |
| `ensembl-vep`                                            | ~20 GB cache on first use                  |
| Illumina EPIC methylation stack                          | methylation projects only                  |
| `[feature.tensorqtl*]` + `[environments]`                | worked example of rung 3, not a dependency |

Uncomment what the project actually loads. Deleting a block is fine too -- the
template is a starting point, not a contract.

## Choosing a profile

`--minimal` is the right default for a new repository whose shape is not yet
known: it solves quickly and `pixi add` grows it one package at a time, so
`pixi.toml` ends up describing what the project really uses. Use the full
profile when the project is another analysis in the same family as the one the
template came from.

## Verify after init

```bash
cd <project>
test -f .pixi/config.toml && cat .pixi/config.toml
grep -n '^name =' pixi.toml
git -C . config --get core.hooksPath
git check-ignore -q .pixi/config.toml && echo "IGNORED -- fix .gitignore" || echo "tracked (correct)"
git add -An .pixi/      # must print: add '.pixi/config.toml'
bash -n tools/format.sh
```

The most common scaffolding failure is a `.gitignore` carrying a bare `.pixi/`
line, which cannot be un-ignored by a later `!` rule because git never descends
into an excluded directory. The block written here uses `.pixi/*` for exactly
that reason.

**Do not add `-v` to `check-ignore` in that test.** With `-v` the exit code
means "some pattern matched", and the negation `!.pixi/config.toml` counts as a
match, so a correct repository exits 0 and the check reads backwards. Without
`-v`, exit 1 means not ignored, which is the pass condition. `git add -An` is
the unambiguous version: it names the file that would be staged.
