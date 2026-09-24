---
name: pixi-env
description: "Scaffold and manage a pixi-based project environment: pixi.toml, the tracked .pixi/config.toml that bioconductor packages require, package.R for GitHub-only R packages, the `pixi run format` and `pixi run lint` tasks and git hooks, and the optional dual-remote lab mirror. Use when starting a new repository, when adding a dependency and choosing between pixi / package.R / a separate pixi environment / tools/opt, when a bioconductor package installs but does not work, when linting R or checking an R script for syntax errors, when a shell or cluster script needs the environment activated, when `pixi shell-hook` fails under set -u, when every batch or cluster task dies within seconds with no output from the program while an interactive shell works fine, or when a repository must publish only a shareable subset to a second remote."
---

# Pixi environment

Every project here is a pixi project: one `pixi.toml` at the repository root owns the R, Python, and CLI toolchain. This skill scaffolds that root from a working template and holds the rules for growing it afterwards.

This skill is portable. The author name and contact that go into the template headers, and any repository-specific track or storage layout, come from the project's `.github/instructions/` bindings.

---

## 1. Scaffold a new project

One command writes the whole root. It copies files only -- it never installs, and it never overwrites an existing file unless `--force` is given.

```bash
bash .github/skills/pixi-env/scripts/init-pixi-project.sh --help
bash .github/skills/pixi-env/scripts/init-pixi-project.sh ~/github/<project> --minimal
bash .github/skills/pixi-env/scripts/init-pixi-project.sh ~/github/<project> --lab-sync
```

| Option | Effect |
| --- | --- |
| `--name <name>` | `workspace.name`; defaults to the target directory basename |
| `--minimal` | toolchain-only `pixi.toml` instead of the full template |
| `--lab-sync` | also install the dual-remote mirror and its pre-push hook |
| `--no-git` | skip `git init` and `core.hooksPath` |
| `--force` | overwrite files that already exist |
| `--dry-run` | report what would be written, touch nothing |

Then, by hand, from the project directory:

```bash
pixi install              # long; run it in tmux, never poll it with sleep
pixi run remote-install   # GitHub-only R packages from package.R
pixi run format           # air + ruff + shfmt + prettier on changed files
pixi run lint             # jarl on changed R files
```

**`pixi install` on the full template takes a long time to solve.** It is a long-running job: launch it in tmux and report the session, do not sit on it. Prefer `--minimal` when the project's shape is not yet known, and grow the manifest with `pixi add`.

`--dry-run` first whenever the target already has files in it. The script reports `skip` for anything that exists, so a re-run is safe, but the dry run is what shows that before it matters.

### What lands

```
pixi.toml            .pixi/config.toml    package.R
air.toml             ruff.toml            .prettierrc
jarl.toml            .prettierignore      .gitignore
tools/format.sh      tools/lint.sh        tools/hooks/pre-commit
tools/opt/README.md  tools/opt/.gitignore
tools/sync-lab.sh    tools/lab-exclude.gitignore   tools/hooks/pre-push  [--lab-sync]
```

An existing `.gitignore` is never rewritten -- only the two-line pixi block is appended, and only when it is absent.

Per-file detail, what the full template comments out and why, and the profile choice: `references/template-inventory.md`.

### Three things to do after the scaffold

1. **Prune or uncomment `pixi.toml`.** The full template ships the general toolkit live and the heavy domain stacks commented out. Uncomment what the project loads; delete what it never will.
2. **Trim `package.R`** to the GitHub packages this project actually needs.
3. **Commit `pixi.toml`, `pixi.lock`, and `.pixi/config.toml` together.**

---

## 2. `.pixi/config.toml` is not optional

```toml
run-post-link-scripts = "insecure"
```

Bioconductor conda packages ship post-link scripts. Pixi does not run them by default, so without this file every `bioconductor-*` package installs cleanly and then fails at load or at first real call. The failure looks like a broken package, not a missing config, which is what makes it expensive.

The file is therefore **tracked**, which needs the `.gitignore` block to be written exactly this way:

```gitignore
.pixi/*
!.pixi/config.toml
```

A bare `.pixi/` line cannot be un-ignored by a later `!` rule -- git never descends into an excluded directory, so the re-include is dead and the config silently never reaches the clone. Check it, do not assume it:

```bash
git check-ignore -q .pixi/config.toml && echo "IGNORED -- fix .gitignore" || echo "tracked (correct)"
git add -An .pixi/     # must print: add '.pixi/config.toml'
```

Do not pass `-v` to `check-ignore` here. With `-v` the exit code reports "a pattern matched", and the negation `!.pixi/config.toml` is a match, so a correctly configured repository still exits 0 and the check reads backwards.

---

## 3. Adding software: take the first rung that works

| Rung | Where | For |
| --- | --- | --- |
| 1 | `pixi add` -> `pixi.toml` | anything on conda-forge or bioconda: R, Python, CLI binaries |
| 2 | `package.R` | an R package that exists only on GitHub |
| 3 | a separate pixi environment | installable, but its versions conflict with the root env |
| 4 | `tools/opt/install-<x>.sh` | on no conda channel; must be built or vendored |

```bash
pixi search r-susier                 # check before assuming it is missing
pixi add r-ggplot2
pixi add --platform linux-64 gcta    # linux-only binary
pixi run remote-install              # after editing package.R
```

**Never resolve a version conflict by downgrading, pinning, or removing a package the root environment already provides.** One working root env serves every track; no single tool is worth breaking it. Escalate to rung 3.

**Never `install.packages()` or `remotes::install_github()` into a live session**, and never install into `$HOME`, a system path, or a module tree. The next environment rebuild drops it and nothing records that it was ever needed.

Rungs 3 and 4 are exceptions, so they get written down: the separate manifest or the installer script is committed, and the `AGENTS.md` of every track that uses the tool carries the exact command that reaches the binary.

Full procedure, the `tools/opt` installer contract and template, and the conda-name rules: `references/install-ladder.md`.

---

## 4. Entering the environment

**Interactive commands, R and Python scripts** -- run from the directory that owns `pixi.toml`:

```bash
pixi run Rscript src/NN-stage/NN-step.R
pixi run <task>
pixi run --manifest-path /abs/path/pixi.toml Rscript ...   # cwd is elsewhere
pixi run --environment <envname> <cmd>                     # a rung-3 env
```

**Shell and cluster scripts** -- activate once for the whole process instead of prefixing every command. Both guards below are load-bearing:

```bash
_pixi_hook=$(pixi shell-hook --manifest-path "${repodir}/pixi.toml") || {
  echo "FATAL: pixi shell-hook failed" >&2
  return 1 2>/dev/null || exit 1
}
if [[ $- == *u* ]]; then had_u=1; set +u; else had_u=0; fi
eval "$(
  printf '%s\n' "${_pixi_hook}" |
    awk '/^for _pixi_f in .*bash-completion/ { skip = 1 }
         skip { if ($0 == "done") skip = 0; next }
         { print }'
)"
if [[ ${had_u} -eq 1 ]]; then set -u; fi
unset _pixi_hook had_u
```

After `shell-hook`, `Rscript`, `python`, and every other pixi-provided binary are on `PATH`. Do not prefix them with `pixi run` again in that script.

Never use conda / mamba / Miniforge, a hard-coded `.../envs/<env>/bin/Rscript`, or a bare `Rscript` from the login shell -- the last one resolves to the system toolchain and silently misses every project package.

### Three traps

- **The hook's bash-completion loop kills non-interactive shells.** `shell-hook` output ends with `for _pixi_f in .../share/bash-completion/completions/*; do ... done`. Several of those completion files use `< <(...)`, which bash 4.4 -- still the system bash on many compute nodes -- cannot parse, and a syntax error inside a sourced file kills a non-interactive shell outright. Symptom: **every batch task exits 2 within seconds, with no output from the real program**, while an interactive login shell activates fine. Forcing the job shell with `-L /bin/bash` does not help; the job shell is already bash. It appears out of nowhere after a pixi env rebuild adds a completion file. Completions are interactive convenience only, so strip the loop with `awk`, as above. Diagnose with a three-line probe job (`echo` / source the config / `echo`), not by reading code -- one queue slot settles it.
- **`pixi shell-hook` breaks under `set -u`.** Activation scripts dereference variables that have no default. Save and restore `-u` around the `eval` as above, rather than a bare `set +u` / `set -u` pair: the bare form turns `-u` **on** for a caller that never had it.
- **conda activation overrides your exported variables.** `activate.d` scripts run after your `export`, so exporting before `pixi run` is discarded. Set it inside the activated environment: `pixi run env VAR=value <cmd>`.

---

## 5. Formatting and linting

Two entry points, same interface, same git-derived file list: `tools/format.sh` rewrites layout, `tools/lint.sh` reports defects. Neither is a substitute for the other.

### Formatting

`tools/format.sh` is the one entry point, wired to `pixi run format`. It routes by extension -- `*.R`/`*.r` to air, `*.py` to ruff, `*.sh`/`*.lsf` and extensionless files with a shell shebang to shfmt, `*.md` to prettier -- and takes its file list from git, so a run only touches what changed.

```bash
pixi run format                    # files differing from HEAD, plus untracked
pixi run format -- --staged        # staged files, then restage them
pixi run format -- --all           # every tracked file
pixi run format -- path/to/file.R
```

`tools/hooks/pre-commit` runs `--staged` on every commit, so what lands in history is already clean. `--staged` deliberately skips partially staged files: rewriting one and restaging it would sweep its unstaged hunks into the commit.

Exclusions live in `.prettierignore`. Width and indent live in `air.toml`, `ruff.toml`, and `.prettierrc`: all three wrap at 80, so R, Python, and Markdown break at the same column; ruff keeps the Python indent of 4 rather than air's R-convention 2. `.prettierrc` also owns `proseWrap`, which decides whether Markdown prose is rewrapped at all -- `"never"` collapses each paragraph and list item onto a single line, `"preserve"` keeps the source line breaks. Check which one is set before editing any `.md` in the repository.

#### Two ways a format run bites back

- **A bare `pixi run format` is not scoped to your change.** The no-argument form takes every file differing from `HEAD` **plus every untracked file**, so it rewrites whatever another window, another session, or a colleague has in progress. Name the paths whenever the change set is not the whole working tree: `pixi run format -- a.md b.R`.
- **A Markdown file is a different file after prettier ran.** Re-read it before building an edit anchor, whoever ran the formatter. Under `proseWrap: "never"` a paragraph that was five source lines becomes one, so every multi-line anchor taken from an earlier read silently stops matching.

`ruff.toml` excludes `*.md`. Ruff formats Python code blocks inside Markdown, so without that line a bare `ruff format .` would fight prettier over every `.md` in the repository. `tools/format.sh` is safe either way -- it hands ruff an explicit `*.py` list, and ruff does not apply exclusions to files named on the command line.

### Linting R

`tools/lint.sh` is the one entry point, wired to `pixi run lint`. It hands the `*.R` / `*.r` files to `jarl check` and takes the same four argument forms as the formatter.

```bash
pixi run lint                      # R files differing from HEAD, plus untracked
pixi run lint -- --staged          # staged R files
pixi run lint -- --all             # every tracked R file
pixi run lint -- src/01-start      # a directory; jarl walks it
```

**This replaces the `Rscript -e 'parse(...)'` syntax rung.** jarl parses before it lints, so a syntax error is reported here with its exact line **and column** and a caret under the offending token, which `parse()` does not give. Run the linter instead of a separate parse pass; never run both.

Exit status is the whole contract, and the three failure codes mean different things:

| Exit | Meaning                                                   |
| ---- | --------------------------------------------------------- |
| 0    | clean                                                     |
| 1    | lint violations found                                     |
| 255  | syntax error -- the file cannot be parsed, fix this first |
| 2    | wrapper failure (no pixi, no jarl, bad option)            |

Rules live in `jarl.toml` at the repository root, discovered by walking up from each linted file, so the same rules apply from any working directory.

#### Two ways a lint run bites back

- **`--all` walks `git ls-files`, which lists tracked directory symlinks.** A repository that symlinks `data/`, `results/`, `logs/`, or `tmp/` to storage outside the tree tracks those links as ordinary entries, and `[[ -d ]]` follows them, so a wrapper that forwards directories to the tool lints the external tree and reports findings from files nobody in the repository owns. `tools/lint.sh` therefore takes a directory **only** when the user named it explicitly; from a git-derived list a dir-symlink falls through `[[ -f ]]` and is skipped. Keep that guard if you touch the wrapper. `tools/format.sh` never accepts a directory at all, so it is immune.
- **Reach jarl through pixi, always.** A `~/.local/bin/jarl` left by an earlier manual install shadows the pixi binary in a login shell and can be an older release with a different rule set and different flags. `pixi run lint` and `pixi run jarl ...` are the only forms; a bare `jarl` proves nothing about what the project's linter does.

`default-exclude = true` already keeps jarl out of `.pixi/` and the symlinked storage roots, so listing them in `exclude` is redundant.

**Jarl has no severity levels.** Every rule violation prints as `warning: <rule>` and is counted as an error in the summary; the only true `error:` is the parse failure that exits 255. "Errors only, no warnings" is a choice of rule category, not a switch:

| Group      | Meaning                                      | In the template |
| ---------- | -------------------------------------------- | --------------- |
| `CORR`     | outright wrong or useless                    | selected        |
| `SUSP`     | most likely wrong or useless                 | selected        |
| `COMM`     | `# jarl-ignore` comments that do nothing     | dropped         |
| `PERF`     | correct, just slower                         | dropped         |
| `READ`     | correct, just less pretty -- air owns layout | dropped         |
| `DPLYR`    | dplyr-specific                               | off by default  |
| `TESTTHAT` | testthat-specific                            | off by default  |

The template also ignores `unused_object`: Jarl analyses one file at a time, so every constant in a sourced `config.R` and every palette in a `color.R` looks unused. Drop that line in a package repository, where the rule is real signal.

`jarl check --fix` rewrites code, so `tools/lint.sh` does not wire it and the template keeps `fixable` inside `select`, so `--fix` is neither a silent no-op nor free to touch everything. Apply fixes by hand, on a clean tree, and re-read the diff:

```bash
pixi run jarl check --fix <path>   # refuses on a dirty tree unless --allow-dirty
```

Suppress a justified violation by widening `ignore` in `jarl.toml`, with a comment naming the rule and why. **Never put a `# jarl-ignore` comment in an R source.** One file carries the whole policy, so a rule is judged once for the project instead of being argued again at each call site, and `grep` on that one file answers "what are we not checking". `COMM` is not selected, so the linter will not catch a stray in-source suppression: this is a review rule, checked with `grep -rn 'jarl-ignore' --include='*.R' .`

---

## 6. Dual-remote lab mirror (`--lab-sync` only)

For a repository that lives in two places: the **personal** remote holds everything, the **lab** remote holds only the shareable subset. A repository that pushes to one remote does not need any of this.

```bash
export LAB_REPO=~/github/<project>-analysis
pixi run sync-lab -- --dry-run     # always first
pixi run sync-lab                  # mirror HEAD, commit, push
```

`tools/lab-exclude.gitignore` is the only place the policy is written -- never edit `.gitignore` inside the lab repository, the next sync overwrites it. The file excludes by category, then excludes `*.md` wholesale, then re-includes an allowlist of exact paths. Last matching rule wins.

**Verify with `--dry-run` after every edit to that file.** A rule that is too broad silently publishes a private document; one that is too narrow silently deletes a file from the lab remote. The dry run prints both lists.

How the mirror works, its safety rails, and the resolution order for the destination: `references/dual-remote-sync.md`.

---

## 7. Verify before reporting success

A zero exit from the init script proves files were copied, nothing more.

```bash
cd <project>
grep -n '^name =' pixi.toml                    # substituted, not PROJECT_NAME
cat .pixi/config.toml                          # run-post-link-scripts present
git check-ignore -q .pixi/config.toml && echo "IGNORED -- fix .gitignore" || echo "tracked (correct)"
git config --get core.hooksPath                # tools/hooks
bash -n tools/format.sh
bash -n tools/lint.sh
pixi info                                      # after pixi install
```

After `pixi install`, one real call is worth more than the exit code:

```bash
pixi run Rscript -e 'cat(R.version.string, "\n")'
pixi run Rscript -e 'library(<a bioconductor package>)'   # proves post-link ran
pixi run jarl --version                                   # proves lint can run
```
