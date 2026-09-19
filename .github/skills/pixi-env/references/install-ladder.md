# Install ladder

Every tool a script calls is installed by the project, so a fresh clone on a
fresh machine rebuilds it. Take the first rung that works, and never skip down
the list to save time.

| Rung | Where                      | For                                                              |
| ---- | -------------------------- | ---------------------------------------------------------------- |
| 1    | `pixi.toml`                | anything on conda-forge or bioconda: R, Python, CLI binaries     |
| 2    | `package.R`                | an R package that exists only on GitHub                          |
| 3    | `tools/<tool>/pixi.toml`   | installable, but its versions conflict with the root environment |
| 4    | `tools/opt/install-<x>.sh` | on no conda channel; must be built or vendored                   |

**Never resolve a version conflict by downgrading, pinning, or removing a
package the root environment already provides.** One working root env serves
every track. Escalate to rung 3 instead.

---

## Rung 1 -- pixi

```bash
pixi add r-ggplot2                 # R package
pixi add --platform linux-64 gcta  # linux-only binary
pixi add "bcftools>=1.23,<2"       # pinned
pixi add --pypi some-package       # PyPI-only Python package
```

`pixi add` edits `pixi.toml` and `pixi.lock` together. Commit both.

Search before assuming a package is missing:

```bash
pixi search r-susier
pixi search bioconductor-minfi
```

The conda name for a CRAN package is `r-<lowercased name>`; for Bioconductor it
is `bioconductor-<lowercased name>`. A dot in the R name survives lowercasing
(`data.table` -> `r-data.table`) and then needs quoting in `pixi.toml`.

### Bioconductor needs `.pixi/config.toml`

```toml
run-post-link-scripts = "insecure"
```

Bioconductor conda packages ship post-link scripts. Pixi does not run them by
default, so without this file those packages install and then fail to work. The
file is tracked: the root `.gitignore` excludes `.pixi/*` and re-includes
`.pixi/config.toml`.

---

## Rung 2 -- `package.R`

For an R package with no conda build, add one line to `package.R` and run the
task. `install_if_missing()` probes with `requireNamespace()` first, so a
re-run never reinstalls a working package.

```r
install_if_missing("chunjie-sam-liu/jutils")
install_if_missing("duckdb/duckdb-r", pkg = "duckdb")   # repo name != package name
```

```bash
pixi run remote-install
```

Never `install.packages()` or `remotes::install_github()` into a live session:
the next environment rebuild silently drops it and nothing records that it was
ever needed.

---

## Rung 3 -- a separate pixi environment

When a tool's constraints would force a downgrade of the root environment, give
it its own feature and environment in the same `pixi.toml` (see the commented
tensorQTL block in the template), or its own manifest under `tools/<tool>/`:

```bash
pixi run --environment <envname> <cmd>
pixi run --manifest-path tools/<tool>/pixi.toml <cmd>
```

Record the conflict that forced the split in a comment above the block, or in
`tools/<tool>/README.md`, and put the exact call command in the `AGENTS.md` of
every track that uses the tool.

---

## Rung 4 -- `tools/opt/`

Only when no conda channel ships the tool. `tools/opt/.gitignore` is an
allowlist: installers and `README.md` are tracked, built artifacts are ignored.

Installer contract:

- named `install-<tool>.sh`, tracked, and idempotent
- writes only inside `tools/opt/`; never `$HOME`, a system path, or a module
- scratch goes to the repository `tmp/` root, never `/tmp`
- verifies the built binary before exiting
- a runtime that **is** on a conda channel still goes in `pixi.toml`; only the
  unavailable piece lives here

```bash
#!/usr/bin/env bash
# @AUTHOR: <Name>
# @CONTACT: <email>
# @DATE: {today}
# @DESCRIPTION: Install <tool> into tools/opt. <tool> is on no conda channel,
#               so it cannot live in pixi.toml; only its runtime does.
# @VERSION: v0.1.0

set -euo pipefail

optdir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
repodir="$(cd "${optdir}/../.." && pwd)"
prefix="${optdir}/<tool>"
tmproot="${repodir}/tmp/<tool>-build"

mkdir -p "${prefix}" "${tmproot}"

# Build or vendor here, calling compilers through the pixi environment:
#   pixi run --manifest-path "${repodir}/pixi.toml" env TMPDIR="${tmproot}" make ...

bin="${prefix}/bin/<tool>"
[[ -x ${bin} ]] || {
  echo "ERROR: expected an executable at ${bin}" >&2
  exit 1
}

echo "==> Verifying"
"${bin}" --version
```

Then add a row to `tools/opt/README.md` and the call command to the consuming
track's `AGENTS.md`.

### Three traps worth knowing

- **conda activation overrides your exported variables.** `activate.d` scripts
  run after your `export`, so exporting before `pixi run` is discarded. Set it
  inside the activated environment instead: `pixi run env VAR=value <cmd>`.
- **The hook's bash-completion loop kills non-interactive shells.** Its output
  ends with a loop that sources every file in
  `share/bash-completion/completions`; several use `< <(...)`, which bash 4.4 --
  still the system bash on many compute nodes -- cannot parse, and a syntax
  error in a sourced file kills a non-interactive shell. Every batch task then
  exits 2 within seconds with no program output, while a login shell is fine.
- **`pixi shell-hook` breaks under `set -u`.** Activation scripts dereference
  variables with no default.

Both hook traps are handled by the same block, and neither guard is optional:

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
