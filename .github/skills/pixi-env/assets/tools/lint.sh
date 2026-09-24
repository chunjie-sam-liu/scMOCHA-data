#!/usr/bin/env bash
# @AUTHOR: Chunjie Liu
# @CONTACT: chunjie.sam.liu.at.gmail.com
# @DATE: 2026-09-22
# @DESCRIPTION: Lint R sources with jarl, using the rules in jarl.toml at the
#   repository root. File lists come from git, so a run only covers what
#   changed. jarl parses before it lints, so this also replaces a separate
#   Rscript parse() pass.
# @VERSION: v0.1.0

set -euo pipefail

# Exit 2, not 1: jarl already uses 1 for "violations found".
die() {
  echo "lint: $*" >&2
  exit 2
}

usage() {
  cat <<'EOF'
Usage: tools/lint.sh [--staged | --all | <path>...]

  (no arguments)  Lint R files that differ from HEAD, plus untracked files.
  --staged        Lint staged R files (pre-commit form).
  --all           Lint every tracked R file.
  <path>...       Lint the given files; a directory is handed to jarl, which
                  walks it itself.
  -h, --help      Show this help.

*.R and *.r -> jarl check. Rules live in jarl.toml at the repository root.
jarl parses before it lints, so a syntax error is reported here with its exact
line and column; no separate parse pass is needed.

Exit status: 0 clean, 1 lint violations, 255 syntax error, 2 wrapper failure.

`jarl check --fix` rewrites code, so it is deliberately not wired here. Run it
by hand on a clean tree: pixi run jarl check --fix <path>
EOF
}

repodir=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)

mode="changed"
paths=()

while [[ $# -gt 0 ]]; do
  case "$1" in
    --staged) mode="staged" ;;
    --all) mode="all" ;;
    -h | --help)
      usage
      exit 0
      ;;
    -*) die "unknown option: $1" ;;
    *)
      mode="explicit"
      paths+=("$1")
      ;;
  esac
  shift
done

cd "${repodir}"

case "${mode}" in
  staged) mapfile -t paths < <(git diff --cached --name-only --diff-filter=ACMR) ;;
  all) mapfile -t paths < <(git ls-files) ;;
  changed)
    mapfile -t paths < <(
      {
        git diff --name-only --diff-filter=ACMR HEAD
        git ls-files --others --exclude-standard
      } | sort -u
    )
    ;;
  explicit) ;;
esac

# A mistyped path must not come back as "nothing to check, all good".
if [[ ${mode} == "explicit" ]]; then
  for f in "${paths[@]}"; do
    [[ -e ${f} ]] || die "no such path: ${f}"
  done
fi

targets=()
for f in "${paths[@]}"; do
  # Only an explicitly named path may be a directory. git also lists the tracked
  # data/results/logs/tmp symlinks, and -d follows them into storage outside the
  # repository.
  if [[ ${mode} == "explicit" && -d ${f} ]]; then
    targets+=("${f}")
    continue
  fi
  [[ -f ${f} ]] || continue
  case "${f##*/}" in
    *.R | *.r) targets+=("${f}") ;;
  esac
done

if ((${#targets[@]} == 0)); then
  echo "lint: no R files to check"
  exit 0
fi

command -v pixi >/dev/null 2>&1 || die "pixi not on PATH; cannot reach jarl"
# Strip the hook's bash-completion loop: those files use `< <(...)`, which bash
# 4.4 cannot parse, and a syntax error in a sourced file kills the shell.
_pixi_hook=$(pixi shell-hook --manifest-path "${repodir}/pixi.toml") || die "pixi shell-hook failed"
set +u
eval "$(
  printf '%s\n' "${_pixi_hook}" |
    awk '/^for _pixi_f in .*bash-completion/ { skip = 1 }
         skip { if ($0 == "done") skip = 0; next }
         { print }'
)"
set -u
unset _pixi_hook

command -v jarl >/dev/null 2>&1 || die "jarl not in the pixi environment; add it with: pixi add jarl"

status=0
jarl check "${targets[@]}" || status=$?

echo "lint: ${#targets[@]} target(s) checked"
exit "${status}"
