#!/usr/bin/env bash
# @AUTHOR: Chunjie Liu
# @CONTACT: chunjie.sam.liu.at.gmail.com
# @DATE: 2026-09-02
# @DESCRIPTION: Format sources in place -- R with air, Python with ruff, shell
#   and LSF scripts with shfmt, Markdown with prettier. File lists come from
#   git, so a run only touches what changed; --staged also restages what it
#   rewrote.
# @VERSION: v0.1.0

set -euo pipefail

die() {
  echo "format: $*" >&2
  exit 1
}

usage() {
  cat <<'EOF'
Usage: tools/format.sh [--staged | --all | <path>...]

  (no arguments)  Format files that differ from HEAD, plus untracked files.
  --staged        Format staged files, then restage them (pre-commit hook).
  --all           Format every tracked file.
  <path>...       Format the given files.
  -h, --help      Show this help.

*.R and *.r     -> air
*.py            -> ruff format
*.sh, *.lsf, and extensionless files with a shell shebang -> shfmt
*.md            -> prettier (exclusions live in .prettierignore)
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

# Rewriting a partially staged file and restaging it would sweep its unstaged
# hunks into the commit, so leave those files untouched.
if [[ ${mode} == "staged" ]]; then
  unstaged=$(git diff --name-only)
  kept=()
  for f in "${paths[@]}"; do
    if [[ -n ${unstaged} ]] && grep -qxF -- "${f}" <<<"${unstaged}"; then
      echo "format: skipped, partially staged: ${f}" >&2
    else
      kept+=("${f}")
    fi
  done
  paths=("${kept[@]}")
fi

has_shell_shebang() {
  local first
  read -r first <"$1" || return 1
  [[ ${first} =~ ^#!.*(^|/|[[:space:]])(ba|da|k|z)?sh([[:space:]]|$) ]]
}

r_files=()
py_files=()
sh_files=()
lsf_files=()
md_files=()

for f in "${paths[@]}"; do
  [[ -f ${f} ]] || continue
  case "${f##*/}" in
    *.R | *.r) r_files+=("${f}") ;;
    *.py) py_files+=("${f}") ;;
    *.sh) sh_files+=("${f}") ;;
    *.lsf) lsf_files+=("${f}") ;;
    *.md) md_files+=("${f}") ;;
    *.*) ;;
    *) has_shell_shebang "${f}" && sh_files+=("${f}") ;;
  esac
done

targets=("${r_files[@]}" "${py_files[@]}" "${sh_files[@]}" "${lsf_files[@]}" "${md_files[@]}")
if ((${#targets[@]} == 0)); then
  echo "format: nothing to format"
  exit 0
fi

command -v pixi >/dev/null 2>&1 || die "pixi not on PATH; cannot reach the formatters"
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

if ((${#r_files[@]})); then
  air format "${r_files[@]}"
fi
# Format only; `ruff check --fix` would rewrite code, not layout.
if ((${#py_files[@]})); then
  ruff format "${py_files[@]}"
fi
if ((${#sh_files[@]})); then
  shfmt -i 2 -ci -w "${sh_files[@]}"
fi
# LSF scripts open with #BSUB directives and carry no shebang, so name the
# dialect instead of letting shfmt guess.
if ((${#lsf_files[@]})); then
  shfmt -i 2 -ci -ln bash -w "${lsf_files[@]}"
fi
if ((${#md_files[@]})); then
  prettier --log-level warn --write "${md_files[@]}"
fi

if [[ ${mode} == "staged" ]]; then
  git add -- "${targets[@]}"
fi

echo "format: ${#targets[@]} file(s) processed"
