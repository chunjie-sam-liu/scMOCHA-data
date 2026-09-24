#!/usr/bin/env bash
# @AUTHOR: Chunjie Liu
# @CONTACT: chunjie.sam.liu.at.gmail.com
# @DATE: 2026-09-15
# @DESCRIPTION: Scaffold a pixi-managed project: pixi.toml, the tracked
#   .pixi/config.toml that bioconductor packages need, package.R, the air /
#   ruff / shfmt / prettier formatter behind `pixi run format`, the jarl
#   linter behind `pixi run lint`, the git hooks, and tools/opt. Copies files
#   only -- nothing is installed and no existing file is overwritten unless
#   --force is given.
# @VERSION: v0.1.0

set -euo pipefail

die() {
  echo "init-pixi-project: $*" >&2
  exit 1
}

usage() {
  cat <<'EOF'
Usage: init-pixi-project.sh [options] [<target-dir>]

  <target-dir>       Project directory (default: the current directory).
                     Created if missing.

  --name <name>      workspace.name in pixi.toml.
                     Default: the basename of <target-dir>.
  --minimal          Scaffold the toolchain-only pixi.toml instead of the full
                     template. Solves in a fraction of the time; grow it with
                     `pixi add`.
  --lab-sync         Also install the dual-remote mirror: tools/sync-lab.sh,
                     tools/lab-exclude.gitignore, and the pre-push hook.
  --no-git           Do not run `git init` and do not set core.hooksPath.
  --force            Overwrite files that already exist.
  --dry-run          Report what would be written; touch nothing.
  -h, --help         Show this help.

Writes, and nothing else:

  pixi.toml  .pixi/config.toml  package.R  air.toml  ruff.toml
  .prettierrc  .prettierignore  jarl.toml
  .gitignore (created, or the pixi block appended to an existing one)
  tools/format.sh  tools/lint.sh  tools/hooks/pre-commit
  tools/opt/{README.md,.gitignore}
  tools/sync-lab.sh  tools/lab-exclude.gitignore  tools/hooks/pre-push  [--lab-sync]

Then, by hand:

  cd <target-dir> && pixi install && pixi run remote-install
EOF
}

skilldir=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
assets="${skilldir}/assets"
[[ -d ${assets} ]] || die "assets directory not found: ${assets}"

name=""
minimal=0
lab_sync=0
use_git=1
force=0
dry_run=0
target=""

# The name lands in pixi.toml through sed, and pixi rejects exotic names anyway.
check_name() {
  [[ $1 =~ ^[A-Za-z0-9][A-Za-z0-9._-]*$ ]] ||
    die "invalid workspace name: $1 (use letters, digits, . _ -)"
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --name)
      [[ $# -ge 2 ]] || die "--name needs a value"
      name="$2"
      # Checked here, before mkdir, so a rejected name leaves no directory.
      check_name "${name}"
      shift
      ;;
    --minimal) minimal=1 ;;
    --lab-sync) lab_sync=1 ;;
    --no-git) use_git=0 ;;
    --force) force=1 ;;
    --dry-run) dry_run=1 ;;
    -h | --help)
      usage
      exit 0
      ;;
    -*) die "unknown option: $1" ;;
    *)
      [[ -z ${target} ]] || die "only one target directory is accepted"
      target="$1"
      ;;
  esac
  shift
done

target="${target:-.}"
if ((dry_run)); then
  [[ -d ${target} ]] || die "--dry-run needs an existing target directory"
else
  mkdir -p "${target}"
fi
target=$(cd "${target}" && pwd)

name="${name:-$(basename "${target}")}"
check_name "${name}"

written=0
skipped=0

# place <asset-relative-path> <target-relative-path> [mode]
place() {
  local src="${assets}/$1" dst="${target}/$2" mode="${3:-}"
  [[ -f ${src} ]] || die "missing asset: ${src}"
  if [[ -e ${dst} ]] && ((! force)); then
    echo "  skip    $2 (exists; --force to overwrite)"
    skipped=$((skipped + 1))
    return 0
  fi
  if ((dry_run)); then
    echo "  would write $2"
    written=$((written + 1))
    return 0
  fi
  mkdir -p "$(dirname "${dst}")"
  cp "${src}" "${dst}"
  if [[ -n ${mode} ]]; then
    chmod "${mode}" "${dst}"
  fi
  echo "  write   $2"
  written=$((written + 1))
}

echo "init-pixi-project: ${target}"
echo "  workspace name: ${name}"
if ((minimal)); then
  echo "  profile: minimal"
else
  echo "  profile: full"
fi
echo

# --- pixi -------------------------------------------------------------------

manifest="pixi.toml"
if ((minimal)); then
  manifest="pixi.minimal.toml"
fi

place "${manifest}" "pixi.toml"
if [[ -f "${target}/pixi.toml" ]] && ((! dry_run)); then
  if grep -q '^name = "PROJECT_NAME"$' "${target}/pixi.toml"; then
    sed -i.bak "s|^name = \"PROJECT_NAME\"\$|name = \"${name}\"|" "${target}/pixi.toml"
    rm -f "${target}/pixi.toml.bak"
    echo "  set     pixi.toml workspace.name = ${name}"
  fi
fi

place "pixi-config.toml" ".pixi/config.toml"
place "package.R" "package.R"

# --- formatter --------------------------------------------------------------

place "air.toml" "air.toml"
place "ruff.toml" "ruff.toml"
place "prettierrc" ".prettierrc"
place "prettierignore" ".prettierignore"
place "tools/format.sh" "tools/format.sh" 755
place "tools/hooks/pre-commit" "tools/hooks/pre-commit" 755

# --- linter -----------------------------------------------------------------

place "jarl.toml" "jarl.toml"
place "tools/lint.sh" "tools/lint.sh" 755

# --- tools/opt --------------------------------------------------------------

place "tools/opt/README.md" "tools/opt/README.md"
place "tools/opt/opt.gitignore" "tools/opt/.gitignore"

# --- dual-remote mirror -----------------------------------------------------

if ((lab_sync)); then
  place "tools/sync-lab.sh" "tools/sync-lab.sh" 755
  place "tools/lab-exclude.gitignore" "tools/lab-exclude.gitignore"
  place "tools/hooks/pre-push" "tools/hooks/pre-push" 755
fi

# --- .gitignore -------------------------------------------------------------

# An existing .gitignore is never rewritten: it may already carry project rules.
# Only the two pixi lines are appended, and only when they are absent.
if [[ -f "${target}/.gitignore" ]]; then
  if grep -qxF '!.pixi/config.toml' "${target}/.gitignore"; then
    echo "  skip    .gitignore (pixi block already present)"
    skipped=$((skipped + 1))
  elif ((dry_run)); then
    echo "  would append the pixi block to .gitignore"
    written=$((written + 1))
  else
    cat >>"${target}/.gitignore" <<'EOF'

# pixi: the built environment is never committed, but .pixi/config.toml is --
# it carries run-post-link-scripts, without which every bioconductor-* package
# installs and then fails to work.
.pixi/*
!.pixi/config.toml
EOF
    echo "  append  .gitignore (pixi block)"
    written=$((written + 1))
  fi
else
  place "gitignore" ".gitignore"
fi

# --- git --------------------------------------------------------------------

if ((use_git)) && ((! dry_run)); then
  if git -C "${target}" rev-parse --git-dir >/dev/null 2>&1; then
    echo "  git     already a repository"
  else
    git -C "${target}" init -q
    echo "  git     init"
  fi
  git -C "${target}" config core.hooksPath tools/hooks
  echo "  git     core.hooksPath = tools/hooks"
fi

# --- report -----------------------------------------------------------------

echo
echo "init-pixi-project: ${written} written, ${skipped} skipped"
if ((dry_run)); then
  exit 0
fi

cat <<EOF

Next, from ${target}:

  pixi install              # build the environment (long; run it in tmux)
  pixi run remote-install   # GitHub-only R packages from package.R
  pixi run format           # air + ruff + shfmt + prettier on changed files
  pixi run lint             # jarl on changed R files

Commit pixi.toml, pixi.lock, and .pixi/config.toml together. A clone without
.pixi/config.toml gets bioconductor packages that install but do not work.
EOF
