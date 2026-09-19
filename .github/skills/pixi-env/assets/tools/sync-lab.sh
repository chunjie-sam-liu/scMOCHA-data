#!/usr/bin/env bash
# @AUTHOR: Chunjie Liu
# @CONTACT: chunjie.sam.liu.at.gmail.com
# @DATE: 2026-08-31
# @DESCRIPTION: Mirror the shareable subset of this repository into the lab
#   repository, then commit and push it there. The shareable set is whatever
#   survives tools/lab-exclude.gitignore, which is appended to the mirrored
#   .gitignore so the lab remote enforces the same policy.
# @VERSION: v0.1.0

set -euo pipefail

die() {
  echo "sync-lab: $*" >&2
  exit 1
}

usage() {
  cat <<'EOF'
Usage: tools/sync-lab.sh [options]

  --rev <commit>    Commit to mirror (default: HEAD).
  -m, --message     Commit message for the lab repository
                    (default: the message of the mirrored commit).
  --lab-repo <dir>  Lab repository to mirror into.
  --no-push         Commit in the lab repository but do not push.
  --dry-run         Report what would change; touch nothing. Implies --no-push.
  -h, --help        Show this help.

The lab repository is resolved in this order:
  --lab-repo <dir>, then $LAB_REPO, then ~/github/<reponame>-analysis.
EOF
}

repodir=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
excludefile="${repodir}/tools/lab-exclude.gitignore"
labdir=${LAB_REPO:-${HOME}/github/$(basename "${repodir}")-analysis}

rev="HEAD"
message=""
do_push=1
dry_run=0

while [[ $# -gt 0 ]]; do
  case "$1" in
    --rev)
      [[ $# -ge 2 ]] || die "--rev needs a value"
      rev="$2"
      shift
      ;;
    -m | --message)
      [[ $# -ge 2 ]] || die "--message needs a value"
      message="$2"
      shift
      ;;
    --lab-repo)
      [[ $# -ge 2 ]] || die "--lab-repo needs a value"
      labdir="$2"
      shift
      ;;
    --no-push) do_push=0 ;;
    --dry-run)
      dry_run=1
      do_push=0
      ;;
    -h | --help)
      usage
      exit 0
      ;;
    *) die "unknown option: $1" ;;
  esac
  shift
done

# --- preconditions ----------------------------------------------------------

[[ -f ${excludefile} ]] || die "missing exclude policy: ${excludefile}"
[[ -e "${labdir}/.git" ]] || die "lab repository not found: ${labdir}"

revsha=$(git -C "${repodir}" rev-parse --verify "${rev}^{commit}" 2>/dev/null) ||
  die "cannot resolve ${rev} in ${repodir}"

git -C "${labdir}" symbolic-ref -q HEAD >/dev/null ||
  die "lab repository is in detached HEAD state: ${labdir}"

if [[ -n $(git -C "${labdir}" status --porcelain) ]]; then
  die "lab repository has uncommitted changes; resolve them first: ${labdir}"
fi

# --- export the mirrored commit --------------------------------------------

tmproot="${tmpdir:-${HOME}/tmp}/sync-lab"
mkdir -p "${tmproot}"
work=$(mktemp -d "${tmproot}/run.XXXXXXXX")
trap 'rm -rf "${work}"' EXIT

tree="${work}/tree"
lists="${work}/lists"
mkdir -p "${tree}" "${lists}"

git -C "${repodir}" archive --format=tar "${revsha}" | tar -x -C "${tree}"

# The lab .gitignore is base + policy, so the remote rejects private files even
# if one ever slips through the copy step.
printf '\n' >>"${tree}/.gitignore"
cat "${excludefile}" >>"${tree}/.gitignore"

# --- decide what is shareable ----------------------------------------------

# Ask git itself which paths the combined .gitignore excludes. A throwaway repo
# on the exported tree keeps this read-only with respect to the lab repository,
# and core.excludesFile=/dev/null keeps any global ignore file out of the answer.
git -C "${tree}" init -q
git -C "${repodir}" ls-tree -r --name-only "${revsha}" >"${lists}/all"
git -C "${tree}" -c core.excludesFile=/dev/null \
  check-ignore --stdin --no-index <"${lists}/all" >"${lists}/ignored" || true
grep -vxF -f "${lists}/ignored" "${lists}/all" >"${lists}/want" || true
rm -rf "${tree}/.git"

[[ -s "${lists}/want" ]] ||
  die "no shareable files after filtering; refusing to wipe ${labdir}"

git -C "${labdir}" ls-files >"${lists}/have"
grep -vxF -f "${lists}/want" "${lists}/have" >"${lists}/gone" || true

# --- apply ------------------------------------------------------------------

if ((dry_run)); then
  echo "sync-lab: would mirror $(wc -l <"${lists}/want") files from ${revsha:0:8} into ${labdir}"
  rsync -a --itemize-changes --dry-run --files-from="${lists}/want" "${tree}/" "${labdir}/"
  if [[ -s "${lists}/gone" ]]; then
    echo "sync-lab: would remove:"
    sed 's/^/  - /' "${lists}/gone"
  fi
  exit 0
fi

rsync -a --files-from="${lists}/want" "${tree}/" "${labdir}/"

if [[ -s "${lists}/gone" ]]; then
  while IFS= read -r path; do
    [[ -n ${path} ]] || continue
    git -C "${labdir}" rm -q --ignore-unmatch -- "${path}"
  done <"${lists}/gone"
fi

git -C "${labdir}" add -A

if git -C "${labdir}" diff --cached --quiet; then
  echo "sync-lab: lab repository already matches ${revsha:0:8}"
else
  if [[ -n ${message} ]]; then
    git -C "${labdir}" commit -q -m "${message}"
  else
    git -C "${repodir}" log -1 --format=%B "${revsha}" |
      git -C "${labdir}" commit -q -F -
  fi
  echo "sync-lab: committed $(git -C "${labdir}" rev-parse --short HEAD) in ${labdir}"
fi

if ((do_push)); then
  git -C "${labdir}" rev-parse --abbrev-ref --symbolic-full-name '@{u}' >/dev/null 2>&1 ||
    die "lab branch has no upstream; run: git -C ${labdir} push -u origin $(git -C "${labdir}" rev-parse --abbrev-ref HEAD)"
  git -C "${labdir}" push
fi
