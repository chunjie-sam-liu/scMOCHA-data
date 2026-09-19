# Dual-remote sync

Some repositories live in two places: the **personal** remote, which holds
everything, and the **lab** remote, which may hold only the shareable subset.
`tools/sync-lab.sh` mirrors one into the other. A repository that pushes to one
remote only does not need any of this -- omit `--lab-sync` at init.

## How it works

1. `git archive` exports the commit being mirrored into a scratch tree under
   the repository's `tmp/` root (`tmpdir` in the env file), falling back to
   `~/tmp/sync-lab/` when the repository has none. The working directory is
   never read, so uncommitted work cannot leak.
2. `tools/lab-exclude.gitignore` is appended to the exported `.gitignore`. The
   lab remote then carries the same policy, so a private file is rejected there
   even if the copy step ever let one through.
3. `git check-ignore` on a throwaway repo over that tree decides which paths are
   shareable. Git answers the question, not a hand-written glob loop.
4. `rsync --files-from` copies the shareable set into the lab repository, files
   that are no longer shareable are `git rm`'d, and the result is committed with
   the original commit message.

Safety rails already in the script: it refuses to run when the lab repository
has uncommitted changes or a detached HEAD, and it refuses to commit when the
shareable set comes out empty rather than wiping the lab repository.

## Setup

```bash
# the lab repository must already exist and be clean
export LAB_REPO=~/github/<project>-analysis          # or pass --lab-repo
git config core.hooksPath tools/hooks                # done by the init script
```

Resolution order for the destination: `--lab-repo <dir>`, then `$LAB_REPO`, then
`~/github/<reponame>-analysis`.

## Use

```bash
pixi run sync-lab -- --dry-run     # always first: shows every file and removal
pixi run sync-lab                  # mirror HEAD, commit, push
bash tools/sync-lab.sh --rev <sha> --no-push
```

With `tools/hooks/pre-push` installed, every push to the personal remote
mirrors the same commit to the lab remote first, and any failure aborts the
push -- the two remotes never drift.

## Editing the policy

`tools/lab-exclude.gitignore` is the only place the policy is written. Never
edit `.gitignore` in the lab repository; the next sync overwrites it.

Last matching rule wins. The file has three blocks:

- **Private by category** -- `.github/`, `AGENTS.md`, editor config, toolchain
  config, `.env*`, `package.R`, and `tools/` itself.
- **All Markdown excluded** -- `*.md` catches paired script docs, plans,
  progress files, and diagram briefs anywhere in the tree, including in
  directories that do not exist yet.
- **Allowlist** -- exact re-includes. Anchor each with a leading slash so it
  names one file:

  ```gitignore
  !README.md
  !/<shared-document>.md
  !/<track>/README.md
  ```

**Verify with `--dry-run` after any edit to this file.** A rule that is too
broad silently publishes a private document; a rule that is too narrow silently
deletes a file from the lab remote. The dry run prints both lists.
