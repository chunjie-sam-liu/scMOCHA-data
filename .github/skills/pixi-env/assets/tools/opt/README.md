# tools/opt

Software that **cannot** be installed by pixi. Everything available on
conda-forge or bioconda belongs in the repository `pixi.toml` instead, and every
R package that exists only on GitHub belongs in `package.R`.

| Path                | What                                        |
| ------------------- | ------------------------------------------- |
| `install-<tool>.sh` | Idempotent installer for one tool. Tracked. |
| `<tool>/`           | Whatever the installer builds. Git-ignored. |

`.gitignore` here is an allowlist: the installers and this README are tracked,
every built artifact is ignored, so the directory rebuilds from scratch on a
new machine.

## Rules

- One installer per tool, named `install-<tool>.sh`, idempotent, and safe to
  re-run.
- The installer writes only inside `tools/opt/`. Never into `$HOME`, a system
  path, or a module tree.
- Scratch goes to the repository `tmp/` root, never `/tmp`.
- A runtime the tool needs that **is** on a conda channel still goes in
  `pixi.toml`; only the unavailable piece lives here.
- Add a row above and record the exact command that reaches the binary in the
  `AGENTS.md` of every track that uses it. A hand-built tool nobody wrote down
  is indistinguishable from a broken environment on the next machine.

```bash
bash tools/opt/install-<tool>.sh    # install or reinstall
```
