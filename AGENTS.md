# Repository Guidelines

## Project Structure & Module Organization
- `src/`: numbered R/Python/Bash pipeline steps (SRA metadata, downloads, variant collection, parsing). Keep numeric prefixes when adding new stages so order is clear.
- `preprocessing/`: dataset-specific prep, QC, and plotting scripts; includes helper shells for muscle alignment and per-study checks.
- `analysis/`, `large-scale/`, `misc/`: exploratory analyses, figures, and one-off experiments; store derived plots in `imgs/` and keep logs in `logs/`.
- `config/` and `scmocha.yaml`: run-time settings and templates; `utils/utils.R` holds shared helpers; `data/` and `bak/` hold local inputs/backups (avoid committing large or sensitive files).

## Build, Test, and Development Commands
- Run R stages directly, e.g. `Rscript src/01-sra-metadata.R --help` or with arguments for new data pulls; mirror numbering when chaining steps.
- Dataset prep examples: `Rscript preprocessing/02-load-variant.R` and `Rscript preprocessing/03-summarize-meta.R`.
- Mixed-language steps: `python3 src/06.2-collect-variants-new.py`, `bash src/06.3-collect-variants-new.sh`. Keep scripts executable (`chmod +x`) and log outputs to `logs/`.
- `task` is available but currently only echoes; prefer direct commands above or add new Taskfile entries if you standardize a workflow.

## Coding Style & Naming Conventions
- R: 2-space indentation, snake_case for objects, and numbered filenames for pipeline order. Favor `data.table` IO (`fwrite`, `fread`), `qs::qsave/qread` for RDS, and keep pipe-heavy chains readable. Use `logger` for progress and `GetoptLong` for CLI args.
- Python: target Python 3.10; keep functions snake_case, modules lowercase, and run `ruff` (optional dev dependency) before committing. Shell scripts should be POSIX-compatible with `set -euo pipefail` where feasible.

## Testing Guidelines
- No formal test suite yet; validate changes by running the affected script on a small sample and comparing outputs in `data/` or derived tables/plots in `imgs/`.
- Keep runs reproducible: pin random seeds when present, persist command-line arguments in commit messages or PR notes, and capture key log lines for reviewer context.

## Commit & Pull Request Guidelines
- Commit messages are short and conventional (`feat: ...`, `fix: ...`, `docs: ...`, `chore: ...`); use imperative tone and scope-specific wording.
- PRs should summarize intent, list commands executed, note data sources or paths touched, and attach representative artifacts (plots/tables). Link related issues or tickets and call out breaking changes or long-running steps.

## Data Handling & Security
- Do not commit large raw datasets or credentials. Use symlinks for external storage when needed (see prior `feat: add symlink...` patterns).
- Strip identifiers before sharing logs; prefer configuration files (`config/`, `*.yaml`) over hardcoded paths or secrets. Ensure output directories exist before writing to avoid leaking to unintended locations.
