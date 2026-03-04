---
name: run-test
description: Agent skill for testing R/Python scripts locally. Use when asked to "test script X", "run script X", "check script X", or "debug script X". Creates a progress log, runs the script, monitors the log output, notes errors and solutions in progress.md, and writes a final report section.
---

# Run-Test Agent Skill

Run and debug a script locally. Two files per session:

| File                    | Purpose                                                    |
| ----------------------- | ---------------------------------------------------------- |
| `{session}.progress.md` | Tasks, agent observations, final report (agent-maintained) |
| `{session}.log.md`      | Raw stdout/stderr from script execution (written by `tee`) |

The agent monitors `log.md` and writes key observations into `progress.md`.

## When to Use

- "test/run/check/fix script X"
- Debug a failing script
- Validate a script before wider use

## Environment

### R Scripts
- **Skill**: Always activate the `jutils` skill when working with R scripts
- **renv**: `renv` is here `/scr1/users/liuc9/tools/miniforge3/envs/renv`
- **Command**: Source conda, activate by absolute path, then run: `source /scr1/users/liuc9/tools/miniforge3/etc/profile.d/conda.sh && conda activate /scr1/users/liuc9/tools/miniforge3/envs/renv && Rscript script.R`
- **Never use** `conda run -n renv` — it may fail to find the environment. Always use `conda activate` with the full absolute path.
- **Note**: Never use `--vanilla` — it skips `.Rprofile` which defines `load_pkg()` and activates renv
- **`load_pkg()`**: Custom helper defined in `~/.Rprofile`. Loads multiple packages at once and installs any that are missing. Use it instead of multiple `library()` calls. Example: `load_pkg(dplyr, ggplot2, data.table)`

### Python Scripts
- **Command**: `uv run script.py`

## Session Manager Commands

```bash
# Create session — prints file paths, run command, and monitor command
python .opencode/skills/run-test/scripts/session.py init path/to/script.R --model {model}

# Append agent observation to progress.md (after reading log.md)
python .opencode/skills/run-test/scripts/session.py log SESSION_ID "message"

# Tick off a task (fuzzy match on task text)
python .opencode/skills/run-test/scripts/session.py check SESSION_ID "run script"

# Write final report section
python .opencode/skills/run-test/scripts/session.py finish SESSION_ID success|failed

# Watch log.md in background and write periodic progress notes to progress.md
python .opencode/skills/run-test/scripts/session.py monitor SESSION_ID [--interval 300] [--timeout 10800]

# List / inspect sessions
python .opencode/skills/run-test/scripts/session.py list
python .opencode/skills/run-test/scripts/session.py status SESSION_ID
```

---

## Long-Running Scripts (>10 minutes)

The Bash tool has a maximum 10-minute timeout. For scripts expected to take longer, use two background tasks running in parallel.

### Workflow

**Step 1: Init** (same as usual — prints both commands)
```bash
python .opencode/skills/run-test/scripts/session.py init src/01-dataset.R --model {model}
```

**Step 2: Start script in background** — use `run_in_background=true` in Bash tool
```bash
source /scr1/users/liuc9/tools/miniforge3/etc/profile.d/conda.sh && conda activate /scr1/users/liuc9/tools/miniforge3/envs/renv && Rscript src/01-dataset.R 2>&1 | tee docs/logs/SESSION.log.md
```

**Step 3: Start monitor in background** — use `run_in_background=true` in Bash tool
```bash
python .opencode/skills/run-test/scripts/session.py monitor SESSION_ID
```

The monitor writes a timestamped snapshot of recent log output to `progress.md` every 5 minutes. It exits automatically when it detects the `SCRIPT_DONE:exit=N` sentinel line written after the script finishes.

**Step 4: Check progress anytime**

Read `{session}.progress.md` using the Read tool — it accumulates observations from the monitor.

**Step 5: On completion**

Both background tasks notify the agent. Then read the final log, verify output, and wrap up:
```bash
python .opencode/skills/run-test/scripts/session.py check SESSION_ID "run script"
python .opencode/skills/run-test/scripts/session.py check SESSION_ID "verify output"
python .opencode/skills/run-test/scripts/session.py finish SESSION_ID success|failed
```

### Monitor options
```bash
python session.py monitor SESSION_ID --interval 120   # check every 2 min instead of 5
python session.py monitor SESSION_ID --timeout 14400  # allow up to 4 hours (default: 3)
```

---

## Workflow

### 1. Initialize

```bash
python .opencode/skills/run-test/scripts/session.py init src/01-dataset.R --model {model}
# Prints:
#   Session ID : 2026-02-25-01-dataset
#   progress   : docs/logs/2026-02-25-01-dataset.progress.md
#   log        : docs/logs/2026-02-25-01-dataset.log.md
#   Run command: source /scr1/users/liuc9/tools/miniforge3/etc/profile.d/conda.sh && conda activate /scr1/users/liuc9/tools/miniforge3/envs/renv && Rscript src/01-dataset.R 2>&1 | tee docs/logs/2026-02-25-01-dataset.log.md
```

### 2. Execute — pipe output to log.md

Use the `Run command` printed by `init` exactly (includes a sentinel line for monitor completion detection):

```bash
{ source /scr1/users/liuc9/tools/miniforge3/etc/profile.d/conda.sh && conda activate /scr1/users/liuc9/tools/miniforge3/envs/renv && Rscript src/01-dataset.R 2>&1; echo "SCRIPT_DONE:exit=$?"; } | tee docs/logs/{session}.log.md
```

This writes all output to `log.md` while also showing it in the terminal.

### 3. Monitor log.md and update progress.md

While the script runs (or after it exits), read `log.md` and write key observations to `progress.md`:

```bash
# Write observation from what you saw in log.md
python .opencode/skills/run-test/scripts/session.py log SESSION_ID "Script loaded 3 samples, processing chr1..."
python .opencode/skills/run-test/scripts/session.py log SESSION_ID "Error on line 42: object 'meta' not found"
```

### 4. On error — stop and mark failed

Log the error cause from `log.md`, then finish as failed:

```bash
python .opencode/skills/run-test/scripts/session.py log SESSION_ID "Error: <cause from log.md>"
python .opencode/skills/run-test/scripts/session.py finish SESSION_ID failed
```

Fill in the Final Report section in `progress.md` explaining what failed and why.

### 5. On success — verify output and finish

```bash
python .opencode/skills/run-test/scripts/session.py check SESSION_ID "run script"
python .opencode/skills/run-test/scripts/session.py check SESSION_ID "verify output"
python .opencode/skills/run-test/scripts/session.py finish SESSION_ID success
```

Fill in the Final Report: duration, output files, any notes.

---

## Agent Checklist

- [ ] Identify script path and type (R/Python)
- [ ] Run `session.py init` — note the session ID, progress path, and log path
- [ ] Run script with `tee` to `{session}.log.md`
- [ ] Monitor `log.md`; write key observations to `progress.md` via `session.py log`
- [ ] On error: log cause, run `finish failed`, stop
- [ ] On success: verify output files, run `finish success`
- [ ] Complete the Final Report section in `progress.md`
