#!/usr/bin/env python3
"""
Run-Test Session Manager

Two files per session:
  {session}.progress.md  — tasks, agent observations, final report
  {session}.log.md       — raw script stdout/stderr (written by tee)

Usage:
    python session.py init path/to/script.R [--model claude]
    python session.py log SESSION_ID "observation from reading log.md"
    python session.py check SESSION_ID "task text"
    python session.py finish SESSION_ID success|failed
    python session.py list
    python session.py status SESSION_ID
"""

from __future__ import annotations

import argparse
import re
import time as _time
from datetime import datetime
from pathlib import Path

SKILL_DIR = Path(__file__).resolve().parent.parent
PROGRESS_TEMPLATE = SKILL_DIR / "templates" / "progress.md"


def project_root() -> Path:
    cur = Path.cwd()
    while cur != cur.parent:
        if (cur / ".git").exists() or (cur / "pyproject.toml").exists():
            return cur
        cur = cur.parent
    return Path.cwd()


def logs_dir() -> Path:
    d = project_root() / "docs" / "logs"
    d.mkdir(parents=True, exist_ok=True)
    return d


def now() -> str:
    return datetime.now().strftime("%Y-%m-%d %H:%M:%S")


def progress_file(session_id: str) -> Path:
    return logs_dir() / f"{session_id}.progress.md"


def log_file(session_id: str) -> Path:
    return logs_dir() / f"{session_id}.log.md"


def make_session_id(script: str) -> str:
    date = datetime.now().strftime("%Y-%m-%d")
    stem = re.sub(r"[^a-zA-Z0-9]+", "-", Path(script).stem).strip("-").lower()
    return f"{date}-{stem}"


CONDA_INIT = "source /scr1/users/liuc9/tools/miniforge3/etc/profile.d/conda.sh"
RENV_PATH = "/scr1/users/liuc9/tools/miniforge3/envs/renv"


def script_type(script: str) -> tuple[str, str]:
    ext = Path(script).suffix.lower()
    if ext in (".r", ".rmd"):
        env_cmd = f"{CONDA_INIT} && conda activate {RENV_PATH} && Rscript"
        return "R", env_cmd
    elif ext == ".py":
        return "Python", "uv run"
    return "Unknown", ""


def cmd_init(script: str, model: str = "claude") -> None:
    sid = make_session_id(script)
    pf = progress_file(sid)
    lf = log_file(sid)

    if pf.exists():
        print(f"Session exists: {sid}")
        print(f"  progress : {pf}")
        print(f"  log      : {lf}")
        return

    stype, env = script_type(script)
    ts = now()

    # Create progress.md
    if PROGRESS_TEMPLATE.exists():
        content = PROGRESS_TEMPLATE.read_text()
        for k, v in {
            "{SCRIPT_NAME}": Path(script).name,
            "{SCRIPT_PATH}": script,
            "{TIMESTAMP}": ts,
            "{SCRIPT_TYPE}": stype,
            "{ENV_CMD}": env,
            "{MODEL}": model,
        }.items():
            content = content.replace(k, v)
    else:
        content = f"""# Script Test: {Path(script).name}

**Script:** `{script}` | **Started:** {ts} | **Type:** {stype} | **Model:** {model}

---

## Tasks

- [ ] Run script
- [ ] Verify output files

---

## Notes

<!-- Agent writes observations here after reading {sid}.log.md -->

---

## Final Report

_Pending._
"""
    pf.write_text(content)

    # Create log.md (empty, will be filled by tee during script execution)
    lf.write_text(f"# Execution Log: {Path(script).name}\n\n**Started:** {ts}\n**Command:** `{env} {script}`\n\n---\n\n")

    run_cmd = f'{{ {env} {script} 2>&1; echo "SCRIPT_DONE:exit=$?"; }} | tee {lf}'
    print(f"Session ID  : {sid}")
    print(f"progress    : {pf}")
    print(f"log         : {lf}")
    print(f"Run command : {run_cmd}")
    print(f"Monitor cmd : python .opencode/skills/run-test/scripts/session.py monitor {sid}")
    print(f"Note        : For scripts >10min, run both commands with run_in_background=true")


def cmd_log(session_id: str, message: str) -> None:
    """Append an agent observation to progress.md (after reading log.md)."""
    pf = progress_file(session_id)
    if not pf.exists():
        print(f"Not found: {pf}")
        return
    ts = now()
    entry = f"\n### {ts}\n\n{message}\n\n---\n"
    content = pf.read_text()
    if "\n## Final Report" in content:
        content = content.replace("\n## Final Report", f"{entry}\n## Final Report")
    else:
        content += entry
    pf.write_text(content)
    print("Logged to progress.md")


def cmd_check(session_id: str, task: str) -> None:
    pf = progress_file(session_id)
    if not pf.exists():
        print(f"Not found: {pf}")
        return
    content = pf.read_text()
    lines = content.splitlines(keepends=True)
    for i, line in enumerate(lines):
        if re.match(r"- \[ \]", line) and task.lower() in line.lower():
            lines[i] = re.sub(r"- \[ \]", "- [x]", line, count=1)
            pf.write_text("".join(lines))
            print(f"Checked: {task}")
            return
    print(f"Task not found (or already checked): {task}")


def cmd_finish(session_id: str, status: str) -> None:
    pf = progress_file(session_id)
    if not pf.exists():
        print(f"Not found: {pf}")
        return
    ts = now()
    icon = "✓ Success" if status.lower() == "success" else "✗ Failed"
    report_section = f"""## Final Report

**Status:** {icon}
**Completed:** {ts}

<!-- Fill in: duration, output files, errors encountered, notes -->
"""
    content = pf.read_text()
    if "## Final Report" in content:
        content = re.sub(r"## Final Report\n[\s\S]*$", report_section, content)
    else:
        content += f"\n---\n\n{report_section}"
    pf.write_text(content)
    print(f"Finished ({icon}): {pf}")


def cmd_monitor(session_id: str, interval: int = 300, timeout: int = 10800) -> None:
    """Watch log.md and write periodic progress notes to progress.md.

    Designed to run as a background task alongside the script. Exits when the
    `SCRIPT_DONE:exit=N` sentinel line is detected in the log (script completed)
    or when the timeout is reached.

    Args:
        interval: Seconds between checks (default 300 = 5 min).
        timeout:  Max seconds to run before giving up (default 10800 = 3 h).
    """
    lf = log_file(session_id)
    start = _time.time()
    last_offset = 0
    check_num = 0

    print(f"Monitor started | session={session_id} interval={interval}s timeout={timeout}s", flush=True)

    while True:
        elapsed = _time.time() - start
        if elapsed > timeout:
            cmd_log(session_id, f"Monitor: timeout after {int(elapsed // 60)}min — stopping without detecting completion")
            print("Monitor: timeout, exiting")
            break

        _time.sleep(interval)
        check_num += 1
        elapsed = _time.time() - start
        elapsed_min = int(elapsed // 60)

        if not lf.exists():
            print("Monitor: log not found yet, waiting...", flush=True)
            continue

        content = lf.read_text(errors="replace")
        new_content = content[last_offset:]
        last_offset = len(content)

        # Detect completion: sentinel line written by run command after script exits
        m = re.search(r"SCRIPT_DONE:exit=(\d+)", new_content)
        if m:
            exit_code = m.group(1)
            status = "success" if exit_code == "0" else f"failed (exit={exit_code})"
            cmd_log(session_id, f"Script completed at +{elapsed_min}min — {status}")
            print("Monitor: completion detected, exiting", flush=True)
            break

        # Snapshot recent meaningful lines (skip blank lines and markdown headers)
        lines = [l.strip() for l in new_content.splitlines() if l.strip() and not l.startswith("#")]
        if lines:
            preview = " | ".join(lines[-5:])[:400]
            cmd_log(session_id, f"Check #{check_num} (+{elapsed_min}min): {preview}")
        else:
            cmd_log(session_id, f"Check #{check_num} (+{elapsed_min}min): no new output")

        print(f"Monitor: check #{check_num} at +{elapsed_min}min", flush=True)


def cmd_list() -> None:
    files = sorted(logs_dir().glob("*.progress.md"), reverse=True)
    if not files:
        print("No sessions found.")
        return
    for f in files:
        sid = f.name.replace(".progress.md", "")
        content = f.read_text()
        done = len(re.findall(r"- \[x\]", content, re.IGNORECASE))
        total = len(re.findall(r"- \[[ x]\]", content, re.IGNORECASE))
        finished = "## Final Report" in content and "_Pending._" not in content
        has_log = log_file(sid).exists()
        state = "done" if finished else "active"
        log_str = " +log" if has_log else ""
        print(f"  {sid}  ({done}/{total} tasks, {state}{log_str})")


def cmd_status(session_id: str) -> None:
    pf = progress_file(session_id)
    lf = log_file(session_id)
    if not pf.exists():
        print(f"Not found: {session_id}")
        return
    content = pf.read_text()
    done = len(re.findall(r"- \[x\]", content, re.IGNORECASE))
    total = len(re.findall(r"- \[[ x]\]", content, re.IGNORECASE))
    finished = "## Final Report" in content and "_Pending._" not in content
    log_size = f"{lf.stat().st_size} bytes" if lf.exists() else "not found"
    print(f"Session  : {session_id}")
    print(f"Tasks    : {done}/{total} done")
    print(f"Finished : {finished}")
    print(f"progress : {pf}")
    print(f"log      : {lf} ({log_size})")


def main() -> None:
    p = argparse.ArgumentParser(description="Run-Test Session Manager")
    sub = p.add_subparsers(dest="cmd")

    s = sub.add_parser("init", help="Create a new session")
    s.add_argument("script", help="Path to script (e.g. src/01-dataset.R)")
    s.add_argument("--model", "-m", default="claude", help="Model name (default: claude)")

    s = sub.add_parser("log", help="Append agent observation to progress.md")
    s.add_argument("session_id")
    s.add_argument("message")

    s = sub.add_parser("check", help="Mark a task complete (fuzzy match)")
    s.add_argument("session_id")
    s.add_argument("task", help="Partial task text to match")

    s = sub.add_parser("finish", help="Write final report section")
    s.add_argument("session_id")
    s.add_argument("status", choices=["success", "failed"])

    s = sub.add_parser("monitor", help="Watch log.md and write periodic progress notes (for long-running scripts)")
    s.add_argument("session_id")
    s.add_argument("--interval", "-i", type=int, default=300, help="Seconds between checks (default: 300)")
    s.add_argument("--timeout", "-t", type=int, default=10800, help="Max seconds to monitor (default: 10800)")

    sub.add_parser("list", help="List all sessions")

    s = sub.add_parser("status", help="Show session status and file paths")
    s.add_argument("session_id")

    args = p.parse_args()
    dispatch = {
        "init": lambda: cmd_init(args.script, args.model),
        "log": lambda: cmd_log(args.session_id, args.message),
        "check": lambda: cmd_check(args.session_id, args.task),
        "finish": lambda: cmd_finish(args.session_id, args.status),
        "monitor": lambda: cmd_monitor(args.session_id, args.interval, args.timeout),
        "list": cmd_list,
        "status": lambda: cmd_status(args.session_id),
    }
    fn = dispatch.get(args.cmd)
    if fn:
        fn()
    else:
        p.print_help()


if __name__ == "__main__":
    main()
