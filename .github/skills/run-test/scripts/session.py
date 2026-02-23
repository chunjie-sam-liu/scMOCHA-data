#!/usr/bin/env python3
"""
Run-Test Session Manager

A utility for managing script testing sessions with AI coding assistants.
Creates and manages log files for tracking script testing progress.

This file is part of the run-test skill and is portable.
Templates are loaded from the same directory.

Usage:
    python .opencode/skills/run-test/session.py init path/to/script.R --model claude
    python .opencode/skills/run-test/session.py list
    python .opencode/skills/run-test/session.py status SESSION_ID
    python .opencode/skills/run-test/session.py append SESSION_ID "Entry text" --phase Execution
    python .opencode/skills/run-test/session.py complete SESSION_ID "Task description"
    python .opencode/skills/run-test/session.py error SESSION_ID "Error message" --type Package
    python .opencode/skills/run-test/session.py finish SESSION_ID --status success
"""

from __future__ import annotations

import argparse
import re
from datetime import datetime
from pathlib import Path

# Skill directory (where this script lives)
SKILL_DIR = Path(__file__).resolve().parent
TEMPLATES_DIR = SKILL_DIR / "templates"


def get_project_root() -> Path:
    """Find project root by looking for pyproject.toml or .git"""
    current = Path(__file__).resolve().parent
    while current != current.parent:
        if (current / "pyproject.toml").exists() or (current / ".git").exists():
            return current
        current = current.parent
    return Path.cwd()


def get_logs_dir() -> Path:
    """Get the logs directory path (always in project docs/logs)"""
    logs_dir = get_project_root() / "docs" / "logs"
    logs_dir.mkdir(parents=True, exist_ok=True)
    return logs_dir


def get_timestamp() -> str:
    """Get current timestamp in ISO format."""
    return datetime.now().strftime("%Y-%m-%d %H:%M:%S")


def get_date() -> str:
    """Get current date in YYYY-MM-DD format."""
    return datetime.now().strftime("%Y-%m-%d")


def sanitize_name(name: str) -> str:
    """Convert script name to safe filename component"""
    name = Path(name).stem
    name = re.sub(r"[^a-zA-Z0-9]+", "-", name)
    name = name.strip("-")
    return name.lower()


def detect_script_type(script_path: str) -> tuple[str, str]:
    """
    Detect script type from extension.
    Returns (type, environment_command)
    """
    ext = Path(script_path).suffix.lower()
    if ext in [".r", ".rmd", ".rscript"]:
        return "R", "conda run -n renv"
    elif ext in [".py", ".pyw"]:
        return "Python", "uv run"
    elif ext in [".sh", ".bash", ".slrm"]:
        return "Bash/SLURM", "bash"
    else:
        return "Unknown", ""


def get_session_id(script_name: str, model: str) -> str:
    """Generate session ID: YYYY-MM-DD-script-name-model"""
    date = get_date()
    safe_name = sanitize_name(script_name)
    safe_model = sanitize_name(model)
    return f"{date}-{safe_name}-{safe_model}"


def load_template(template_name: str) -> str | None:
    """Load a template file from the skill's templates directory."""
    template_path = TEMPLATES_DIR / f"{template_name}.md"
    if template_path.exists():
        return template_path.read_text()
    return None


def init_session(script_name: str, model: str = "claude") -> dict:
    """Initialize a new testing session with log files"""
    logs_dir = get_logs_dir()

    session_id = get_session_id(script_name, model)
    timestamp = get_timestamp()
    date = get_date()
    script_type, env_cmd = detect_script_type(script_name)
    workdir = str(get_project_root())

    files = {
        "todo": logs_dir / f"{session_id}.todo.md",
        "progress": logs_dir / f"{session_id}.progress.md",
        "report": logs_dir / f"{session_id}.report.md",
    }

    # Check if session already exists
    if files["todo"].exists():
        print(f"Session already exists: {session_id}")
        print(f"  Todo: {files['todo']}")
        print(f"  Progress: {files['progress']}")
        return {"session_id": session_id, "files": files, "status": "exists"}

    # Replacements for templates
    script_path = str(Path(script_name).resolve()) if Path(script_name).exists() else script_name
    replacements = {
        "{SCRIPT_NAME}": Path(script_name).name,
        "{SCRIPT_PATH}": script_path,
        "{DATE}": date,
        "{TIMESTAMP}": timestamp,
        "{MODEL}": model,
        "{WORK_DIR}": workdir,
        "{SCRIPT_TYPE}": script_type,
        "{ENV_CMD}": env_cmd,
        # Legacy placeholders
        "{R/Python}": script_type,
    }

    # Create todo.md
    todo_template = load_template("todo")
    if todo_template:
        content = todo_template
        for placeholder, value in replacements.items():
            content = content.replace(placeholder, value)
        files["todo"].write_text(content)
    else:
        files["todo"].write_text(f"""# Script Testing: {Path(script_name).name}

- **Date**: {date}
- **Model**: {model}
- **Script**: {script_name}
- **Type**: {script_type}
- **Environment**: {env_cmd}
- **Working Directory**: {workdir}

---

## Pre-flight Checks

- [ ] Identify script type ({script_type})
- [ ] Check environment ({env_cmd})
- [ ] Check dependencies installed
- [ ] Validate input files exist
- [ ] Check `.env` configuration
- [ ] Verify output directory exists

## Execution Tasks

- [ ] Run script (Attempt 1)
- [ ] Capture output and logs
- [ ] Check for errors
- [ ] Verify output files generated

## Post-execution

- [ ] Document execution time
- [ ] Record memory usage (if available)
- [ ] Generate final report
- [ ] Provide SLURM recommendations

---

## Errors to Fix

<!-- Add errors here as they occur -->

---

## Notes

<!-- Additional context, observations, decisions -->

---

## Session Summary

- **Total Attempts**: 0
- **Errors Fixed**: 0
- **Final Status**: Pending
""")
    print(f"Created: {files['todo']}")

    # Create progress.md
    progress_template = load_template("progress")
    if progress_template:
        content = progress_template
        for placeholder, value in replacements.items():
            content = content.replace(placeholder, value)
        # Add initialization entry
        init_entry = f"""
### [{timestamp}] Session Initialized

- Created log files
- Script type: {script_type}
- Environment: {env_cmd}
- Beginning pre-flight checks

---
"""
        content = content.replace(
            "<!-- Entries are added chronologically as the session progresses -->",
            f"<!-- Entries are added chronologically as the session progresses -->\n{init_entry}",
        )
        files["progress"].write_text(content)
    else:
        files["progress"].write_text(f"""# Progress Log: {Path(script_name).name}

## Session Info

| Field             | Value           |
| ----------------- | --------------- |
| Started           | {timestamp}     |
| Script            | {script_name}   |
| Type              | {script_type}   |
| Environment       | {env_cmd}       |
| Working Directory | {workdir}       |
| Model             | {model}         |

---

## Execution Log

### [{timestamp}] Session Initialized

- Created log files
- Script type: {script_type}
- Environment: {env_cmd}
- Beginning pre-flight checks

---
""")
    print(f"Created: {files['progress']}")
    print(f"Report will be: {files['report']}")

    return {"session_id": session_id, "files": files, "status": "created", "script_type": script_type, "env_cmd": env_cmd}


def list_sessions() -> list:
    """List all testing sessions"""
    logs_dir = get_logs_dir()
    if not logs_dir.exists():
        print("No sessions found (logs directory doesn't exist)")
        return []

    sessions = []
    for todo_file in logs_dir.glob("*.todo.md"):
        if todo_file.name.startswith("TEMPLATE"):
            continue
        session_id = todo_file.stem.replace(".todo", "")
        has_report = (logs_dir / f"{session_id}.report.md").exists()
        sessions.append({
            "session_id": session_id,
            "todo": todo_file,
            "progress": logs_dir / f"{session_id}.progress.md",
            "report": logs_dir / f"{session_id}.report.md",
            "has_report": has_report,
        })

    if not sessions:
        print("No sessions found.")
        return []

    sessions = sorted(sessions, key=lambda x: x["session_id"], reverse=True)
    print(f"Found {len(sessions)} session(s):\n")
    for s in sessions:
        status_info = get_session_status(s["session_id"])
        status_str = "Complete" if s["has_report"] else "In Progress"
        progress = status_info.get("progress", "N/A")
        print(f"  {s['session_id']}")
        print(f"    Status: {status_str} | Tasks: {progress}")
    print()

    return sessions


def get_session_status(session_id: str) -> dict:
    """Get status of a specific session"""
    logs_dir = get_logs_dir()
    files = {
        "todo": logs_dir / f"{session_id}.todo.md",
        "progress": logs_dir / f"{session_id}.progress.md",
        "report": logs_dir / f"{session_id}.report.md",
    }

    if not files["todo"].exists():
        return {"status": "not_found", "session_id": session_id}

    todo_content = files["todo"].read_text()

    unchecked = len(re.findall(r"- \[ \]", todo_content))
    checked = len(re.findall(r"- \[x\]", todo_content, re.IGNORECASE))
    total_tasks = checked + unchecked

    attempts_match = re.search(r"\*\*Total Attempts\*\*:\s*(\d+)", todo_content)
    errors_match = re.search(r"\*\*Errors Fixed\*\*:\s*(\d+)", todo_content)
    final_status_match = re.search(r"\*\*Final Status\*\*:\s*(\w+)", todo_content)

    result = {
        "status": "active",
        "session_id": session_id,
        "files": files,
        "total_tasks": total_tasks,
        "completed_tasks": checked,
        "pending_tasks": unchecked,
        "progress": f"{checked}/{total_tasks}" if total_tasks > 0 else "N/A",
        "has_report": files["report"].exists(),
    }

    if attempts_match:
        result["attempts"] = int(attempts_match.group(1))
    if errors_match:
        result["errors_fixed"] = int(errors_match.group(1))
    if final_status_match:
        result["final_status"] = final_status_match.group(1)

    return result


def append_to_progress(session_id: str, entry: str, phase: str = None) -> bool:
    """Append an entry to the progress log."""
    logs_dir = get_logs_dir()
    progress_file = logs_dir / f"{session_id}.progress.md"

    if not progress_file.exists():
        print(f"Error: Progress file not found: {progress_file}")
        return False

    timestamp = get_timestamp()
    header = f"### [{timestamp}]"
    if phase:
        header += f" Phase: {phase}"

    full_entry = f"\n{header}\n\n{entry}\n\n---\n"

    with open(progress_file, "a") as f:
        f.write(full_entry)

    print("Appended entry to progress log")
    return True


def update_todo_task(session_id: str, task: str, completed: bool = True) -> bool:
    """Update a task in the todo list."""
    logs_dir = get_logs_dir()
    todo_file = logs_dir / f"{session_id}.todo.md"

    if not todo_file.exists():
        print(f"Error: Todo file not found: {todo_file}")
        return False

    content = todo_file.read_text()

    if completed:
        pattern = re.compile(r"- \[ \] " + re.escape(task))
        if pattern.search(content):
            content = pattern.sub(f"- [x] {task}", content)
            todo_file.write_text(content)
            print(f"Marked complete: {task}")
            return True
        else:
            print(f"Task not found (or already complete): {task}")
            return False
    else:
        pattern = re.compile(r"- \[x\] " + re.escape(task), re.IGNORECASE)
        if pattern.search(content):
            content = pattern.sub(f"- [ ] {task}", content)
            todo_file.write_text(content)
            print(f"Marked pending: {task}")
            return True
        else:
            print(f"Task not found (or already pending): {task}")
            return False


def add_error_to_todo(session_id: str, error_msg: str, error_type: str = "Unknown") -> bool:
    """Add an error to the Errors to Fix section."""
    logs_dir = get_logs_dir()
    todo_file = logs_dir / f"{session_id}.todo.md"

    if not todo_file.exists():
        print(f"Error: Todo file not found: {todo_file}")
        return False

    content = todo_file.read_text()
    timestamp = get_timestamp()

    error_entry = f"\n- [ ] **[{error_type}]** {error_msg} _(added {timestamp})_"

    if "## Errors to Fix" in content:
        content = content.replace(
            "## Errors to Fix\n",
            f"## Errors to Fix\n{error_entry}\n",
        )
        todo_file.write_text(content)
        print(f"Added error to todo: [{error_type}] {error_msg[:50]}...")
        return True
    else:
        print("Error: Could not find 'Errors to Fix' section")
        return False


def update_session_summary(session_id: str, attempts: int = None, errors_fixed: int = None, final_status: str = None) -> bool:
    """Update the session summary in todo.md."""
    logs_dir = get_logs_dir()
    todo_file = logs_dir / f"{session_id}.todo.md"

    if not todo_file.exists():
        print(f"Error: Todo file not found: {todo_file}")
        return False

    content = todo_file.read_text()

    if attempts is not None:
        content = re.sub(
            r"\*\*Total Attempts\*\*:\s*\d+",
            f"**Total Attempts**: {attempts}",
            content
        )

    if errors_fixed is not None:
        content = re.sub(
            r"\*\*Errors Fixed\*\*:\s*\d+",
            f"**Errors Fixed**: {errors_fixed}",
            content
        )

    if final_status is not None:
        content = re.sub(
            r"\*\*Final Status\*\*:\s*\w+",
            f"**Final Status**: {final_status}",
            content
        )

    todo_file.write_text(content)
    print("Updated session summary")
    return True


def generate_report(session_id: str, status: str = "Success") -> bool:
    """Generate the final report from the session data."""
    logs_dir = get_logs_dir()
    todo_file = logs_dir / f"{session_id}.todo.md"
    report_file = logs_dir / f"{session_id}.report.md"

    if not todo_file.exists():
        print(f"Error: Todo file not found: {todo_file}")
        return False

    session_status = get_session_status(session_id)
    timestamp = get_timestamp()

    # Parse session_id for script name and model
    parts = session_id.split("-")
    if len(parts) >= 4:
        # Format: YYYY-MM-DD-scriptname-model
        script_name = "-".join(parts[3:-1])
        model = parts[-1]
    else:
        script_name = session_id
        model = "unknown"

    report_template = load_template("report")

    if report_template:
        content = report_template
        content = content.replace("{SCRIPT_NAME}", script_name)
        content = content.replace("{TIMESTAMP}", timestamp)
        content = content.replace("{SCRIPT_PATH}", session_id)
        content = content.replace("{MODEL}", model)
        content = content.replace("✓ Success / ✗ Failed", f"{'✓ Success' if status.lower() == 'success' else '✗ Failed'}")
    else:
        content = f"""# Test Report: {session_id}

**Generated:** {timestamp}
**Status:** {'✓ Success' if status.lower() == 'success' else '✗ Failed'}

## Summary

- **Tasks Completed**: {session_status.get('completed_tasks', 'N/A')}/{session_status.get('total_tasks', 'N/A')}
- **Attempts**: {session_status.get('attempts', 'N/A')}
- **Errors Fixed**: {session_status.get('errors_fixed', 'N/A')}

---

*Report generated by run-test skill*
"""

    report_file.write_text(content)
    print(f"Generated report: {report_file}")

    update_session_summary(session_id, final_status=status)

    return True


def main():
    parser = argparse.ArgumentParser(
        description="Run-Test Session Manager",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
    # Initialize a new session
    python .opencode/skills/run-test/session.py init eqtl/01-eqtl-to-duckdb.R --model claude

    # List all sessions
    python .opencode/skills/run-test/session.py list

    # Check session status
    python .opencode/skills/run-test/session.py status 2026-02-03-01-eqtl-to-duckdb-claude

    # Append to progress log
    python .opencode/skills/run-test/session.py append SESSION_ID "Script started" --phase Execution

    # Mark task as complete
    python .opencode/skills/run-test/session.py complete SESSION_ID "Run script (Attempt 1)"

    # Add an error
    python .opencode/skills/run-test/session.py error SESSION_ID "Missing package dplyr" --type Package

    # Generate final report
    python .opencode/skills/run-test/session.py finish SESSION_ID --status success

Environment Commands:
    R scripts:      conda run -n renv Rscript --vanilla script.R
    Python scripts: uv run script.py
        """,
    )

    subparsers = parser.add_subparsers(dest="command", help="Commands")

    # init command
    init_parser = subparsers.add_parser("init", help="Initialize a new testing session")
    init_parser.add_argument("script", help="Path to the script to test")
    init_parser.add_argument("--model", "-m", default="claude", help="Model name (default: claude)")

    # status command
    status_parser = subparsers.add_parser("status", help="Get session status")
    status_parser.add_argument("session_id", help="Session ID to check")

    # list command
    subparsers.add_parser("list", help="List all sessions")

    # append command
    append_parser = subparsers.add_parser("append", help="Append entry to progress log")
    append_parser.add_argument("session_id", help="Session ID")
    append_parser.add_argument("entry", help="Entry text to append")
    append_parser.add_argument("--phase", "-p", help="Phase name (e.g., Execution, Pre-flight)")

    # complete command
    complete_parser = subparsers.add_parser("complete", help="Mark a task as complete")
    complete_parser.add_argument("session_id", help="Session ID")
    complete_parser.add_argument("task", help="Task text to mark complete")

    # uncomplete command
    uncomplete_parser = subparsers.add_parser("uncomplete", help="Mark a task as pending")
    uncomplete_parser.add_argument("session_id", help="Session ID")
    uncomplete_parser.add_argument("task", help="Task text to mark pending")

    # error command
    error_parser = subparsers.add_parser("error", help="Add an error to the todo list")
    error_parser.add_argument("session_id", help="Session ID")
    error_parser.add_argument("error_msg", help="Error message")
    error_parser.add_argument("--type", "-t", default="Unknown", help="Error type (Package, File, Memory, etc.)")

    # summary command
    summary_parser = subparsers.add_parser("summary", help="Update session summary")
    summary_parser.add_argument("session_id", help="Session ID")
    summary_parser.add_argument("--attempts", "-a", type=int, help="Total attempts")
    summary_parser.add_argument("--errors", "-e", type=int, help="Errors fixed")
    summary_parser.add_argument("--status", "-s", help="Final status")

    # finish command
    finish_parser = subparsers.add_parser("finish", help="Generate final report and finish session")
    finish_parser.add_argument("session_id", help="Session ID")
    finish_parser.add_argument("--status", "-s", default="Success", help="Final status (Success/Failed)")

    args = parser.parse_args()

    if args.command == "init":
        result = init_session(args.script, args.model)
        print(f"\nSession ID: {result['session_id']}")
        print(f"Status: {result['status']}")
        if result['status'] == 'created':
            print(f"Script Type: {result.get('script_type', 'Unknown')}")
            print(f"Environment: {result.get('env_cmd', 'N/A')}")

    elif args.command == "status":
        result = get_session_status(args.session_id)
        print(f"Session: {result['session_id']}")
        print(f"Status: {result['status']}")
        if result["status"] == "active":
            print(f"Progress: {result['progress']}")
            if "attempts" in result:
                print(f"Attempts: {result.get('attempts', 'N/A')}")
            if "final_status" in result:
                print(f"Final Status: {result.get('final_status', 'N/A')}")
            print(f"Has Report: {result.get('has_report', False)}")

    elif args.command == "list":
        list_sessions()

    elif args.command == "append":
        append_to_progress(args.session_id, args.entry, args.phase)

    elif args.command == "complete":
        update_todo_task(args.session_id, args.task, completed=True)

    elif args.command == "uncomplete":
        update_todo_task(args.session_id, args.task, completed=False)

    elif args.command == "error":
        add_error_to_todo(args.session_id, args.error_msg, args.type)

    elif args.command == "summary":
        update_session_summary(
            args.session_id,
            attempts=args.attempts,
            errors_fixed=args.errors,
            final_status=args.status
        )

    elif args.command == "finish":
        generate_report(args.session_id, args.status)

    else:
        parser.print_help()


if __name__ == "__main__":
    main()
