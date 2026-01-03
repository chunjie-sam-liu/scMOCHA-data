---
name: runAndMonitorScript
description: Run a script, monitor execution, identify and fix issues, and provide progress summaries
argument-hint: The script file to run and monitor
---
# Run and Monitor Script Execution

1. **Understand the script**: Analyze the provided script file to understand its purpose and structure.

2. **Set up environment**: Use the appropriate environment activation (e.g., `conda activate renv` for R scripts).

3. **Run the script**: Execute the script in the background, redirecting output to a temporary log file.

4. **Monitor execution**: Periodically check the process status, log file growth, and progress indicators.

5. **Identify issues**: Scan the log for errors, warnings, or unexpected behavior.

6. **Fix issues automatically**: If issues are found, analyze the root cause and apply fixes to the code.

7. **Continue monitoring**: Resume execution after fixes and ensure smooth operation.

8. **Provide summaries**: Generate comprehensive summaries of progress, metrics, and completion status.

**Key monitoring points:**
- Process status (running/stopped)
- Log file size and line count
- Progress metrics (variants processed, files created, etc.)
- Error/warning counts
- Output file generation

**Automatic fixes to consider:**
- Path corrections (relative to absolute)
- Directory creation
- Logic fixes for integration steps
- Memory management improvements

**Summary format:**
- Process status
- Progress metrics
- Output files created
- Issues encountered and resolved
- Completion status
