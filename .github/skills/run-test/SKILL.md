---
name: run-test
description: Agent skill for testing R/Python scripts before large-scale SLURM submission. Use when asked to "test script X", "run script X", "check script X", or "debug script X". This skill monitors execution, logs progress, fixes errors in a loop, and generates a final report with SLURM recommendations.
---

# Run-Test Agent Skill

This skill provides a structured workflow for testing bioinformatics scripts (R/Python) before submitting them for large-scale batch processing on SLURM clusters.

## When to Use This Skill

Activate this skill when the user:
- Says "test script X" or "run script X"
- Says "check script X" (pre-flight only)
- Says "fix script X" (error analysis)
- Wants to validate a script before SLURM submission
- Needs to debug a failing script with logging

## Environment Configuration

### R Scripts
- **load_pkg()**: Defined in ~/.Rprofile
- **Environment**: Use conda environment `renv`
- **Command**: `conda run -n renv Rscript --vanilla script.R`
- **Activation**: `conda activate renv`

### Python Scripts
- **Tool**: Use `uv` for Python execution
- **Command**: `uv run script.py`
- **No virtual environment activation needed** - uv handles it

## Skill File Structure

This skill is self-contained in `skills/run-test/`:

```
skills/run-test/
├── SKILL.md                    # This file
├── scripts/session.py          # Session manager utility
└── templates/
    ├── todo.md                 # Task tracking template
    ├── progress.md             # Execution log template
    └── report.md               # Final report template
```

Session logs are created in `docs/logs/` with naming:
```
docs/logs/{YYYY-MM-DD}-{script-name}-{model}.todo.md
docs/logs/{YYYY-MM-DD}-{script-name}-{model}.progress.md
docs/logs/{YYYY-MM-DD}-{script-name}-{model}.report.md
```

---

## Quick Reference

### Session Manager Commands

```bash
# Initialize a new session
python skills/run-test/scripts/session.py init path/to/script.R --model {model}

# List all sessions
python skills/run-test/scripts/session.py list

# Check session status
python skills/run-test/scripts/session.py status SESSION_ID

# Append to progress log
python skills/run-test/scripts/session.py append SESSION_ID "Entry text" --phase Execution

# Mark task complete
python skills/run-test/scripts/session.py complete SESSION_ID "Task description"

# Add error to fix
python skills/run-test/scripts/session.py error SESSION_ID "Error message" --type Package

# Generate final report
python skills/run-test/scripts/session.py finish SESSION_ID --status success
```

### User Commands

| User Says        | Agent Action                                                  |
| ---------------- | ------------------------------------------------------------- |
| "test script X"  | Full workflow: init → preflight → execute → fix loop → report |
| "check script X" | Pre-flight checks only, no execution                          |
| "fix script X"   | Resume from last error, attempt fixes                         |
| "status"         | Show current session status from todo.md                      |
| "report"         | Generate final report from current session                    |

---

## Workflow Phases

### Phase 1: Initialize Session

1. **Identify the script** - Get the full path and script type (R/Python)
2. **Generate session ID** - Format: `{date}-{script-basename}-{model}`
3. **Create log files** - Initialize todo.md and progress.md from templates
4. **Log initialization** - Record session start in progress.md

```bash
# Using session manager
python skills/run-test/scripts/session.py init path/to/script.R --model {model}
```

### Phase 2: Pre-flight Checks

Before running the script, verify:

#### For R Scripts (conda renv):
1. **Conda environment**: `conda env list | grep renv`
2. **R installation**: `conda run -n renv R --version`
3. **Package dependencies**: Parse script for `library()` and `load_pkg()` calls
4. **Check packages installed**: `conda run -n renv Rscript -e "installed.packages()[,'Package']"`
5. **Environment files**: Check for `.env` in script directory or project root
6. **Input files**: Look for `import()`, `read_*()`, or file path arguments
7. **Output directories**: Verify output paths exist or can be created

```bash
# Check conda renv environment
conda env list | grep renv

# Check R packages in renv
conda run -n renv Rscript -e "if (!require('dplyr')) stop('Package not found')"

# Check jutils
conda run -n renv Rscript -e "library(jutils)"
```

#### For Python Scripts (uv):
1. **uv installation**: `which uv` or `uv --version`
2. **Check pyproject.toml**: Verify dependencies are defined
3. **Sync environment**: `uv sync` if needed
4. **Package dependencies**: Parse imports or check pyproject.toml
5. **Input files**: Look for file path arguments or hardcoded paths
6. **Output directories**: Verify output paths exist

```bash
# Check uv
uv --version

# Check Python packages
uv pip list

# Test import
uv run python -c "import pandas; print(pandas.__version__)"
```

### Phase 3: Execute Script

1. **Construct the command**:
   - R: `conda run -n renv Rscript --vanilla path/to/script.R [args]`
   - Python: `uv run path/to/script.py [args]`

2. **Run with output capture**:
   ```bash
   # R script
   { time conda run -n renv Rscript --vanilla script.R 2>&1; } 2>&1 | tee output.log

   # Python script
   { time uv run script.py 2>&1; } 2>&1 | tee output.log
   ```

3. **Monitor execution**:
   - Log start timestamp
   - Capture stdout/stderr
   - Log completion timestamp
   - Record exit code

4. **Update progress.md** with execution details

### Phase 4: Error Handling Loop

If the script fails, enter the error-fix-retry loop:

```
┌─────────────────────────────────────────┐
│           Execute Script                │
└─────────────────┬───────────────────────┘
                  │
                  ▼
        ┌─────────────────┐
        │  Check Exit     │
        │     Code        │
        └────────┬────────┘
                 │
       ┌─────────┴─────────┐
       │                   │
       ▼                   ▼
   Exit = 0           Exit != 0
   (Success)           (Error)
       │                   │
       ▼                   ▼
   ┌───────┐      ┌───────────────┐
   │ Done! │      │ Analyze Error │
   └───────┘      └───────┬───────┘
                          │
                          ▼
                 ┌─────────────────┐
                 │  Classify Error │
                 │  Type & Cause   │
                 └────────┬────────┘
                          │
                          ▼
                 ┌─────────────────┐
                 │   Apply Fix     │
                 │  (if possible)  │
                 └────────┬────────┘
                          │
                          ▼
                 ┌─────────────────┐
                 │ Log Fix in      │
                 │ progress.md     │
                 └────────┬────────┘
                          │
                          ▼
                 ┌─────────────────┐
                 │  Increment      │
                 │  Attempt Count  │
                 └────────┬────────┘
                          │
           ┌──────────────┴──────────────┐
           │                             │
           ▼                             ▼
    Attempts < MAX              Attempts >= MAX
           │                             │
           ▼                             ▼
    ┌──────────────┐            ┌─────────────────┐
    │ Re-execute   │            │ Report Failure  │
    │   Script     │──────┐     │ (Need User Help)│
    └──────────────┘      │     └─────────────────┘
           │              │
           └──────────────┘
```

#### Error Classification

| Error Type         | Examples                               | Auto-Fix Possible      |
| ------------------ | -------------------------------------- | ---------------------- |
| **Package/Import** | Missing R package, Python import error | Sometimes (install)    |
| **File Not Found** | Input file missing, wrong path         | No (need user input)   |
| **Permission**     | Cannot write to directory              | Sometimes (mkdir)      |
| **Memory**         | Cannot allocate vector of size X       | No (need SLURM config) |
| **Syntax**         | Parse error, invalid code              | Yes (edit script)      |
| **Data**           | Column not found, type mismatch        | Sometimes              |
| **Environment**    | Missing env variable                   | Sometimes (dotenv)     |
| **Timeout**        | Script takes too long                  | No (need user input)   |

#### Fix Application

1. **Log the error** in progress.md with full traceback
2. **Analyze root cause** - 3-5 potential causes
3. **Propose fix** - Explain what needs to change
4. **Apply fix** if auto-fixable:
   - Install missing package (R: `conda run -n renv Rscript -e "install.packages('pkg')"`)
   - Install missing package (Python: `uv add pkg`)
   - Create output directory
   - Fix syntax error in script
   - Set environment variable
5. **Document changes** in progress.md
6. **Update todo.md** - Mark error as fixed or blocked

### Phase 5: Verification

After successful execution:

1. **Check output files exist** with expected sizes
2. **Validate output format** (if possible)
3. **Record resource usage**:
   - Execution time
   - Peak memory (from `/usr/bin/time -v` if available)
   - CPU usage

### Phase 6: Generate Report

Create the final report in `{session}.report.md`:

1. **Executive Summary**
   - Final status (Success/Failed)
   - Total duration
   - Attempt count
   - Errors encountered and fixed

2. **Execution History Table**
   - Each attempt with timestamp, status, duration

3. **Errors & Resolutions**
   - Each error with type, message, cause, fix

4. **Output Files**
   - List with sizes and verification status

5. **Resource Profile**
   - Observed metrics
   - SLURM recommendations

6. **Large-Scale Readiness Checklist**
   - All criteria for batch submission

---

## SLURM Recommendations

Based on test execution, recommend SLURM parameters:

### Memory Estimation

| Observed Peak | Recommended `--mem` |
| ------------- | ------------------- |
| < 1 GB        | 2G                  |
| 1-4 GB        | 8G                  |
| 4-16 GB       | 32G                 |
| 16-64 GB      | 64G                 |
| > 64 GB       | 128G or more        |

### Time Estimation

| Observed Time | Recommended `--time` |
| ------------- | -------------------- |
| < 1 min       | 00:10:00             |
| 1-10 min      | 00:30:00             |
| 10-60 min     | 02:00:00             |
| 1-4 hours     | 08:00:00             |
| > 4 hours     | 24:00:00 or more     |

### CPU Recommendations

- Single-threaded R/Python: `--cpus-per-task=1`
- Parallel processing (mclapply, pbmclapply): `--cpus-per-task=8`
- Heavy computation: `--cpus-per-task=16`

### Example SLURM Script (R)

```bash
#!/bin/bash
#SBATCH --job-name={script_name}
#SBATCH --time={recommended_time}
#SBATCH --mem={recommended_memory}
#SBATCH --cpus-per-task={cpus}
#SBATCH --output=logs/%x.%j.out
#SBATCH --error=logs/%x.%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user={user_email}

# Activate conda environment
source $(conda info --base)/etc/profile.d/conda.sh
conda activate renv

# Set working directory
cd $SLURM_SUBMIT_DIR

# Run R script
Rscript --vanilla {script_path}
```

### Example SLURM Script (Python)

```bash
#!/bin/bash
#SBATCH --job-name={script_name}
#SBATCH --time={recommended_time}
#SBATCH --mem={recommended_memory}
#SBATCH --cpus-per-task={cpus}
#SBATCH --output=logs/%x.%j.out
#SBATCH --error=logs/%x.%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user={user_email}

# Set working directory
cd $SLURM_SUBMIT_DIR

# Run Python script with uv
uv run {script_path}
```

---

## Error Recovery Strategies

### Package Not Found (R - conda renv)

```bash
# Error: there is no package called 'dplyr'

# Fix: Install in conda renv environment
conda run -n renv Rscript -e "install.packages('dplyr', repos='https://cloud.r-project.org')"

# Or install via conda
conda install -n renv -c conda-forge r-dplyr
```

### Module Import Error (Python - uv)

```bash
# Error: ModuleNotFoundError: No module named 'pandas'

# Fix: Add with uv
uv add pandas

# Then sync
uv sync
```

### File Not Found

```
# Error: cannot open file 'data/input.csv': No such file or directory

# Analysis:
# 1. Check if file exists: ls -la data/
# 2. Check working directory: pwd
# 3. Check relative vs absolute path

# Cannot auto-fix - need user to provide correct path or file
```

### Memory Error (R)

```r
# Error: cannot allocate vector of size 2.5 Gb

# Cannot auto-fix for test run
# Recommendation: Run with SLURM --mem=32G or higher
# Or: Process data in chunks
```

### Permission Denied

```bash
# Error: cannot create directory '/output/results': Permission denied

# Fix: Create in accessible location
mkdir -p ~/project/output/results
# Update script path accordingly
```

---

## Integration with jutils

When testing R scripts that use jutils:

1. **Check jutils is installed**: `conda run -n renv Rscript -e "library(jutils)"`
2. **Load .env automatically**: jutils `dotenv()` loads environment
3. **Use import/export**: Check for `import()` and `export()` calls
4. **DuckDB connections**: Watch for `db_conn()` and ensure cleanup
5. **Parallel processing**: `pbmclapply` uses 8 cores by default

---

## Example Session

### User Request
```
test script eqtl/01-eqtl-to-duckdb.R
```

### Agent Response

1. **Initialize Session**
   ```
   Creating session: 2026-02-03-01-eqtl-to-duckdb-claude
   - docs/logs/2026-02-03-01-eqtl-to-duckdb-claude.todo.md
   - docs/logs/2026-02-03-01-eqtl-to-duckdb-claude.progress.md
   ```

2. **Pre-flight Checks**
   ```
   ✓ Script type: R
   ✓ Conda env renv: available
   ✓ R version: 4.3.0 (via conda run -n renv R --version)
   ✓ Required packages: jutils, dplyr, arrow (all installed)
   ✓ .env file found: eqtl/.env
   ✓ Input files: /data/eqtl/*.parquet (exist)
   ✓ Output directory: /data/duckdb/ (exists)
   ```

3. **Execute Script**
   ```
   Running: conda run -n renv Rscript --vanilla eqtl/01-eqtl-to-duckdb.R
   [Execution output...]

   Status: Error (exit code 1)
   ```

4. **Error Analysis**
   ```
   Error: Column 'rsid' not found in data

   Analysis:
   - Input file has 'rs_id' instead of 'rsid'
   - This is a column naming mismatch

   Fix: Update script to use 'rs_id' column name
   ```

5. **Apply Fix & Retry**
   ```
   Applied fix: Changed line 45 from `rsid` to `rs_id`
   Re-running script...

   Running: conda run -n renv Rscript --vanilla eqtl/01-eqtl-to-duckdb.R
   Status: Success (exit code 0)
   ```

6. **Generate Report**
   ```
   Final Report: docs/logs/2026-02-03-01-eqtl-to-duckdb-claude.report.md

   Summary:
   - Status: Success
   - Attempts: 2
   - Duration: 3m 45s
   - Memory: ~2.1 GB

   SLURM Recommendation:
   --time=00:15:00 --mem=8G --cpus-per-task=1
   ```

---

## Checklist for Agent

When testing a script:

- [ ] Identify script path and type (R/Python)
- [ ] Create session log files using session.py
- [ ] For R: Verify conda renv environment is available
- [ ] For Python: Verify uv is available
- [ ] Check all dependencies are available
- [ ] Verify input files exist and are readable
- [ ] Check environment variables are set
- [ ] Ensure output directory exists
- [ ] Run script with correct command (conda run -n renv / uv run)
- [ ] If error: analyze, fix if possible, retry
- [ ] If success: verify output files
- [ ] Record execution time and resource usage
- [ ] Generate final report with SLURM recommendations
- [ ] Mark session complete in todo.md
