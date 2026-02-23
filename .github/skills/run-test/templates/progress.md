# Progress Log: {SCRIPT_NAME}

## Session Info

| Field             | Value           |
| ----------------- | --------------- |
| Started           | {TIMESTAMP}     |
| Script            | {SCRIPT_PATH}   |
| Type              | {SCRIPT_TYPE}   |
| Environment       | {ENV_CMD}       |
| Working Directory | {WORK_DIR}      |
| Model             | {MODEL}         |

---

## Execution Log

<!-- Entries are added chronologically as the session progresses -->

### [{TIMESTAMP}] Session Initialized

- Created log files
- Beginning pre-flight checks

---

### [{TIMESTAMP}] Phase: Pre-flight

**Environment Check:**
```
{environment check output}
```

**Dependencies Check:**
```
{dependency check output}
```

**Input Files Check:**
- [ ] File 1: {path} - {status}
- [ ] File 2: {path} - {status}

**Environment Variables Check:**
- `.env` file: {found/not found}
- Key variables: {list}

---

### [{TIMESTAMP}] Phase: Execution - Attempt {N}

**Command:**
```bash
{command}
```

**Status:** Running / Complete / Error

**Output:**
```
{stdout output}
```

**Errors (if any):**
```
{stderr output}
```

---

### [{TIMESTAMP}] Error Analysis

**Error Type:** {type}
**Root Cause:** {analysis}
**Proposed Fix:** {fix description}

---

### [{TIMESTAMP}] Fix Applied

**Issue:** {description}
**Solution:** {what was changed}
**Files Modified:**
- {file1}
- {file2}

---

### [{TIMESTAMP}] Phase: Execution - Attempt {N+1}

**Command:**
```bash
{command}
```

**Status:** Success ✓

**Output:**
```
{stdout output}
```

---

### [{TIMESTAMP}] Phase: Verification

**Output Files Generated:**
| File    | Size   | Status |
| ------- | ------ | ------ |
| {file1} | {size} | ✓      |
| {file2} | {size} | ✓      |

**Validation:**
- [ ] Output format correct
- [ ] Data integrity verified
- [ ] No warnings in log

---

### [{TIMESTAMP}] Session Complete

**Summary:**
- Total Duration: {duration}
- Attempts: {N}
- Status: Success/Failed
