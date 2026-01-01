# scMOCHA-data Agent Coding Guidelines

**Version:** 1.0.0
**Last Updated:** 2024-12-12
**Maintainer:** Chun-Jie Liu (chunjie.sam.liu.at.gmail.com)

---

## Quick Reference

### When to Use Which Script Type

| Script Type               | Naming Pattern            | Use Case                                              | Location                 |
| ------------------------- | ------------------------- | ----------------------------------------------------- | ------------------------ |
| **Numbered Analysis**     | `NN-descriptive-name.R`   | Sequential analysis workflows, figures for manuscript | `analysis/`              |
| **Sub-numbered Analysis** | `NN.N-descriptive-name.R` | Related sub-analyses within a numbered workflow       | `analysis/`              |
| **Uppercase Utility**     | `UPPERCASE_NAME.{R\|py}`  | Data processing utilities, database operations        | `analysis/`              |
| **Processing**            | `NN-action-description.R` | Data preprocessing, ETL pipelines                     | `preprocessing/`, `src/` |
| **Python Processing**     | `descriptive-name.py`     | Batch processing, parallel data operations            | `src/`, `analysis/`      |

### Essential Libraries Checklist

**R Scripts:**
- ✓ `load_pkg(jutils)` - Always first (includes magrittr and utilities)
- ✓ `load_pkg(ggplot2, data.table, dplyr)` - Core data manipulation and plotting
- ✓ `load_pkg(GetoptLong)` - For CLI arguments
- ✓ `load_pkg(logger)` - For logging
- ✓ Use `=` for assignment in args section, not `<-`

**Python Scripts:**
- ✓ `from pathlib import Path` - For file paths
- ✓ `import polars as pl` - For DataFrames (preferred over pandas)
- ✓ `import typer` - For CLI interface
- ✓ `from rich import print` - For rich console output
- ✓ `from rich.logging import RichHandler` - For logging

---

## 1. File Naming Conventions

### 1.1 Numbered Analysis Scripts

**Pattern:** `NN-descriptive-name.R` or `NN.NN-descriptive-name.R`

**Examples:**
- `00-colors.R` - Define color schemes
- `01-dataset-celltype-stats.R` - Main analysis script
- `15.1-Variant-celltype-specific-filter-plot.R` - Sub-analysis
- `15.2-Variant-celltype-specific-findpeaks.R` - Related sub-analysis

**Rules:**
- Use two-digit zero-padded numbers: `00`, `01`, `02`, ..., `99`
- Sub-numbering uses decimal notation: `15.1`, `15.2`, `15.3`
- Use hyphens to separate words: `variant-diff`, not `variant_diff` or `variantdiff`
- Capitalize significant words: `Disease-AD-variant-diff`, not `disease-ad-variant-diff`
- Keep names descriptive but concise (< 50 chars)

### 1.2 Uppercase Utility Scripts

**Pattern:** `UPPERCASE_DESCRIPTION.{R|py}`

**Examples:**
- `BARCODE_CELLTYPE.py` - Extract barcode-celltype mappings
- `COUNT2DUCKDB.py` - Convert count data to DuckDB
- `DUCKDB.R` - DuckDB database operations
- `EXPRCELL.py` - Expression-cell data processing

**Rules:**
- ALL UPPERCASE for entire filename (excluding extension)
- Use underscores to separate words: `ALL_VARIANT_CELLS_AF.py`
- Used for reusable utilities and data processing tools
- Can exist in any folder (`analysis/`, `src/`, `preprocessing/`)

### 1.3 Processing Scripts

**Pattern:** `NN-action-description.R` or `descriptive-name.R`

**Examples:**
- `01-sra-metadata.R` - Processing SRA metadata
- `02-load-variant.R` - Load variant data
- `plot-meta-variant.R` - Plot metadata and variants
- `somatic.R` - Somatic variant analysis

**Rules:**
- Numbered if part of sequential pipeline
- Action verbs preferred: `load-`, `plot-`, `collect-`, `parse-`
- Descriptive nouns okay for standalone scripts: `somatic.R`, `muscle.sh`

### 1.4 Python Processing Scripts

**Pattern:** `descriptive-name.py` or `NN-descriptive-name.py`

**Examples:**
- `scMOCHA.collectvariant.py` - Collect variants from scMOCHA
- `gds_xml2json.py` - Convert GDS XML to JSON
- `biosample_runinfo2csv.py` - Convert biosample to CSV

**Rules:**
- Use lowercase with underscores for internal separators
- Can use dots for namespace-like organization: `scMOCHA.collectvariant.py`
- Numbered if part of pipeline: `01-gather-files.py`, `02-update-high-quality-variant.py`

---

## 2. Script Headers and Metadata

### 2.1 R Script Header Template

**Mandatory header for ALL R scripts:**

```r
#!/usr/bin/env Rscript
# Metainfo ----------------------------------------------------------------
# @AUTHOR: Chun-Jie Liu
# @CONTACT: chunjie.sam.liu.at.gmail.com
# @DATE: YYYY-MM-DD HH:MM:SS
# @DESCRIPTION: Brief description of script purpose
# @VERSION: v0.0.1
```

**Alternative with --vanilla flag:**

```r
#!/usr/bin/env Rscript --vanilla
# Metainfo ----------------------------------------------------------------
# @AUTHOR: Chun-Jie Liu
# @CONTACT: chunjie.sam.liu.at.gmail.com
# @DATE: YYYY-MM-DD HH:MM:SS
# @DESCRIPTION: Brief description of script purpose
# @VERSION: v0.0.1
```

**Rules:**
- Always include shebang line: `#!/usr/bin/env Rscript`
- Use `--vanilla` flag for scripts that should avoid loading `.Rprofile`
- `@AUTHOR` is always "Chun-Jie Liu"
- `@CONTACT` is always "chunjie.sam.liu.at.gmail.com"
- `@DATE` format: `YYYY-MM-DD HH:MM:SS` (24-hour format)
- `@DESCRIPTION` should be one concise sentence
- `@VERSION` starts at `v0.0.1`

### 2.2 Python Script Header Template

**Mandatory header for ALL Python scripts:**

```python
#!/usr/bin/env python
# -*- coding:utf-8 -*-
# @AUTHOR: Chun-Jie Liu
# @CONTACT: chunjie.sam.liu.at.gmail.com
# @DATE: YYYY-MM-DD HH:MM:SS
# @DESCRIPTION: Brief description of script purpose
# @VERSION: v0.0.1
```

**Rules:**
- Always include shebang: `#!/usr/bin/env python`
- Always include encoding declaration: `# -*- coding:utf-8 -*-`
- Same metadata fields as R scripts
- Maintain consistent formatting with R header style

---

## 3. Code Organization

### 3.1 R Script Section Structure

**Standard section order (use all or subset as needed):**

```r
#!/usr/bin/env Rscript
# Metainfo ----------------------------------------------------------------
# @AUTHOR: Chun-Jie Liu
# @CONTACT: chunjie.sam.liu.at.gmail.com
# @DATE: 2024-12-12 10:00:00
# @DESCRIPTION: Template for analysis script
# @VERSION: v0.0.1

# Library -----------------------------------------------------------------
load_pkg(jutils)
load_pkg(ggplot2, patchwork, data.table, GetoptLong, logger)

# args --------------------------------------------------------------------
# s: string, i: integer, f: float, !: boolean
# @: array, %: hash
GetoptLong.options(help_style = "two-column")
VERSION = "v0.0.1"

verbose = FALSE
param = "default_value"

GetoptLong(
  "param=s",
  "parameter description",
  "verbose!",
  "print messages"
)

# src ---------------------------------------------------------------------
source("00-colors.R")  # If needed for color schemes

# header ------------------------------------------------------------------
log_threshold(TRACE)
log_layout(layout_glue_colors)

# function ----------------------------------------------------------------
fn_process_data <- function(.data) {
  # Function implementation
}

# load data ---------------------------------------------------------------
input_data <- import("/path/to/data.qs")

# body --------------------------------------------------------------------
# Main analysis logic here

# save --------------------------------------------------------------------
export(results, "output.qs")
```

**Section Header Rules:**
- Use format: `# Section Name -------...` with dashes extending to ~70 characters
- Keep sections in consistent order (Library → args → src → header → function → load data → body → save)
- Omit sections not needed (e.g., if no CLI args, skip `args` section)
- Use lowercase for section names: `# load data`, not `# Load Data`

### 3.2 Python Script Section Structure

**Standard structure:**

```python
#!/usr/bin/env python
# -*- coding:utf-8 -*-
# @AUTHOR: Chun-Jie Liu
# @CONTACT: chunjie.sam.liu.at.gmail.com
# @DATE: 2024-12-12 10:00:00
# @DESCRIPTION: Template for Python processing script
# @VERSION: v0.0.1

# Standard library imports
import os
import subprocess
from pathlib import Path
from typing import Annotated, Optional

# Third-party imports
import polars as pl
import typer
from rich import print
from rich.logging import RichHandler

# Constants and configuration
BASEDIR = Path("/liulab/chunjie/data/scMOCHA")
OUTDIR = BASEDIR / "output"
CONFIG = {
    "param1": "value1",
    "param2": "value2",
}

# Logging setup
import logging
logging.basicConfig(
    level=logging.INFO,
    format="%(message)s",
    datefmt="[%X]",
    handlers=[RichHandler(rich_tracebacks=True)]
)
logger = logging.getLogger("scMOCHA")

# Class definitions
class DataProcessor:
    def __init__(self, config: dict):
        self.config = config

    def process(self):
        pass

# Function definitions
def process_file(filepath: Path) -> pl.DataFrame:
    """Process a single file."""
    pass

# Main execution
def main(
    input_file: Annotated[str, typer.Argument(help="Input file path")],
    verbose: Annotated[bool, typer.Option("--verbose", "-v")] = False,
):
    """Main entry point."""
    logger.info(f"Processing {input_file}...")
    # Implementation

if __name__ == "__main__":
    typer.run(main)
```

**Structure Rules:**
- Group imports: standard library → third-party → local
- Use blank lines to separate import groups
- Constants in UPPERCASE
- Logging setup before class/function definitions
- Main logic in `if __name__ == "__main__":` block

---

## 4. Data I/O Patterns

### 4.1 R Data I/O

**Reading Data:**

```r
# Preferred: qs format (fast, compressed)
data <- import("/path/to/file.qs")

# Alternative: fst format
data <- import("/path/to/file.fst")

# Legacy: RDS compressed
data <- readr::read_rds("file.rds.gz")

# CSV files
data <- data.table::fread("file.csv")

# DuckDB databases
conn <- DBI::dbConnect(duckdb::duckdb(), "database.duckdb")
data <- dplyr::tbl(conn, "table_name") |> as.data.table()
DBI::dbDisconnect(conn, shutdown = TRUE)
```

**Writing Data:**

```r
# Preferred: qs format
export(data, "output.qs")

# Alternative: RDS compressed
readr::write_rds(data, "output.rds.gz", compress = "gz")

# CSV output
data.table::fwrite(data, "output.csv")

# DuckDB
conn <- DBI::dbConnect(duckdb::duckdb(), "database.duckdb")
DBI::dbWriteTable(conn, "table_name", data)
DBI::dbDisconnect(conn, shutdown = TRUE)
```

### 4.2 Python Data I/O

**Reading Data:**

```python
# Preferred: Polars CSV
df = pl.read_csv("file.csv")

# Polars Parquet
df = pl.read_parquet("file.parquet")

# DuckDB
import duckdb
conn = duckdb.connect("database.duckdb")
df = conn.execute("SELECT * FROM table").pl()
conn.close()
```

**Writing Data:**

```python
# Polars CSV
df.write_csv("output.csv")

# Polars Parquet
df.write_parquet("output.parquet")

# DuckDB
conn = duckdb.connect("database.duckdb")
conn.execute("CREATE TABLE table_name AS SELECT * FROM df")
conn.close()
```

### 4.3 Path Patterns

**Standard base paths:**

```r
# R
basedir <- "/liulab/chunjie/data/scMOCHA"
datadir <- fs::path(basedir, "data")
outdir <- fs::path(basedir, "analysis/zzz/MANUSCRIPTFIGURES")

# Using glue for dynamic paths
filepath <- glue::glue("{basedir}/data/{gseid}/{srrid}/variant.qs")
```

```python
# Python
BASEDIR = Path("/liulab/chunjie/data/scMOCHA")
DATADIR = BASEDIR / "data"
OUTDIR = BASEDIR / "analysis" / "zzz" / "MANUSCRIPTFIGURES"

# Dynamic paths
filepath = BASEDIR / "data" / gseid / srrid / "variant.qs"
```

**Directory structure:**
- `/liulab/chunjie/data/scMOCHA/` - Base directory
- `data/` - Raw and processed data
- `analysis/zzz/` - Analysis outputs
- `analysis/zzz/MANUSCRIPTFIGURES/` - Final publication figures

**File format preferences:**
1. **R:** `.qs` > `.fst` > `.rds.gz` > `.csv`
2. **Python:** `.parquet` > `.csv`
3. **Shared:** DuckDB (`.duckdb`) for large datasets

---

## 5. Common Libraries and Tools

### 5.1 R Core Libraries (Load in Order)

```r
# Use load_pkg() from jutils - loads packages quietly
load_pkg(jutils)  # ALWAYS FIRST - includes magrittr and utility functions

# Load additional packages as comma-separated arguments
load_pkg(
  ggplot2,      # Core plotting
  patchwork,    # Plot composition
  data.table,   # Fast data operations
  dplyr,        # Data manipulation
  tidyr,        # Data tidying
  glue,         # String interpolation
  fs,           # File system operations
  GetoptLong,   # CLI arguments
  logger        # Logging
)

# Or load packages individually if preferred
load_pkg(ggplot2, data.table, dplyr)
load_pkg(DBI, duckdb)  # Database packages
```

**Rules:**
- `load_pkg(jutils)` ALWAYS first - provides pipe operator and utilities
- Use `load_pkg()` instead of `library()` - automatically suppresses startup messages
- Can load multiple packages in single call: `load_pkg(pkg1, pkg2, pkg3)`
- Group related libraries together for readability
- Only load what you need - don't copy-paste full template if unused

### 5.2 R Specialized Libraries

```r
# Single-cell analysis
library(Seurat)
library(SeuratObject)

# Parallel processing
library(parallel)
library(furrr)

# Statistical analysis
library(stats)
library(ggpubr)

# Functional programming
library(purrr)

# String operations
library(stringr)

# Advanced plotting
library(ggrepel)
library(ggsci)
library(ComplexHeatmap)
library(circlize)
```

### 5.3 Python Core Libraries

```python
# Standard library (always at top)
import os
import subprocess
from pathlib import Path
from typing import Annotated, Optional
from concurrent.futures import ProcessPoolExecutor

# Data processing (preferred: Polars over Pandas)
import polars as pl
# import pandas as pd  # Only if Polars not suitable

# CLI and console
import typer
from rich import print
from rich.logging import RichHandler
from rich.progress import Progress, SpinnerColumn, BarColumn

# Database
import duckdb

# Logging
import logging
```

### 5.4 Python Specialized Libraries

```python
# Single-cell analysis
import scanpy as sc
import anndata as ad

# Numerical computing
import numpy as np

# Plotting (if needed, prefer R for complex plots)
import matplotlib.pyplot as plt
import seaborn as sns
```

---

## 6. Plotting and Visualization

### 6.1 Color Schemes

**Always source color definitions from:**

```r
source("00-colors.R")  # Contains all project color schemes
```

**Standard color palettes defined:**

```r
# Disease colors
color_disease <- c(
  "Alzheimer's Disease" = "#BC3C29FF",
  "COVID-19" = "#E18727FF",
  "Healthy" = "#0072B5FF",
  "Unknown" = "grey50"
)

# Cell type colors
color_celltype <- c(
  "B" = "#66C2A5FF",
  "CD4 T" = "#FC8D62FF",
  "CD8 T" = "#8DA0CBFF",
  "NK" = "#E78AC3FF",
  "Monocyte" = "#A6D854FF",
  # ... etc
)

# Dataset colors
color_dataset <- c(...)

# Sex colors
color_sex <- c(
  "Female" = "#E41A1C",
  "Male" = "#377EB8",
  "Unknown" = "grey50"
)
```

**Using color schemes:**

```r
# Convert to factor with color order
data <- data |>
  dplyr::mutate(
    disease = factor(disease, levels = names(color_disease))
  )

# Apply in ggplot
ggplot(data, aes(x = x, y = y, color = disease)) +
  geom_point() +
  scale_color_manual(values = color_disease)
```

### 6.2 Output Directory Standards

```r
# Manuscript figures
outdir <- "/liulab/chunjie/data/scMOCHA/analysis/zzz/MANUSCRIPTFIGURES"

# Ensure directory exists
fs::dir_create(outdir)

# Save plots with descriptive names
ggsave(
  filename = fs::path(outdir, "01-celltype-distribution.pdf"),
  plot = p,
  width = 8,
  height = 6
)
```

### 6.3 Python Visualization (Scanpy)

```python
import scanpy as sc

# Set output directory
sc.settings.figdir = OUTDIR

# Save plots with scanpy
sc.pl.umap(adata, color="celltype", save=f"_{gseid}_umap.pdf")
```

---

## 7. Documentation Style

### 7.1 Comments

```r
# ! Important note or warning
# TODO: Future improvement or fix needed
# Debug code retained for reference (but commented out)

# Short descriptive comment above code block
result <- complex_operation()

# Inline comments for clarification
data |>
  dplyr::filter(n > 10) |>  # Keep only well-represented groups
  dplyr::mutate(...)
```

### 7.2 Function Documentation

```r
# Function: fn_process_data
# Description: Process variant data with filtering and normalization
# Parameters:
#   .data - input data.frame
#   threshold - numeric threshold for filtering
# Returns: processed data.frame
fn_process_data <- function(.data, threshold = 0.05) {
  # Implementation
}
```

### 7.3 Section Dividers

```r
# Section Name ------------------------------------------------------------
# ^ Use exactly this format: lowercase, dashes to ~70 chars

# Subsection can use fewer dashes if needed
# subsection ------------------------------
```

### 7.4 Commented Code Retention

**Keep commented-out code for:**
- Alternative approaches tried
- Debugging steps that might be needed later
- Example usage
- Logger examples

```r
# Alternative approach using base R
# result <- lapply(data, function(x) { ... })

# Debug output
# print(head(data))
# log_info("Processing {nrow(data)} rows")
```

---

## 8. Variable Naming Conventions

### 8.1 R Naming Rules

**Use snake_case for all variables:**

```r
# Good
gse_dataset_metadata_full <- ...
all_hetero_af <- ...
variant_cell_count <- ...

# Bad
gseDatasetMetadataFull <- ...  # camelCase
AllHeteroAf <- ...             # PascalCase
variant.cell.count <- ...      # dot.notation
```

**Common prefixes:**
- `.x`, `.y` - Lambda/anonymous function parameters in `purrr`
- `.d`, `.m`, `.fit` - Temporary variables in pipelines
- `fn_` - Custom function names

**Boolean flags:**

```r
verbose <- FALSE
dir_exists <- fs::file_exists(path)
is_valid <- check_validity()
```

**Standard variable names:**

```r
basedir <- "/path/to/base"
datadir <- fs::path(basedir, "data")
outdir <- fs::path(basedir, "output")
infile <- "input.qs"
outfile <- "output.qs"
```

### 8.2 Python Naming Rules

**UPPERCASE for constants:**

```python
BASEDIR = Path("/liulab/chunjie/data/scMOCHA")
DATADIR = BASEDIR / "data"
OUTDIR = BASEDIR / "output"
MAX_WORKERS = 20
CONFIG = {...}
```

**snake_case for variables and functions:**

```python
variant_data = pl.read_csv("variants.csv")
cell_count = len(cells)

def process_variant_data(data: pl.DataFrame) -> pl.DataFrame:
    return data
```

**PascalCase for classes:**

```python
class VariantProcessor:
    def __init__(self, config: dict):
        self.config = config
```

**Private methods (single underscore):**

```python
def _internal_helper():
    pass
```

---

## 9. Error Handling and Logging

### 9.1 R Logging Setup

```r
# In header section
log_threshold(TRACE)
log_layout(layout_glue_colors)

# Usage throughout script (often commented but present)
# log_info('Starting analysis...')
# log_debug('Processing {nrow(data)} rows')
# log_success('Analysis complete!')
# log_warn('Low sample size detected: n = {n}')
# log_error('Failed to process: {error_msg}')
```

### 9.2 R Error Handling

```r
# tryCatch for operations that might fail
result <- tryCatch(
  {
    risky_operation()
  },
  error = function(e) {
    log_error("Operation failed: {e$message}")
    return(NULL)  # Return sensible default
  }
)

# Conditional return for validation
if (sum(.n_disease < 10) > 0) {
  log_warn("Insufficient samples in some groups")
  return(NULL)
}

# Existence checks
if (!fs::file_exists(infile)) {
  log_error("Input file not found: {infile}")
  stop("Missing required input file")
}
```

### 9.3 Python Logging Setup

```python
import logging
from rich.logging import RichHandler

logging.basicConfig(
    level=logging.INFO,
    format="%(message)s",
    datefmt="[%X]",
    handlers=[RichHandler(rich_tracebacks=True)]
)
logger = logging.getLogger("scMOCHA")

# Usage
logger.info(f"Processing {len(items)} items...")
logger.warning(f"Low quality detected: {qc_score}")
logger.error(f"Failed to process {item}: {e}")
```

### 9.4 Python Error Handling

```python
# Try-except with logging
try:
    result = process_data(data)
except Exception as e:
    logger.error(f"Processing failed: {e}")
    return None

# Validation with early return
if not filepath.exists():
    logger.error(f"File not found: {filepath}")
    raise FileNotFoundError(f"Missing required file: {filepath}")

# Context managers for cleanup
with duckdb.connect("database.duckdb") as conn:
    result = conn.execute("SELECT ...").fetchall()
# Connection automatically closed
```

---

## 10. Workflow Patterns

### 10.1 Parallel Processing

**R - mclapply:**

```r
library(parallel)

results <- parallel::mclapply(
  X = data_list,
  FUN = function(.x) {
    # Processing logic
    process_item(.x)
  },
  mc.cores = 20  # Typically 20 cores
)
```

**R - purrr with functional style:**

```r
data |>
  dplyr::mutate(
    result = purrr::map(.x = column, .f = \(.x) {
      # Lambda function processing
      process(.x)
    })
  )
```

**Python - ProcessPoolExecutor:**

```python
from concurrent.futures import ProcessPoolExecutor

def process_row(row: dict):
    # Processing logic
    return result

with ProcessPoolExecutor(max_workers=20) as executor:
    results = list(executor.map(process_row, data.iter_rows(named=True)))
```

### 10.2 Nested Data Structures (R tidyr)

```r
# Nest-map-unnest pattern
data |>
  tidyr::nest(.by = c(groupvar), .key = "nested_data") |>
  dplyr::mutate(
    result = purrr::map(.x = nested_data, .f = function(.x) {
      # Process each nested group
      analyze(.x)
    })
  ) |>
  tidyr::unnest(result)
```

### 10.3 Command-line Arguments

**R - GetoptLong:**

```r
# args section
# s: string, i: integer, f: float, !: boolean
# @: array, %: hash
GetoptLong.options(help_style = "two-column")
VERSION = "v0.0.1"

# Set default values using = (not <-)
verbose = FALSE
gseid = "GSE123456"
basedir = "/liulab/chunjie/data/scMOCHA"
ncores = 20

# Direct argument specification (no spec string)
GetoptLong(
  "gseid=s",
  "GSE accession ID",
  "basedir=s",
  "base directory path",
  "ncores=i",
  "number of cores",
  "verbose!",
  "print detailed messages"
)
```

**Python - Typer:**

```python
import typer
from typing import Annotated

def main(
    gseid: Annotated[str, typer.Argument(help="GSE accession ID")],
    basedir: Annotated[str, typer.Option("--basedir", "-b")] = "/default/path",
    verbose: Annotated[bool, typer.Option("--verbose", "-v")] = False,
    ncores: Annotated[int, typer.Option("--ncores", "-n")] = 20,
):
    """Process scMOCHA data for given GSE ID."""
    if verbose:
        logger.setLevel(logging.DEBUG)

    # Implementation

if __name__ == "__main__":
    typer.run(main)
```

### 10.4 DuckDB Patterns

**R - DuckDB operations:**

```r
# Connect
conn <- DBI::dbConnect(duckdb::duckdb(), "database.duckdb")

# Query with dplyr
result <- dplyr::tbl(conn, "table_name") |>
  dplyr::filter(variant_type == "heteroplasmic") |>
  dplyr::select(position, ref, alt, vaf) |>
  as.data.table()

# Write table
DBI::dbWriteTable(conn, "new_table", data, overwrite = TRUE)

# Always disconnect with shutdown
DBI::dbDisconnect(conn, shutdown = TRUE)
```

**Python - DuckDB with context manager:**

```python
class DuckDBManager:
    def __init__(self, dbfile: str):
        self._connection = duckdb.connect(dbfile)

    def __enter__(self):
        return self._connection

    def __exit__(self, exc_type, exc_val, exc_tb):
        self._connection.close()

# Usage
with DuckDBManager("database.duckdb") as conn:
    df = conn.execute("SELECT * FROM variants WHERE vaf > 0.01").pl()
```

### 10.5 Data Processing Pipelines (R)

**Standard pipeline pattern:**

```r
# Load -> Filter -> Transform -> Nest -> Process -> Unnest -> Export
data <- import("input.qs") |>
  dplyr::filter(!is.na(variant)) |>
  dplyr::mutate(
    vaf_pct = vaf * 100,
    category = case_when(
      vaf < 0.01 ~ "low",
      vaf < 0.1 ~ "medium",
      TRUE ~ "high"
    )
  ) |>
  tidyr::nest(.by = c(gseid, celltype), .key = "nested") |>
  dplyr::mutate(
    stats = purrr::map(.x = nested, .f = \(.x) {
      compute_statistics(.x)
    })
  ) |>
  tidyr::unnest(stats)

# Save result
export(data, "output.qs")
```

### 10.6 Batch Processing (Python)

```python
# Read metadata
SRR = pl.read_csv("gse_srrid_srrdir.csv")

# Define processing function
def process_sample(row: dict):
    gseid = row["gseid"]
    srrid = row["srrid"]

    cmd = [
        "Rscript",
        "process-sample.R",
        "-g", gseid,
        "-s", srrid
    ]

    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        logger.error(f"Failed {gseid}/{srrid}: {result.stderr}")
    return result.returncode

# Parallel execution
with ProcessPoolExecutor(max_workers=5) as executor:
    results = list(executor.map(
        process_sample,
        SRR.iter_rows(named=True)
    ))
```

---

## 11. Common Patterns and Idioms

### 11.1 Custom Functions (R)

**Naming:**
- Prefix all custom functions with `fn_`
- Use descriptive names: `fn_parse_log`, `fn_ks_test`, `fn_cor_test_variant`

**Lambda notation (modern R):**

```r
# Use \(.x) syntax for anonymous functions
data |>
  purrr::map(\(.x) {
    process(.x)
  })

# Instead of older function(.x) syntax
data |>
  purrr::map(function(.x) {
    process(.x)
  })
```

### 11.2 Data Validation (R)

```r
# Sample size checks
if (sum(.n_disease < 10) > 0) {
  log_warn("Some groups have < 10 samples")
  return(NULL)
}

# File existence
if (!fs::file_exists(infile)) {
  log_error("Missing input: {infile}")
  stop("Required file not found")
}

# Directory creation
if (!fs::dir_exists(outdir)) {
  fs::dir_create(outdir, recurse = TRUE)
}

# Data quality checks
data <- data |>
  dplyr::filter(
    !is.na(variant),
    vaf > 0,
    vaf <= 1
  )
```

### 11.3 Rich Progress Indicators (Python)

```python
from rich.progress import (
    Progress,
    SpinnerColumn,
    BarColumn,
    TimeElapsedColumn,
    TextColumn
)

with Progress(
    SpinnerColumn(),
    TextColumn("[progress.description]{task.description}"),
    BarColumn(),
    TextColumn("[progress.percentage]{task.percentage:>3.0f}%"),
    TimeElapsedColumn(),
) as progress:
    task = progress.add_task("[cyan]Processing samples...", total=len(samples))

    for sample in samples:
        process_sample(sample)
        progress.update(task, advance=1)
```

### 11.4 Factor Ordering by Color Scheme (R)

```r
# Order factors by color palette to ensure consistent plotting
data <- data |>
  dplyr::mutate(
    disease = factor(disease, levels = names(color_disease)),
    celltype = factor(celltype, levels = names(color_celltype))
  )
```

---

## 12. Version Control and Maintenance

### 12.1 Script Versioning

**Version number format:** `v[MAJOR].[MINOR].[PATCH]`

**Increment rules:**
- `PATCH` (v0.0.1 → v0.0.2): Bug fixes, typos, minor adjustments
- `MINOR` (v0.1.0 → v0.2.0): New features, additional analyses
- `MAJOR` (v1.0.0 → v2.0.0): Breaking changes, complete rewrites

**Update @VERSION and @DATE when:**
- Making substantial changes to logic
- Adding new features or analyses
- Fixing significant bugs
- Before committing to version control

### 12.2 Commented Code Management

**Retain commented code for:**
- Alternative approaches that were tried
- Debugging snippets that might be needed
- Example usage patterns
- Performance comparison alternatives

**Remove commented code when:**
- It's been superseded by better approach for >3 months
- It no longer applies to current data structure
- It causes confusion about current implementation

### 12.3 Dependency Management

**R packages:**
- Use `renv` for project-specific package management
- Document required packages in script header if unusual
- Keep `library()` calls in script (don't rely on `.Rprofile`)

**Python packages:**
- Use `uv` for package management (not pip)
- Keep `requirements.txt` or `pyproject.toml` updated
- Pin versions for critical dependencies

---

## 13. Pre-Commit Checklist

Before committing new or modified scripts, verify:

### R Scripts
- [ ] Header present with all metadata fields
- [ ] `load_pkg(jutils)` first, then other packages with `load_pkg()`
- [ ] Section headers in standard order
- [ ] Variables use snake_case
- [ ] Args use `=` assignment (not `<-`)
- [ ] GetoptLong direct call (no spec string)
- [ ] Custom functions prefixed with `fn_`
- [ ] Paths use `fs::path()` or `glue::glue()`
- [ ] Color schemes sourced from `00-colors.R` (if plotting)
- [ ] Output directory exists or is created
- [ ] Logging setup (even if commented)
- [ ] DuckDB connections properly closed
- [ ] File saved as `.qs` format (preferred for data)

### Python Scripts
- [ ] Header present with encoding and metadata
- [ ] Imports grouped: standard → third-party → local
- [ ] Constants in UPPERCASE
- [ ] Polars preferred over Pandas
- [ ] Typer used for CLI arguments
- [ ] Rich logging configured
- [ ] Context managers for file/DB operations
- [ ] Type hints on function parameters
- [ ] `if __name__ == "__main__":` block present
- [ ] Main function uses `typer.run()`

### Both
- [ ] File naming follows conventions
- [ ] No hardcoded personal paths (use relative or config)
- [ ] Comments explain "why", not "what"
- [ ] No sensitive information (passwords, tokens)
- [ ] Version number updated if substantial changes

---

## 14. Common Pitfalls to Avoid

### 14.1 R-Specific

❌ **Don't:**
```r
# Loading with library() instead of load_pkg()
library(magrittr)
library(ggplot2)  # Verbose, one at a time

# Using <- in args section
verbose <- FALSE
gseid <- "GSE123456"

# Using spec string for GetoptLong
spec <- "Usage: ...\nOptions:\n..."
GetoptLong(spec)

# Not ordering factors by color scheme
ggplot(aes(color = disease)) +
  scale_color_manual(values = color_disease)
# Mismatch if disease levels != color_disease names

# Forgetting to shutdown DuckDB
DBI::dbDisconnect(conn)  # Leaves files locked

# Using paste() for paths
path <- paste(basedir, "data", "file.txt", sep = "/")
```

✅ **Do:**
```r
# Use load_pkg() for all packages
load_pkg(jutils)  # Always first
load_pkg(ggplot2, data.table, dplyr)  # Load multiple at once

# Use = in args section
verbose = FALSE
gseid = "GSE123456"

# Direct GetoptLong() call
GetoptLong(
  "gseid=s",
  "GSE accession ID",
  "verbose!",
  "print messages"
)

# Order factors to match color scheme
data |>
  dplyr::mutate(disease = factor(disease, levels = names(color_disease)))

# Always shutdown DuckDB
DBI::dbDisconnect(conn, shutdown = TRUE)

# Use fs::path() or glue::glue()
path <- fs::path(basedir, "data", "file.txt")
```

### 14.2 Python-Specific

❌ **Don't:**
```python
# Using pandas when Polars would work
import pandas as pd
df = pd.read_csv("large_file.csv")  # Slower

# String concatenation for paths
path = basedir + "/data/" + filename

# Not using type hints
def process(data):
    return result

# Bare except clauses
try:
    process()
except:  # Too broad
    pass
```

✅ **Do:**
```python
# Prefer Polars for data processing
import polars as pl
df = pl.read_csv("large_file.csv")  # Faster

# Use pathlib
from pathlib import Path
path = basedir / "data" / filename

# Use type hints
def process(data: pl.DataFrame) -> pl.DataFrame:
    return result

# Specific exception handling
try:
    process()
except FileNotFoundError as e:
    logger.error(f"File not found: {e}")
```

### 14.3 General

❌ **Don't:**
- Hardcode absolute paths specific to your machine
- Mix tabs and spaces (use spaces only)
- Leave print statements for debugging in production code
- Use ambiguous variable names (`data`, `df`, `temp`, `x`)
- Skip error handling on file operations
- Forget to create output directories before writing

✅ **Do:**
- Use relative paths or configurable base directories
- Configure editor for 2-space indentation (R) or 4-space (Python)
- Use logger instead of print (can be disabled)
- Use descriptive names (`variant_data`, `celltype_counts`)
- Wrap file operations in try-catch with logging
- Check/create directories: `fs::dir_create()` or `Path.mkdir(parents=True, exist_ok=True)`

---

## 15. Contributing to These Guidelines

### 15.1 When to Update Guidelines

Update this document when you:
- Discover a new pattern used consistently across multiple scripts
- Identify a better practice that should become standard
- Find inconsistencies in current guidelines
- Add new tools or libraries to standard stack
- Deprecate old patterns in favor of new ones

### 15.2 How to Propose Changes

1. Document the current pattern vs. proposed pattern
2. Provide examples from existing codebase
3. Explain rationale (performance, readability, maintainability)
4. Update version number and date
5. Submit for review

### 15.3 Guideline Versioning

**Version format:** `[MAJOR].[MINOR].[PATCH]`

- `MAJOR`: Significant paradigm shifts (e.g., switching from pandas to polars)
- `MINOR`: New sections, substantial additions (e.g., adding new library standards)
- `PATCH`: Clarifications, typo fixes, minor adjustments

**Current version:** 1.0.0 (2024-12-12)

---

## 16. Quick Start Templates

For quick access to complete script templates, see:
- [`.github/templates/template-analysis-numbered.R`](.github/templates/template-analysis-numbered.R)
- [`.github/templates/template-utility-uppercase.R`](.github/templates/template-utility-uppercase.R)
- [`.github/templates/template-processing.R`](.github/templates/template-processing.R)
- [`.github/templates/template-processing.py`](.github/templates/template-processing.py)

---

## Appendix A: Common Use Cases

### A.1 Creating a New Analysis Script (Numbered)

1. **Choose number**: Next available in sequence or sub-number if related
2. **Copy template**: Use [template-analysis-numbered.R](.github/templates/template-analysis-numbered.R)
3. **Update header**: Date, description, version
4. **Load required libraries**: Core + specialized as needed
5. **Source colors**: `source("00-colors.R")` if plotting
6. **Set paths**: basedir, datadir, outdir
7. **Implement analysis**: Load → Process → Visualize → Save
8. **Test run**: Verify outputs and plots
9. **Commit**: With descriptive message

### A.2 Creating a Data Processing Utility (Python)

1. **Choose name**: UPPERCASE_DESCRIPTIVE.py
2. **Copy template**: Use [template-processing.py](.github/templates/template-processing.py)
3. **Update header**: Date, description, version
4. **Define constants**: BASEDIR, DATADIR, etc.
5. **Set up logging**: Configure Rich logger
6. **Implement logic**: Process → Transform → Save
7. **Add CLI**: Typer arguments and options
8. **Test**: Run with sample data
9. **Parallelize if needed**: Use ProcessPoolExecutor for large datasets
10. **Commit**: With usage example in commit message

### A.3 Adding a Preprocessing Step

1. **Determine position**: Sequential number in pipeline
2. **Copy template**: Use [template-processing.R](.github/templates/template-processing.R)
3. **Update header**: Date, description, version
4. **Check dependencies**: What input files are needed?
5. **Implement processing**: Read → Clean → Transform → Write
6. **Log progress**: Use logger for transparency
7. **Handle errors**: Graceful failures with informative messages
8. **Test pipeline**: Run full sequence from previous step
9. **Update documentation**: Add step to pipeline docs
10. **Commit**: With pipeline position noted

---

## Appendix B: Library Quick Reference

### B.1 When to Use Which R Library

| Task                   | Library         | Usage                                     |
| ---------------------- | --------------- | ----------------------------------------- |
| Read/write data        | `rio`           | `import()`, `export()`                    |
| Fast data ops          | `data.table`    | `fread()`, `fwrite()`, `[.data.table`     |
| Data manipulation      | `dplyr`         | `filter()`, `mutate()`, `select()`        |
| Data reshaping         | `tidyr`         | `nest()`, `unnest()`, `pivot_*()`         |
| Plotting               | `ggplot2`       | `ggplot()`, `geom_*()`                    |
| Combine plots          | `patchwork`     | `+`, `                                    | `, `/` operators |
| String interpolation   | `glue`          | `glue()`, `glue_data()`                   |
| File operations        | `fs`            | `path()`, `dir_create()`, `file_exists()` |
| Functional programming | `purrr`         | `map()`, `map2()`, `walk()`               |
| Parallel processing    | `parallel`      | `mclapply()`, `mcmapply()`                |
| Database               | `DBI`, `duckdb` | `dbConnect()`, `dbWriteTable()`           |
| CLI arguments          | `GetoptLong`    | `GetoptLong()`                            |
| Logging                | `logger`        | `log_info()`, `log_error()`               |

### B.2 When to Use Which Python Library

| Task                | Library              | Usage                                             |
| ------------------- | -------------------- | ------------------------------------------------- |
| DataFrames          | `polars`             | `pl.read_csv()`, `pl.DataFrame()`                 |
| Legacy DataFrames   | `pandas`             | `pd.read_csv()` (use only if Polars insufficient) |
| File paths          | `pathlib`            | `Path()`, `/` operator                            |
| CLI interface       | `typer`              | `typer.run()`, `typer.Argument()`                 |
| Console output      | `rich`               | `print()`, `Progress()`                           |
| Logging             | `logging` + `rich`   | `RichHandler()`, `logger.info()`                  |
| Parallel processing | `concurrent.futures` | `ProcessPoolExecutor()`                           |
| Database            | `duckdb`             | `duckdb.connect()`                                |
| Single-cell         | `scanpy`             | `sc.read_*()`, `sc.pl.*()`                        |
| Numerical           | `numpy`              | `np.array()`, `np.mean()`                         |

---

## Appendix C: File Format Decision Tree

**Should I use `.qs`, `.rds.gz`, `.fst`, `.csv`, `.parquet`, or `.duckdb`?**

```
Is data > 1GB?
├─ YES → Use DuckDB (.duckdb) for queryable database
└─ NO
    ├─ Need to share with Python?
    │   └─ Use Parquet (.parquet) or CSV (.csv)
    └─ R-only data?
        ├─ Need fast read/write?
        │   └─ Use qs (.qs)
        ├─ Need column-oriented access?
        │   └─ Use fst (.fst)
        └─ Need maximum compatibility?
            └─ Use RDS compressed (.rds.gz)
```

**File format characteristics:**

| Format     | Speed | Compression | R Support | Python Support | Size   |
| ---------- | ----- | ----------- | --------- | -------------- | ------ |
| `.qs`      | ⚡⚡⚡   | ⭐⭐⭐         | ✅         | ❌              | Small  |
| `.fst`     | ⚡⚡⚡   | ⭐⭐          | ✅         | ⚠️ Limited      | Medium |
| `.rds.gz`  | ⚡⚡    | ⭐⭐⭐         | ✅         | ❌              | Small  |
| `.parquet` | ⚡⚡    | ⭐⭐⭐         | ✅         | ✅              | Small  |
| `.csv`     | ⚡     | ❌           | ✅         | ✅              | Large  |
| `.duckdb`  | ⚡⚡    | ⭐⭐          | ✅         | ✅              | Medium |

**Recommendation hierarchy:**
1. **R-only, < 1GB:** `.qs`
2. **R-only, > 1GB:** `.duckdb`
3. **R + Python, < 1GB:** `.parquet`
4. **R + Python, > 1GB:** `.duckdb`
5. **Human-readable needed:** `.csv`

---

**End of Guidelines v1.0.0**

For questions or suggestions, contact: chunjie.sam.liu.at.gmail.com
