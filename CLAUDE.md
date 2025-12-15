# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

**scMOCHA-data** is a single-cell mitochondrial variant analysis pipeline for large-scale scRNA-seq datasets. The project processes data from multiple GEO datasets to identify, classify, and analyze mitochondrial DNA variants at the single-cell level. Key focus areas include heteroplasmic and homoplasmic variants, somatic mutations, and cell-type-specific variant patterns.

**Primary Technologies:**
- R (analysis, visualization, statistics)
- Python (data processing, batch operations, DuckDB)
- WDL (Workflow Description Language for batch processing)
- Conda environment (`scmocha.yaml`)

## Repository Structure

```
├── src/              # Data acquisition and initial processing (01-06)
├── preprocessing/    # Data loading, summarization, quality control
├── analysis/         # Numbered analysis scripts (00-50) and utilities
├── fisher/           # Fisher test variant analyses
├── config/           # Reference data (GTF, FASTA, color schemes)
├── data/             # Symlink to /mnt/isilon/.../raw (raw data)
├── large-scale/      # Symlink to shared project directory
└── .github/          # Coding guidelines and templates
```

**Key Entry Points:**
- `src/01-sra-metadata.R` → SRA data acquisition
- `preprocessing/02-load-variant.R` → Load variant calls
- `analysis/00-colors.R` → Project color schemes (source this for plots)
- `analysis/NN-*.R` → Numbered analysis workflows

## Common Commands

**Environment Setup:**
```bash
# Activate conda environment
conda activate scmocha

# Or use renv for R-only work
conda activate renv
```

**Running Analysis Scripts:**
```bash
# R scripts with CLI arguments (use GetoptLong)
Rscript analysis/01-dataset-celltype-stats.R --gseid=GSE226602 --verbose

# Python utilities (use typer)
python analysis/COUNT2DUCKDB.py input.csv --verbose

# Parallel processing (typically 20 cores)
Rscript analysis/03-individual-variants-correlation-with-age-disease.R --ncores=20
```

**Common Workflows:**
```bash
# 1. Acquire SRA metadata
Rscript src/01-sra-metadata.R --gseid=GSE123456

# 2. Download and process FASTQ
Rscript src/02-sra-download-dump.R --gseid=GSE123456

# 3. Load variants
Rscript preprocessing/02-load-variant.R --gseid=GSE123456

# 4. Run numbered analysis
Rscript analysis/01-dataset-celltype-stats.R --gseid=GSE123456
```

## Architecture and Patterns

### Data Flow Pipeline

1. **Data Acquisition** (`src/01-04`): SRA metadata → FASTQ download → scMOCHA variant calling
2. **Preprocessing** (`preprocessing/`): Load variant calls → summarize metadata → quality control
3. **Analysis** (`analysis/00-50`): Numbered sequential analyses for manuscript figures
4. **Utilities** (`analysis/UPPERCASE_NAME.*`): Reusable data processing tools

### Variant Classification System

Mitochondrial variants are classified hierarchically (see `src/06.4-somatic.md`):

```
All Variants
├─ Unreliable (excluded): mis-alignment, RNA editing, <10 cells
└─ Reliable
   ├─ Haplogroup-defining (ethnicity-related)
   ├─ Homoplasmic (AF ≥ 0.95 across cell types)
   └─ Heteroplasmic (0.05 ≤ AF < 0.95)
      └─ Somatic (cell-type-specific heteroplasmic subset)
```

**Key Thresholds:**
- `CUTOFF_MIN_READS`: 10
- `CUTOFF_MIN_CELLS`: 10 per cell type
- `CUTOFF_HETEROPLASMIC`: 0.05
- `CUTOFF_HOMOPLASMIC`: 0.95
- `CUTOFF_SOMATIC_MIN_N_CELLS`: 10

### File Naming Conventions

**R Analysis Scripts:**
- `NN-descriptive-name.R` (e.g., `01-dataset-celltype-stats.R`)
- `NN.N-descriptive-name.R` for sub-analyses (e.g., `15.2-Variant-celltype-specific-findpeaks.R`)

**Utility Scripts:**
- `UPPERCASE_NAME.{R|py}` (e.g., `COUNT2DUCKDB.py`, `BARCODE_CELLTYPE.py`)

**Processing Scripts:**
- `NN-action-description.R` (e.g., `02-load-variant.R`)

### Standard Script Structure

**R Script Template:**
```r
#!/usr/bin/env Rscript
# Metainfo ----------------------------------------------------------------
# @AUTHOR: Chun-Jie Liu
# @CONTACT: chunjie.sam.liu.at.gmail.com
# @DATE: YYYY-MM-DD HH:MM:SS
# @DESCRIPTION: Brief description
# @VERSION: v0.0.1

# Library -----------------------------------------------------------------
load_pkg(jutils)  # ALWAYS FIRST
load_pkg(ggplot2, data.table, dplyr)

# args --------------------------------------------------------------------
GetoptLong.options(help_style = "two-column")
gseid = "GSE123456"  # Use = not <-
verbose = FALSE

GetoptLong(
  "gseid=s", "GSE accession ID",
  "verbose!", "print messages"
)

# src ---------------------------------------------------------------------
source("00-colors.R")  # For visualization

# header ------------------------------------------------------------------
log_threshold(TRACE)
log_layout(layout_glue_colors)

# function ----------------------------------------------------------------
# body --------------------------------------------------------------------
# save --------------------------------------------------------------------
```

**Python Script Template:**
```python
#!/usr/bin/env python
# -*- coding:utf-8 -*-
# @AUTHOR: Chun-Jie Liu
# @CONTACT: chunjie.sam.liu.at.gmail.com
# @DATE: YYYY-MM-DD HH:MM:SS
# @DESCRIPTION: Brief description
# @VERSION: v0.0.1

from pathlib import Path
import polars as pl  # Prefer over pandas
import typer
from rich import print

BASEDIR = Path("/liulab/chunjie/data/scMOCHA")
```

### Data I/O Standards

**R Preferences:**
```r
# Read: prefer qs > fst > rds.gz > csv
data <- import("file.qs")  # jutils::import

# Write
export(data, "output.qs")  # jutils::export

# Paths
filepath <- fs::path(basedir, "data", "file.qs")

# DuckDB
conn <- DBI::dbConnect(duckdb::duckdb(), "db.duckdb")
data <- dplyr::tbl(conn, "table") |> as.data.table()
DBI::dbDisconnect(conn, shutdown = TRUE)  # Always shutdown
```

**Python Preferences:**
```python
# Read: prefer Polars over Pandas
df = pl.read_csv("file.csv")

# Paths
from pathlib import Path
filepath = BASEDIR / "data" / "file.csv"

# DuckDB
import duckdb
conn = duckdb.connect("db.duckdb")
df = conn.execute("SELECT * FROM table").pl()
conn.close()
```

### Critical Coding Rules

**R-Specific:**
- `load_pkg(jutils)` ALWAYS first (provides pipe, utilities)
- Use `load_pkg(pkg1, pkg2)` NOT `library()`
- Args section: use `=` NOT `<-` for assignments
- GetoptLong: direct call, no spec string
- Factor ordering: match color schemes (`factor(x, levels = names(color_x))`)
- DuckDB: ALWAYS `shutdown = TRUE` on disconnect

**Python-Specific:**
- Prefer `polars` over `pandas`
- Use `typer` for CLI
- Use `rich` for console output
- Constants in UPPERCASE
- Type hints on functions

**Both:**
- Follow header template exactly
- Use templates from `.github/templates/`
- Section headers: `# section name -------...` (lowercase, ~70 chars)
- Paths: `fs::path()` (R) or `pathlib.Path` (Python)
- Logging: configure but often commented for production

### Key Libraries

**R Core Stack:**
- `jutils`: utilities, magrittr pipe (always first)
- `data.table`: fast operations
- `dplyr`, `tidyr`: manipulation
- `ggplot2`, `patchwork`: visualization
- `GetoptLong`: CLI args
- `logger`: logging
- `DBI`, `duckdb`: database

**R Specialized:**
- `Seurat`: single-cell
- `ComplexHeatmap`: heatmaps
- `parallel`: mclapply (typically 20 cores)

**Python Stack:**
- `polars`: DataFrames
- `typer`: CLI
- `rich`: console output
- `duckdb`: database
- `scanpy`: single-cell (if needed)

### Color Schemes

**Always source:** `source("00-colors.R")` in analysis scripts

**Available palettes:**
- `color_disease`: Alzheimer's, COVID-19, Healthy
- `color_celltype`: B, CD4 T, CD8 T, NK, Monocyte, etc.
- `color_chemistry`: SC5P-PE, SC5P-R2, SC3Pv3, SC3Pv2
- `color_gender`: Female, Male, Unknown

**Usage:**
```r
data |>
  dplyr::mutate(disease = factor(disease, levels = names(color_disease)))

ggplot(aes(color = disease)) +
  scale_color_manual(values = color_disease)
```

### Output Conventions

**Standard paths:**
```r
basedir <- "/liulab/chunjie/data/scMOCHA"
outdir <- fs::path(basedir, "analysis/zzz/MANUSCRIPTFIGURES")
```

**Create directories:**
```r
fs::dir_create(outdir)
```

**Save plots:**
```r
ggsave(
  filename = fs::path(outdir, "01-figure-name.pdf"),
  plot = p,
  width = 8,
  height = 6
)
```

## Package Management

**R:** Use `conda activate renv` for R environments

**Python:** Use `uv` NOT pip:
```bash
uv add package-name
uv remove package-name
uv sync
```

## Reference Files

**Configuration data:** `config/`
- `rCRS.MT.fasta`: Reference mitochondrial genome
- `Homo_sapiens.GRCh38.107.gtf.id_name_length_genetype.fst`: Gene annotations
- `Mito-Genome-Loci-MitoMAP-Foswiki.fst`: Mitochondrial genome features
- `metacolors.json`: Project color schemes

**Templates:** `.github/templates/`
- `template-analysis-numbered.R`
- `template-utility-uppercase.R`
- `template-processing.R`
- `template-processing.py`

**Complete guidelines:** `.github/scmocha-agent-coding-guidelines.md` (comprehensive 1500+ line reference)

## Data Locations

**Symlinks (already configured):**
- `data/` → `/mnt/isilon/u01_project/large-scale/liuc9/raw`
- `large-scale/` → `/mnt/isilon/u01_project/large-scale`

**Base directory:** `/liulab/chunjie/data/scMOCHA`

**Output directory:** `analysis/zzz/MANUSCRIPTFIGURES`

## Pre-Flight Checklist

Before committing R scripts:
- [ ] Header with @AUTHOR, @CONTACT, @DATE, @DESCRIPTION, @VERSION
- [ ] `load_pkg(jutils)` first
- [ ] Args use `=` not `<-`
- [ ] Paths use `fs::path()`
- [ ] Source `00-colors.R` if plotting
- [ ] DuckDB shutdown with `shutdown = TRUE`

Before committing Python scripts:
- [ ] Header with encoding and metadata
- [ ] `polars` preferred over `pandas`
- [ ] `typer` for CLI
- [ ] Constants in UPPERCASE
- [ ] Type hints on functions

## Common Patterns

**Parallel processing (R):**
```r
library(parallel)
results <- mclapply(
  X = data_list,
  FUN = function(.x) process(.x),
  mc.cores = 20
)
```

**Nest-map-unnest (R):**
```r
data |>
  tidyr::nest(.by = c(groupvar), .key = "nested_data") |>
  dplyr::mutate(
    result = purrr::map(.x = nested_data, .f = \(.x) analyze(.x))
  ) |>
  tidyr::unnest(result)
```

**Batch processing (Python):**
```python
from concurrent.futures import ProcessPoolExecutor

with ProcessPoolExecutor(max_workers=20) as executor:
    results = list(executor.map(process_row, data.iter_rows(named=True)))
```

## Important Notes

- **Never** use `<-` in args sections (R)
- **Always** load `load_pkg(jutils)` first (R)
- **Never** forget `shutdown = TRUE` for DuckDB (R)
- **Prefer** `polars` over `pandas` (Python)
- **Always** order factors by color scheme levels
- **Never** hardcode personal paths
- **Use** templates from `.github/templates/`
- **Check** `.github/scmocha-agent-coding-guidelines.md` for comprehensive standards

## WDL Workflow

For batch processing scMOCHA variant calling:
- **Workflow:** `scMOCHA.batch.wdl`
- **Inputs:** `scMOCHA.batch.inputs.json`
- Default: 20 cores, 50GB memory, Docker image `chunjiesamliu/scmocha`
