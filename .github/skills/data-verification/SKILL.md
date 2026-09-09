---
name: data-verification
description: "Verify data schemas before writing code and verify generated files before reporting success. Use when writing parse, merge, join, filter, or rename logic; when column names, sample IDs, or key formats are involved; when inspecting TSV/CSV/Parquet/Excel/RDS/qs2/VCF/PLINK/compressed inputs; or when confirming that an output file is fresh, non-empty, and correct after a script or job finishes. Covers cheapest-possible header checks, join-key format checks, row-count checks, staleness checks, and how to tell a genuine empty result from a silent failure."
---

# Data verification

Two failure modes dominate data work: assuming a schema that does not exist,
and trusting an exit code that does not mean success. Both are cheap to
prevent and expensive to discover late.

---

## 1. Verify the schema before writing code

**Never assume column names, types, or key formats.** Inspect the real file
first, using the cheapest check that answers the question.

| Format            | Check                                                                                                                                                                                                                               |
| ----------------- | ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| TSV / CSV         | `head -1 file.tsv`, `head -3 file.tsv`, `awk -F'\t' 'NR==1{print; exit}' file.tsv`                                                                                                                                                  |
| Compressed text   | `zcat file.tsv.gz \| head -3`                                                                                                                                                                                                       |
| Parquet           | `duckdb -c "describe select * from 'file.parquet'"`, or `pixi run Rscript -e 'print(arrow::open_dataset(f)$schema)'`; for rows, `dplyr::collect(head(arrow::open_dataset(f), 3))`                                                   |
| Excel             | `pixi run Rscript -e 'readxl::excel_sheets(f)'`, then `readxl::read_excel(f, n_max = 3)`                                                                                                                                            |
| RDS               | `pixi run Rscript -e 'str(readRDS(f), max.level=1)'`                                                                                                                                                                                |
| qs2               | `pixi run Rscript -e 'str(qs2::qs_read(f), max.level=1)'`                                                                                                                                                                           |
| VCF               | `bcftools view -h file.vcf.gz \| tail -1`; the same for `.vcf.bgz`. Also `bcftools view -h f \| grep -m1 '##reference'` for the build, and `bcftools query -l f \| head -3` for the sample-ID format                                |
| PLINK             | `head -2 file.fam`, `head -2 file.bim`, `head -1 file.psam`                                                                                                                                                                         |
| IDAT              | `pixi run Rscript -e 'x <- illuminaio::readIDAT(f); str(x[c("ChipType","nSNPsRead")])'`; one file, never a directory                                                                                                                |
| Matrix (beta / M) | shape and keys, not values: `ncol` equals the expected array count, `colnames` are the row-unique sample key, `nrow` equals the probe universe; `pixi run Rscript -e 'm <- qs2::qs_read(f); str(dimnames(m), max.level=1); dim(m)'` |
| Database          | list tables, then `describe`/`PRAGMA table_info` on the target                                                                                                                                                                      |

Rules:

- Never hard-code a column name before seeing it in the actual header.
- `arrow::read_parquet(f, n_max = 3)` is not a cheap peek: `n_max` is silently
  ignored and the whole file is read. Use `open_dataset()` for the schema.
  `readRDS` and `qs_read` also load the entire object, so on a large file check
  the header or schema instead.
- Before any join, check the key format on **both** sides. The same identifier
  is often written differently by different tools (`FID_IID` vs `IID`,
  `chr1:100:A:G` vs `1_100_A_G`, with or without a `chr` prefix).
- PLINK-style headers start with `#`; strip it before use
  (`setnames(d, sub("^#", "", names(d)))`).
- Confirm real filenames with `ls` rather than assuming a tool's extension.
- Check row counts (`wc -l`) so a "small" file is not silently truncated.

## 2. Verify the output after writing it

A zero exit code does not guarantee a correct file. Before reporting success:

- **Fresh, not stale** — confirm this run wrote the file, not a previous one.
  Compare mtimes with `stat -c '%Y %n' <input> <output>`; the output must be
  the newer one. Use `ls -l` or `stat -c '%y %n'` for a human-readable check.
- **Not empty** — `wc -l`, `wc -c`, `ls -l`, or a quick `head`. Zero rows or
  zero bytes is a red flag.
- **Empty must be explained** — decide whether it is a genuine empty result
  (a filter legitimately matched nothing) or an empty caused by a crash, wrong
  path, or upstream failure. Read the logs and exit status to tell them apart.
  Never report an empty file as success without confirming the cause.
- **Text outputs** — inspect the first few lines.
- **Tabular outputs** — verify row counts and representative values, not just
  file existence.
- **Batch / array jobs** — count the produced files against the expected count
  (`ls <pattern> | wc -l`). A DONE status can mask a swallowed failure.

## 3. Report the actual check

When reporting back, state the exact command that was run and its actual
result — a row count, an exit code, the head of the output. "Should work" is
not verification. If a check was skipped, say so explicitly.
