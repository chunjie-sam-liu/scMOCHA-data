---
name: excel-export
description: Write R code that exports .xlsx deliverables with openxlsx2. Use when a script must produce an Excel workbook, add or restyle a sheet, format P-values / betas / counts / percentages, build a multi-sheet funnel or lookup workbook, add a summary sheet with several stacked tables, freeze panes, autofilter, color column-group headers, or append a sheet to a workbook that already exists. Covers the two export tiers (jutils::export vs openxlsx2), wb_dims, the copy-on-write vs R6-chain gotcha, the number-format contract, the header and block-fill palette, sheet naming, migrating openxlsx v1 code, and the verification step before reporting success.
---

# Excel Export

An `.xlsx` is a deliverable a collaborator opens by hand. It must be readable
without reformatting: P-values in scientific notation, betas at fixed decimals,
counts with thousand separators, a frozen header row, and an autofilter.

**Use `openxlsx2`, not `openxlsx`.** Older scripts may still call the v1 API;
do not copy them. The mapping table near the bottom converts them.

This skill is portable and names no project value. The track color file and its
`color_xlsx_*` object names come from the repository's `.github/instructions/`
bindings or `AGENTS.md`.

## Pick the tier first

| Tier                  | Use when                                                                          | Tool                                                |
| --------------------- | --------------------------------------------------------------------------------- | --------------------------------------------------- |
| 1. Plain dump         | One table, internal use, no P columns, no styling needed                          | `export(df, "out.xlsx")` or `writexl::write_xlsx()` |
| 2. Styled deliverable | Anything shared, anything with P-values, anything multi-sheet, any summary tables | `openxlsx2`                                         |

```r
# Tier 1 - single table, no P columns
export(dt, as.character(outdir / "duplicates.xlsx"))

# Tier 1 - named list becomes a multi-sheet book, jutils auto-styles headers
export(list(cars = mtcars, flowers = iris), as.character(outdir / "two.xlsx"))
```

Tier 1 loses control of number formats. The moment a column holds a P-value,
move to tier 2 - `export()` renders `1.2e-09` as `0.00`.

## The one gotcha that bites first

`openxlsx2` offers two calling conventions for the same operations, and they
have **different mutation semantics**:

```r
# Wrapper form: copy-on-write. You MUST reassign or the change is lost.
wb <- wb_add_data(wb, sheet = nm, x = dt)
wb <- wb |> wb_add_worksheet(nm) |> wb_add_data(sheet = nm, x = dt)

# R6 chain form: mutates wb in place. No reassignment.
wb$add_worksheet(nm)$add_data(sheet = nm, x = dt)
```

Silently dropping a change because a `wb_*()` result was not reassigned is the
most common openxlsx2 bug.

**Prefer the R6 form for anything that uses helper functions.** A helper that
styles a region can then call `wb$add_fill(...)` and return nothing, instead of
threading `wb` in and out of every call. Use the wrapper/pipe form only for
short linear scripts.

## Tier-2 contract

1. **Number formats are assigned by column-name list, not by column position.**
   Declare `P_COLS`, `NUM_COLS`, `INT_COLS` once, then intersect them with each
   sheet's names. A sheet that lacks a column simply gets no style for it.
2. **Every data sheet gets `with_filter = TRUE`,
   `wb_freeze_pane(first_active_row = 2, first_active_col = 2)`, and
   `wb_set_col_widths(widths = "auto")`.**
3. **Sheet names are zero-padded and ordered**: `00_Summary`, `01_All_...`,
   `02_..._validated`. Add `00_Summary` first, or fix the order at the end with
   `wb_set_order()`.
4. **Every color comes from the track `color.R`** — header fills, block fills,
   note greys, cohort colors alike — never a raw hex typed into the export
   script. Load `palette` before adding one.
5. **A workbook ships with its own definitions table.** Reviewers must not have
   to ask what a column means.

### Number formats

`wb_add_numfmt(wb, sheet, dims, numfmt)` takes the OOXML format code directly.

| Meaning                                   | `numfmt`     |
| ----------------------------------------- | ------------ |
| Every P-value, Q-value, heterogeneity P   | `"0.00E+00"` |
| beta, SE, AF, MAF, call rate, correlation | `"0.0000"`   |
| N, MAC, carriers, POS, any count          | `"#,##0"`    |
| A proportion already stored as 0-1        | `"0.0%"`     |

Store a percentage as the fraction and let `0.0%` render it. Never multiply by
100 and append a `%` character - that turns a number into text.

Alignment is a separate call:
`wb_add_cell_style(dims = ..., horizontal = "right")`.

### Palette

Workbook chrome is a color like any other: it lives in the track `color.R`
under `color_xlsx_*` names, and the export script references those names. If
the names are not there yet, add them there first with the values below.

| Element               | Fill                    | Font                               |
| --------------------- | ----------------------- | ---------------------------------- |
| Default header row    | `color_xlsx_hdr`        | `color_xlsx_white`, bold, centered |
| Identity column block | `color_xlsx_hdr`        | `color_xlsx_white`                 |
| Cohort column block   | that cohort's own color | `color_xlsx_white`                 |
| Summary table title   | -                       | project accent, bold, size 12      |
| Summary table note    | -                       | `color_xlsx_note`, italic, size 9  |
| Total row             | -                       | bold, top border `color_xlsx_hdr`  |

Canonical values to seed `color.R` with: header fill `#4D4D4D`, note grey
`#7F7F7F`, header font `#FFFFFF`. **The object names above are this skill's
placeholders, not a naming mandate** - a project that already names its color
objects differently keeps its own convention, and the repository bindings say
which names are in force. Read the track's `color.R` before adding to it rather
than seeding a second vocabulary beside the existing one.

The export script loads them the same way a figure script does, from the track
it lives in:

```r
source("src/color.R")   # the color file of the track this script lives in
```

`wb_color(hex = color_xlsx_hdr)` accepts 6-digit hex with the leading `#` and
returns ARGB. `wb_color("white")` also works. A color file entry carrying an
alpha suffix (`"#374E55FF"`) must be trimmed to 6 digits first.

## `wb_dims()` is how you address cells

There is no `rows =` / `cols =` on the styling functions - everything takes
`dims`, an A1-notation range. Build it with `wb_dims()`, never by hand.

```r
wb_dims(rows = 5, cols = 1)                # "A5"
wb_dims(rows = 5:8, cols = seq_len(4))     # "A5:D8"
wb_dims(x = dt)                            # "A1:D4"  header + data
wb_dims(x = dt, select = "col_names")      # "A1:D1"  header only
wb_dims(x = dt, select = "data")           # "A2:D4"  body only
wb_dims(x = dt, cols = "P")                # "C2:C4"  one column's body
wb_dims(x = dt, from_row = 5)              # "A5:D8"  same table, moved down
wb_dims(x = dt, from_row = 5, cols = "P")  # "C6:C8"
```

`x = dt` plus `from_row` is what makes stacked summary blocks safe: the column
lookup is by name and the row offset is applied for you.

## Minimal working example

```r
P_COLS <- c("DISC_P", "REP_P", "REF_P")
NUM_COLS <- c("DISC_BETA", "DISC_SE", "REP_BETA", "REF_MAF")
INT_COLS <- c("POS", "DISC_N", "REP_N", "MAC")

HDR_FILL <- color_xlsx_hdr     # both from the track color.R
WHITE <- color_xlsx_white

wb <- openxlsx2::wb_workbook()

# R6 methods mutate wb in place, so these helpers need no reassignment.
fmt_hdr <- function(sheet, dims, fill = HDR_FILL) {
  wb$add_fill(sheet = sheet, dims = dims, color = wb_color(hex = fill))
  wb$add_font(
    sheet = sheet, dims = dims,
    color = wb_color(hex = WHITE), bold = TRUE
  )
  wb$add_cell_style(sheet = sheet, dims = dims, horizontal = "center")
  wb$add_border(
    sheet = sheet, dims = dims,
    top_color = wb_color(hex = WHITE), bottom_color = wb_color(hex = WHITE),
    left_color = wb_color(hex = WHITE), right_color = wb_color(hex = WHITE)
  )
  invisible(NULL)
}

fmt_num <- function(sheet, d, cols, numfmt, from_row = 1) {
  for (cn in intersect(cols, names(d))) {
    dd <- wb_dims(x = d, cols = cn, from_row = from_row)
    wb$add_numfmt(sheet = sheet, dims = dd, numfmt = numfmt)
    wb$add_cell_style(sheet = sheet, dims = dd, horizontal = "right")
  }
  invisible(NULL)
}

for (nm in names(sheets)) {
  d <- sheets[[nm]]
  wb$add_worksheet(nm)
  wb$add_data(sheet = nm, x = d, with_filter = TRUE)

  fmt_hdr(nm, wb_dims(x = d, select = "col_names"))
  if (nrow(d) > 0) {
    fmt_num(nm, d, P_COLS, "0.00E+00")
    fmt_num(nm, d, NUM_COLS, "0.0000")
    fmt_num(nm, d, INT_COLS, "#,##0")
  }

  wb$freeze_pane(sheet = nm, first_active_row = 2, first_active_col = 2)
  wb$set_col_widths(sheet = nm, cols = seq_along(d), widths = "auto")
}

out_xlsx <- as.character(outdir / "priority.xlsx")
wb$save(out_xlsx, overwrite = TRUE)
log_info("Workbook written: {out_xlsx}")
for (nm in names(sheets)) log_info("  {nm}: {nrow(sheets[[nm]])} rows")
```

Unlike openxlsx v1 there are no reusable `createStyle()` objects in the
high-level API. The equivalent of "declare the styles once" is "declare these
two helpers once" - which is why the R6 form matters.

## The summary sheet

A styled workbook opens on `00_Summary`: several small titled tables stacked
down one sheet, each with a one-line italic note explaining how to read it.
Build it with one `add_block()` helper that returns the next free row, so the
tables cannot collide. Full helper in
[references/openxlsx2_recipes.md](references/openxlsx2_recipes.md).

```r
r <- 1L
r <- add_block(r, "Table 1. Prioritization funnel (overall)", tbl_funnel,
               note = "One row = one (run, lead variant) pair.",
               pct = "Pct_of_all_leads")
r <- add_block(r, "Table 2. Funnel by outcome", tbl_outcome, total_row = TRUE)
r <- add_block(r, "Definitions", tbl_defs)
```

Rules for the summary sheet:

- Number the tables (`Table 1.`, `Table 2.`) and reference the data sheet each
  one summarizes.
- Any table with a `Total` row must have that total equal the rows the table
  covers. Build category levels with a `lvl_of()`-style helper that appends
  unexpected labels instead of dropping them, or the total will silently
  disagree with the sheet.
- The last block is always `Definitions`.
- Column 1 is wide (`widths = 60`), the rest fixed (`widths = 18`);
  `widths = "auto"` looks wrong on a sheet of stacked tables.

## Column-group blocks

When a wide sheet mixes cohorts or analysis stages, define the layout as a
named list of blocks mapping output name -> source column, then color each
block's header. This is what makes a 50-column lookup sheet legible.

```r
BLOCKS <- list(
  Identity = c(Variant_ID = "ID", CHR = "CHR", POS = "POS"),
  Discovery = c(DISC_N = "n_disc", DISC_BETA = "b_disc", DISC_P = "p_disc"),
  Reference = c(REF_N = "n_ref", REF_BETA = "b_ref", REF_P = "p_ref")
)
# Every fill comes from the track color.R; no hex is typed here. The names are
# this skill's placeholders; the repository bindings say which names are in force.
BLOCK_FILL <- c(
  Identity = color_xlsx_hdr,
  Discovery = color_cohort[["Discovery"]],
  Reference = color_cohort[["Reference"]]
)

col_map <- unlist(BLOCKS, use.names = FALSE)
names(col_map) <- unlist(lapply(BLOCKS, names), use.names = FALSE)
wide <- full[, ..col_map]
data.table::setnames(wide, names(col_map))

# One fmt_hdr call per block; wb_dims resolves the names to a contiguous range.
for (blk in names(BLOCKS)) {
  fmt_hdr(
    nm,
    wb_dims(x = wide, cols = names(BLOCKS[[blk]]), select = "col_names"),
    fill = BLOCK_FILL[[blk]]
  )
}
```

The block list doubles as the rename map and the column order, so renaming a
source column is a one-line edit in one place.

## Appending to an existing workbook

```r
wb <- openxlsx2::wb_load(as.character(wb_in))
if (SHEET_NEW %in% wb_get_sheet_names(wb)) {
  wb$remove_worksheet(SHEET_NEW)
}
wb$add_worksheet(SHEET_NEW)
wb$add_data(sheet = SHEET_NEW, x = new_tbl, with_filter = TRUE)
```

Removing before adding makes the script idempotent - a re-run replaces the
sheet instead of erroring on a duplicate name. Note `wb_get_sheet_names()`
takes the loaded workbook object, not a file path.

Write the result to a new versioned filename unless the user asked to update in
place - see `data-result-layout`.

## Migrating openxlsx v1 code

| openxlsx (v1)                                      | openxlsx2                                                            |
| -------------------------------------------------- | -------------------------------------------------------------------- |
| `createWorkbook()`                                 | `wb_workbook()`                                                      |
| `addWorksheet(wb, nm)`                             | `wb$add_worksheet(nm)`                                               |
| `writeData(wb, nm, dt, withFilter = TRUE)`         | `wb$add_data(sheet = nm, x = dt, with_filter = TRUE)`                |
| `writeData(wb, nm, x, startRow = r, startCol = 1)` | `wb$add_data(sheet = nm, x = x, dims = wb_dims(rows = r, cols = 1))` |
| `createStyle(numFmt = "0.00E+00")` + `addStyle`    | `wb$add_numfmt(dims = ..., numfmt = "0.00E+00")`                     |
| `createStyle(fgFill = ...)` + `addStyle`           | `wb$add_fill(dims = ..., color = wb_color(hex = ...))`               |
| `createStyle(fontColour =, textDecoration =)`      | `wb$add_font(dims = ..., color = ..., bold = TRUE)`                  |
| `createStyle(halign = "center")`                   | `wb$add_cell_style(dims = ..., horizontal = "center")`               |
| `createStyle(border = ...)`                        | `wb$add_border(dims = ..., top_border = "thin", ...)`                |
| `addStyle(..., rows, cols, gridExpand = TRUE)`     | one call with `dims = wb_dims(rows = ..., cols = ...)`               |
| `addStyle(..., stack = TRUE)`                      | not needed; each `add_*` touches only its own style attribute        |
| `freezePane(wb, nm, firstActiveRow = 2)`           | `wb$freeze_pane(sheet = nm, first_active_row = 2)`                   |
| `setColWidths(wb, nm, cols, widths = "auto")`      | `wb$set_col_widths(sheet = nm, cols = ..., widths = "auto")`         |
| `saveWorkbook(wb, f, overwrite = TRUE)`            | `wb$save(f, overwrite = TRUE)`                                       |
| `loadWorkbook(f)`                                  | `wb_load(f)`                                                         |
| `getSheetNames(f)`                                 | `wb_get_sheet_names(wb)` - takes the workbook, not a path            |
| `removeWorksheet(wb, nm)`                          | `wb$remove_worksheet(nm)`                                            |
| `read.xlsx(f, sheet = nm)`                         | `wb_to_df(f, sheet = nm)` (alias `read_xlsx()`)                      |

Two behavioural differences worth repeating: v1 `addStyle()` replaced the whole
cell style unless `stack = TRUE`, whereas each openxlsx2 `add_*` touches only
its own attribute; and v1 mutated `wb` in place from any call, whereas
openxlsx2 only does so through the `wb$method()` form.

## Before reporting success

A zero exit code is not proof the workbook is right. Run this:

```bash
pixi run Rscript -e '
  suppressMessages(library(openxlsx2))
  f <- "<path>/book.xlsx"
  wb <- wb_load(f)
  s <- unname(wb_get_sheet_names(wb))
  print(s)
  for (n in s) cat(n, nrow(wb_to_df(f, sheet = n)), "rows\n")
  print(head(wb_to_df(f, sheet = s[2], apply_numfmts = TRUE), 3))
'
```

`apply_numfmts = TRUE` renders cells the way Excel will, so it is the only
check that actually proves the P column is scientific and the percent column is
a percent.

Confirm, and state in the report:

- sheet names and order match the intended list, `00_Summary` first;
- each sheet's row count matches the in-R count that was logged;
- the P column reads back as `1.00E-09`, not `0.00`;
- an empty sheet is genuinely empty (explain why), not a failed filter;
- the file mtime is from this run.

## References

Both reference files are self-contained: every block builds its own synthetic
input and runs as written, so a recipe can be executed before it is adapted.

- [references/worked_examples.md](references/worked_examples.md) - one complete
  workbook end to end: block column layout, funnel stage masks, summary sheet
  of five stacked tables, four data sheets with colored block headers, a sheet
  appended after saving, and the readback that verifies it.
- [references/openxlsx2_recipes.md](references/openxlsx2_recipes.md) - the
  individual pieces on their own: style helpers, `add_block()`, funnel-count
  helpers, column widths, sheet order, text-column helpers, `add_data_table()`.

Upstream docs: `vignette("openxlsx2")` and `vignette("openxlsx2_style_manual")`,
mirrored at <https://janmarvin.github.io/openxlsx2/>.

Related skills: `jutils` for `export()`, `palette` before introducing any new
color, `data-result-layout` to choose the output path.
