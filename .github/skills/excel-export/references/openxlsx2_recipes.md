# openxlsx2 recipes

Concrete helpers for styled workbooks. All code here was executed against
`openxlsx2` 1.28 before being written down. Copy and adapt; do not re-derive.

Everything assumes the R6 form, where `wb$method()` mutates `wb` in place:

```r
suppressMessages(library(openxlsx2))
wb <- wb_workbook()
```

If you use the `wb_*()` wrapper form instead, every helper below must take `wb`
as an argument and return it, because wrappers are copy-on-write.

---

## Style helpers

`openxlsx2` has no reusable `createStyle()` object in the high-level API.
Styling is applied per cell range through `add_fill` / `add_font` /
`add_cell_style` / `add_border` / `add_numfmt`. Declare these two helpers once
at the top of the script and the rest of the code stays short.

```r
HDR_FILL <- "#4D4D4D"     # in real code: color_xlsx_hdr, from the track color.R
ACCENT   <- "#3C5488"     # in real code: the project accent from color.R
NOTE_GREY <- "#7F7F7F"    # in real code: color_xlsx_note
WHITE    <- "#FFFFFF"     # in real code: color_xlsx_white

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
```

`fmt_num()` takes the data frame, not just the sheet name, because `wb_dims()`
resolves the column name against `x = d`. That is what keeps the formatting
tied to column names rather than positions.

Number format codes:

```r
P_FMT   <- "0.00E+00"   # P, Q, heterogeneity P
NUM_FMT <- "0.0000"     # beta, SE, AF, MAF, call rate
INT_FMT <- "#,##0"      # N, MAC, carriers, POS
PCT_FMT <- "0.0%"       # a fraction in 0-1
PCT1_FMT <- "0.0"       # already scaled 0-100, no % sign (e.g. power)
```

---

## `add_block()` - stacked titled tables on one sheet

Returns the next free row, so blocks chain and cannot overlap. Leaves three
blank rows between blocks.

```r
SUM <- "00_Summary"
wb$add_worksheet(SUM)

add_block <- function(row, title, d, note = NULL, total_row = FALSE,
                      pct = NULL, sci = NULL) {
  cell <- wb_dims(rows = row, cols = 1)
  wb$add_data(sheet = SUM, x = title, dims = cell)
  wb$add_font(
    sheet = SUM, dims = cell,
    bold = TRUE, size = 12, color = wb_color(hex = ACCENT)
  )
  row <- row + 1L

  if (!is.null(note)) {
    cell <- wb_dims(rows = row, cols = 1)
    wb$add_data(sheet = SUM, x = note, dims = cell)
    wb$add_font(
      sheet = SUM, dims = cell,
      italic = TRUE, size = 9, color = wb_color(hex = NOTE_GREY)
    )
    row <- row + 1L
  }

  wb$add_data(sheet = SUM, x = d, dims = wb_dims(x = d, from_row = row))
  fmt_hdr(SUM, wb_dims(x = d, from_row = row, select = "col_names"))

  if (nrow(d) > 0) {
    num <- names(d)[vapply(d, is.numeric, logical(1))]
    fmt_num(SUM, d, setdiff(num, c(pct, sci)), INT_FMT, from_row = row)
    fmt_num(SUM, d, pct, PCT_FMT, from_row = row)
    fmt_num(SUM, d, sci, P_FMT, from_row = row)

    if (total_row) {
      tr <- wb_dims(rows = row + nrow(d), cols = seq_along(d))
      wb$add_font(sheet = SUM, dims = tr, bold = TRUE)
      wb$add_border(
        sheet = SUM, dims = tr,
        top_border = "thin", top_color = wb_color(hex = HDR_FILL),
        bottom_border = NULL, left_border = NULL, right_border = NULL
      )
    }
  }

  row + nrow(d) + 3L
}
```

`setdiff(num, c(pct, sci))` means a new count column is formatted correctly
without touching the helper, while the named percent and scientific columns
keep their own format. Setting the unwanted borders to `NULL` is how you get a
top-only rule; the defaults draw all four sides.

Usage and sheet finish:

```r
r <- 1L
r <- add_block(r, "Table 1. Prioritization funnel (overall)", tbl_funnel,
               note = "One row = one (run, lead variant) pair.",
               pct = "Pct_of_all_leads")
r <- add_block(r, "Table 2. Funnel by outcome", tbl_outcome, total_row = TRUE)
r <- add_block(r, "Definitions", tbl_defs)

wb$set_col_widths(sheet = SUM, cols = 1, widths = 60)
wb$set_col_widths(sheet = SUM, cols = 2:12, widths = 18)
wb$freeze_pane(sheet = SUM, first_active_row = 2)
```

---

## Funnel-count helpers

A funnel workbook counts the same categories at each filtering stage. Define
the stage masks once; every table is then a function of them, so no table can
disagree with another.

```r
STAGE_FLAGS <- list(
  All_leads   = rep(TRUE, nrow(wide)),
  Replicated  = keep_rep,
  Ref_not_sig = keep_rep & keep_ref,
  Priority    = keep_rep & keep_ref & keep_dir
)
STAGE_COLS <- names(STAGE_FLAGS)

# Keeps the expected order but never silently drops an unexpected label, so a
# Total row always equals the number of rows the table covers.
lvl_of <- function(values, expected) {
  c(expected, setdiff(sort(unique(values)), expected))
}

count_vec <- function(values, keyname, levels) {
  base <- data.table::data.table(K = levels)
  data.table::setnames(base, "K", keyname)
  for (nm in STAGE_COLS) {
    tmp <- data.table::data.table(K = values[STAGE_FLAGS[[nm]]])[, .N, by = K]
    data.table::setnames(tmp, c(keyname, nm))
    base <- merge(base, tmp, by = keyname, all.x = TRUE, sort = FALSE)
  }
  for (nm in STAGE_COLS) {
    data.table::set(base, which(is.na(base[[nm]])), nm, 0L)
  }
  base[]
}

with_total <- function(dt) {
  keycol <- names(dt)[1L]
  tot <- data.table::as.data.table(as.list(
    dt[, vapply(.SD, \(x) as.integer(sum(x)), integer(1)), .SDcols = STAGE_COLS]
  ))
  tot[, (keycol) := "Total"]
  rbind(dt, tot[, c(keycol, STAGE_COLS), with = FALSE])
}

n_uniq <- function(col) {
  vapply(STAGE_FLAGS, \(k) data.table::uniqueN(wide[[col]][k]), integer(1))
}
```

Usage:

```r
tbl_trait <- with_total(
  count_vec(wide$Trait, "Trait", lvl_of(wide$Trait, TRAIT_LEVELS))
)
```

---

## Data sheets with per-block colored headers

```r
FUNNEL_SHEETS <- names(sheets)[1:4]

for (nm in names(sheets)) {
  d <- sheets[[nm]]
  wb$add_worksheet(nm)
  wb$add_data(sheet = nm, x = d, with_filter = TRUE)

  if (nm %in% FUNNEL_SHEETS) {
    for (blk in names(BLOCKS)) {
      fmt_hdr(
        nm,
        wb_dims(x = d, cols = names(BLOCKS[[blk]]), select = "col_names"),
        fill = BLOCK_FILL[[blk]]
      )
    }
  } else {
    fmt_hdr(nm, wb_dims(x = d, select = "col_names"))
  }

  if (nrow(d) > 0) {
    fmt_num(nm, d, P_COLS, P_FMT)
    fmt_num(nm, d, NUM_COLS, NUM_FMT)
    fmt_num(nm, d, INT_COLS, INT_FMT)
    wb$add_cell_style(
      sheet = nm, dims = wb_dims(x = d, select = "data"), vertical = "center"
    )
  }

  wb$freeze_pane(sheet = nm, first_active_row = 2, first_active_col = 2)
  wb$set_col_widths(sheet = nm, cols = seq_along(d), widths = "auto")
}
```

Apply the body-wide `vertical = "center"` before the per-column number formats
if you prefer; unlike openxlsx v1 there is no `stack =` argument, because each
`add_*` call only touches its own style attribute and leaves the rest alone.

---

## Column widths

`widths = "auto"` is right for data sheets, but a long variant ID column
stretches the sheet past the screen. Pin those:

```r
wb$set_col_widths(sheet = nm, cols = seq_along(d), widths = "auto")
for (idcol in c("Variant_ID", "ID")) {
  j <- match(idcol, names(d))
  if (!is.na(j)) wb$set_col_widths(sheet = nm, cols = j, widths = 22)
}
```

---

## Sheet order

`00_Summary` must be the first tab. Either add it before the data sheets, or
fix the order at the end:

```r
wb$set_order(c(SUM, names(sheets)))
```

---

## Text columns that must stay text

Formatted intervals, ratios with CIs, and yes/no flags are built as character
columns in R, never as Excel formulas.

```r
yesno <- function(x) data.table::fifelse(x %in% TRUE, "Yes", "No")

lab <- function(x, empty = "(not reported)") {
  x <- as.character(x)
  data.table::fifelse(is.na(x) | x == "", empty, x)
}

or_ci <- function(b, se) {
  data.table::fifelse(
    is.finite(b) & is.finite(se),
    sprintf("%.2f (%.2f-%.2f)", exp(b), exp(b - 1.96 * se), exp(b + 1.96 * se)),
    NA_character_
  )
}
```

`x %in% TRUE` is intentional: it maps `NA` to `"No"` instead of `NA`, which is
what a reviewer expects in a flag column.

---

## Appending a sheet to an existing workbook

```r
wb <- wb_load(as.character(wb_in))
for (nm in c(SHEET_NEW, SHEET_NEW_SUMMARY)) {
  if (nm %in% wb_get_sheet_names(wb)) wb$remove_worksheet(nm)
}
wb$add_worksheet(SHEET_NEW)
wb$add_data(sheet = SHEET_NEW, x = new_tbl, with_filter = TRUE)
fmt_hdr(SHEET_NEW, wb_dims(x = new_tbl, select = "col_names"))
wb$freeze_pane(sheet = SHEET_NEW, first_active_row = 2, first_active_col = 2)
wb$save(out_xlsx, overwrite = TRUE)
```

Removing before adding makes the script idempotent - a re-run replaces the
sheet instead of erroring on a duplicate name.

---

## Built-in table styling as a shortcut

For a simple sheet that needs no per-block header colors,
`wb_add_data_table()` gives banded rows, a filter, and a header style in one
call:

```r
wb$add_data_table(
  sheet = nm, x = d,
  table_style = "TableStyleLight9", with_filter = TRUE
)
```

Number formats still have to be applied afterwards with `fmt_num()`. Do not use
this for the funnel sheets - the banding fights with the block header colors.

---

## Saving

```r
out_xlsx <- as.character(outdir / "priority-v2.xlsx")
wb$save(out_xlsx, overwrite = TRUE)

log_info("Workbook written: {out_xlsx}")
for (nm in names(sheets)) log_info("  {nm}: {nrow(sheets[[nm]])} rows")
```

`wb$save()` needs a character path; an `fs_path` from `outdir / "x.xlsx"`
must be wrapped in `as.character()`.

Log the per-sheet row counts. They are the numbers to check the written file
against afterwards, with
`wb_to_df(f, sheet = n, apply_numfmts = TRUE)`.
