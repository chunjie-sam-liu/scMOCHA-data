# Worked openxlsx2 example

One complete workbook, built from synthetic data, end to end: four funnel data
sheets with colored column-group headers, a summary sheet of five stacked
tables, a sheet appended after the file was already saved, and the readback
that verifies it. Every block below was executed before being written down.

Run the sections in order; they share state.

Colors are written as literal hexes below only so the workbook builds
standalone. In real code every one of them comes from the track `color.R`.
-> `palette`

---

## 0. Setup and synthetic input

```r
suppressMessages({
  library(jutils)
  load_pkg(openxlsx2)
})

set.seed(1)
n <- 120
TRAITS <- c("EF", "GLS", "LAVI", "CM")

full <- data.table::data.table(
  trait   = sample(TRAITS, n, TRUE),
  ID      = sprintf("chr%d:%d:A:T", sample(1:22, n, TRUE),
                    sample(1e6:2.4e8, n, TRUE)),
  gene    = paste0("GENE", seq_len(n)),
  n_disc  = sample(500:700, n, TRUE),
  b_disc  = rnorm(n, 0.3, 0.2),
  se_disc = runif(n, 0.05, 0.15),
  p_disc  = 10^-runif(n, 5, 10),
  n_rep   = sample(100:130, n, TRUE),
  b_rep   = rnorm(n, 0.2, 0.3),
  se_rep  = runif(n, 0.08, 0.25),
  p_rep   = runif(n, 1e-4, 1),
  maf_ref = c(rep(0, 20), runif(n - 20, 0, 0.4)),
  p_ref   = runif(n, 1e-6, 1)
)
full[, run := paste0(ifelse(trait == "CM", "cm.", "echo."), tolower(trait))]
full[, CHR := as.integer(sub("^chr([0-9]+):.*$", "\\1", ID))]
full[, POS := as.integer(sub("^chr[0-9]+:([0-9]+):.*$", "\\1", ID))]
full[, tier := data.table::fcase(
  p_disc < 5e-8, "genome_wide",
  p_disc < 1e-6, "strong",
  default = "suggestive"
)]
full[, ref_status := data.table::fcase(
  maf_ref < 0.01, "absent",
  p_ref > 0.05, "null",
  sign(b_disc) == sign(b_rep), "shared",
  default = "opposite"
)]
```

---

## 1. Column layout as blocks

The block list is simultaneously the column order, the rename map, and the
header-color map. Renaming a source column is then a one-line edit in one
place.

```r
ACCENT    <- "#3C5488"
HDR_FILL  <- "#4D4D4D"
NOTE_GREY <- "#7F7F7F"
WHITE     <- "#FFFFFF"
COHORT_COLORS <- c(Discovery = "#E64B35", Replication = "#F39B7F",
                   Reference = "#374E55")

BLOCKS <- list(
  Identity = c(Trait = "trait", Run = "run", Tier = "tier",
               Variant_ID = "ID", CHR = "CHR", POS = "POS", Gene = "gene"),
  Discovery = c(DISC_N = "n_disc", DISC_BETA = "b_disc",
                DISC_SE = "se_disc", DISC_P = "p_disc"),
  Replication = c(REP_N = "n_rep", REP_BETA = "b_rep",
                  REP_SE = "se_rep", REP_P = "p_rep"),
  Reference = c(REF_MAF = "maf_ref", REF_P = "p_ref", REF_status = "ref_status")
)
BLOCK_FILL <- c(
  Identity    = HDR_FILL,
  Discovery   = COHORT_COLORS[["Discovery"]],
  Replication = COHORT_COLORS[["Replication"]],
  Reference   = COHORT_COLORS[["Reference"]]
)

col_map <- unlist(BLOCKS, use.names = FALSE)
names(col_map) <- unlist(lapply(BLOCKS, names), use.names = FALSE)
wide <- full[, ..col_map]
data.table::setnames(wide, names(col_map))
data.table::setorder(wide, DISC_P)

P_COLS   <- c("DISC_P", "REP_P", "REF_P")
NUM_COLS <- c("DISC_BETA", "DISC_SE", "REP_BETA", "REP_SE", "REF_MAF")
INT_COLS <- c("POS", "DISC_N", "REP_N")
```

---

## 2. Funnel stages and summary tables

Define the stage masks once. Every summary table is then a function of them, so
no two tables can disagree.

```r
keep_rep <- wide$REP_P < 0.05
keep_ref <- wide$REF_status %in% c("absent", "null")
keep_dir <- sign(wide$DISC_BETA) == sign(wide$REP_BETA)

sheets <- list(
  `01_All_leads`            = wide,
  `02_Replicated`           = wide[keep_rep],
  `03_Not_sig_in_reference` = wide[keep_rep & keep_ref],
  `04_Priority`             = wide[keep_rep & keep_ref & keep_dir]
)
FUNNEL_SHEETS <- names(sheets)

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
  for (nm in STAGE_COLS) data.table::set(base, which(is.na(base[[nm]])), nm, 0L)
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

tbl_funnel <- data.table::data.table(
  Step = c("01 All discovery leads (P < 1e-5)",
           "02 Replicated (P < 0.05)",
           "03 + Not significant in reference",
           "04 + Same effect direction (PRIORITY)"),
  Sheet = FUNNEL_SHEETS,
  N_rows = vapply(STAGE_FLAGS, sum, integer(1)),
  N_unique_variants = n_uniq("Variant_ID"),
  Pct_of_all_leads = vapply(STAGE_FLAGS, sum, integer(1)) / nrow(wide)
)
tbl_trait <- with_total(count_vec(wide$Trait, "Trait",
                                  lvl_of(wide$Trait, TRAITS)))
tbl_tier <- with_total(count_vec(
  wide$Tier, "Tier",
  lvl_of(wide$Tier, c("genome_wide", "strong", "suggestive"))
))
tbl_ref <- with_total(count_vec(
  wide$REF_status, "Reference_status",
  lvl_of(wide$REF_status, c("absent", "null", "shared", "opposite"))
))
tbl_defs <- data.table::data.table(
  Column = c("DISC_P", "REP_P", "REF_MAF", "REF_status"),
  Definition = c("discovery P", "replication P at the exact lead",
                 "reference-cohort MAF", "absent / null / shared / opposite")
)
```

`Pct_of_all_leads` is stored as a fraction. The `0.0%` number format renders it
as a percentage; multiplying by 100 here would turn it into text.

---

## 3. Style helpers

`openxlsx2` has no reusable `createStyle()` object in the high-level API, so
the equivalent of "declare the styles once" is "declare these two helpers
once". Both rely on the R6 form mutating `wb` in place.

```r
wb <- wb_workbook()
SUM <- "00_Summary"
wb$add_worksheet(SUM)          # added first so it is the leftmost tab

fmt_hdr <- function(sheet, dims, fill = HDR_FILL) {
  wb$add_fill(sheet = sheet, dims = dims, color = wb_color(hex = fill))
  wb$add_font(sheet = sheet, dims = dims,
              color = wb_color(hex = WHITE), bold = TRUE)
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

---

## 4. The summary sheet

`add_block()` returns the next free row, so blocks chain and cannot overlap.

```r
add_block <- function(row, title, d, note = NULL, total_row = FALSE,
                      pct = NULL, sci = NULL) {
  cell <- wb_dims(rows = row, cols = 1)
  wb$add_data(sheet = SUM, x = title, dims = cell)
  wb$add_font(sheet = SUM, dims = cell, bold = TRUE, size = 12,
              color = wb_color(hex = ACCENT))
  row <- row + 1L

  if (!is.null(note)) {
    cell <- wb_dims(rows = row, cols = 1)
    wb$add_data(sheet = SUM, x = note, dims = cell)
    wb$add_font(sheet = SUM, dims = cell, italic = TRUE, size = 9,
                color = wb_color(hex = NOTE_GREY))
    row <- row + 1L
  }

  wb$add_data(sheet = SUM, x = d, dims = wb_dims(x = d, from_row = row))
  fmt_hdr(SUM, wb_dims(x = d, from_row = row, select = "col_names"))

  if (nrow(d) > 0) {
    num <- names(d)[vapply(d, is.numeric, logical(1))]
    fmt_num(SUM, d, setdiff(num, c(pct, sci)), "#,##0", from_row = row)
    fmt_num(SUM, d, pct, "0.0%", from_row = row)
    fmt_num(SUM, d, sci, "0.00E+00", from_row = row)

    if (total_row) {
      tr <- wb_dims(rows = row + nrow(d), cols = seq_along(d))
      wb$add_font(sheet = SUM, dims = tr, bold = TRUE)
      wb$add_border(sheet = SUM, dims = tr,
                    top_border = "thin", top_color = wb_color(hex = HDR_FILL),
                    bottom_border = NULL, left_border = NULL,
                    right_border = NULL)
    }
  }

  row + nrow(d) + 3L
}

r <- 1L
r <- add_block(r, "Table 1. Prioritization funnel (overall)", tbl_funnel,
               note = "One row = one (run, lead variant) pair.",
               pct = "Pct_of_all_leads")
r <- add_block(r, "Table 2. Funnel by trait", tbl_trait, total_row = TRUE)
r <- add_block(r, "Table 3. Funnel by discovery significance tier", tbl_tier,
               note = "genome_wide P < 5e-8; strong P < 1e-6; suggestive P < 1e-5.",
               total_row = TRUE)
r <- add_block(r, "Table 4. Reference-cohort status", tbl_ref, total_row = TRUE)
r <- add_block(r, "Definitions", tbl_defs)

wb$set_col_widths(sheet = SUM, cols = 1, widths = 60)
wb$set_col_widths(sheet = SUM, cols = 2:12, widths = 18)
wb$freeze_pane(sheet = SUM, first_active_row = 2)
```

Setting the unwanted borders to `NULL` is how the total row gets a top-only
rule; the defaults draw all four sides.

---

## 5. Data sheets with per-block colored headers

```r
for (nm in names(sheets)) {
  d <- sheets[[nm]]
  wb$add_worksheet(nm)
  wb$add_data(sheet = nm, x = d, with_filter = TRUE)

  if (nm %in% FUNNEL_SHEETS) {
    for (blk in names(BLOCKS)) {
      fmt_hdr(nm,
              wb_dims(x = d, cols = names(BLOCKS[[blk]]), select = "col_names"),
              fill = BLOCK_FILL[[blk]])
    }
  } else {
    fmt_hdr(nm, wb_dims(x = d, select = "col_names"))
  }

  if (nrow(d) > 0) {
    fmt_num(nm, d, P_COLS, "0.00E+00")
    fmt_num(nm, d, NUM_COLS, "0.0000")
    fmt_num(nm, d, INT_COLS, "#,##0")
  }

  wb$freeze_pane(sheet = nm, first_active_row = 2, first_active_col = 2)
  wb$set_col_widths(sheet = nm, cols = seq_along(d), widths = "auto")
  j <- match("Variant_ID", names(d))
  if (!is.na(j)) wb$set_col_widths(sheet = nm, cols = j, widths = 22)
}

out_xlsx <- "priority.xlsx"
wb$save(out_xlsx, overwrite = TRUE)
log_info("Workbook written: {out_xlsx}")
for (nm in names(sheets)) log_info("  {nm}: {nrow(sheets[[nm]])} rows")
```

`wb_dims(x = d, cols = names(BLOCKS[[blk]]), select = "col_names")` resolves a
block's column names to one contiguous header range, so each block gets its
fill in a single call. This requires the block columns to be adjacent, which
building `wide` from `BLOCKS` guarantees.

Log output from this run:

```
Workbook written: priority.xlsx
  01_All_leads: 120 rows
  02_Replicated: 6 rows
  03_Not_sig_in_reference: 6 rows
  04_Priority: 5 rows
```

---

## 6. Appending a sheet after the file was saved

Load, drop the sheet if it is already there, re-add. Never rebuild the book.

```r
power_tbl <- data.table::data.table(
  Run = unique(wide$Run),
  N_leads = as.integer(table(wide$Run)[unique(wide$Run)]),
  Power_80_pct = runif(length(unique(wide$Run)), 20, 95)
)

wb2 <- wb_load(out_xlsx)
if ("05_Power" %in% wb_get_sheet_names(wb2)) wb2$remove_worksheet("05_Power")
wb2$add_worksheet("05_Power")
wb2$add_data(sheet = "05_Power", x = power_tbl, with_filter = TRUE)
wb2$add_fill(sheet = "05_Power",
             dims = wb_dims(x = power_tbl, select = "col_names"),
             color = wb_color(hex = HDR_FILL))
wb2$add_font(sheet = "05_Power",
             dims = wb_dims(x = power_tbl, select = "col_names"),
             color = wb_color(hex = WHITE), bold = TRUE)
# already on a 0-100 scale, so "0.0" and not "0.0%"
wb2$add_numfmt(sheet = "05_Power",
               dims = wb_dims(x = power_tbl, cols = "Power_80_pct"),
               numfmt = "0.0")
wb2$freeze_pane(sheet = "05_Power", first_active_row = 2, first_active_col = 2)
wb2$save(out_xlsx, overwrite = TRUE)
```

Removing before adding makes the script idempotent - a re-run replaces the
sheet instead of erroring on a duplicate name. `wb_get_sheet_names()` takes the
loaded workbook object, not a file path.

---

## 7. Verify before reporting success

```r
wbv <- wb_load(out_xlsx)
s <- unname(wb_get_sheet_names(wbv))
cat("sheets:", paste(s, collapse = ", "), "\n")
for (nx in s) cat(" ", nx, nrow(wb_to_df(out_xlsx, sheet = nx)), "rows\n")

print(head(
  wb_to_df(out_xlsx, sheet = "01_All_leads", apply_numfmts = TRUE)[
    , c("Variant_ID", "DISC_N", "DISC_BETA", "DISC_P", "REF_MAF")
  ], 3
))
```

Actual output of this run:

```
sheets: 00_Summary, 01_All_leads, 02_Replicated, 03_Not_sig_in_reference, 04_Priority, 05_Power
  00_Summary 41 rows
  01_All_leads 120 rows
  02_Replicated 6 rows
  03_Not_sig_in_reference 6 rows
  04_Priority 5 rows
  05_Power 4 rows

           Variant_ID DISC_N DISC_BETA   DISC_P REF_MAF
2   chr9:68304626:A:T    598    0.3669 1.28E-10  0.3169
3  chr18:48441924:A:T    615    0.4830 1.67E-10  0.0297
4 chr20:159790242:A:T    631   -0.0016 1.71E-10  0.2093
```

`apply_numfmts = TRUE` renders cells the way Excel will. It is the only check
that proves `DISC_P` is scientific rather than a rounded `0.00`, and it is what
turns "the script exited 0" into "the workbook is correct".

The summary sheet reads back the same way, with the percentage column rendered:

```
                                         A                       B      C                 D                E
1 Table 1. Prioritization funnel (overall)                    <NA>   <NA>              <NA>             <NA>
2  One row = one (run, lead variant) pair.                    <NA>   <NA>              <NA>             <NA>
3                                     Step                   Sheet N_rows N_unique_variants Pct_of_all_leads
4        01 All discovery leads (P < 1e-5)            01_All_leads    120               120           100.0%
5                 02 Replicated (P < 0.05)           02_Replicated      6                 6             5.0%
6        03 + Not significant in reference 03_Not_sig_in_reference      6                 6             5.0%
7    04 + Same effect direction (PRIORITY)             04_Priority      5                 5             4.2%
```
