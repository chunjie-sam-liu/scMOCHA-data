#!/usr/bin/env Rscript
# Metainfo ----------------------------------------------------------------
# @AUTHOR: Chun-Jie Liu
# @CONTACT: chunjie.sam.liu.at.gmail.com
# @DATE: 2025-11-30 23:32:49
# @DESCRIPTION: this script is used for ...

# Library -----------------------------------------------------------------

load_pkg(jutils)

# args --------------------------------------------------------------------

# s: string, i: integer, f: float, !: boolean, @: array, %: list
GetoptLong.options(help_style = "two-column")
VERSION = "v0.0.1"

# default: default value specified here.

verbose = TRUE

GetoptLong("verbose!", "print messages")


logger::log_threshold(logger::TRACE)
logger::log_layout(logger::layout_glue_colors)
log_success("Logger is configured.")
# header ------------------------------------------------------------------

# load data ---------------------------------------------------------------

# load conn ---------------------------------------------------------------

# src ---------------------------------------------------------------------

# function ----------------------------------------------------------------

# body --------------------------------------------------------------------

# footer ------------------------------------------------------------------

# save image --------------------------------------------------------------
m <- "/scr1/users/liuc9/tmp/gse_data_variant_heteroplasmic_fisher.qs"


d <- import(
  m,
  nthreads = 8
)

export(
  d,
  file = "/scr1/users/liuc9/tmp/gse_data_variant_heteroplasmic_fisher.qs",
)


export(
  d,
  file = "/scr1/users/liuc9/tmp/gse_data_variant_heteroplasmic_fisher.new.qs",
  preset = "fast",
  nthreads = 4,
)

a1 <- import(
  "/scr1/users/liuc9/tmp/gse_data_variant_heteroplasmic_fisher.qs",
  nthreads = 4
)
a2 <- import(
  "/scr1/users/liuc9/tmp/gse_data_variant_heteroplasmic_fisher.new.qs",
  nthreads = 4
)


wb <- wb_workbook()$add_worksheet(
  "cj"
)$add_data(x = mtcars)$add_fill(
  dims = wb_dims(x = mtcars, rows = 1:5), # only 1st 5 rows of x data
  color = wb_color("yellow")
)$add_fill(
  dims = wb_dims(x = mtcars, select = "col_names"), # only column names
  color = wb_color("cyan2")
)$add_fill(
  dims = wb_dims(x = mtcars, cols = 2:3), # entire data
  color = wb_color("red")
)

p <- ggplot(mtcars, aes(x = mpg, fill = as.factor(gear))) +
  ggtitle("Distribution of Gas Mileage") +
  geom_density(alpha = 0.5)

print(p)
wb$add_worksheet("add_plot")$add_plot(
  width = 5,
  height = 3.5,
  file_type = "png",
  units = "in"
)

Sys.setenv(R_ZIPCMD = Sys.which("zip"))

wb$save("mtcars.xlsx")

wb_dims(x = mtcars)
