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

tmp <- tempfile(fileext = ".xml")
rvg::dml_xlsx(file = tmp, fonts = list(sans = "Bradley Hand"))
print(p)
dev.off()

# Add rvg to the workbook
wb$add_worksheet("add_drawing")$add_drawing(xml = tmp)$add_drawing(
  xml = tmp,
  dims = NULL
)

library(mschart) # mschart >= 0.4 for openxlsx2 support

## create chart from mschart object (this creates new input data)
mylc <- ms_linechart(
  data = browser_ts,
  x = "date",
  y = "freq",
  group = "browser"
)

wb$add_worksheet("add_mschart")$add_mschart(dims = "A10:G25", graph = mylc)


## create chart referencing worksheet cells as input
# write data starting at B2
wb$add_worksheet("add_mschart - wb_data")$add_data(
  x = mtcars,
  dims = "B2"
)$add_data(x = data.frame(name = rownames(mtcars)), dims = "A2")

# create wb_data object this will tell this mschart
# from this PR to create a file corresponding to openxlsx2
dat <- wb_data(wb, dims = "A2:G10")

# create a few mscharts
scatter_plot <- ms_scatterchart(
  data = dat,
  x = "mpg",
  y = c("disp", "hp")
)

bar_plot <- ms_barchart(
  data = dat,
  x = "name",
  y = c("disp", "hp")
)

area_plot <- ms_areachart(
  data = dat,
  x = "name",
  y = c("disp", "hp")
)

line_plot <- ms_linechart(
  data = dat,
  x = "name",
  y = c("disp", "hp"),
  labels = c("disp", "hp")
)

# add the charts to the data
wb$add_mschart(dims = "F4:L20", graph = scatter_plot)$add_mschart(
  dims = "F21:L37",
  graph = bar_plot
)$add_mschart(dims = "M4:S20", graph = area_plot)$add_mschart(
  dims = "M21:S37",
  graph = line_plot
)

# add chartsheet
wb$add_chartsheet()$add_mschart(graph = scatter_plot)

Sys.setenv(R_ZIPCMD = Sys.which("zip"))

wb$save("mtcars.xlsx")

wb_dims(x = mtcars)
