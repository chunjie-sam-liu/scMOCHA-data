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
