#!/usr/bin/env Rscript
# Metainfo ----------------------------------------------------------------
# @AUTHOR: Chun-Jie Liu
# @CONTACT: chunjie.sam.liu.at.gmail.com
# @DATE: 2026-03-25 13:27:43
# @DESCRIPTION: this script is used for ...

# Reproducibility ----------------------------------------------------------
set.seed(1)
# Library -----------------------------------------------------------------

suppressMessages({
  load_pkg(jutils)
})

# Args --------------------------------------------------------------------

# s: string, i: integer, f: float, !: boolean, @: array, %: list
VERSION = "v0.0.1"

GetoptLong.options(help_style = "two-column")

# default: default value specified here.

nthread = 8
GetoptLong(
  "nthread=i",
  "Number of threads to use",
  "verbose",
  "Enable verbose logging"
)


# Logger ------------------------------------------------------------------

log_layout(layout_glue_colors)

if (isTRUE(verbose)) {
  log_threshold(TRACE)
  log_info("Verbose mode enabled")
} else {
  log_threshold(INFO)
}


# Load data ---------------------------------------------------------------
load_pkg(jutils)
dotenv()
suppressMessages({
  conflicted::conflict_prefer("filter", "dplyr")
})


# Source ---------------------------------------------------------------------

# Conn ---------------------------------------------------------------

# Function ----------------------------------------------------------------

# Main --------------------------------------------------------------------
variant_disease_reported <- c(
  "961T>C",
  "13271T>C",
  "4175G>A",
  "14831G>A",
  "9804G>A",
  "1382A>C",
  "9025G>A",
  "8362T>G",
  "3667T>G",
  "7065G>A"
)
variant_pathogenic <- c(
  "3572G>G",
  "3727T>C",
  "3728C>T",
  "5343T>C",
  "6967G>A",
  "7583T>G",
  "1111T>G",
  "14520C>A",
  "15666T>C"
)

variant_vus <- c("1670A>G", "2636G>A", "3734A>G")

thevariants <- c(
  "4864C>T",
  "6190G>A",
  "6667C>A",
  "6668C>G",
  "6997T>G",
  "7847G>A",
  "8361G>T",
  "8849T>C",
  "11372G>A",
  "11708A>T"
)
# Save  --------------------------------------------------------------

# Session info -------------------------------------------------------------
if (isTRUE(verbose)) {
  sessionInfo()
}
