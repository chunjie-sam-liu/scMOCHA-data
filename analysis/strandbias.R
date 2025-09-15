#!/usr/bin/env Rscript
# Metainfo ----------------------------------------------------------------
# @AUTHOR: Chun-Jie Liu
# @CONTACT: chunjie.sam.liu.at.gmail.com
# @DATE: 2025-09-12 14:00:19
# @DESCRIPTION: filename
# @VERSION: v0.0.1

# Library -----------------------------------------------------------------

suppressPackageStartupMessages(library(magrittr))
library(ggplot2)
library(patchwork)
library(prismatic)
library(paletteer)
library(data.table)
#library(rlang)
library(glue)
library(parallel)
library(GetoptLong)
library(logger)
library(scales)

# args --------------------------------------------------------------------

# s: string, i: integer, f: float, !: boolean
# @: array
# %: hash
# default: default value specified here.
verbose <- FALSE
spec <- "
Usage: Rscript foorbar.R [options]
Options:

<verbose!> Print messages
"

GetoptLong.options(help_style = "two-column")
GetoptLong(spec, template_control = list(opt_width = 21))

# header ------------------------------------------------------------------

# future: :plan(future: :multisession, workers = 10)

# load data ---------------------------------------------------------------

# load conn ---------------------------------------------------------------

# src ---------------------------------------------------------------------

# function ----------------------------------------------------------------

# body --------------------------------------------------------------------
# rf <- 170
# rr <- 1315
# af <- 7
# ar <- 499

rf <- 167
rr <- 1296
af <- 50
ar <- 465

table <- matrix(c(rf, rr, af, ar), nrow = 2, byrow = T)

colnames(table) <- c("Forward", "Reverse")
rownames(table) <- c("Ref", "Alt")
table
result <- fisher.test(table)
result

strand_ratio <- max(af, ar) / (af + ar)
strand_ratio

# footer ------------------------------------------------------------------

# future: :plan(future: :sequential)

# save image --------------------------------------------------------------
