#!/usr/bin/env Rscript
# Metainfo ----------------------------------------------------------------
# @AUTHOR: Chun-Jie Liu
# @CONTACT: chunjie.sam.liu.at.gmail.com
# @DATE: 2025-11-05 14:39:20
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
library(fs)

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

fn_load_count <- function(thepath, type = c("cluster", "cell")) {
  type <- match.arg(type)

  pattern <- if (type == "cluster") {
    "*cluster.*.txt.gz*"
  } else {
    "*cell.*.txt.gz*"
  }

  tibble::tibble(
    path = list.files(
      thepath,
      pattern,
      full.names = T
    )
  ) |>
    dplyr::filter(!grepl("coverage", x = path)) |>
    dplyr::mutate(
      d = parallel::mclapply(
        path,
        data.table::fread,
        mc.cores = 4
      )
    ) |>
    dplyr::mutate(n = basename(path)) |>
    dplyr::mutate(n = gsub(paste0(type, ".|.txt.gz"), "", n)) |>
    dplyr::select(n, d) |>
    tidyr::unnest(cols = d) |>
    as.data.table() |>
    dplyr::mutate(nv = V3 + V4) |>
    dplyr::select(
      gt = n,
      pos = V1,
      group = V2,
      fw = V3,
      rv = V4,
      nv
    ) -> cluster_n

  fasta <- Biostrings::readDNAStringSet(
    "/home/liuc9/github/scMOCHA/fasta/rCRS.chrM.fasta"
  )

  fasta$chrM |>
    as.data.table() |>
    tibble::rownames_to_column(var = "pos") |>
    dplyr::rename(ref = x) |>
    dplyr::mutate(posref = glue::glue("{pos}{ref}")) |>
    dplyr::mutate(pos = as.integer(pos)) -> fasta_df
  # data.table::fwrite(
  #   fasta_df,
  #   file = "/home/liuc9/github/scMOCHA-data/config/rCRS.MT.fasta.csv",
  #   sep = ","
  # )

  cluster_n |>
    dtplyr::lazy_dt() |>
    dplyr::left_join(fasta_df, by = "pos") |>
    dplyr::mutate(gt = factor(gt, levels = c("A", "G", "C", "T"))) |>
    as.data.table() -> cluster_n_temp

  cluster_n_temp[, ratio := nv / sum(nv), by = .(group, pos)]

  cluster_n_temp |>
    dplyr::mutate(
      label = glue::glue(
        "total = {nv} \n forward = {fw} \n reverse = {rv} \n ratio = {round(ratio, 3) * 100}%"
      )
    ) -> cluster_n_forplot

  cluster_n_forplot
}

# body --------------------------------------------------------------------
thepath <- "/mnt/isilon/u01_project/large-scale/ting/raw/GSE233304/Running_results/GSM8244396"


cellcount <- fn_load_count(thepath, type = "cell")
clustercount <- fn_load_count(thepath, type = "cluster")
barcode <- data.table::fread(
  "/mnt/isilon/u01_project/large-scale/ting/raw/GSE233304/Running_results/GSM8244396/barcode_cluster.tsv",
  col.names = c("barcode", "cj", "celltype")
)

cellcount |>
  dplyr::filter(pos == 3173) |>
  dplyr::select(gt, group, fw, rv, nv) |>
  dplyr::left_join(barcode, by = c("group" = "barcode")) |>
  dplyr::select(-cj) |>
  # dplyr::group_by(gt, celltype) |>
  dplyr::group_by(gt) |>
  dplyr::summarise(nv = sum(nv)) |>
  print(n = Inf)

clustercount |>
  dplyr::filter(pos == 3173, group == "CD4_T")

# footer ------------------------------------------------------------------
library(Rsamtools)

#
#
# ? bam --------------------------------------------------------------------
#
#
bam <- "/mnt/isilon/u01_project/large-scale/ting/raw/GSE235050/cromwell-executions/scMOCHABatch/66a50c5e-b7f0-4d1c-9c2f-9dec97c3a7f7/call-scMOCHA/shard-0/sub.scMOCHA/bf8e6c39-e9e6-4a4c-828e-af5935932eaf/call-call_mt_variants/execution/cell/temp/barcoded_bams/barcodes.1.sort.bam"

# which <- GRanges("MT", IRanges(3173, 3173))
what <- c(
  "qname",
  "flag",
  "rname",
  "strand",
  "pos",
  "qwidth",
  "mapq",
  "cigar",
  "seq",
  "qual"
)
tags <- c(
  "NH",
  "HI",
  "AS",
  "nM",
  "RE",
  "CR",
  "CY",
  "CB",
  "UR",
  "UY",
  "UB",
  "MU"
)
param <- ScanBamParam(tag = tags, what = scanBamWhat())
reads <- scanBam(bam, param = param)

library(Rsamtools)
reads[[1]] |>
  as.data.frame() |>
  as.data.table() -> reads_dt
target_pos <- 3173

mtchr <- "MT"
proper_pair <- TRUE # same as "True" in Python
NHmax <- 1
NMmax <- 4

reads_dt |>
  dplyr::filter(
    # 1. Chromosome / contig
    rname == mtchr
  ) |>
  dplyr::filter(
    # 2. Tag filters: NH <= NHmax, NM/nM <= NMmax
    (is.na(tag.NH) | as.numeric(tag.NH) <= NHmax) &
      (is.na(tag.nM) |
        as.numeric(tag.nM) <= NMmax |
        as.numeric(tag.nM) <= NMmax)
  ) |>
  dplyr::filter(
    # 3. Proper pair if requested
    bitwAnd(flag, 0x2) != 0
  )

# reads_dt |>
#   dplyr::mutate(
#     read_end      = pos + nchar(seq) - 1,
#     covers_target = pos <= target_pos & read_end >= target_pos,
#     offset        = target_pos - pos + 1,
#   ) |>
#   dplyr::filter(
#   covers_target, offset > 0, offset <= nchar(seq)
# ) |>
# dplyr::mutate(
#   base      = substr(seq, offset, offset),
#     base_qual = utf8ToInt(substr(qual, offset, offset)) - 33L
# ) |>
#   dplyr::filter(!is.na(tag.CB)) |>
#   dplyr::filter(!is.na(tag.CJ)) |>
#   dplyr::count(base)

# future: :plan(future: :sequential)

# save image --------------------------------------------------------------
gse_data <- import(
  "/home/liuc9/github/scMOCHA-data/analysis/zzz/clean-data/gse_data.qs"
)
