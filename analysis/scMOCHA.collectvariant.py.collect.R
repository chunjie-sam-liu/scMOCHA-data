#!/usr/bin/env Rscript
# Metainfo ----------------------------------------------------------------
# @AUTHOR: Chun-Jie Liu
# @CONTACT: chunjie.sam.liu.at.gmail.com
# @DATE: 2025-06-10 17:24:07
# @DESCRIPTION: filename
# @VERSION: v0.0.1

# Library -----------------------------------------------------------------

suppressPackageStartupMessages(library(magrittr))
library(ggplot2)
library(patchwork)
library(prismatic)
library(paletteer)
library(data.table)
# library(rlang)
library(GetoptLong)
library(logger)
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

# src ---------------------------------------------------------------------

# header ------------------------------------------------------------------

# future: :plan(future: :multisession, workers = 10)

# load data ---------------------------------------------------------------
srr_filename <- "/home/liuc9/github/scMOCHA-data/analysis/zzz/clean-data/gse_srrid_srrdir.csv"
srr <- import(srr_filename)
# function ----------------------------------------------------------------
fn_strandbias <- function(d) {
  d |>
    # head(20) |>
    tidyr::replace_na(list(
      AF = 0,
      AR = 0,
      CF = 0,
      CR = 0,
      GF = 0,
      GR = 0,
      TF = 0,
      TR = 0
    )) |>
    dplyr::mutate(
      strandbias = parallel::mcmapply(
        variant = variant,
        variant_type = variant_type,
        AF = AF,
        AR = AR,
        CF = CF,
        CR = CR,
        GF = GF,
        GR = GR,
        TF = TF,
        TR = TR,
        FUN = function(variant, variant_type, AF, AR, CF, CR, GF, GR, TF, TR) {
          if (variant_type != "colorful") {
            return(
              tibble::tibble(
                pvalue = NA_real_,
                strand_ratio = NA_real_
              )
            )
          }
          tryCatch(
            expr = {
              ref <- gsub("\\d*|>.*", "", variant)
              alt <- gsub(".*>", "", variant)

              # rf <- ifelse(is.na(get(paste0(ref, "F"))), 0, get(paste0(ref, "F")))
              # rr <- ifelse(is.na(get(paste0(ref, "R"))), 0, get(paste0(ref, "R")))
              # af <- ifelse(is.na(get(paste0(alt, "F"))), 0, get(paste0(alt, "F")))
              # ar <- ifelse(is.na(get(paste0(alt, "R"))), 0, get(paste0(alt, "R")))
              rf <- get(paste0(ref, "F"))
              rr <- get(paste0(ref, "R"))
              af <- get(paste0(alt, "F"))
              ar <- get(paste0(alt, "R"))

              table <- matrix(c(rf, rr, af, ar), nrow = 2, byrow = T)
              colnames(table) <- c("Forward", "Reverse")
              rownames(table) <- c("Ref", "Alt")
              result <- fisher.test(table)
              strand_ratio <- max(af, ar) / (af + ar)
              return(
                tibble::tibble(
                  pvalue = result$p.value,
                  strand_ratio = strand_ratio
                )
              )
            },
            error = \(e) {
              return(
                tibble::tibble(
                  pvalue = NA_real_,
                  strand_ratio = NA_real_
                )
              )
            }
          )
        },
        mc.cores = 20,
        SIMPLIFY = FALSE
      )
    ) |>
    tidyr::unnest(cols = strandbias) -> d_strandbias

  d_strandbias
}

fn_strandbias_sequential <- function(d) {
  d |>
    # head(20) |>
    tidyr::replace_na(list(
      AF = 0,
      AR = 0,
      CF = 0,
      CR = 0,
      GF = 0,
      GR = 0,
      TF = 0,
      TR = 0
    )) |>
    dplyr::mutate(
      strandbias = purrr::pmap(
        list(
          variant = variant,
          variant_type = variant_type,
          AF = AF,
          AR = AR,
          CF = CF,
          CR = CR,
          GF = GF,
          GR = GR,
          TF = TF,
          TR = TR
        ),
        function(variant, variant_type, AF, AR, CF, CR, GF, GR, TF, TR) {
          if (variant_type != "colorful") {
            return(
              tibble::tibble(
                pvalue = NA_real_,
                strand_ratio = NA_real_
              )
            )
          }
          tryCatch(
            expr = {
              ref <- gsub("\\d*|>.*", "", variant)
              alt <- gsub(".*>", "", variant)

              rf <- get(paste0(ref, "F"))
              rr <- get(paste0(ref, "R"))
              af <- get(paste0(alt, "F"))
              ar <- get(paste0(alt, "R"))

              table <- matrix(c(rf, rr, af, ar), nrow = 2, byrow = T)
              colnames(table) <- c("Forward", "Reverse")
              rownames(table) <- c("Ref", "Alt")
              result <- fisher.test(table)
              strand_ratio <- max(af, ar) / (af + ar)
              return(
                tibble::tibble(
                  pvalue = result$p.value,
                  strand_ratio = strand_ratio
                )
              )
            },
            error = \(e) {
              return(
                tibble::tibble(
                  pvalue = NA_real_,
                  strand_ratio = NA_real_
                )
              )
            }
          )
        }
      )
    ) |>
    tidyr::unnest(cols = strandbias) -> d_strandbias

  d_strandbias
}
# body --------------------------------------------------------------------
srr |>
  # head(2) |>
  dplyr::mutate(
    load = parallel::mclapply(
      srrdir,
      function(srrdir) {
        import(
          path(
            srrdir,
            "variant_info_from_heatmap.qs"
          )
        ) -> d
        dd <- fn_strandbias_sequential(d)
        rm(d)
        gc()
        dd
      },
      mc.cores = 30
    )
  ) -> srr_load


srr_load |>
  dplyr::filter(
    purrr::map_lgl(load, ~ all(class(.x) == "try-error"))
  ) |>
  nrow()

srr_load |>
  dplyr::mutate(
    load = purrr::map2(
      .x = srrdir,
      .y = load,
      .f = function(srrdir, load) {
        # srrdir <- srr_load$srrdir[[457]]
        # load <- srr_load$load[[457]]
        if (all(class(load) == "try-error")) {
          import(
            path(
              srrdir,
              "variant_info_from_heatmap.qs"
            )
          ) -> d
          dd <- fn_strandbias(d)
          rm(d)
          gc()
          return(dd)
        } else {
          return(load)
        }
      }
    )
  ) -> srr_load_


# srr_load |>
#   dplyr::mutate

srr_load_ |>
  dplyr::select(-srrdir) |>
  tidyr::unnest(cols = load) |>
  dplyr::rename(
    AFO = AF,
    ARE = AR,
    CFO = CF,
    CRE = CR,
    GFO = GF,
    GRE = GR,
    TFO = TF,
    TRE = TR,
    fisher_test_pvalue = pvalue,
    alt_strand_ratio = strand_ratio
  ) -> srr_load_unnest

export(
  srr_load_unnest,
  file = "/home/liuc9/github/scMOCHA-data/analysis/zzz/clean-data/all_variant_cell.fishertest.csv",
  format = "both",
)
export(
  srr_load_unnest |> as.data.table(),
  file = "/home/liuc9/github/scMOCHA-data/analysis/zzz/clean-data/all_variant_cell.fishertest.qs",
)


# srr_load_unnest <- import("/home/liuc9/github/scMOCHA-data/analysis/zzz/clean-data/all_variant_cell.qs")

conn <- DBI::dbConnect(
  duckdb::duckdb(),
  dbdir = "/mnt/isilon/u01_project/large-scale/liuc9/raw/zzz/clean-data/all_hetero_af.cell.duckdb.1.2.1" |>
    glue::glue()
)
DBI::dbListTables(conn)
DBI::dbWriteTable(
  conn,
  "allvariants_cell_fishertest",
  srr_load_unnest,
  temporary = FALSE,
  overwrite = TRUE
)

DBI::dbDisconnect(conn, shutdown = TRUE)

#
#
# ? update variant_type --------------------------------------------------------------------
#
#
conn <- DBI::dbConnect(
  duckdb::duckdb(),
  dbdir = "/mnt/isilon/u01_project/large-scale/liuc9/raw/zzz/clean-data/all_hetero_af.cell.duckdb.1.2.1" |>
    glue::glue()
)

allvariants_cell_fishertest <- dplyr::tbl(
  conn,
  "allvariants_cell_fishertest"
) |>
  data.table::as.data.table()

allvariants_cell_fishertest |>
  # head(2000000) |>
  dplyr::mutate(
    variant_type_fisher_test = parallel::mcmapply(
      variant = variant,
      variant_type = variant_type,
      AFO = AFO,
      ARE = ARE,
      CFO = CFO,
      CRE = CRE,
      GFO = GFO,
      GRE = GRE,
      TFO = TFO,
      TRE = TRE,
      fisher_test_pvalue = fisher_test_pvalue,
      alt_strand_ratio = alt_strand_ratio,
      FUN = function(
        variant,
        variant_type,
        AFO,
        ARE,
        CFO,
        CRE,
        GFO,
        GRE,
        TFO,
        TRE,
        fisher_test_pvalue,
        alt_strand_ratio
      ) {
        if (variant_type != "colorful") {
          return(variant_type)
        }
        ref <- gsub("\\d*|>.*", "", variant)
        alt <- gsub(".*>", "", variant)

        reff <- switch(ref, A = AFO, C = CFO, G = GFO, T = TFO)
        refr <- switch(ref, A = ARE, C = CRE, G = GRE, T = TRE)
        altf <- switch(alt, A = AFO, C = CFO, G = GFO, T = TFO)
        altr <- switch(alt, A = ARE, C = CRE, G = GRE, T = TRE)

        if (fisher_test_pvalue < 0.05) {
          return("black")
        } else {
          if (alt_strand_ratio < 0.1 | alt_strand_ratio > 0.9) {
            return("black")
          } else {
            if (altf >= 2 & altr >= 2) {
              return("colorful")
            } else {
              return("black")
            }
          }
        }
      },
      mc.cores = 50,
      SIMPLIFY = TRUE
    )
  ) -> allvariants_cell_fishertest_varianttype


DBI::dbWriteTable(
  conn,
  "allvariants_cell_fishertest",
  allvariants_cell_fishertest_varianttype,
  temporary = FALSE,
  overwrite = TRUE
)
DBI::dbDisconnect(conn, shutdown = TRUE)

# footer ------------------------------------------------------------------

# future: :plan(future: :sequential)
