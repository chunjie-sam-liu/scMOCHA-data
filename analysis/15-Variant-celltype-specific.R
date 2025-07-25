#!/usr/bin/env Rscript
# Metainfo ----------------------------------------------------------------
# @AUTHOR: Chun-Jie Liu
# @CONTACT: chunjie.sam.liu.at.gmail.com
# @DATE: 2025-06-05 16:28:44
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
cleandatadir <- "/home/liuc9/github/scMOCHA-data/data/zzz/clean-data"

all_heteroplasmic_af <- import(
  file.path(cleandatadir, "all_hetero_af.cluster.fst"),
) |>
  dplyr::select(-num_variants) |>
  dplyr::rename(celltype = barcode)


# function ----------------------------------------------------------------
fn_ks_test <- function(.gseid_srrid) {
  # .gseid_srrid <- "GSE226602_GSM7080017"
  .filename <- "/home/liuc9/github/scMOCHA-data/analysis/zzz/clean-data/all_hetero_af.cell/all_hetero_af.cell.{.gseid_srrid}.qs" |>
    glue::glue()
  d <- import(.filename)
  d |>
    tidyr::pivot_longer(
      cols = -c(celltype, barcode),
      names_to = "variant",
      values_to = "af"
    ) |>
    tidyr::nest(
      .by = "variant",
      .key = "celltype_af"
    ) -> d_nest

  d_nest |>
    # head(100) |>
    dplyr::mutate(
      ks_test = parallel::mclapply(
        X = celltype_af,
        FUN = function(.m) {
          .m |>
            dplyr::filter(af > 0) |>
            dplyr::filter(celltype != "other") -> .mm
          tryCatch(
            {
              kruskal.test(
                af ~ celltype,
                data = .mm
              ) -> .fit
              .fit |>
                broom::tidy()
            },
            error = \(e) {
              return(tibble::tibble(
                statistic = NA_real_,
                p.value = NA_real_,
                method = "Kruskal-Wallis test",
                parameter = NA_real_
              ))
            }
          )
        },
        mc.cores = 20
      )
    ) |>
    tidyr::unnest(ks_test) |>
    dplyr::arrange(p.value) -> d_ks

  .filename_out <- "/home/liuc9/github/scMOCHA-data/analysis/zzz/db/all_hetero_af.cell.ks_test/{.gseid_srrid}.ks_test.qs" |>
    glue::glue()
  export(
    d_ks,
    .filename_out,
  )
}

# body --------------------------------------------------------------------
#
# fn_ks_test("GSE155673_GSM4712885")

all_heteroplasmic_af |>
  dplyr::select(gseid, srrid) |>
  dplyr::distinct() |>
  dplyr::mutate(
    gseid_srrid = paste(gseid, srrid, sep = "_")
  ) -> gseid_srrid
# gseid_srrid |>
#   dplyr::mutate(
#     a = purrr::map(
#       gseid_srrid,
#       fn_ks_test
#     )
#   )

# gseid_srrid |>
#   # head(10) |>
#   dplyr::mutate(
#     load = purrr::map(
#       gseid_srrid,
#       \(.x) {
#         import(
#           file.path(
#             "/home/liuc9/github/scMOCHA-data/analysis/zzz/db/all_hetero_af.cell.ks_test",
#             glue::glue("{.x}.ks_test.qs")
#           )
#         )
#       }
#     )
#   ) |>
#   tidyr::unnest(cols = load) |>
#   dplyr::arrange(p.value) -> gseid_srrid_ks_load

# export(
#   gseid_srrid_ks_load,
#   file.path(
#     "/home/liuc9/github/scMOCHA-data/analysis/zzz/db/all_hetero_af.cell.ks_test",
#     "a_gseid_srrid_ks_load.qs"
#   )
# )

gseid_srrid |>
  # head(10) |>
  dplyr::mutate(
    load = parallel::mclapply(
      X = gseid_srrid,
      FUN = \(.x) {
        import(
          file.path(
            "/home/liuc9/github/scMOCHA-data/analysis/zzz/db/all_hetero_af.cell.ks_test",
            glue::glue("{.x}.ks_test.qs")
          )
        ) |>
          dplyr::select(-celltype_af, -parameter, -method)
      },
      mc.cores = 20
    )
  ) |>
  tidyr::unnest(cols = load) |>
  dplyr::arrange(p.value) -> gseid_srrid_ks_load


gseid_srrid_ks_load
export(
  gseid_srrid_ks_load,
  file.path(
    "/home/liuc9/github/scMOCHA-data/analysis/zzz/db/all_hetero_af.cell.ks_test/a_gseid_srrid_ks_load.nocellaf.qs"
  )
)


# ? for duckdb --------------------------------------------------------------------

b_gseid_srrid_ks_load_p0.05_s25 <- import(
  file.path(
    "/home/liuc9/github/scMOCHA-data/analysis/zzz/db/all_hetero_af.cell.ks_test/b_gseid_srrid_ks_load_p0.05_s25.qs"
  )
)

b_gseid_srrid_ks_load_p0.05_s25 |>
  tidyr::unnest(cols = celltype_af) -> b_gseid_srrid_ks_load_p0.05_s25_unnest

v <- packageVersion("duckdb")
conn <- DBI::dbConnect(
  duckdb::duckdb(),
  dbdir = "/home/liuc9/github/scMOCHA-data/analysis/zzz/db/all_hetero_af.cell.ks_test/b_gseid_srrid_ks_load_p0.05_s25.duckdb.{v}" |>
    glue::glue()
)
DBI::dbWriteTable(
  conn,
  "gseid_srrid_ks_load_p0.05_s25_unnest",
  b_gseid_srrid_ks_load_p0.05_s25_unnest,
  temporary = FALSE,
  overwrite = TRUE
)
DBI::dbDisconnect(conn, shutdown = TRUE)


# conn <- DBI::dbConnect(
#   duckdb::duckdb(),
#   dbdir = "/home/liuc9/github/scMOCHA-data/analysis/zzz/db/all_hetero_af.cell.ks_test/b_gseid_srrid_ks_load_p0.05_s25.duckdb.1.2.1"
# )
# dplyr::tbl(
#   conn,
#   "gseid_srrid_ks_load_p0.05_s25_unnest"
# )
# DBI::dbDisconnect(conn, shutdown = TRUE)

# ? don't run below --------------------------------------------------------------------

source("./analysis/00-colors.R")

gseid_srrid_ks_load |>
  dplyr::filter(p.value < 0.05) |>
  ggplot(
    aes(x = statistic)
  ) +
  geom_histogram(
    aes(y = after_stat(density)),
    bins = 100,
    fill = "grey50",
    color = "black",
    alpha = 0.5
  ) +
  geom_vline(
    xintercept = 50,
    linetype = "dashed",
    color = "red"
  )


# gseid_srrid_ks_load |>
#   dplyr::filter(
#     p.value < 0.05,
#     statistic > 100
#   ) |>
#   dplyr::filter(
#     variant %in% variants_df$short
#   )

gseid_srrid_ks_load |>
  dplyr::filter(
    p.value < 0.05,
    statistic > 100
  )


thevariant <- "7833T>C"
thevariant <- "4175G>A"
thevariant <- "2101C>A"
thevariant <- "3664G>A"
thevariant <- "3243A>G"

gseid_srrid_ks_load |>
  dplyr::filter(variant == thevariant) |>
  dplyr::filter(
    p.value < 0.05,
    # statistic > 100
  ) |>
  dplyr::slice(1) |>
  tidyr::unnest(cols = celltype_af) |>
  dplyr::filter(af > 0) |>
  dplyr::mutate(
    celltype = gsub(
      "_",
      " ",
      celltype
    )
  ) |>
  dplyr::mutate(
    celltype = factor(celltype, levels = names(color_celltype) |> rev())
  ) |>
  ggplot(aes(
    x = af,
    y = celltype,
    fill = celltype
  )) +
  ggjoy::geom_joy(
    # scale = 0.9,
    # alpha = 0.8,
    rel_min_height = 0.01,
    size = 0.1
  ) +
  scale_fill_manual(
    values = color_celltype,
    na.value = "grey50"
  ) +
  ggjoy::theme_joy()


#
#
#
#
#
#
#
d <- import(
  "/home/liuc9/github/scMOCHA-data/analysis/zzz/clean-data/all_hetero_af.cell.ks_test/GSE226602_GSM7080017.ks_test.qs"
)

d -> d_nest

d_nest |>
  # head(100) |>
  dplyr::mutate(
    ks_test = parallel::mclapply(
      X = celltype_af,
      FUN = function(.m) {
        .m |>
          dplyr::filter(af > 0) |>
          dplyr::filter(celltype != "other") -> .mm
        tryCatch(
          {
            kruskal.test(
              af ~ celltype,
              data = .mm
            ) -> .fit
            .fit |>
              broom::tidy()
          },
          error = \(e) {
            return(tibble::tibble(
              statistic = NA_real_,
              p.value = NA_real_,
              method = "Kruskal-Wallis test",
              parameter = NA_real_
            ))
          }
        )
      },
      mc.cores = 10
    )
  ) |>
  tidyr::unnest(ks_test) |>
  dplyr::filter(p.value < 0.05) |>
  dplyr::arrange(p.value) -> d_ks

d_ks$celltype_af[[1]] -> .x
source("/home/liuc9/github/scMOCHA-data/analysis/00-colors.R")

.x |>
  dplyr::filter(af > 0) |>
  dplyr::mutate(
    celltype = gsub(
      "_",
      " ",
      celltype
    )
  ) |>
  dplyr::mutate(
    celltype = factor(celltype, levels = names(color_celltype)),
  ) -> .xx

.xx |>
  ggplot(aes(
    x = af,
    y = celltype,
    fill = celltype
  )) +
  ggjoy::geom_joy() +
  scale_fill_manual(
    values = color_celltype,
    name = "Cell type"
  )


FSA::dunnTest(
  af ~ celltype,
  data = .xx,
  method = "bh"
) -> .dunn

.dunn$res |>
  dplyr::filter(P.adj < 0.05) |>
  dplyr::arrange(P.adj)

# footer ------------------------------------------------------------------

# future: :plan(future: :sequential)

# save image --------------------------------------------------------------
