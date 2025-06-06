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
  .filename <- "/home/liuc9/github/scMOCHA-data/analysis/zzz/clean-data/all_hetero_af.cell/all_hetero_af.cell.{.gseid_srrid}.qs" |> glue::glue()
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
    ) ->
  d_nest

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
            error = \(e){
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
    dplyr::arrange(p.value) ->
  d_ks

  .filename_out <- "/home/liuc9/github/scMOCHA-data/analysis/zzz/clean-data/all_hetero_af.cell.ks_test/{.gseid_srrid}.ks_test.qs" |> glue::glue()
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
  ) ->
gseid_srrid
# gseid_srrid |>
#   dplyr::mutate(
#     a = purrr::map(
#       gseid_srrid,
#       fn_ks_test
#     )
#   )

gseid_srrid |>
  # head(10) |>
  dplyr::mutate(
    load = purrr::map(
      gseid_srrid,
      \(.x) {
        import(
          file.path(
            "/home/liuc9/github/scMOCHA-data/analysis/zzz/clean-data/all_hetero_af.cell.ks_test",
            glue::glue("{.x}.ks_test.qs")
          )
        )
      }
    )
  ) |>
  tidyr::unnest(cols = load) |>
  dplyr::arrange(p.value) ->
gseid_srrid_ks_load

export(
  gseid_srrid_ks_load,
  file.path(
    "/home/liuc9/github/scMOCHA-data/analysis/zzz/clean-data/all_hetero_af.cell.ks_test",
    "a_gseid_srrid_ks_load.qs"
  )
)

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
  )

variants_df <- data.frame(
  short = c(
    "3571insC", "3664G>A", "3916G>A", "4063G>A", "4142G>A",
    "5910G>A", "5949G>A", "6124T>C", "6253T>C", "6340C>T",
    "6567C>T", "6663A>G", "6924G>T", "8932C>T", "10398A>G",
    "11778G>A", "11872insC", "12425delA", "13802C>T", "13937delAC",
    "14429delG", "15342insT", "3243A>G", "8344A>G", "8993T>G"
  ),
  detailed = c(
    "OC m.3571insC", "PC, OC m.3664G>A", "OC m.3916G>A", "OC m.4063G>A", "PC m.4142G>A",
    "PC m.5910G>A", "PC m.5949G>A", "PC m.6124T>C", "PC m.6253T>C", "PC 6340C>T",
    "OC m.6567C>T", "PC m.6663A>G", "PC m.6924G>T", "PC m.8932C>T", "PC m.10398A>G",
    "LHON m.11778G>A", "OC m.11872insC", "OC m.12425delA", "PC m.13802C>T", "OC m.13937delAC",
    "OC m.14429delG", "OC m.15342insT", "MELAS syndrome m.3243A>G", "MERRF syndrome m.8344A>G", "NARP m.8993T>G"
  ),
  note = c(
    "", "", "", "", "",
    "", "", "", "", "",
    "", "", "", "", "",
    "LHON", "", "", "", "",
    "", "", "MELAS syndrome", "MERRF syndrome", "NARP"
  ),
  stringsAsFactors = FALSE
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
  ) |>
  dplyr::count(variant) |>
  dplyr::arrange(desc(n)) |>
  print(n = 20)


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
d <- import("/home/liuc9/github/scMOCHA-data/analysis/zzz/clean-data/all_hetero_af.cell.ks_test/GSE226602_GSM7080017.ks_test.qs")

d ->
d_nest

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
          error = \(e){
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
  dplyr::arrange(p.value) ->
d_ks

d_ks$celltype_af[[1]] -> .x
source("/home/liuc9/github/scMOCHA-data/analysis/00-colors.R")

.x |>
  dplyr::filter(af > 0) |>
  dplyr::mutate(celltype = gsub(
    "_",
    " ",
    celltype
  )) |>
  dplyr::mutate(
    celltype = factor(celltype, levels = names(color_celltype)),
  ) ->
.xx

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
