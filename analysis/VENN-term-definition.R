#!/usr/bin/env Rscript
# Metainfo ----------------------------------------------------------------
# @AUTHOR: Chun-Jie Liu
# @CONTACT: chunjie.sam.liu.at.gmail.com
# @DATE: 2025-11-29 18:34:50
# @DESCRIPTION: filename
# @VERSION: v0.0.1

# Library -----------------------------------------------------------------

load_pkg(
  ggplot2,
  patchwork,
  prismatic,
  paletteer,
  data.table,
  glue,
  parallel,
  GetoptLong,
  scales,
  fs,
  jutils,
  logger
)

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

# load data ---------------------------------------------------------------

gse_data <- import(
  "/home/liuc9/github/scMOCHA-data/analysis/zzz/clean-data/gse_data.qs"
)


somatic_variants <- import(
  file.path(
    "/home/liuc9/github/scMOCHA-data/analysis/zzz/clean-data/",
    "real_somatic_variant_celltype.qs"
  )
)
# somatic_variants$variant |>
#   unique() |>
#   length()
# load conn ---------------------------------------------------------------

conn <- DBI::dbConnect(
  duckdb::duckdb(),
  "/home/liuc9/github/scMOCHA-data/analysis/zzz/clean-data/all_hetero_af.cell.duckdb.1.2.1"
)


DBI::dbListTables(conn)


# src ---------------------------------------------------------------------

# function ----------------------------------------------------------------

# body --------------------------------------------------------------------
gse_data |>
  dplyr::select(gseid, srrid, anno) |>
  tidyr::unnest(cols = c(anno)) |>
  dplyr::mutate(variant = glue::glue("{Position}{Ref}>{Alt}")) |>
  dplyr::select(variant, Haplogroup) |>
  dplyr::distinct() |>
  dplyr::filter(Haplogroup != "") |>
  dplyr::filter(!is.na(Haplogroup)) |>
  tidyr::nest(.by = variant) -> dt_gse_anno


tbl_meta <- dplyr::tbl(conn, "meta")

tbl_allvariants <- dplyr::tbl(
  conn,
  "allvariants"
)

dt_meta <- as.data.table(tbl_meta)
dt_allvariants <- as.data.table(tbl_allvariants)

dt_allvariants |>
  dplyr::filter(
    issomatic %in% c("homoplasmic", "heteroplasmic")
  ) |>
  dplyr::mutate(
    ishaplo = ifelse(
      variant %in% dt_gse_anno$variant,
      TRUE,
      FALSE
    )
  ) |>
  dplyr::mutate(
    isrealsomatic = ifelse(
      variant %in% somatic_variants$variant,
      TRUE,
      FALSE
    )
  ) -> dt_allvariants_


v_ethnicity <- dt_allvariants_ |>
  dplyr::filter(ishaplo) |>
  dplyr::pull(variant) |>
  unique()
v_homo <- dt_allvariants_ |>
  dplyr::filter(issomatic == "homoplasmic") |>
  dplyr::pull(variant) |>
  unique()
v_hete <- dt_allvariants_ |>
  dplyr::filter(issomatic == "heteroplasmic") |>
  dplyr::pull(variant) |>
  unique()
v_somatic <- dt_allvariants_ |>
  dplyr::filter(isrealsomatic) |>
  dplyr::pull(variant) |>
  unique()

library(ggvenn)
ggvenn(
  data = list(
    Ethnicity = v_ethnicity,
    Homoplasmic = v_homo,
    Heteroplasmic = v_hete,
    Somatic = v_somatic
  )
)

library(eulerr)
fruits[, 1:3] |> class()

set.seed(1)
plot(
  euler(
    c(
      Ethnicity = length(v_ethnicity),
      Homoplasmic = length(v_homo),
      Heteroplasmic = length(v_hete),
      Somatic = length(v_somatic)
    ),
    shape = "ellipse"
  ),
  quantities = TRUE
)

dt_allvariants_ |>
  dplyr::filter(ish)

dt_allvariants_ |>
  dplyr::mutate(
    highaf = ifelse(af > .95, TRUE, FALSE),
  ) |>
  dplyr::mutate(
    is_hete = ifelse(
      !highaf,
      TRUE,
      FALSE
    )
  ) |>
  dplyr::mutate(
    is_homo = ifelse(
      highaf,
      TRUE,
      FALSE
    )
  ) |>
  dplyr::select(
    Ethnicity = ishaplo,
    Homoplasmic = is_homo,
    Heteroplasmic = is_hete,
    Somatic = isrealsomatic
  ) -> dt_allvariants_euler

plot(
  euler(
    dt_allvariants_euler[, c(
      "Ethnicity",
      "Homoplasmic",
      "Heteroplasmic",
      "Somatic"
    )],
    # shape = "ellipse",
    control = list(extraopt = FALSE)
  ),
  quantities = list(type = c("counts"), font = 3),
  labels = list(fontfamily = "serif"),
  edges = list(lty = 3),
  fills = c("#BEBADAFF", "#8DD3C7FF", "#FFFFB3FF", "red")
)
RColorBrewer::brewer.pal(10, "Set3") |> color()

# footer ------------------------------------------------------------------

# save image --------------------------------------------------------------
