#!/usr/bin/env Rscript
# Metainfo ----------------------------------------------------------------
# @AUTHOR: Chun-Jie Liu
# @CONTACT: chunjie.sam.liu.at.gmail.com
# @DATE: 2025-06-11 11:48:13
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
library(glue)
library(parallel)
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
conn <- DBI::dbConnect(
  duckdb::duckdb(),
  dbdir = "/home/liuc9/github/scMOCHA-data/analysis/zzz/clean-data/all_variant_cell.duckdb.1.2.1"
)

# function ----------------------------------------------------------------


# body --------------------------------------------------------------------
dplyr::tbl(conn, "all_variant_cell") |>
  dplyr::filter(
    variant_in_cell_cluster == "cell"
  ) |>
  dplyr::select(
    gseid, srrid, variant
  ) |>
  dplyr::distinct() |>
  as.data.table() ->
gseid_srrid_variant

all_variant_cell_table <- dplyr::tbl(conn, "all_variant_cell")

gseid_srrid_variant |>
  # head(100) |>
  dplyr::mutate(
    co = parallel::mcmapply(
      .x = srrid,
      .y = variant,
      FUN = \(.x, .y) {
        # .x <- "GSM4762179"
        # .y <- "11251A>G"

        log_trace(
          glue::glue(
            "Processing variant {.y} for srrid {.x}"
          )
        )
        all_variant_cell_table |>
          dplyr::filter(
            variant == .y,
            srrid == .x,
            variant_in_cell_cluster == "cell"
          ) |>
          dplyr::select(
            barcode, af, depth, variant_type, celltype
          ) |>
          as.data.table() ->
        .d
        .d |>
          dplyr::group_by(celltype) |>
          dplyr::summarise(sum_depth = sum(depth, na.rm = TRUE), mean_depth = mean(depth, na.rm = T)) ->
        .dd
        log_trace("has data in database ", nrow(.d))
        .d |>
          dplyr::count(
            celltype, variant_type
          ) |>
          dplyr::left_join(
            .dd,
            by = "celltype"
          )
      },
      mc.cores = 20,
      SIMPLIFY = FALSE
    )
  ) ->
gseid_srrid_variant_co


gseid_srrid_variant_co |>
  tidyr::unnest(cols = co) |>
  tidyr::pivot_wider(
    names_from = variant_type,
    values_from = c(n),
  ) |>
  tidyr::nest(
    .by = c(gseid, srrid, variant),
    .key = "variant_celltype"
  ) ->
gseid_srrid_variant_celltype



gseid_srrid_variant_celltype |>
  dplyr::mutate(
    n_colorful = parallel::mcmapply(
      .x = variant_celltype,
      FUN = \(.x) {
        .x |>
          tidyr::pivot_longer(
            cols = -c(celltype, sum_depth, mean_depth),
            names_to = "group",
            values_to = "n",
          ) |>
          dplyr::mutate(
            n = ifelse(
              n >= 4,
              n,
              NA_real_
            )
          ) |>
          dplyr::filter(
            !is.na(n)
          ) |>
          dplyr::count(group) |>
          tidyr::pivot_wider(
            names_from = group,
            values_from = n,
            names_prefix = "n_"
          )
      },
      mc.cores = 20,
      SIMPLIFY = FALSE
    )
  ) |>
  tidyr::unnest(n_colorful) ->
gseid_srrid_variant_celltype_n


gseid_srrid_variant_celltype_n |>
  dplyr::filter(
    !is.na(n_black),
    n_black == 8,
    n_colorful < 2
  ) |>
  # dplyr::slice(6) |>
  dplyr::filter(
    srrid == "GSM7080031"
  ) |>
  tidyr::unnest(cols = variant_celltype)



# ? real somatic mutation --------------------------------------------------------------------

ALLVARIANTS <- import(file.path(
  "/home/liuc9/github/scMOCHA-data/analysis/zzz/clean-data/", "all_variant.qs"
)) |>
  dplyr::filter(
    issomatic == "heteroplasmic"
  )

META <- import("/home/liuc9/github/scMOCHA-data/analysis/zzz/clean-data/gse_dataset_metadata_full.qs") |>
  dplyr::select(gseid, srrid, Age_new, Age_group)


gseid_srrid_variant_celltype_n |>
  dplyr::filter(
    # srrid == "GSM7080031"
    variant %in% ALLVARIANTS$variant
  ) |>
  dplyr::filter(
    n_black >= 6,
    n_colorful < 6
  ) ->
somatic_variants

export(
  somatic_variants,
  file = file.path(
    "/home/liuc9/github/scMOCHA-data/analysis/zzz/clean-data/",
    "real_somatic_variant_celltype.qs"
  )
)

META |>
  dplyr::left_join(
    somatic_variants |>
      dplyr::count(gseid, srrid)
  ) |>
  dplyr::filter(
    Age_group != "Unknown"
  ) |>
  dplyr::mutate(
    n = ifelse(
      is.na(n),
      0,
      n
    )
  ) |>
  ggpubr::ggscatter(
    x = "Age_new", y = "n",
    color = "black", shape = 20, size = 3, # Points color, shape and size
    add = "loess", # Add regressin line
    add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
    conf.int = TRUE, # Add confidence interval
    cor.coef = TRUE, # Add correlation coefficient. see ?stat_cor
    cor.coeff.args = list(method = "pearson", label.x = 3, label.sep = "\n"),
    jitter = 0.2,
    xlab = "Age (years)", ylab = "Number of somatic variants",
  ) ->
p_real_somatic_variants_age

p_real_somatic_variants_age |>
  ggplot2::ggsave(
    filename = file.path(
      "/home/liuc9/github/scMOCHA-data/analysis/zzz/plot-real-somatic-variant",
      "real_somatic_variants_age.pdf"
    ),
    width = 8, height = 6
  )



# ? somatic variant example --------------------------------------------------------------------

somatic_variants |>
  dplyr::filter(
    n_colorful <= 2
  ) |>
  dplyr::slice(2) |>
  tidyr::unnest(cols = variant_celltype)

# thevariant <- "7757G>A"
# thesrrid <- "GSM7437874"

thevariant <- "6967G>A"
thesrrid <- "GSM7080026"

all_variant_cell_table |>
  dplyr::filter(
    srrid == thesrrid,
    variant == thevariant
  ) |>
  dplyr::collect() |>
  dplyr::mutate(
    variant_type = dplyr::case_match(
      variant_type,
      "colorful" ~ "red",
      "black" ~ "darkblue",
      "white" ~ "white",
      "grey" ~ "gray",
      NA ~ "white"
    )
  ) |>
  dplyr::mutate(
    variant_type = factor(
      variant_type,
      levels = c("red", "darkblue", "gray", "white")
    )
  ) |>
  dplyr::arrange(
    variant_type,
    -af
  ) ->
forplot_



forplot_ |>
  dplyr::mutate(
    barcode = factor(
      barcode,
      levels = forplot_$barcode
    )
  ) ->
forplot
source("analysis/00-colors.R")

thetheme <- theme(
  panel.background = element_blank(),
  panel.grid = element_blank(),
  axis.ticks = element_blank(),
  axis.text = element_blank(),
  axis.title.x = element_blank(),
)


forplot |>
  dplyr::mutate(
    celltype = gsub(
      "_",
      " ",
      celltype
    )
  ) |>
  dplyr::mutate(
    celltype = factor(
      celltype,
      names(color_celltype)
    )
  ) |>
  ggplot(aes(
    x = barcode,
    y = 1,
    fill = celltype
  )) +
  geom_col() +
  scale_fill_manual(
    name = "Cell Type",
    values = color_celltype,
  ) +
  thetheme +
  labs(
    y = "Cell Type",
  ) ->
p1_celltype

forplot |>
  dplyr::mutate(
    af = ifelse(
      af < 0.01,
      NA_real_,
      af
    )
  ) |>
  ggplot(aes(
    x = barcode,
    y = af,
    fill = af
  )) +
  geom_col() +
  scale_fill_gradient(
    name = "Allele Frequency",
    high = "red",
    low = "white"
  ) +
  scale_y_continuous(
    expand = expansion(mult = c(0, 0)),
  ) +
  theme(
    panel.background = element_blank(),
    panel.grid = element_blank(),
    axis.line = element_line(colour = "black"),
    axis.ticks.x = element_blank(),
    axis.text.x = element_blank(),
    axis.title.x = element_blank(),
  ) +
  labs(
    y = "Allele Frequency",
  ) ->
p2_af



forplot |>
  dplyr::mutate(
    variant_type = as.character(variant_type),
  ) |>
  ggplot(aes(
    x = barcode,
    y = 1,
    fill = variant_type
  )) +
  geom_col() +
  scale_fill_identity(
    guide = "legend",
    name = "Variant cell",
    breaks = c("red", "darkblue", "gray", "white"),
    labels = c("Heteroplasmy", "Suficcient reads", "No sufficient reads", "No reads")
  ) +
  thetheme +
  labs(
    y = "Variant cells",
  ) ->
p3_variant_cells


forplot |>
  dplyr::mutate(
    depth = log2(depth + 1) # log2 transform to reduce skewness
  ) |>
  ggplot(aes(
    x = barcode,
    y = depth,
    fill = depth
  )) +
  geom_col() +
  scale_fill_gradient(
    name = "log2(depth + 1)",
    high = "gold",
    low = "white"
  ) +
  scale_y_continuous(
    expand = expansion(mult = c(0, 0)),
  ) +
  theme(
    panel.background = element_blank(),
    panel.grid = element_blank(),
    axis.line = element_line(colour = "black"),
    axis.ticks.x = element_blank(),
    axis.text.x = element_blank(),
    axis.title.x = element_blank(),
  ) +
  labs(
    y = "Log2(Depth + 1)",
  ) ->
p4_depth

wrap_plots(
  p2_af,
  p4_depth,
  p3_variant_cells,
  p1_celltype,
  ncol = 1,
  heights = c(15, 15, 10, 10),
  guides = "collect"
) +
  plot_annotation(
    title =
      glue::glue(
        "Variant {thevariant} in {thesrrid}"
      ),
    theme = theme(
      plot.title = element_text(hjust = 0.5, size = 16, face = "bold")
    )
  ) ->
p_all
p_all

ggsave(
  p_all,
  filename = file.path(
    "/home/liuc9/github/scMOCHA-data/analysis/zzz/plot-real-somatic-variant",
    glue::glue("somatic_variant_{thevariant}_{thesrrid}.pdf")
  ),
  width = 12, height = 8
)
# footer ------------------------------------------------------------------

# future: :plan(future: :sequential)

# save image --------------------------------------------------------------
