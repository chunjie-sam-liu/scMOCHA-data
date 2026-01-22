#!/usr/bin/env Rscript
# Metainfo ----------------------------------------------------------------
# @AUTHOR: Chun-Jie Liu
# @CONTACT: chunjie.sam.liu.at.gmail.com
# @DATE: 2026-01-21 22:09:20
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


# src ---------------------------------------------------------------------

# header ------------------------------------------------------------------

# load data ---------------------------------------------------------------
load_pkg(jutils)
dotenv()
color_celltype <- c(
  "homo hetero" = "blue",
  "Bimodal" = "red",
  "Mono" = "#A6D854FF",
  "B" = "#66C2A5FF",
  "T cell" = "#E5C494FF",
  "CD8 T" = "#8DA0CBFF",
  "NK" = "#FFD92FFF",
  "other" = "black"
)


thevariants <- tibble(
  variant = c(
    "709G>A", # homo hetero
    "14905G>A", # Bimodal
    "8138A>G", # Bimodal
    "2011G>A", # Bimodal
    "7751T>C", # Bimodal
    "7609T>C", # Bimodal
    "4813T>C", # Mono
    "7159T>C", # B cell
    "7833T>C", # T cell
    "10500G>A", # T cell
    "10097A>G", # T cell
    "8005T>C", # CD8 T
    "7850G>A", # CD8 T
    "9033A>G", # NK
    "7757G>A", # NK
    "9390A>G", # NK
    "6374T>C", # NK
    "10236A>G", # NK
    "1474G>A", # NK
    "9609T>C", # NK
    "2636G>A", # NK
    "15612G>A", # NK
    "2343G>A", # NK
    "7837T>C", # NK
    "6928T>C", # NK
    "2666T>C" # NK
  ),
  type = c(
    "homo hetero", # homo hetero
    "Bimodal",
    "Bimodal",
    "Bimodal",
    "Bimodal",
    "Bimodal", # Bimodal
    "Mono", # Mono
    "B", # B cell
    "T cell",
    "T cell",
    "T cell", # T cell
    "CD8 T",
    "CD8 T", # CD8 T
    "NK",
    "NK",
    "NK",
    "NK",
    "NK",
    "NK",
    "NK",
    "NK",
    "NK",
    "NK",
    "NK",
    "NK",
    "NK" # NK
  )
)

ALLVARIANTS_TEST <- import(
  outdir / "VARIANT-KRUSKAL-WALLIS-TEST.xlsx"
)

ALLVARIANTS_TEST_SIG <- ALLVARIANTS_TEST |>
  dplyr::filter(p.value < 0.05, statistic > 20)


# load conn ---------------------------------------------------------------

# function ----------------------------------------------------------------

# body --------------------------------------------------------------------

source(
  path(
    Sys.getenv("HIGHRESDIR"),
    "plot_celltype_specific_variant.R"
  )
)
source(
  path(
    Sys.getenv("HIGHRESDIR"),
    "plot_individual_proportion.R"
  )
)

thevariants |>
  left_join(ALLVARIANTS_TEST_SIG, by = "variant") |>
  mutate(
    label = gsub(".>.", "", "*{variant}*{srrid}*" |> glue())
  ) |>
  pull(label) |>
  paste0(collapse = "; ")
print(n = Inf)

# footer ------------------------------------------------------------------

# save image --------------------------------------------------------------
