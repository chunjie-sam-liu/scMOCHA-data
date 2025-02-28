#!/usr/bin/env Rscript
# Metainfo ----------------------------------------------------------------
# @AUTHOR: Chun-Jie Liu
# @CONTACT: chunjie.sam.liu.at.gmail.com
# @DATE: 2025-02-26 15:37:33
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
log_threshold(TRACE)
log_layout(layout_glue_colors)

# future: :plan(future: :multisession, workers = 10)

# function ----------------------------------------------------------------

fn_plot_pie <- function(.d, .colors = NULL) {
  .d |>
    dplyr::select(group = 1, n) |>
    dplyr::arrange(-n) |>
    dplyr::mutate(csum = rev(cumsum(rev(n)))) %>%
    dplyr::mutate(pos = n / 2 + dplyr::lead(csum, 1)) %>%
    dplyr::mutate(pos = dplyr::if_else(is.na(pos), n / 2, pos)) %>%
    dplyr::mutate(percentage = n / sum(n)) |>
    dplyr::mutate(group = factor(group, levels = group)) ->
  .dd

  .scalefill <- if (is.null(.colors)) {
    ggsci::scale_fill_aaas(
      name = NULL
    )
  } else {
    scale_fill_manual(
      name = NULL,
      values = .colors
    )
  }
  .scalecolor <- if (is.null(.colors)) {
    ggsci::scale_color_aaas(
      name = NULL
    )
  } else {
    scale_color_manual(
      name = NULL,
      values = .colors
    )
  }

  .dd |>
    ggplot(aes(
      x = "",
      y = n,
    )) +
    geom_bar(
      aes(fill = group),
      stat = "identity",
      width = 1,
      color = "white",
      show.legend = FALSE
    ) +
    .scalefill +
    ggrepel::geom_label_repel(
      aes(
        y = pos,
        label = glue::glue("{group}\n{n} ({scales::percent(percentage)})"),
        color = group,
      ),
      size = 6,
      nudge_x = 1,
      nudge_y = 0,
      show.legend = FALSE,
      max.overlaps = Inf,
    ) +
    .scalecolor +
    coord_polar(theta = "y", start = 0) +
    theme_void() +
    theme(
      plot.title = element_text(
        hjust = 0.5,
        size = 22,
      ),
      # legend.position = "none"
    )
}



# load data ---------------------------------------------------------------
project_source_sra <- readr::read_rds("/home/liuc9/github/scMOCHA-data/data/scfoundation/out/project_source_sra.rds.gz")

# project_source_sra |> dplyr::glimpse()

# body --------------------------------------------------------------------
project_source_sra |>
  dplyr::filter(proj_source == "GEO") |>
  dplyr::select(proj_ID, source_name) |>
  dplyr::distinct() ->
project_source_sra_proj_ID_source_name


project_source_sra_proj_ID_source_name |>
  dplyr::count(source_name) |>
  dplyr::arrange(-n) |>
  dplyr::filter(
    grepl(
      pattern = "pbmc|blood",
      x = source_name,
      ignore.case = TRUE
    )
  )


project_source_sra_proj_ID_source_name |>
  dplyr::count(source_name) |>
  dplyr::arrange(-n) |>
  tibble::rownames_to_column(var = "idx") |>
  dplyr::select(-n) |>
  data.table::fwrite("/home/liuc9/github/scMOCHA-data/data/scfoundation/out/project_source_sra_proj_ID_source_name.tissue_count.tsv", sep = "\t")


project_source_sra |>
  dplyr::select(-sample_attribute_new) |>
  dplyr::arrange(sample_attribute) |>
  DT::datatable(options = list(server = TRUE))


# project_source_sra |>
#   dplyr::select(-sample_attribute_new) |>
#   dplyr::filter(!is.na(source_name)) |>
#   dplyr::select(proj_ID, samp_ID, library_construction_protocol, sample_attribute) |>
#   dplyr::distinct() |>
#   data.table::fwrite(
#     "/home/liuc9/github/scMOCHA-data/data/scfoundation/out/project_source_sra_sample_attribute.tsv",
#     sep = "\t"
#   )


# project_source_sra_proj_ID_source_name |>
#   dplyr::filter(!is.na(source_name)) |>
#   dplyr::distinct() |>
#   data.table::fwrite(
#     "/home/liuc9/github/scMOCHA-data/data/scfoundation/out/project_source_sra_tissue_grouped.tsv",
#     sep = "\t"
#   )


tissue_grouped <- data.table::fread("/home/liuc9/github/scMOCHA-data/data/scfoundation/out/project_source_sra_tissue_grouped.tsv") |>
  dplyr::mutate(
    tissue_group = stringr::str_to_title(tissue_group)
  ) |>
  dplyr::mutate(
    tissue_group = ifelse(
      grepl("Pbmc", tissue_group, ignore.case = TRUE),
      "PBMC",
      tissue_group
    )
  ) |>
  dplyr::select(proj_ID, tissue_group) |>
  dplyr::distinct()

tissue_grouped |>
  dplyr::count(proj_ID) |>
  dplyr::filter(n > 1) |>
  dplyr::arrange(-n)

project_source_sra |>
  dplyr::left_join(
    tissue_grouped,
    by = c("proj_ID")
  ) |>
  dplyr::filter(!is.na(tissue_group)) ->
project_source_sra_tissue_grouped

readr::write_rds(
  project_source_sra_tissue_grouped,
  "/home/liuc9/github/scMOCHA-data/data/scfoundation/out/project_source_sra_tissue_grouped.rds.gz"
)

project_source_sra_tissue_grouped |>
  dplyr::select(proj_ID, tissue_group) |>
  dplyr::distinct() |>
  dplyr::count(tissue_group) |>
  fn_plot_pie() ->
tissue_grouped_pie

ggsave(
  filename = "/home/liuc9/github/scMOCHA-data/data/scfoundation/out/project_source_tissue_group_pie.pdf",
  plot = tissue_grouped_pie,
  width = 7,
  height = 6
)



# ! pbmc --------------------------------------------------------------------

project_source_sra_tissue_grouped |>
  dplyr::filter(tissue_group == "PBMC") |>
  dplyr::select(proj_ID, samp_ID, tissue_group) |>
  dplyr::distinct() ->
project_source_sra_tissue_grouped_pbmc


processed_gseid <- c(
  "GSE155673",
  "GSE157344",
  "GSE149689",
  "GSE171555",
  "GSE155223",
  "GSE163668",
  "GSE226602",
  "GSE166992",
  "GSE181279",
  # "WT",
  # plus ting,
  "GSE161354", # done
  "GSE175524", # done
  "GSE206283", # done
  "GSE226598", # done
  "GSE235050", # done
  "GSE261140", # some errors
  "GSE279945",
  # run by ting
  "GSE214865",
  "GSE220189",
  "GSE233844",
  "GSE175499",
  # scfoundation
  "GSE149313",
  "GSE154386",
  "GSE159117",
  "GSE162117",
  "GSE188632"
)

ggvenn::ggvenn(
  data = list(
    "Processed PBMC dataset" = processed_gseid,
    "GEO PBMC dataset" = project_source_sra_tissue_grouped_pbmc$proj_ID |> unique()
  ),
  show_percentage = FALSE,
  fill_color = ggsci::pal_aaas()(2) |> rev(),
  fill_alpha = 0.9,
  stroke_size = 0.5,
  set_name_size = 6,
  text_size = 8
) ->
p_venn_pbmc
ggsave(
  filename = "/home/liuc9/github/scMOCHA-data/data/scfoundation/out/project_source_tissue_group_pbmc_venn.pdf",
  plot = p_venn_pbmc,
  device = "pdf",
  width = 7,
  height = 6
)
pcc <- readr::read_tsv(file = "https://raw.githubusercontent.com/chunjie-sam-liu/chunjie-sam-liu.life/master/public/data/pcc.tsv")

project_source_sra_tissue_grouped_pbmc |>
  dplyr::filter(!proj_ID %in% processed_gseid) |> # nrow()
  dplyr::count(proj_ID) |>
  dplyr::arrange(-n) |>
  dplyr::mutate(
    proj_ID = ifelse(n < 20, "Others", proj_ID)
  ) |>
  dplyr::group_by(proj_ID) |>
  dplyr::summarize(n = sum(n)) |>
  dplyr::slice_head(n = 10) |>
  fn_plot_pie() ->
tissue_grouped_pbmc_geo_pie

ggsave(
  filename = "/home/liuc9/github/scMOCHA-data/data/scfoundation/out/project_source_tissue_group_pbmc_geo_pie.pdf",
  plot = tissue_grouped_pbmc_geo_pie,
  width = 7,
  height = 6
)

# bone marrow --------------------------------------------------------------------

project_source_sra_tissue_grouped |>
  dplyr::filter(tissue_group == "Bone Marrow") |>
  dplyr::select(proj_ID, samp_ID, tissue_group) |>
  dplyr::distinct() ->
project_source_sra_tissue_grouped_bone_marrow

project_source_sra_tissue_grouped_bone_marrow |>
  dplyr::count(proj_ID) |>
  dplyr::arrange(-n) |>
  # dplyr::mutate(
  #   proj_ID = ifelse(n < 20, "Others", proj_ID)
  # ) |>
  # dplyr::group_by(proj_ID) |>
  # dplyr::summarize(n = sum(n)) |>
  # dplyr::slice_head(n = 10) |>
  fn_plot_pie() ->
tissue_grouped_bone_marrow_geo_pie

ggsave(
  filename = "/home/liuc9/github/scMOCHA-data/data/scfoundation/out/project_source_tissue_group_bone_marrow_geo_pie.pdf",
  plot = tissue_grouped_bone_marrow_geo_pie,
  width = 7,
  height = 6
)



# footer ------------------------------------------------------------------

# future: :plan(future: :sequential)

# save image --------------------------------------------------------------
