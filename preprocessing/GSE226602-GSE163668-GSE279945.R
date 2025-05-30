#!/usr/bin/env Rscript
# Metainfo ----------------------------------------------------------------
# @AUTHOR: Chun-Jie Liu
# @CONTACT: chunjie.sam.liu.at.gmail.com
# @DATE: 2025-01-22 15:31:06
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


# load data ---------------------------------------------------------------
gseid_list <- c("GSE226602", "GSE163668", "GSE279945", "GSE162117")

tibble::tibble(
  gseid = c("GSE226602", "GSE163668", "GSE279945", "GSE162117"),
  chem = c("SC5P-PE", "SC5P-R2", "SC3Pv3", "SC3Pv2")
) -> gseid_list

basedir <- "/home/liuc9/github/scMOCHA-data/data"
outdir <- "/home/liuc9/github/scMOCHA-data/data/out_variant_check"
# body --------------------------------------------------------------------



# ! load data --------------------------------------------------------------------


gseid_list |>
  dplyr::mutate(
    cell_ratio_variant = purrr::map(
      gseid,
      ~ {
        data.table::fread(
          file.path(
            basedir,
            .x,
            "out",
            glue::glue("{.x}.cell_ratio_and_variant_clean.csv")
          )
        )
      }
    )
  ) |>
  dplyr::mutate(
    anno = purrr::map(
      gseid,
      ~ {
        readr::read_rds(
          file.path(
            basedir,
            .x,
            "out",
            glue::glue("{.x}.scmocha.out.rds.gz")
          )
        )
      }
    )
  ) ->
gseid_list_anno



# ! merge data --------------------------------------------------------------------



dplyr::inner_join(
  gseid_list_anno |>
    dplyr::select(-chem) |>
    tidyr::unnest(cols = cell_ratio_variant) |>
    dplyr::select(-anno),
  gseid_list_anno |>
    dplyr::select(-chem) |>
    tidyr::unnest(cols = anno),
  by = c("gseid", "srrid")
) ->
gseid_list_anno_merged


gseid_list_anno_merged |> dplyr::glimpse()
selected_srrid <- c("GSM4995425", "GSM4995448", "GSM4933442", "GSM7080044", "GSM8583898")

# selected_srrid <- gseid_list_anno_merged$srrid |>
#   unique() |>
#   as.character()

gseid_list_anno_merged |>
  dplyr::mutate(
    label = ifelse(srrid %in% selected_srrid, srrid, "")
  ) |>
  ggplot(aes(
    x = `Depth read mean`,
    y = `# of somatic variants`,
    color = chemistry,
    label = label
  )) +
  geom_point() +
  ggrepel::geom_text_repel(
    show.legend = FALSE,
  ) +
  theme_minimal() -> p
p

ggsave(
  file.path(outdir, "p_plotly_select_samples.pdf"),
  p,
  width = 8,
  height = 5
)

plotly::ggplotly(p) -> p_plotly_select_samples

p_plotly_select_samples

reticulate::py_run_string("import sys")
htmlwidgets::saveWidget(
  p_plotly_select_samples,
  file.path(outdir, "p_plotly_select_samples.html"),
  selfcontained = TRUE
)




# ! selected samples --------------------------------------------------------------------

gseid_list_anno_merged |>
  dplyr::filter(srrid %in% selected_srrid) |>
  dplyr::select(gseid, srrid, chemistry, `# of somatic variants`, `# of variants`, srrdir, somatic_variant) |>
  dplyr::mutate(
    sv = purrr::map(
      somatic_variant,
      ~ {
        .x$somatic
      }
    )
  ) -> gseid_list_anno_merged_selected





gseid_list_anno_merged_selected$sv[[1]]
gseid_list_anno_merged_selected$sv[[2]]


gseid_list_anno_merged_selected$srrdir

gseid_list_anno_merged_selected$srrid[[2]]
gseid_list_anno_merged_selected |>
  # dplyr::select(gseid, srrid, chemistry) |>
  dplyr::mutate(label = glue::glue("{srrid}-{chemistry}")) |>
  dplyr::select(label, sv) ->
gseid_list_anno_merged_selected_



gseid_list_anno_merged_selected_ |>
  tidyr::spread(key = label, value = sv) |>
  as.list() |>
  purrr::map(
    ~ {
      .x[[1]]
    }
  ) ->
variant_list

library(ggVennDiagram)
# variant_list <- list(
#   "GSE226602-GSM7080044-SC5P-PE" = gseid_list_anno_merged_selected$sv[[1]],
#   "GSE163668-GSM4995425-SC5P-R2" = gseid_list_anno_merged_selected$sv[[2]],
#   "GSE163668-GSM4995448-SC5P-R2" = gseid_list_anno_merged_selected$sv[[3]],
#   "GSE163668-GSM4995459-SC5P-R2" = gseid_list_anno_merged_selected$sv[[4]],
#   "GSE279945-GSM8583898-SC3Pv3" = gseid_list_anno_merged_selected$sv[[5]]
# )

variant_list_df <- variant_list |>
  ggVennDiagram::Venn() |>
  ggVennDiagram::process_data()

ggplot() +
  geom_path(aes(X, Y, color = id, group = id),
    data = ggVennDiagram::venn_setedge(variant_list_df),
    show.legend = FALSE
  ) +
  ggsci::scale_color_npg() +
  geom_text(aes(X, Y, label = name),
    data = ggVennDiagram::venn_setlabel(variant_list_df)
  ) +
  geom_label(aes(X, Y, label = count),
    data = ggVennDiagram::venn_regionlabel(variant_list_df)
  ) +
  coord_equal() +
  theme_void() ->
p_venn
p_venn

ggsave(
  file.path(outdir, "p_venn.pdf"),
  p_venn,
  width = 8,
  height = 5
)

variant_list


# stats --------------------------------------------------------------------

gseid_list_anno_merged_selected |>
  dplyr::mutate(
    stats = purrr::map2(
      .x = gseid,
      .y = srrid,
      ~ {
        data.table::fread(
          file.path(
            basedir,
            .x,
            "final",
            .y,
            "cell.variant_stats.tsv.gz"
          )
        )
      }
    )
  ) ->
gseid_list_anno_merged_selected_stats


gseid_list_anno_merged_selected_stats$sv[[1]]
gseid_list_anno_merged_selected_stats$stats[[1]] |>
  dplyr::mutate(
    vmr_log = log10(vmr)
  ) |>
  # dplyr::filter(
  #   variant %in% gseid_list_anno_merged_selected_stats$sv[[1]]
  # ) |>
  ggplot(aes(x = strand_correlation, y = vmr_log)) +
  geom_point() +
  geom_vline(
    xintercept = 0.65,
    linetype = 20,
    color = "red"
  ) +
  geom_hline(
    yintercept = log10(0.01),
    linetype = 20,
    color = "red"
  )


gseid_list_anno_merged_selected_stats$stats[[1]] |>
  dplyr::filter(
    strand_correlation == 1
  )


# footer ------------------------------------------------------------------

# future: :plan(future: :sequential)

# save image --------------------------------------------------------------
