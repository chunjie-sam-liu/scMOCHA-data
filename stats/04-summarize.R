#!/usr/bin/env Rscript --vanilla
# Metainfo ----------------------------------------------------------------

# @AUTHOR: Chun-Jie Liu
# @CONTACT: chunjie.sam.liu.at.gmail.com
# @DATE: Fri Sep 27 15:01:38 2024
# @DESCRIPTION: filename

# Library -----------------------------------------------------------------

suppressPackageStartupMessages(library(magrittr))
library(ggplot2)
library(patchwork)
library(prismatic)
library(paletteer)
library(data.table)
#library(rlang)
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

# future::plan(future::multisession, workers = 10)

# function ----------------------------------------------------------------


# load data ---------------------------------------------------------------
basedir <- "/mnt/isilon/u01_project/large-scale/liuc9/raw"

gseidlist <- c(
  "GSE163668",
  "GSE149689",
  "GSE155223",
  "GSE155673",
  "GSE157344",
  "GSE166992",
  "GSE171555"
)

# body --------------------------------------------------------------------

tibble::tibble(
  gse = gseidlist
) |>
  dplyr::mutate(
    meta = purrr::map(
      .x = gse,
      .f = \(.gse) {
        data.table::fread(
          file.path(
            basedir, .gse, "out",
            "{.gse}.meta.csv" |> glue::glue()
          )
        )
      }
    )
  ) ->
  metadata

metadata |>
  dplyr::mutate(
    meta = purrr::map(
      .x = meta,
      .f = \(.meta) {
        if(!"age" %in% colnames(.meta)) {
          age = NA
        }
        if(!"gender" %in% colnames(.meta)) {
          gender = NA
        }
        if(!"sex" %in% colnames(.meta)) {
          sex = NA
        }
        if(!"race" %in% colnames(.meta)) {
          race = NA
        }
        .meta |>
          dplyr::mutate(
            age = age,
            gender = gender,
            sex = sex,
            race
          )
      }
    )
  ) |>
  tidyr::unnest(cols = meta) ->
  meta_un

meta_un$disease |> unique() |> sort()
meta_un$gender |> unique()
meta_un$age |> unique()
meta_un$sex |> unique()
meta_un$race |> unique()

meta_un |>
  dplyr::mutate(
    dis = ifelse(
      disease %in% c("", "CONTROL", "healthy", "Healthy", "Healthy control", "NEGATIVE"),
      "Health",
      "Covid-19"
    )
  ) |>
  dplyr::select(gse, srrid, dis, age, gender) ->
  meta_un_meta


meta_un_meta |>
  dplyr::group_by(gse, dis) |>
  dplyr::count() |>
  dplyr::ungroup() |>
  tidyr::spread(key = dis, value = n) |>
  dplyr::slice(match(gseidlist, gse)) |>
  dplyr::mutate(
    a = glue::glue("{`Covid-19`},{Health}")
  ) |>
  dplyr::select(a)

meta_un_meta |>
  dplyr::mutate(has_age = ifelse(is.na(age), "no", "yes")) |>
  dplyr::select(gse, has_age) |>
  dplyr::distinct() |>
  dplyr::slice(match(gseidlist, gse))

meta_un_meta |>
  dplyr::mutate(has_age = ifelse(is.na(gender), "no", "yes")) |>
  dplyr::select(gse, has_age) |>
  dplyr::distinct() |>
  dplyr::slice(match(gseidlist, gse))


metadata |>
  dplyr::mutate(
    anno = purrr::map(
      .x = gse,
      .f = \(.gse) {
        file.path(
          basedir, .gse, "out",
          "{.gse}.scmocha.out.rds" |> glue::glue()
        ) |>
          readr::read_rds()
      }
    )
  ) ->
  variants

variants |>
  dplyr::select(gse, anno) |>
  tidyr::unnest(anno) |>
  dplyr::inner_join(
    meta_un_meta |>
      dplyr::select(-gse),
    by = "srrid"
  ) ->
  metadata_anno

metadata_anno |>
  dplyr::filter(!gse %in% c("GSE157344", "GSE171555")) |>
  dplyr::arrange(dis, age, gender) |>
  dplyr::select(srrid, dis, celltype_ratio) |>
  # dplyr::arrange(disease) |>
  dplyr::mutate(color = dplyr::case_match(
    dis,
    "Covid-19" ~ggsci::pal_jama()(4)[[1]],
    "Health" ~ ggsci::pal_jama()(4)[[2]]
  )) |>
  dplyr::arrange(dplyr::desc(dplyr::row_number())) ->
  for_ratio_plot


for_ratio_plot |>
  tidyr::unnest(cols = celltype_ratio) |>
  # dplyr::filter(!grepl("cluster", celltype)) |>
  dplyr::mutate(celltype = factor(
    celltype
  )) |>
  ggplot(aes(
    x = `ratio`,
    y = `srrid`,
    fill = celltype
  )) +
  geom_col() +
  # scale_fill_manual(
  #   # values = paletteer::paletteer_d(
  #   #   palette = "ggsci::springfield_simpsons",
  #   #   direction = -1
  #   # ),
  #   # limits = paste("cluster", c(1:11, 13), sep = "_"),
  #   values = ggsci::scale_color_aaas()[8],
  #   name = "Cell type"
  # ) +
  ggsci::scale_fill_nejm(name = "Cell type") +
  scale_x_continuous(
    expand = expansion(mult = 0, add = 0),
    labels = scales::percent_format()
  ) +
  scale_y_discrete(
    limits = for_ratio_plot$`srrid`
  ) +
  theme(
    panel.background = element_blank(),
    panel.grid = element_blank(),
    axis.text = element_text(color = "black", size = 12, face = "bold"),
    axis.text.y = element_text(
      color = for_ratio_plot$color
    ),
    axis.ticks.y = element_blank(),
    axis.line.x = element_line(color = "black", size = 0.5),
    axis.title = element_text(color = "black", size = 14, face = "bold"),
    axis.title.y = element_blank(),
    legend.position = "right"
  ) +
  labs(x = "Cell ratio") ->
  p_cellratio;p_cellratio


ggsave(
  filename = "Cell_ratio.pdf",
  plot = p_cellratio,
  device = "pdf",
  width = 18,
  height = 10,
  path = "/home/liuc9/github/scMOCHA-data/data/out"
)


fn_plot_gene <- function() {
  mt_exons_df <- "/home/liuc9/github/scMOCHA/fasta/mt_exons.df.rds.gz"


  gtf_gene_df <-
    readr::read_rds(
      file = mt_exons_df
    )
  library(gggenes)
  ggplot(gtf_gene_df, aes(xmin = start, xmax = end, y = seqnames)) +
    # geom_gene_arrow() +
    geom_gene_arrow(
      aes(
        fill = gene_biotype
      ),
      arrowhead_height = unit(3, "mm"), arrowhead_width = unit(1, "mm")
    ) +
    scale_fill_brewer(
      palette = "Set1",
      name = "Gene type",
      labels = c("MT rRNA", "MT tRNA", "Protein coding")
    ) +
    ggrepel::geom_text_repel(
      aes(x = (start + end) / 2, label = gene_name, color = gene_biotype),
      # fill = "white",
      # nudge_x =1,
      # nudge_y = -0.1,
      size = 3,
      show.legend = F,
      max.overlaps = Inf,
    ) +
    scale_color_brewer(palette = "Set1") +
    scale_x_continuous(
      limits = c(0, 17000),
      breaks = seq(0, 17000, 2000),
      expand = expansion(mult = c(0, 0.03)),
    ) +
    scale_y_discrete(
      expand = expansion(mult = c(0, 0), add = c(0, 0))
    ) +
    theme_genes() +
    theme(
      legend.position = "bottom",
      axis.title = element_blank(),
      axis.text.y = element_blank(),
      axis.text.x = element_text(size = 14),
      legend.text = element_text(size = 14)
    ) ->
    pg;pg
}


metadata_anno |>
  dplyr::filter(!gse %in% c("GSE157344", "GSE171555")) |>
  dplyr::arrange(dis, age, gender) |>
  dplyr::select(srrid, dis, depth) |>
  # dplyr::arrange(disease) |>
  dplyr::mutate(color = dplyr::case_match(
    dis,
    "Covid-19" ~ggsci::pal_jama()(4)[[1]],
    "Health" ~ ggsci::pal_jama()(4)[[2]]
  )) |>
  dplyr::arrange(dplyr::desc(dplyr::row_number())) ->
  for_depth_plot


for_depth_plot |>
  tidyr::unnest(cols = depth) |>
  ggplot(aes(x=pos, y = depth, fill = `color`)) +
  geom_bar(stat = "identity") +
  scale_x_continuous(
    expand = expansion(mult = c(0.01, 0)),
    limits = c(1, 17000),
    breaks = seq(0, 17000, 2000),
    labels = seq(0, 17000, 2000)
  ) +
  scale_y_continuous(
    expand = c(0.01, 0),
    # limits = c(0, 520000),
    label = scales::label_number()
  ) +
  scale_fill_identity(
    name = "Sample"
  ) +
  theme(
    plot.margin = margin(t = 0, b = 0, unit = "cm"),
    panel.background = element_blank(),
    panel.grid = element_blank(),
    axis.line.y.left = element_line(color = "black"),
    # axis.line.x.bottom = element_line(color = "black"),
    axis.ticks.x = element_blank(),
    axis.text.x = element_blank(),
    axis.line.x = element_blank(),
    axis.title.x = element_blank(),
    legend.position = c(0.8, 0.5),
    legend.key = element_blank(),
    axis.title.y = element_text(size = 16, color = "black"),
    axis.text.y = element_text( color = "black"),
    legend.text = element_text(
      size = 14,
      color = "black"
    ),
    legend.title = element_text(
      size = 16,
      colour = "black"
    ),
    strip.background = element_blank(),
    strip.text = element_text(
      size = 8,
      color = "black",
      face = "bold"
    )
  ) +
  facet_wrap(
    facets = ~dis,
    ncol = 1,
    strip.position = "right"
  ) +
  labs(y = "Depth") ->
  p_mt_depth;p_mt_depth


wrap_plots(
  p_mt_depth,
  fn_plot_gene(),
  ncol = 1,
  heights = c(0.9, 0.1)
) ->
  p_depth;p_depth

ggsave(
  filename = "Sample_depth_merge.pdf",
  plo = p_depth,
  device = "pdf",
  width = 15,
  height = 8,
  path = "/home/liuc9/github/scMOCHA-data/data/out"
)

# footer ------------------------------------------------------------------

# future::plan(future::sequential)

# save image --------------------------------------------------------------
