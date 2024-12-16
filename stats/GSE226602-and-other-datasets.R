#!/usr/bin/env Rscript --vanilla
# Metainfo ----------------------------------------------------------------
# @AUTHOR: Chun-Jie Liu
# @CONTACT: chunjie.sam.liu.at.gmail.com
# @DATE: `r date()`
# @DESCRIPTION: filename

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
basedir <- "/home/liuc9/github/scMOCHA-data/data"

outdir <- file.path(basedir, "out")

gseids <- c("GSE149689", "GSE155223", "GSE155673", "GSE157344", "GSE163668", "GSE166992", "GSE171555", "GSE181279", "GSE226602")

gseids_meta <- tibble::tibble(
  GSE_ID = c("GSE163668", "GSE149689", "GSE155223", "GSE155673", "GSE157344", "GSE166992", "GSE171555", "GSE226602", "GSE181279"),
  samples = c(38, 20, 18, 12, 33, 9, 48, 50, 5),
  Disease = c("COVID-19", "COVID-19", "COVID-19", "COVID-19", "COVID-19", "COVID-19", "COVID-19", "AD", "AD"),
  Source = c("PBMC", "PBMC", "PBMC", "PBMC", "PBMC", "PBMC", "PBMC", "PBMC", "PBMC"),
  Chemistry = c("SC5P-R2", "SC3Pv3", "SC5P-R2", "SC3Pv3", "SC3Pv3", "SC5P-PE", "SC5P-R2", "SC5P-PE", "SC5P-PE"),
  Publication = c("Nature, 2021", "Exp Mol Med, 2022", "Cell Rep, 2023", "Science, 2020", "Nat Commun, 2021", "Cell Rep, 2021", "Med, 2021", "Neuron, 2024", "Front Immunol., 2021")
)

# body --------------------------------------------------------------------

# load gse cell ratio and variant data

tibble::tibble(
  gseid = gseids
) |>
  dplyr::mutate(
    cell_ratio_variant = purrr::map(
      .x = gseid,
      .f = \(.gseid) {
        data.table::fread(
          file.path(basedir, .gseid, "out", glue::glue("{.gseid}.cell_ratio_and_variant_clean.csv"))
        )
      }
    )
  ) |>
  dplyr::mutate(
    anno = purrr::map(
      .x = gseid,
      .f = \(.gseid) {
        readr::read_rds(
          file.path(basedir, .gseid, "out", glue::glue("{.gseid}.scmocha.out.rds.gz"))
        )
      }
    )
  ) |>
  dplyr::left_join(
    gseids_meta,
    by = c("gseid" = "GSE_ID")
  ) ->
gse_cell_ratio_variant_meta

# save gse cell ratio and variant data ------------------------------------
gse_cell_ratio_variant_meta |>
  dplyr::mutate(
    `Avg # of somatic variants` = purrr::map_dbl(
      cell_ratio_variant,
      ~ mean(.x$`# of somatic variants`)
    )
  ) ->
gse_cell_ratio_variant_meta_xlsx

gse_cell_ratio_variant_meta_xlsx |>
  dplyr::select(-cell_ratio_variant, -anno) |>
  dplyr::arrange(`Avg # of somatic variants`) |>
  writexl::write_xlsx(
    path = file.path(outdir, "gses_cell_ratio_variant_meta.xlsx")
  )

gse_cell_ratio_variant_meta_xlsx |>
  dplyr::arrange(dplyr::desc(`Avg # of somatic variants`)) |>
  dplyr::mutate(
    label = glue::glue("{gseid} ({round(`Avg # of somatic variants`, 2)})")
  ) |>
  dplyr::mutate(
    label = factor(label, levels = label)
  ) ->
gseid_ranked

# plot average depth and number of somatic variants ------------------------
gse_cell_ratio_variant_meta |>
  dplyr::left_join(
    gseid_ranked |> dplyr::select(gseid, label),
    by = "gseid"
  ) |>
  tidyr::unnest(cols = cell_ratio_variant) |>
  dplyr::mutate(
    Chemistry = factor(Chemistry, levels = c("SC3Pv3", "SC5P-R2", "SC5P-PE") |> rev()),
  ) ->
forplot

cor.test(~ `Depth mean` + `# of somatic variants`, data = forplot) |> broom::tidy() -> cor_test_all
cor.test(~ `Depth mean` + `# of somatic variants`, data = forplot, subset = gseid != "GSE181279") |> broom::tidy() -> cor_test_250k

forplot |>
  ggplot(aes(
    x = `Depth mean`,
    y = `# of somatic variants`,
  )) +
  geom_point(aes(
    shape = Chemistry,
    color = label
  )) +
  geom_smooth(
    method = "loess", se = FALSE, color = "black",
    linetype = 21,
  ) +
  geom_hline(
    yintercept = 10,
    linetype = 21,
    color = "red"
  ) +
  ggsci::scale_color_aaas(
    name = "GSE ID",
  ) +
  scale_x_continuous(
    labels = scales::label_number(),
    limits = c(0, 60000),
    breaks = seq(0, 60000, 10000),
  ) +
  scale_y_continuous(
    labels = scales::label_number(),
    limits = c(0, 80),
    breaks = seq(0, 80, 10),
  ) +
  theme(
    panel.background = element_blank(),
    panel.grid = element_blank(),
    axis.line = element_line(size = 0.5, color = "black"),
    axis.title = element_text(
      size = 16,
      color = "black",
      face = "bold"
    ),
    axis.text.y = element_text(
      size = 14,
      color = "black"
    ),
    plot.title = element_text(
      hjust = 0.5,
      color = "black",
      size = 16,
      face = "bold"
    )
  ) +
  labs(
    x = "Average Depth",
    y = "Number of Somatic Variants",
    title = glue::glue("All samples test, Pearson's r = {round(cor_test_all$estimate, 2)}, p-value = {scales::pvalue(cor_test_all$p.value)}"),
    subtitle = glue::glue("Exclude GSE181279 samples, Pearson's r = {round(cor_test_250k$estimate, 2)}, p-value = {scales::pvalue(cor_test_250k$p.value)}")
  ) ->
p_somatic_variant

ggsave(
  filename = file.path(outdir, "somatic_variant.pdf"),
  plot = p_somatic_variant,
  width = 10,
  height = 6,
  dpi = 300
)

# SC5P-PE samples ---------------------------------------------------------
gse_cell_ratio_variant_meta |>
  dplyr::filter(Chemistry == "SC5P-PE") |>
  dplyr::select(-cell_ratio_variant) ->
sc5ppe_anno

sc5ppe_anno |>
  dplyr::mutate(
    somatic_variants = purrr::map(
      .x = anno,
      .f = \(.anno) {
        .anno |>
          dplyr::select(srrid, somatic_variant) |>
          dplyr::mutate(
            somatic_variant = purrr::map(
              .x = somatic_variant,
              .f = \(.x) {
                tibble::tibble(
                  variant = .x$somatic
                )
              }
            )
          ) ->
        .anno_somatic

        .anno_somatic$somatic_variant |>
          purrr::map(~ .x$variant) |>
          purrr::reduce(union) ->
        union_variants

        tibble::tibble(
          anno_somatic = list(.anno_somatic),
          union_variants = list(union_variants)
        )
      }
    )
  ) |>
  tidyr::unnest(cols = somatic_variants) ->
sc5ppe_anno_somatic

sc5ppe_anno_somatic

# venn diagram dataset ------------------------------------------------------------


ggvenn::ggvenn(
  data = list(
    "GSE166992 (n=9)" = sc5ppe_anno_somatic$union_variants[[1]],
    "GSE181279 (n=50)" = sc5ppe_anno_somatic$union_variants[[2]],
    "GSE226602 (n=5)" = sc5ppe_anno_somatic$union_variants[[3]]
  ),
  fill_color = ggsci::pal_npg()(3),
) ->
p_venn_sc5ppe

ggsave(
  filename = file.path(outdir, "venn_gse_sc5ppe.pdf"),
  plot = p_venn_sc5ppe,
  width = 7,
  height = 5,
  dpi = 300
)

# venn diagram individual -------------------------------------------------------
sc5ppe_anno_somatic |>
  tidyr::unnest(cols = anno_somatic) |>
  dplyr::select(srrid, somatic_variant) |>
  dplyr::mutate(somatic_variant = purrr::map(somatic_variant, ~ .x$variant)) |>
  dplyr::pull(somatic_variant) |>
  purrr::reduce(intersect)


sc5ppe_anno_somatic |>
  dplyr::select(gseid, anno) |>
  tidyr::unnest(cols = anno) |>
  dplyr::select(gseid, srrid, coverage, haplo_variant, haplo_violin, somatic_variant) ->
sc5ppe_anno_somatic_detail


thevariant <- "3173G>A"
# thevariant <- "1670A>G"

sc5ppe_anno_somatic_detail |>
  dplyr::select(gseid, srrid, haplo_violin) |>
  tidyr::unnest(cols = haplo_violin) |>
  dplyr::mutate(
    gseid = factor(gseid, levels = c("GSE166992", "GSE226602", "GSE181279"))
  ) |>
  dplyr::filter(variant == thevariant) ->
sel_variant

pcc <- readr::read_tsv(file = "https://raw.githubusercontent.com/chunjie-sam-liu/chunjie-sam-liu.life/master/public/data/pcc.tsv") |>
  dplyr::arrange(cancer_types)


sel_variant |>
  dplyr::filter(cluster == "B") |>
  dplyr::group_by(gseid, srrid) |>
  dplyr::summarise(
    mean_af = mean(af, na.rm = TRUE)
  ) |>
  dplyr::ungroup() |>
  dplyr::arrange(gseid, mean_af) ->
sel_variant_ranked

sel_variant |>
  dplyr::mutate(
    srrid = factor(srrid, levels = sel_variant_ranked$srrid)
  ) |>
  # dplyr::filter(cluster == "B") |>
  ggplot() +
  ggh4x::facet_wrap2(
    ~cluster,
    ncol = 1,
    strip.position = "right",
    strip = ggh4x::strip_themed(
      background_y = ggh4x::elem_list_rect(
        fill = pcc$color
      ),
      text_y = ggh4x::elem_list_text(
        colour = "white",
        face = c("bold")
      ),
      by_layer_y = FALSE,
    )
  ) +
  ggbeeswarm::geom_quasirandom(
    aes(
      x = srrid,
      y = af,
      color = af
    ),
    size = 1,
    dodge.width = .75,
    alpha = .5,
    varwidth = TRUE
  ) +
  scale_color_gradient2(
    name = "AF",
    low = "white",
    mid = "red",
    high = "#3B0049",
    midpoint = 0.5,
  ) +
  scale_y_continuous(
    expand = c(0.01, 0),
    limits = c(0, 1),
  ) +
  theme(
    plot.margin = margin(t = 0, b = 0, unit = "cm"),
    panel.background = element_blank(),
    panel.grid = element_blank(),
    axis.line.y.left = element_line(color = "black"),
    axis.ticks.x = element_blank(),
    axis.text.x = element_blank(),
    axis.line.x = element_blank(),
    axis.title.x = element_blank(),
    legend.position = "right",
    legend.key = element_blank(),
    axis.title.y = element_text(color = "black"),
    axis.text.y = element_text(color = "black"),
  ) +
  labs(y = "AF") ->
p_af

sel_variant |>
  dplyr::mutate(
    srrid = factor(srrid, levels = sel_variant_ranked$srrid)
  ) |>
  dplyr::filter(cluster == "B") |>
  dplyr::group_by(gseid) |>
  dplyr::mutate(
    mid_srrid = srrid[ceiling(dplyr::n() / 2)]
  ) |>
  ggplot(aes(
    x = srrid,
    y = 1
  )) +
  geom_tile(
    aes(
      fill = gseid
    )
  ) +
  geom_text(
    aes(
      y = 1,
      label = ifelse(srrid == mid_srrid, as.character(gseid), "")
    ),
  ) +
  ggsci::scale_fill_jco() +
  theme(
    panel.background = element_blank(),
    panel.grid = element_blank(),
    axis.line = element_blank(),
    axis.title = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(), legend.position = "none"
  ) ->
p_tile

ggsave(
  filename = file.path(outdir, "haplo_violin_{thevariant}.pdf" |> glue::glue()),
  plot =
    wrap_plots(
      p_af,
      p_tile,
      ncol = 1,
      heights = c(30, 1)
    ),
  width = 24,
  height = 12,
  dpi = 300
)
# cell type ratio ---------------------------------------------------------

sc5ppe_anno_somatic |>
  dplyr::select(gseid, anno) |>
  tidyr::unnest(cols = anno) |>
  dplyr::select(gseid, srrid, celltype_ratio) |>
  tidyr::unnest(cols = celltype_ratio) ->
sc5ppe_anno_celltype



sc5ppe_anno_celltype |>
  dplyr::mutate(
    srrid = factor(srrid, levels = sel_variant_ranked$srrid)
  ) |>
  ggplot(aes(
    x = srrid,
    y = n,
  )) +
  geom_col(
    aes(
      fill = celltype
    ),
    position = "stack"
  ) +
  scale_fill_manual(
    name = "Cell Type",
    values = pcc$color
  ) +
  scale_y_continuous(
    expand = expansion(add = c(0.005, 0.005)),
  ) +
  theme(
    plot.margin = margin(t = 0, b = 0, unit = "cm"),
    panel.background = element_blank(),
    panel.grid = element_blank(),
    axis.line.y.left = element_line(color = "black"),
    axis.ticks.x = element_blank(),
    axis.text.x = element_blank(),
    axis.line.x = element_blank(),
    axis.title.x = element_blank(),
    legend.position = "right",
    legend.key = element_blank(),
    axis.title.y = element_text(color = "black"),
    axis.text.y = element_text(color = "black"),
  ) +
  labs(y = "# of cells") ->
p_celltype_count

sc5ppe_anno_celltype |>
  dplyr::mutate(
    srrid = factor(srrid, levels = sel_variant_ranked$srrid)
  ) |>
  ggplot(aes(
    x = srrid,
    y = ratio,
  )) +
  geom_col(
    aes(
      fill = celltype
    ),
    position = "stack"
  ) +
  scale_fill_manual(
    name = "Cell Type",
    values = pcc$color
  ) +
  scale_y_continuous(
    expand = expansion(add = c(0.005, 0.005)),
  ) +
  theme(
    plot.margin = margin(t = 0, b = 0, unit = "cm"),
    panel.background = element_blank(),
    panel.grid = element_blank(),
    axis.line.y.left = element_line(color = "black"),
    axis.ticks.x = element_blank(),
    axis.text.x = element_blank(),
    axis.line.x = element_blank(),
    axis.title.x = element_blank(),
    legend.position = "right",
    legend.key = element_blank(),
    axis.title.y = element_text(color = "black"),
    axis.text.y = element_text(color = "black"),
  ) +
  labs(y = "cell ratio") ->
p_celltype_ratio


sc5ppe_anno_somatic |>
  dplyr::select(gseid, anno) |>
  tidyr::unnest(cols = anno) |>
  dplyr::mutate(
    mean_depth = purrr::map_dbl(
      .x = depth,
      .f = \(.depth) {
        mean(.depth$depth, na.rm = TRUE)
      }
    )
  ) |>
  dplyr::select(srrid, mean_depth) |>
  dplyr::mutate(
    srrid = factor(srrid, levels = sel_variant_ranked$srrid)
  ) |>
  ggplot(aes(
    x = srrid,
    y = mean_depth
  )) +
  geom_col(
    fill = "grey"
  ) +
  scale_y_continuous(
    expand = expansion(add = c(0.005, 0.005)),
  ) +
  theme(
    plot.margin = margin(t = 0, b = 0, unit = "cm"),
    panel.background = element_blank(),
    panel.grid = element_blank(),
    axis.line.y.left = element_line(color = "black"),
    axis.ticks.x = element_blank(),
    axis.text.x = element_blank(),
    axis.line.x = element_blank(),
    axis.title.x = element_blank(),
    legend.position = "right",
    legend.key = element_blank(),
    axis.title.y = element_text(color = "black"),
    axis.text.y = element_text(color = "black"),
  ) +
  labs(y = "Average Depth") ->
p_depth_srrid

ggsave(
  filename = file.path(outdir, "celltype_ratio.pdf" |> glue::glue()),
  plot = wrap_plots(
    p_depth_srrid,
    p_celltype_ratio,
    p_celltype_count,
    p_tile,
    ncol = 1,
    heights = c(15, 15, 15, 1),
    guides = "collect"
  ),
  width = 24,
  height = 13,
  dpi = 300
)
# footer ------------------------------------------------------------------

# future: :plan(future: :sequential)

# save image --------------------------------------------------------------
