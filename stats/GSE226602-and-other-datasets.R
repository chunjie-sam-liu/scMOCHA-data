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
fn_plot_mtdna <- function() {
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
      breaks = seq(0, 17000, 1000),
      labels = seq(0, 17000, 1000),
      expand = expansion(mult = c(0, 0.01)),
    ) +
    scale_y_discrete(
      expand = expansion(mult = c(0, 0), add = c(0, 0))
    ) +
    # theme_genes() +
    theme(
      legend.position = "bottom",
      axis.title = element_blank(),
      axis.text.y = element_blank(),
      # axis.text.x = element_text(size = 14),
      # legend.text = element_text(size = 14),
      panel.background = element_blank(),
      panel.grid = element_blank(),
      axis.ticks.y = element_blank(),
      axis.line.x = element_line(color = "black"),
      axis.text.x = element_text(
        vjust = -1,
      ),
    ) +
    coord_cartesian(xlim = c(0, 17000)) ->
  pg
  pg
}


# load data ---------------------------------------------------------------
basedir <- "/home/liuc9/github/scMOCHA-data/data"

outdir <- file.path(basedir, "out_new")
dir.create(outdir, showWarnings = FALSE)

gseids <- c(
  "GSE155673",
  "GSE157344",
  "GSE149689",
  "GSE171555",
  "GSE155223",
  "GSE163668",
  "GSE226602",
  "GSE166992",
  "GSE181279",
  "WT",
  # plus ting,
  "GSE161354", # done
  "GSE175524", # done
  "GSE206283", # done
  "GSE226598", # done
  "GSE235050", # done
  "GSE261140", # some errors
  "GSE279945"
)

gseids_meta_raw <- tibble::tibble(
  GSE_ID = c("GSE163668", "GSE149689", "GSE155223", "GSE155673", "GSE157344", "GSE166992", "GSE171555", "GSE226602", "GSE181279", "WT"),
  samples = c(38, 20, 18, 12, 33, 9, 48, 50, 5, 1),
  Disease = c("COVID-19", "COVID-19", "COVID-19", "COVID-19", "COVID-19", "COVID-19", "COVID-19", "AD", "AD", "-"),
  Source = c("PBMC", "PBMC", "PBMC", "PBMC", "PBMC", "PBMC", "PBMC", "PBMC", "PBMC", "Mixed"),
  Chemistry = c("SC5P-R2", "SC3Pv3", "SC5P-R2", "SC3Pv3", "SC3Pv3", "SC5P-PE", "SC5P-R2", "SC5P-PE", "SC5P-PE", "SC3Pv3"),
  Publication = c("Nature, 2021", "Exp Mol Med, 2022", "Cell Rep, 2023", "Science, 2020", "Nat Commun, 2021", "Cell Rep, 2021", "Med, 2021", "Neuron, 2024", "Front Immunol., 2021", "-")
)
gseids_meta_ting <- tibble::tibble(
  GSE_ID = c(
    "GSE161354", # done
    "GSE175524", # done
    "GSE206283", # done
    "GSE226598", # done
    "GSE235050", # done
    "GSE261140", # some errors
    "GSE279945"
  ),
) |>
  dplyr::mutate(
    samples = purrr::map(
      GSE_ID,
      .f = \(.x) {
        data.table::fread(
          file.path(basedir, .x, "out", glue::glue("{.x}.cell_ratio_and_variant_clean.csv"))
        ) ->
        .d
        tibble::tibble(
          samples = nrow(.d),
          Disease = "-",
          Source = "PBMC",
          Chemistry = unique(.d$Chemistry)[[1]],
          Publication = "-"
        )
      }
    ),
  ) |>
  tidyr::unnest(cols = samples)

gseids_meta <- dplyr::bind_rows(
  gseids_meta_raw,
  gseids_meta_ting
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
        ) |>
          dplyr::select(-Chemistry)
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
  ) ->
gse_data_loaded


gse_data_loaded |>
  dplyr::select(-anno) |>
  tidyr::unnest(cols = cell_ratio_variant) |>
  dplyr::group_by(gseid) |>
  dplyr::summarise(
    `Avg. mutation` = mean(`# of variants`, na.rm = TRUE),
    `Avg. somatic mutation` = mean(`# of somatic variants`, na.rm = TRUE),
    `Avg. total reads` = mean(`Total reads`, na.rm = TRUE),
    `Avg. mapped reads` = mean(`Depth read mean`, na.rm = TRUE),
    `Avg. call depth` = mean(`Depth mean`, na.rm = TRUE),
  ) |>
  dplyr::left_join(
    # gseids_meta,
    gseids_meta,
    by = c("gseid" = "GSE_ID")
  ) |>
  dplyr::slice(match(gseids, gseid)) |>
  writexl::write_xlsx(
    path = file.path(outdir, "gses_meta_read.xlsx")
  )

gse_data_loaded |>
  dplyr::left_join(
    gseids_meta,
    by = c("gseid" = "GSE_ID")
  ) ->
gse_cell_ratio_variant_meta

# save gse cell ratio and variant data ------------------------------------
gse_cell_ratio_variant_meta |>
  # dplyr::mutate(
  #   `Avg # of variants` = purrr::map_dbl(
  #     cell_ratio_variant,
  #     ~ mean(.x$`# of variants`)
  #   )
  # ) |>
  # dplyr::mutate(
  #   `Avg # of somatic variants` = purrr::map_dbl(
  #     cell_ratio_variant,
  #     ~ mean(.x$`# of somatic variants`)
  #   )
  # ) ->
  tidyr::unnest(cols = cell_ratio_variant) ->
gse_cell_ratio_variant_meta_xlsx

gse_cell_ratio_variant_meta_xlsx |>
  dplyr::select(-anno) |>
  dplyr::arrange(`# of somatic variants`) |>
  writexl::write_xlsx(
    path = file.path(outdir, "gses_cell_ratio_variant_meta.xlsx")
  )

gse_cell_ratio_variant_meta_xlsx |>
  dplyr::group_by(gseid, Chemistry) |>
  dplyr::summarise(
    `Avg # of somatic variants` = mean(`# of somatic variants`, na.rm = TRUE),
    `# of samples` = dplyr::n()
  ) |>
  dplyr::ungroup() |>
  dplyr::arrange(dplyr::desc(`Avg # of somatic variants`)) |>
  dplyr::mutate(
    label = glue::glue("{gseid} ({Chemistry}, {`# of samples`}, {round(`Avg # of somatic variants`, 2)})")
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

pcc <- readr::read_tsv(file = "https://raw.githubusercontent.com/chunjie-sam-liu/chunjie-sam-liu.life/master/public/data/pcc.tsv") |>
  dplyr::arrange(cancer_types)

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
  # ggsci::scale_color_aaas(
  #   name = "GSE ID",
  # ) +
  scale_color_manual(
    name = "GSE ID",
    values = pcc$color
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
  width = 13,
  height = 7,
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


# ggvenn::ggvenn(
#   data = list(
#     "GSE166992 (n=9)" = sc5ppe_anno_somatic$union_variants[[1]],
#     "GSE181279 (n=50)" = sc5ppe_anno_somatic$union_variants[[2]],
#     "GSE226602 (n=5)" = sc5ppe_anno_somatic$union_variants[[3]],
#     "GSE161354 (n=8)" = sc5ppe_anno_somatic$union_variants[[4]],
#     "GSE235050 (n=12)" = sc5ppe_anno_somatic$union_variants[[5]]
#   ),
#   fill_color = ggsci::pal_npg()(5),
# ) ->
# p_venn_sc5ppe

# ggsave(
#   filename = file.path(outdir, "venn_gse_sc5ppe.pdf"),
#   plot = p_venn_sc5ppe,
#   width = 7,
#   height = 5,
#   dpi = 300
# )
library(ggVennDiagram)
variant_list <- list(
  "GSE166992 (n=9)" = sc5ppe_anno_somatic$union_variants[[1]],
  "GSE181279 (n=50)" = sc5ppe_anno_somatic$union_variants[[2]],
  "GSE226602 (n=5)" = sc5ppe_anno_somatic$union_variants[[3]],
  "GSE161354 (n=8)" = sc5ppe_anno_somatic$union_variants[[4]],
  "GSE235050 (n=12)" = sc5ppe_anno_somatic$union_variants[[5]]
)
variant_list_df <- variant_list |>
  ggVennDiagram::Venn() |>
  ggVennDiagram::process_data()

ggplot() +
  # geom_polygon(
  #   aes(X, Y, fill = name, group = id),
  #   data = ggVennDiagram::venn_regionedge(variant_list_df),
  #   show.legend = FALSE
  # ) +
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
    gseid = factor(gseid, levels = c("GSE166992", "GSE226602", "GSE181279", "GSE161354", "GSE235050"))
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

# depth 2173G>A -----------------------------------------------------------
thevariants <- c("1670A>G", "3173G>A")

theposes <- thevariants |>
  purrr::map(~ gsub(pattern = "[>|AGCT]", "", x = .)) |>
  purrr::map_int(as.integer)

gse_cell_ratio_variant_meta |>
  dplyr::select(gseid, anno, Chemistry) |>
  tidyr::unnest(cols = anno) |>
  dplyr::select(gseid, srrid, depth, Chemistry) |>
  tidyr::unnest(cols = depth) ->
all_gseid_depth

all_gseid_depth |>
  dplyr::filter(pos %in% theposes) ->
theposes_depth

theposes_depth |>
  dplyr::filter(pos == theposes[2]) |>
  dplyr::group_by(gseid) |>
  dplyr::summarise(
    mean_depth = mean(depth, na.rm = TRUE)
  ) |>
  dplyr::arrange(dplyr::desc(mean_depth)) ->
theposes_depth_ranked

theposes_depth |>
  dplyr::filter(pos == theposes[2]) |>
  dplyr::mutate(
    gseid = factor(gseid, levels = theposes_depth_ranked$gseid),
    Chemistry = factor(Chemistry, levels = c("SC3Pv3", "SC5P-R2", "SC5P-PE") |> rev())
  ) ->
theposes_depth_forplot

theposes_depth_forplot |>
  ggplot() +
  geom_violin(
    aes(
      x = gseid,
      y = depth,
      fill = Chemistry
    ),
    alpha = 0.5,
    # size = 1,
    color = NA,
    show.legend = FALSE
  ) +
  scale_fill_manual(
    name = "Chemistry",
    values = c(viridis::viridis_pal(option = "D")(3), "red")
  ) +
  ggbeeswarm::geom_quasirandom(
    aes(
      x = gseid,
      y = depth,
      color = Chemistry
    ),
    size = 1,
    dodge.width = .75,
    alpha = .5,
    varwidth = TRUE
  ) +
  scale_color_manual(
    name = "Chemistry",
    values = c(viridis::viridis_pal(option = "D")(3), "red")
  ) +
  scale_x_discrete(
    # expand = expansion(mult = c(0.001, 0.001))
  ) +
  scale_y_continuous(
    labels = scales::label_number(),
    limits = c(0, 40000),
    breaks = seq(0, 40000, 10000),
    # expand = expansion(add = c(0.1, 0.1))
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
    ),
    axis.title.x = element_blank(),
    legend.position = "inside",
    legend.position.inside = c(0.8, 0.7)
  ) +
  labs(
    x = "GSE ID",
    y = "Depth",
    title = glue::glue("Depth of {theposes[2]}"),
  ) ->
p_depth_3173

ggsave(
  filename = file.path(outdir, "depth_3173.pdf"),
  plot = p_depth_3173,
  width = 13,
  height = 5,
  dpi = 300
)

{
  # all position in sc5ppe -------------------------------------------------

  sc5ppe_anno_somatic$union_variants |>
    purrr::reduce(union) |>
    sort() ->
  all_sc5ppe_variants

  all_sc5ppe_poss <- all_sc5ppe_variants |>
    purrr::map(~ gsub(pattern = "[>|AGCT]", "", x = .)) |>
    purrr::map_int(as.integer) |>
    sort()

  all_gseid_depth |>
    dplyr::filter(pos %in% all_sc5ppe_poss) |>
    dplyr::mutate(
      Chemistry = factor(Chemistry, levels = c("SC3Pv3 Mixed", "SC3Pv3", "SC5P-R2", "SC5P-PE") |> rev())
    ) ->
  all_sc5ppe_depth

  all_sc5ppe_depth |>
    dplyr::group_by(Chemistry, srrid) |>
    dplyr::summarise(
      depth = mean(depth, na.rm = TRUE)
    ) |>
    dplyr::arrange(Chemistry, -depth) |>
    dplyr::ungroup() ->
  all_sc5ppe_depth_ranked


  all_sc5ppe_depth_ranked |>
    dplyr::mutate(
      srrid = factor(srrid, levels = all_sc5ppe_depth_ranked$srrid),
    ) |>
    ggplot(aes(
      x = 1,
      y = srrid
    )) +
    geom_tile(
      aes(
        fill = Chemistry
      )
    ) +
    scale_fill_manual(
      name = "Chemistry",
      values = c(viridis::viridis_pal(option = "D")(3), "red")
    ) +
    theme(
      panel.background = element_blank(),
      panel.grid = element_blank(),
      axis.line = element_blank(),
      axis.title = element_blank(),
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      legend.position = "none",
      plot.margin = margin(0, 0, 0, 0, unit = "cm")
    ) ->
  p_tile_all_srrid

  all_sc5ppe_depth |>
    dplyr::mutate(
      srrid = factor(srrid, levels = all_sc5ppe_depth_ranked$srrid),
      depth_log2 = log2(depth + 1)
    ) |>
    dplyr::mutate(
      posc = as.character(pos)
    ) |>
    ggplot(aes(
      x = posc,
      y = srrid,
      fill = depth_log2
    )) +
    geom_tile() +
    scale_fill_gradient2(
      name = "log2(Depth + 1)",
      low = "white",
      mid = "red",
      high = "#3B0049",
      midpoint = 0.5
    ) +
    theme(
      panel.background = element_blank(),
      panel.grid = element_blank(),
      axis.line = element_blank(),
      axis.title = element_blank(),
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      legend.position = "right",
      legend.key = element_blank(),
      plot.margin = margin(0, 0, 0, 0, unit = "cm")
    ) ->
  p_tile_all_depth

  ggsave(
    filename = file.path(outdir, "depth_all_sc5ppe_position.pdf"),
    plot = wrap_plots(
      p_tile_all_srrid,
      p_tile_all_depth,
      ncol = 2,
      widths = c(1, 40)
    ),
    width = 14,
    height = 5,
    dpi = 300
  )
}

# total metrics -----------------------------------------------------------

gse_cell_ratio_variant_meta |>
  dplyr::select(-cell_ratio_variant) |>
  tidyr::unnest(cols = anno) |>
  dplyr::select(gseid, srrid, metrics, Chemistry) |>
  tidyr::unnest(cols = metrics) |>
  dplyr::select(gseid, srrid, `Number of Reads`, Chemistry) |>
  dplyr::mutate(
    gseid = factor(gseid, levels = theposes_depth_ranked$gseid),
    Chemistry = factor(Chemistry, levels = c("SC3Pv3 Mixed", "SC3Pv3", "SC5P-R2", "SC5P-PE") |> rev())
  ) |>
  ggplot() +
  geom_violin(
    aes(
      x = gseid,
      y = `Number of Reads`,
      fill = Chemistry
    ),
    alpha = 0.5,
    # size = 1,
    color = NA,
    show.legend = FALSE
  ) +
  scale_fill_manual(
    name = "Chemistry",
    values = c(viridis::viridis_pal(option = "D")(3), "red")
  ) +
  ggbeeswarm::geom_quasirandom(
    aes(
      x = gseid,
      y = `Number of Reads`,
      color = Chemistry
    ),
    size = 1,
    dodge.width = .75,
    alpha = .5,
    varwidth = TRUE
  ) +
  scale_color_manual(
    name = "Chemistry",
    values = c(viridis::viridis_pal(option = "D")(3), "red")
  ) +
  scale_x_discrete(
    # expand = expansion(mult = c(0.001, 0.001))
  ) +
  scale_y_log10(
    # labels = scales::label_log(),
    # limits = c(0, 40000),
    # breaks = seq(0, 40000, 10000),
    expand = expansion(add = c(0.1, 0.1))
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
    ),
    axis.title.x = element_blank(),
    legend.position = "inside",
    legend.position.inside = c(0.7, 0.7)
  ) +
  labs(
    x = "GSE ID",
    y = "Total reads",
  ) ->
p_total_reads

ggsave(
  filename = file.path(outdir, "depth_total_reads.pdf"),
  plot = p_total_reads,
  width = 13,
  height = 5,
  dpi = 300
)

gse_cell_ratio_variant_meta |>
  dplyr::select(-cell_ratio_variant) |>
  tidyr::unnest(cols = anno) |>
  dplyr::select(gseid, srrid, metrics, Chemistry) |>
  tidyr::unnest(cols = metrics) |>
  dplyr::group_by(gseid) |>
  dplyr::summarise(avg_totla_reads = mean(`Number of Reads`, na.rm = TRUE)) ->
avg_totla_reads

# all position in three chemistry -------------------------------------------------

all_gseid_depth |>
  dplyr::group_by(
    Chemistry, pos
  ) |>
  dplyr::summarise(
    depth = mean(depth, na.rm = TRUE)
  ) |>
  dplyr::ungroup() |>
  dplyr::mutate(
    Chemistry = factor(Chemistry, levels = c("SC3Pv3", "SC3Pv3", "SC5P-R2", "SC5P-PE") |> rev())
  ) ->
all_gseid_depth_forplot

all_gseid_depth_forplot |>
  ggplot() +
  ggh4x::facet_wrap2(
    ~Chemistry,
    ncol = 1,
    strip.position = "right",
    scales = "free_y",
    strip = ggh4x::strip_themed(
      background_y = ggh4x::elem_list_rect(
        fill = c(viridis::viridis_pal(option = "D")(3), "red")
      ),
      text_y = ggh4x::elem_list_text(
        colour = "white",
        face = c("bold")
      ),
      by_layer_y = FALSE,
    ),
  ) +
  geom_col(
    aes(
      x = pos,
      y = depth,
      fill = Chemistry
    ),
  ) +
  scale_fill_manual(
    name = "Chemistry",
    values = c(viridis::viridis_pal(option = "D")(3), "red")
  ) +
  scale_x_continuous(
    limits = c(0, 17000),
    breaks = seq(0, 17000, 1000),
    labels = seq(0, 17000, 1000),
    expand = expansion(mult = c(0, 0.01)),
  ) +
  scale_y_continuous(
    expand = c(0.01, 0),
    label = scales::label_number(),
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
    legend.position = "none",
    legend.key = element_blank(),
    axis.title.y = element_text(color = "black"),
    axis.text.y = element_text(color = "black"),
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
  coord_cartesian(xlim = c(0, 17000)) +
  labs(y = "Depth") ->
p_depth_all

p_mtdna <- fn_plot_mtdna()

ggsave(
  filename = file.path(outdir, "depth_all_position.free_y.pdf"),
  plot = wrap_plots(
    p_depth_all,
    p_mtdna,
    ncol = 1,
    heights = c(15, 1)
  ),
  width = 17,
  height = 9,
  dpi = 300
)

all_gseid_depth |>
  dplyr::group_by(gseid, srrid) |>
  dplyr::summarise(
    mean_depth = mean(depth, na.rm = TRUE)
  ) |>
  dplyr::summarise(
    avg_mito_reads = mean(mean_depth, na.rm = TRUE)
  ) ->
avg_mito_reads

# all position in three chemistry separated dataset -------------------------------------------------

all_gseid_depth |>
  dplyr::group_by(
    Chemistry, gseid, pos
  ) |>
  dplyr::summarise(
    depth = mean(depth, na.rm = TRUE)
  ) |>
  dplyr::ungroup() |>
  dplyr::mutate(
    gseid = factor(gseid, levels = theposes_depth_ranked$gseid),
    Chemistry = factor(Chemistry, levels = c("SC3Pv3", "SC3Pv3", "SC5P-R2", "SC5P-PE") |> rev())
  ) ->
all_gseid_depth_forplot_separated

all_gseid_depth_forplot_separated |>
  ggplot() +
  ggh4x::facet_wrap2(
    ~gseid,
    ncol = 1,
    strip.position = "right",
    # scales = "free_y",
    strip = ggh4x::strip_themed(
      background_y = ggh4x::elem_list_rect(
        fill = c(rep(viridis::viridis_pal(option = "D")(3), each = 3), "red")
      ),
      text_y = ggh4x::elem_list_text(
        colour = "white",
        face = c("bold")
      ),
      by_layer_y = FALSE,
    ),
  ) +
  geom_col(
    aes(
      x = pos,
      y = depth,
      fill = Chemistry
    ),
  ) +
  scale_fill_manual(
    name = "Chemistry",
    values = c(viridis::viridis_pal(option = "D")(3), "red")
  ) +
  scale_x_continuous(
    limits = c(0, 17000),
    breaks = seq(0, 17000, 1000),
    labels = seq(0, 17000, 1000),
    expand = expansion(mult = c(0, 0.01)),
  ) +
  scale_y_continuous(
    expand = c(0.01, 0),
    label = scales::label_number(),
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
    legend.position = "none",
    legend.key = element_blank(),
    axis.title.y = element_text(color = "black"),
    axis.text.y = element_text(color = "black"),
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
  coord_cartesian(xlim = c(0, 17000)) +
  labs(y = "Depth") ->
p_depth_all_separated

p_mtdna <- fn_plot_mtdna()

ggsave(
  filename = file.path(outdir, "depth_all_position.separated.pdf"),
  plot = wrap_plots(
    p_depth_all_separated,
    p_mtdna,
    ncol = 1,
    heights = c(15, 1)
  ),
  width = 17,
  height = 9,
  dpi = 300
)

# add mixed four cell line reads -------------------------------------------------
mixed_four_cellline_coverage <- data.table::fread("/mnt/isilon/u01_project/large-scale/liuc9/raw/out/mixed_four_cellline_cluster.coverage.txt.gz", col.names = c("pos", "celltype", "count")) |>
  dplyr::group_by(pos) |>
  dplyr::summarise(depth = sum(count)) |>
  dplyr::mutate(
    Chemistry = "SC3Pv3 Mixed",
    gseid = "Mixed cellline"
  ) |>
  dplyr::select(Chemistry, gseid, pos, depth)


all_gseid_depth_forplot_separated |>
  # dplyr::bind_rows(mixed_four_cellline_coverage) |>
  dplyr::mutate(
    gseid = factor(gseid, levels = c(theposes_depth_ranked$gseid, "Mixed cellline")),
    Chemistry = factor(Chemistry, levels = c("SC3Pv3 Mixed", "SC3Pv3", "SC5P-R2", "SC5P-PE") |> rev())
  ) ->
all_gseid_depth_forplot_separated_mixed

all_gseid_depth_forplot_separated_mixed |>
  ggplot() +
  ggh4x::facet_wrap2(
    ~gseid,
    ncol = 1,
    strip.position = "right",
    scales = "free_y",
    strip = ggh4x::strip_themed(
      background_y = ggh4x::elem_list_rect(
        fill = c(rep(viridis::viridis_pal(option = "D")(3), each = 3), "red")
      ),
      text_y = ggh4x::elem_list_text(
        colour = "white",
        face = c("bold")
      ),
      by_layer_y = FALSE,
    ),
  ) +
  geom_col(
    aes(
      x = pos,
      y = depth,
      fill = Chemistry
    ),
  ) +
  # scale_fill_viridis_d(
  #   name = "Chemistry",
  #   option = "D",
  # ) +
  scale_fill_manual(
    name = "Chemistry",
    values = c(viridis::viridis_pal(option = "D")(3), "red")
  ) +
  scale_x_continuous(
    limits = c(0, 17000),
    breaks = seq(0, 17000, 1000),
    labels = seq(0, 17000, 1000),
    expand = expansion(mult = c(0, 0.01)),
  ) +
  scale_y_continuous(
    expand = c(0.01, 0),
    label = scales::label_number(),
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
    legend.position = "none",
    legend.key = element_blank(),
    axis.title.y = element_text(color = "black"),
    axis.text.y = element_text(color = "black"),
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
  coord_cartesian(xlim = c(0, 17000)) +
  labs(y = "Depth") ->
p_depth_all_separated_mixed

ggsave(
  filename = file.path(outdir, "depth_all_position.separated.mixed_cellline.free_y.pdf"),
  plot = wrap_plots(
    p_depth_all_separated_mixed,
    p_mtdna,
    ncol = 1,
    heights = c(15, 1)
  ),
  width = 24,
  height = 12,
  dpi = 300
)


# gse_cell_ratio_variant_meta |>
#   dplyr::select(-cell_ratio_variant) |>
#   tidyr::unnest(cols = anno) |>
#   dplyr::select(gseid, srrid, metrics, Chemistry) |>
#   tidyr::unnest(cols = metrics) |>
#   dplyr::select(gseid, srrid, `Number of Reads`, Chemistry) |>
#   dplyr::bind_rows(
#     tibble::tibble(
#       gseid = "Mixed cellline",
#       srrid = "Mixed cellline",
#       `Number of Reads` = 561027285,
#       Chemistry = "SC3Pv3 Mixed"
#     )
#   ) |>
#   dplyr::mutate(
#     gseid = factor(gseid, levels = c(theposes_depth_ranked$gseid, "Mixed cellline")),
#     Chemistry = factor(Chemistry, levels = c("SC3Pv3 Mixed", "SC3Pv3", "SC5P-R2", "SC5P-PE") |> rev())
#   ) |>
#   ggplot() +
#   geom_violin(
#     aes(
#       x = gseid,
#       y = `Number of Reads`,
#       fill = Chemistry
#     ),
#     alpha = 0.5,
#     # size = 1,
#     color = NA,
#     show.legend = FALSE
#   ) +
#   scale_fill_manual(
#     name = "Chemistry",
#     values = c(viridis::viridis_pal(option = "D")(3), "red")
#   ) +
#   ggbeeswarm::geom_quasirandom(
#     aes(
#       x = gseid,
#       y = `Number of Reads`,
#       color = Chemistry
#     ),
#     size = 1,
#     dodge.width = .75,
#     alpha = .5,
#     varwidth = TRUE
#   ) +
#   scale_color_manual(
#     name = "Chemistry",
#     values = c(viridis::viridis_pal(option = "D")(3), "red")
#   ) +
#   scale_x_discrete(
#     # expand = expansion(mult = c(0.001, 0.001))
#   ) +
#   scale_y_log10(
#     # labels = scales::label_log(),
#     # limits = c(0, 40000),
#     # breaks = seq(0, 40000, 10000),
#     expand = expansion(add = c(0.1, 0.1))
#   ) +
#   theme(
#     panel.background = element_blank(),
#     panel.grid = element_blank(),
#     axis.line = element_line(size = 0.5, color = "black"),
#     axis.title = element_text(
#       size = 16,
#       color = "black",
#       face = "bold"
#     ),
#     axis.text.y = element_text(
#       size = 14,
#       color = "black"
#     ),
#     plot.title = element_text(
#       hjust = 0.5,
#       color = "black",
#       size = 16,
#       face = "bold"
#     ),
#     axis.title.x = element_blank(),
#     legend.position = "inside",
#     legend.position.inside = c(0.7, 0.7)
#   ) +
#   labs(
#     x = "GSE ID",
#     y = "Total reads",
#   ) ->
# p_total_reads_mixed

# ggsave(
#   filename = file.path(outdir, "depth_total_reads_mixed_cellline.pdf"),
#   plot = p_total_reads_mixed,
#   width = 13,
#   height = 5,
#   dpi = 300
# )


# 10x kit version or seq strategy matters -----------------------------------
gse_cell_ratio_variant_meta_xlsx |>
  dplyr::select(-cell_ratio_variant, -anno) |>
  dplyr::left_join(avg_mito_reads, by = "gseid") |>
  dplyr::left_join(avg_totla_reads, by = "gseid") |>
  dplyr::relocate(Publication, .after = avg_totla_reads) |>
  dplyr::mutate(
    Chemistry = factor(Chemistry, levels = c("SC3Pv3", "SC5P-R2", "SC5P-PE"))
  ) |>
  dplyr::arrange(Chemistry, `Avg # of somatic variants`) ->
gse_cell_ratio_variant_meta_xlsx_final

writexl::write_xlsx(
  x = gse_cell_ratio_variant_meta_xlsx_final,
  path = file.path(outdir, "10x-kit-version-matters.xlsx")
)


# footer ------------------------------------------------------------------

# future: :plan(future: :sequential)

# save image --------------------------------------------------------------
outdir <- "/home/liuc9/github/scMOCHA-data/data/out"
load(file = file.path(outdir, "GSE226602-and-other-datasets.rda"))
save.image(
  file = file.path(outdir, "GSE226602-and-other-datasets.rda")
)
