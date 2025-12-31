#!/usr/bin/env Rscript
# Metainfo ----------------------------------------------------------------
# @AUTHOR: Chun-Jie Liu
# @CONTACT: chunjie.sam.liu.at.gmail.com
# @DATE: 2025-01-20 16:07:57
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
      arrowhead_height = unit(3, "mm"),
      arrowhead_width = unit(1, "mm")
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
    coord_cartesian(xlim = c(0, 17000)) -> pg
  pg
}


# load data ---------------------------------------------------------------

basedir <- "/home/liuc9/github/scMOCHA-data/data/scfoundation"

outdir <- file.path(basedir, "out")
dir.create(outdir, showWarnings = FALSE)

gseids <- c(
  "GSE140881", # error, sra decompressed, only one R2 file
  "GSE142595", # error, sra decompressed, it R1 -> I1 error, corrected, on the run. Seq-Well, not 10x
  "GSE149313",
  "GSE154386",
  "GSE159117",
  "GSE162117",
  "GSE167825", # not run, not dumped, dumping, on the run
  "GSE179566", # not run, sra decompressed, only one R2 file
  "GSE188632",
  "GSE192391" # not run, not dumped, checking
)


gseids_meta_1 <- tibble::tibble(
  GSE_ID = gseids,
) |>
  dplyr::mutate(
    samples = purrr::map(
      GSE_ID,
      .f = \(.x) {
        .filename <- file.path(
          basedir,
          .x,
          "out",
          glue::glue("{.x}.cell_ratio_and_variant_clean.csv")
        )
        if (!file.exists(.filename)) {
          logger::log_error("File not found: {.filename}")
          return(NULL)
        }
        data.table::fread(.filename) -> .d
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

gseids_meta <- gseids_meta_1
# body --------------------------------------------------------------------

tibble::tibble(
  gseid = gseids
) |>
  dplyr::mutate(
    cell_ratio_variant = purrr::map(
      .x = gseid,
      .f = \(.gseid) {
        .filename <- file.path(
          basedir,
          .gseid,
          "out",
          glue::glue("{.gseid}.cell_ratio_and_variant_clean.csv")
        )
        if (!file.exists(.filename)) {
          logger::log_error("File not found: {.filename}")
          return(NULL)
        }
        data.table::fread(.filename) |>
          dplyr::select(-Chemistry)
      }
    )
  ) |>
  dplyr::filter(purrr::map_lgl(cell_ratio_variant, ~ !is.null(.x))) |>
  dplyr::mutate(
    anno = purrr::map(
      .x = gseid,
      .f = \(.gseid) {
        readr::read_rds(
          file.path(
            basedir,
            .gseid,
            "out",
            glue::glue("{.gseid}.scmocha.out.rds.gz")
          )
        )
      }
    )
  ) -> gse_data_loaded

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
  ) -> gse_cell_ratio_variant_meta


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
  tidyr::unnest(cols = cell_ratio_variant) -> gse_cell_ratio_variant_meta_xlsx

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
    label = glue::glue(
      "{gseid} ({Chemistry}, {`# of samples`}, {round(`Avg # of somatic variants`, 2)})"
    )
  ) |>
  dplyr::mutate(
    label = factor(label, levels = label)
  ) -> gseid_ranked

# plot average depth and number of somatic variants ------------------------
gse_cell_ratio_variant_meta |>
  dplyr::left_join(
    gseid_ranked |> dplyr::select(gseid, label),
    by = "gseid"
  ) |>
  tidyr::unnest(cols = cell_ratio_variant) |>
  dplyr::mutate(
    Chemistry = factor(
      Chemistry,
      levels = c("SC3Pv2", "SC3Pv3", "SC5P-R2", "SC5P-PE") |> rev()
    ),
  ) -> forplot

cor.test(~ `Depth mean` + `# of somatic variants`, data = forplot) |>
  broom::tidy() -> cor_test_all
cor.test(
  ~ `Depth mean` + `# of somatic variants`,
  data = forplot,
  subset = gseid != "GSE181279"
) |>
  broom::tidy() -> cor_test_250k

pcc <- readr::read_tsv(
  file = "https://raw.githubusercontent.com/chunjie-sam-liu/chunjie-sam-liu.life/master/public/data/pcc.tsv"
) |>
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
    method = "loess",
    se = FALSE,
    color = "black",
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
    # limits = c(0, 60000),
    # breaks = seq(0, 60000, 10000),
  ) +
  scale_y_continuous(
    labels = scales::label_number(),
    # limits = c(0, 80),
    # breaks = seq(0, 80, 10),
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
    title = glue::glue(
      "All samples test, Pearson's r = {round(cor_test_all$estimate, 2)}, p-value = {scales::pvalue(cor_test_all$p.value)}"
    ),
    # subtitle = glue::glue("Exclude GSE181279 samples, Pearson's r = {round(cor_test_250k$estimate, 2)}, p-value = {scales::pvalue(cor_test_250k$p.value)}")
  ) -> p_somatic_variant

ggsave(
  filename = file.path(outdir, "somatic_variant.pdf"),
  plot = p_somatic_variant,
  width = 13,
  height = 7,
  dpi = 300
)


# ! depth --------------------------------------------------------------------

gse_cell_ratio_variant_meta |>
  dplyr::select(gseid, anno, Chemistry) |>
  tidyr::unnest(cols = anno) |>
  dplyr::select(gseid, srrid, depth, Chemistry) |>
  tidyr::unnest(cols = depth) -> all_gseid_depth

all_gseid_depth |>
  dplyr::group_by(
    Chemistry,
    pos
  ) |>
  dplyr::summarise(
    depth = mean(depth, na.rm = TRUE)
  ) |>
  dplyr::ungroup() |>
  dplyr::mutate(
    Chemistry = factor(
      Chemistry,
      levels = c("SC3Pv2", "SC3Pv3", "SC5P-R2", "SC5P-PE") |> rev()
    )
  ) -> all_gseid_depth_forplot


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
    legend.position = "top",
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
  labs(y = "Depth") -> p_depth_all

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

# all position in three chemistry separated dataset -------------------------------------------------

thevariants <- c("1670A>G", "3173G>A")

theposes <- thevariants |>
  purrr::map(~ gsub(pattern = "[>|AGCT]", "", x = .)) |>
  purrr::map_int(as.integer)

all_gseid_depth |>
  dplyr::filter(pos %in% theposes) -> theposes_depth

theposes_depth |>
  dplyr::filter(pos == theposes[2]) |>
  dplyr::group_by(gseid, Chemistry) |>
  dplyr::summarise(
    mean_depth = mean(depth, na.rm = TRUE)
  ) |>
  dplyr::ungroup() |>
  dplyr::mutate(
    Chemistry = factor(
      Chemistry,
      levels = c("SC3Pv2", "SC3Pv3", "SC5P-R2", "SC5P-PE") |> rev()
    )
  ) |>
  dplyr::arrange(Chemistry, dplyr::desc(mean_depth)) -> theposes_depth_ranked


all_gseid_depth |>
  dplyr::group_by(
    Chemistry,
    gseid,
    pos
  ) |>
  dplyr::summarise(
    depth = mean(depth, na.rm = TRUE)
  ) |>
  dplyr::ungroup() |>
  dplyr::mutate(
    gseid = factor(gseid, levels = theposes_depth_ranked$gseid),
    Chemistry = factor(
      Chemistry,
      levels = c("SC3Pv2", "SC3Pv3", "SC5P-R2", "SC5P-PE") |> rev()
    )
  ) -> all_gseid_depth_forplot_separated

chem_color <- tibble::tibble(
  color = c(viridis::viridis_pal(option = "D")(4)),
  Chemistry = c("SC3Pv2", "SC3Pv3", "SC5P-R2", "SC5P-PE") |> rev()
)

theposes_depth_ranked |>
  dplyr::left_join(
    chem_color,
    by = "Chemistry"
  ) -> theposes_depth_ranked_chem_color

all_gseid_depth_forplot_separated |>
  ggplot() +
  ggh4x::facet_wrap2(
    ~gseid,
    ncol = 1,
    strip.position = "right",
    scales = "free_y",
    strip = ggh4x::strip_themed(
      background_y = ggh4x::elem_list_rect(
        # fill = c(rep(viridis::viridis_pal(option = "D")(3), each = 3), "red")
        fill = theposes_depth_ranked_chem_color$color
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
    values = theposes_depth_ranked_chem_color$color |> unique()
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
    legend.position = "top",
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
  labs(y = "Depth") -> p_depth_all_separated

p_mtdna <- fn_plot_mtdna()

ggsave(
  filename = file.path(outdir, "depth_all_position.separated.free_y.pdf"),
  plot = wrap_plots(
    p_depth_all_separated,
    p_mtdna,
    ncol = 1,
    heights = c(15, 1)
  ),
  width = 23,
  height = 12,
  dpi = 300
)

# footer ------------------------------------------------------------------

# future: :plan(future: :sequential)

# save image --------------------------------------------------------------
