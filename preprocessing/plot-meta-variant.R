#!/usr/bin/env Rscript
# Metainfo ----------------------------------------------------------------

# @AUTHOR: Chun-Jie Liu
# @CONTACT: chunjie.sam.liu.at.gmail.com
# @DATE: Thu Aug 22 13:26:42 2024
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

# log_info('Starting the script...')
# log_debug('This is the second log line')
# log_trace('Note that the 2nd line is being placed right after the 1st one.')
# log_success('Doing pretty well so far!')
# log_warn('But beware, as some errors might come :/')
# log_error('This is a problem')
# log_debug('Note that getting an error is usually bad')
# log_error('This is another problem')
# log_fatal('The last problem')

# future::plan(future::multisession, workers = 10)

# function ----------------------------------------------------------------


# load data ---------------------------------------------------------------

basedir <- "/home/liuc9/github/scMOCHA-data/data"
# basedir <- "/home/liuc9/github/scMOCHA/06-bigdata"
gseid <- "GSE226602"
datadir <- file.path(
  basedir, gseid
)

outdir <- file.path(
  datadir,
  "out"
)


# body --------------------------------------------------------------------

pheno <- data.table::fread(
  file = file.path(
    datadir,
    "{gseid}.pheno.select.csv" |> glue::glue()
  )
)

variant <- readr::read_rds(
  file.path(
    outdir,
    "{gseid}.scmocha.out.rds" |> glue::glue()
  )
)


# stat --------------------------------------------------------------------


# variant |>
#   dplyr::mutate(
#     nmut = purrr::map_int(
#       .x = anno,
#       .f = \(.x) {
#         if (is.null(.x)) {
#           return(NA_integer_)
#         }
#         nrow(.x)
#       }
#     )
#   ) |>
#   dplyr::mutate(
#     haplogroup = purrr::map2(
#       .x = anno,
#       .y = srrid,
#       .f = \(.x, .y) {
#         message(.y)
#         if (is.null(.x)) {
#           return(
#             tibble::tibble(
#               Haplogroup = NA_character_,
#               Verbose_haplogroup = NA_character_
#             )
#           )
#         }
#         .x |>
#           dplyr::select(Haplogroup, Verbose_haplogroup) |>
#           dplyr::filter(!is.na(Haplogroup)) |>
#           dplyr::filter(Haplogroup != "") |>
#           dplyr::distinct() |>
#           dplyr::mutate_all(.funs = as.character) ->
#         .xx

#         if (nrow(.xx) == 0) {
#           tibble::tibble(
#             Haplogroup = NA_character_,
#             Verbose_haplogroup = NA_character_
#           )
#         } else {
#           .xx
#         }
#       }
#     )
#   ) |>
#   tidyr::unnest(cols = haplogroup) |>
#   dplyr::inner_join(
#     pheno,
#     by = "srrid"
#   ) |>
#   tidyr::unnest(
#     cols = cell_stats
#   ) |>
#   dplyr::mutate(
#     ratio = round(`number of cells after filtering` / `estimated number of cells`, 2)
#   ) ->
# metadata_anno

variant |>
  dplyr::inner_join(
    pheno,
    by = "srrid"
  ) ->
metadata_anno

# metadata_anno |>
#   dplyr::arrange(disease, age, gender) |>
#   dplyr::select(
#     srrid, age, gender,
#     disease, genotype,
#     `Median UMI/cell` = `median UMI counts per cell`,
#     `Median genes/cell` = `median genes per cell`,
#     `# of cells` = `estimated number of cells`,
#     `# cells after filter` = `number of cells after filtering`,
#     `Cell ratio` = ratio,
#     `# of variants` = nmut,
#     Haplogroup = Haplogroup,
#     Haplogroup_v = Verbose_haplogroup
#   ) ->
# metadata_clean

# metadata_clean |>
#   writexl::write_xlsx(
#     path = file.path(
#       outdir,
#       "/metadata_clean.xlsx"
#     )
#   )


# cell ratio --------------------------------------------------------------

metadata_anno |>
  dplyr::arrange(disease, age, gender) |>
  dplyr::select(srrid, disease, celltype_ratio) |>
  # dplyr::arrange(disease) |>
  dplyr::mutate(color = dplyr::case_match(
    disease,
    "Alzheimers Disease" ~ ggsci::pal_jama()(4)[[1]],
    "Healthy Control" ~ ggsci::pal_jama()(4)[[2]]
  )) |>
  dplyr::arrange(dplyr::desc(dplyr::row_number())) ->
for_ratio_plot

for_ratio_plot |>
  tidyr::unnest(cols = celltype_ratio) |>
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
p_cellratio


ggsave(
  filename = "Cell_ratio.pdf",
  plo = p_cellratio,
  device = "pdf",
  width = 18,
  height = 10,
  path = outdir
)


# depth -------------------------------------------------------------------


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
  pg
  pg
}

metadata_anno |>
  dplyr::arrange(disease, age, gender) |>
  dplyr::select(srrid, disease, depth) |>
  # dplyr::arrange(disease) |>
  dplyr::mutate(color = dplyr::case_match(
    disease,
    "Alzheimers Disease" ~ ggsci::pal_jama()(4)[[1]],
    "Healthy Control" ~ ggsci::pal_jama()(4)[[2]]
  )) |>
  dplyr::arrange(dplyr::desc(dplyr::row_number())) ->
for_depth_plot

for_depth_plot |>
  tidyr::unnest(cols = depth) |>
  ggplot(aes(x = pos, y = depth, fill = `color`)) +
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
  facet_wrap(
    facets = ~disease,
    ncol = 1,
    strip.position = "right"
  ) +
  labs(y = "Depth") ->
p_mt_depth
p_mt_depth


wrap_plots(
  p_mt_depth,
  fn_plot_gene(),
  ncol = 1,
  heights = c(0.9, 0.1)
) ->
p_depth
p_depth

ggsave(
  filename = "Sample_depth_merge.pdf",
  plo = p_depth,
  device = "pdf",
  width = 15,
  height = 8,
  path = outdir
)


# all factor correlations -------------------------------------------------

metadata_anno |>
  dplyr::arrange(disease, age, gender) |>
  # dplyr::arrange(disease) |>
  dplyr::mutate(color = dplyr::case_match(
    disease,
    "Alzheimers Disease" ~ ggsci::pal_jama()(4)[[1]],
    "Healthy Control" ~ ggsci::pal_jama()(4)[[2]]
  )) |>
  dplyr::arrange(dplyr::desc(dplyr::row_number())) |>
  dplyr::glimpse() |>
  dplyr::select(
    srrid,
    genotype,
    disease,
    age,
    Sex = gender,
    nmut,
    `median UMI counts per cell`,
    `number of cells after filtering`,
    depth
  ) |>
  dplyr::mutate(
    n_na = purrr::map(
      .x = depth,
      .f = \(.d) {
        .d |>
          dplyr::summarise(
            dep_s = sum(depth),
            dep_mea = mean(depth),
            dep_med = median(depth)
          )
      }
    )
  ) |>
  dplyr::select(-depth) |>
  tidyr::unnest(cols = n_na) |>
  dplyr::mutate(
    Sex = factor(Sex),
    genotype = factor(genotype),
    disease = factor(disease)
  ) ->
metadata_anno_depth_dep

correlation::correlation(
  metadata_anno_depth_dep |>
    dplyr::select(-dep_s, -dep_mea),
  p_adjust = "none"
) |>
  summary(redundant = TRUE) ->
cor_summr

plot(cor_summr)

cor_summr |>
  as.data.frame() |>
  tidyr::pivot_longer(
    cols = -Parameter,
    names_to = "var2",
    values_to = "pval"
  ) |>
  dplyr::filter(!is.na(pval)) |>
  dplyr::filter(Parameter != var2) |>
  ggplot(aes(
    x = Parameter,
    y = var2,
    fill = pval
  )) +
  geom_tile() +
  geom_text(aes(label = round(pval, 2))) +
  scale_fill_gradient2(
    name = "R",
    low = "blue",
    mid = "white",
    high = "red",
    midpoint = 0,
    # space = "Lab",
    # na.value = "grey50",
    guide = "colourbar",
    # aesthetics = "colour"
  ) +
  scale_x_discrete(
    limits = c("nmut", "dep_med", "number of cells after filtering", "median UMI counts per cell", "age"),
    labels = c("# variants", "median depth", "# cells", "median UMI/cell", "Age") |> stringr::str_to_sentence()
  ) +
  scale_y_discrete(
    limits = c("nmut", "dep_med", "number of cells after filtering", "median UMI counts per cell", "age") |> rev(),
    labels = c("# variants", "median depth", "# cells", "median UMI/cell", "Age") |> rev() |> stringr::str_to_sentence()
  ) +
  theme(
    panel.background = element_blank(),
    panel.grid = element_blank(),
    axis.ticks = element_blank(),
    axis.title = element_blank(),
    axis.text = element_text(
      color = "black",
      size = 14
    )
  ) ->
p_cor
p_cor

ggsave(
  filename = "All-factor-correlations.pdf",
  plo = p_cor,
  device = "pdf",
  width = 9,
  height = 6,
  path = outdir
)


# sex ---------------------------------------------------------------------


t.test(nmut ~ Sex, data = metadata_anno_depth_dep) |>
  broom::tidy()
metadata_anno_depth_dep |>
  ggplot(aes(
    x = Sex,
    y = nmut,
  )) +
  # geom_violin() +
  geom_boxplot(
    aes(fill = disease),
    width = 0.5,
    show.legend = T
  ) +
  geom_point(position = position_jitter(width = 0.3)) +
  theme_bw() +
  ggsci::scale_fill_aaas(
    name = "Disease"
  ) +
  scale_x_discrete(
    limits = c("male", "female"),
    labels = c("Male", "Female")
  ) +
  theme(
    axis.title = element_text(size = 16, face = "bold", colour = "black"),
    axis.text.x = element_text(size = 14, face = "bold", colour = "black"),
    axis.title.x = element_blank(),
    legend.position = "top"
  ) +
  labs(
    x = "Gender",
    y = "# of variants"
  ) ->
gender_cor_plot
gender_cor_plot


ggsave(
  filename = "All-factor-correlations-sex.pdf",
  plo = gender_cor_plot,
  device = "pdf",
  width = 4,
  height = 4,
  path = outdir
)
# Age ---------------------------------------------------------------------


cor.test(
  formula = ~ nmut + age,
  data = metadata_anno_depth_dep
) ->
cta

cor.test(
  formula = ~ nmut + age,
  data = metadata_anno_depth_dep |>
    dplyr::filter(disease == "Healthy Control")
) ->
cta_mci

cor.test(
  formula = ~ nmut + age,
  data = metadata_anno_depth_dep |>
    dplyr::filter(disease == "Alzheimers Disease")
) ->
cta_ad

yhight <- 32
xwidth <- 50

metadata_anno_depth_dep |>
  # dplyr::filter(nmut >10) |>
  dplyr::mutate(dia = disease) |>
  dplyr::mutate(Age = age) |>
  dplyr::mutate(
    label = glue::glue(
      "N variants = {nmut}\n Median depth = {dep_med}\n Gender = {Sex}"
    )
  ) |>
  dplyr::mutate(
    dia = factor(dia, levels = c("Healthy Control", "Alzheimers Disease"))
  ) |>
  ggplot(aes(
    x = Age,
    y = nmut
  )) +
  geom_point(aes(color = dia), show.legend = FALSE) +
  geom_smooth(method = "loess", se = FALSE, color = "black", linetype = 21) +
  geom_smooth(aes(color = dia), method = "glm", se = FALSE) +
  ggrepel::geom_text_repel(
    aes(label = label),
    # box.padding = 0.5,
    max.overlaps = 10,
    # max.overlaps = Inf
    size = 3,
    min.segment.length = 0,
    seed = 42,
    box.padding = 0.5
  ) +
  ggsci::scale_color_jama(
    name = "Disease type"
  ) +
  annotate(
    geom = "segment",
    x = xwidth,
    y = yhight,
    xend = xwidth + 1,
    yend = yhight,
    linetype = 21,
    colour = "black",
    linewidth = 1
  ) +
  annotate(
    geom = "text",
    x = xwidth + 3,
    y = yhight,
    size = 5,
    label = latex2exp::TeX(glue::glue("$\\rho$={round(cta$estimate, 2)}, $P$={round(cta$p.value,3)}")),
    fontface = "bold",
  ) +
  annotate(
    geom = "segment",
    x = xwidth,
    y = yhight - 2,
    xend = xwidth + 1,
    yend = yhight - 2,
    linetype = 1,
    colour = ggsci::pal_jama()(2)[1],
    linewidth = 1
  ) +
  annotate(
    geom = "text",
    x = xwidth + 3,
    y = yhight - 2,
    size = 5,
    label = latex2exp::TeX(glue::glue("$\\rho$={round(cta_mci$estimate, 2)}, $P$={round(cta_mci$p.value,3)}")),
    fontface = "bold",
    color = ggsci::pal_jama()(2)[1],
  ) +
  annotate(
    geom = "segment",
    x = xwidth,
    y = yhight - 4,
    xend = xwidth + 1,
    yend = yhight - 4,
    linetype = 1,
    colour = ggsci::pal_jama()(2)[2],
    linewidth = 1
  ) +
  annotate(
    geom = "text",
    x = xwidth + 3,
    y = yhight - 4,
    size = 5,
    label = latex2exp::TeX(glue::glue("$\\rho$={round(cta_ad$estimate, 2)}, $P$={round(cta_ad$p.value,3)}")),
    fontface = "bold",
    color = ggsci::pal_jama()(2)[2]
  ) +
  theme_bw() +
  theme(
    # panel.grid = element_blank(),
    axis.text = element_text(size = 14, colour = "black"),
    axis.title = element_text(size = 16, face = "bold", colour = "black"),
    legend.position = "top"
  ) +
  labs(
    x = "Age",
    y = "# of variants"
  ) ->
p_linear_1
p_linear_1

ggsave(
  filename = "All-factor-correlations-linear-age-nvariant.pdf",
  plo = p_linear_1,
  device = "pdf",
  width = 15,
  height = 8,
  path = outdir
)


# Inter section -----------------------------------------------------------

metadata_anno |>
  dplyr::mutate(
    variant = purrr::map2(
      .x = anno,
      .y = srrdir,
      .f = function(.x, .y) {
        if (is.na(.y)) {
          return(NULL)
        }
        .x |>
          dplyr::mutate(
            variant = glue::glue("{Position}{Ref}>{Alt}")
          ) |>
          dplyr::select(variant)
      }
    )
  ) ->
metadata_anno_depth_variant

metadata_anno_depth_variant |>
  dplyr::mutate(color = dplyr::case_match(
    disease,
    "Alzheimers Disease" ~ ggsci::pal_jama()(4)[[1]],
    "Healthy Control" ~ ggsci::pal_jama()(4)[[2]]
  )) |>
  dplyr::select(srrid, source_name = disease, variant, color) |>
  dplyr::filter(!purrr::map_lgl(.x = variant, .f = is.null)) ->
for_variant

fn_upset_plot <- function(.x) {
  # .x <- "nCoV_PBMC(severe)"
  library(ggupset)
  for_variant |>
    dplyr::filter(source_name == .x) ->
  d

  d |>
    tidyr::unnest(cols = variant) |>
    dplyr::select(-source_name) |>
    dplyr::group_by(variant) |>
    tidyr::nest() |>
    dplyr::ungroup() |>
    dplyr::mutate(
      srrid = purrr::map(
        .x = data,
        .f = function(.x) {
          .x |> dplyr::pull(srrid)
        }
      )
    ) |>
    dplyr::select(-data) ->
  dd


  dd |>
    ggplot(aes(x = srrid)) +
    geom_bar(width = 0.6, fill = d$color[1]) +
    geom_text(
      stat = "count",
      aes(label = after_stat(count)),
      vjust = -0.5,
      color = "black",
      size = 3,
      fontface = "bold"
    ) +
    scale_x_upset(order_by = "degree") +
    scale_y_continuous(
      expand = expansion(mult = c(0, 0.1), add = 0)
    ) +
    theme_combmatrix(
      combmatrix.label.make_space = TRUE,
      combmatrix.panel.point.color.fill = d$color[1],
      combmatrix.panel.line.size = 0,
      combmatrix.label.text = element_text(
        # size = 12,
        color = "black",
        face = "bold"
      ),
      combmatrix.label.extra_spacing = 1,
      combmatrix.panel.striped_background.color.one = "white",
      combmatrix.panel.striped_background.color.two = "grey",
    ) +
    labs(
      y = "# of Variants",
      x = "",
      title = .x
    ) +
    theme(
      panel.background = element_blank(),
      panel.grid = element_blank(),
      axis.line = element_line(size = 0.5, color = "black"),
      axis.title.y = element_text(
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
    ) ->
  .p_up

  ggsave(
    plot = .p_up,
    filename = "upset-{.x}.pdf" |> glue::glue(),
    path = outdir,
    width = 15,
    height = 7,
    device = "pdf"
  )

  dd |>
    dplyr::mutate(n = purrr::map_int(
      .x = srrid,
      .f = length
    )) |>
    dplyr::mutate(
      sharing = purrr::map_chr(
        .x = srrid,
        .f = paste0,
        collapse = ","
      )
    ) |>
    dplyr::arrange(n) |>
    dplyr::select(-srrid) ->
  .v

  list(
    v = .v,
    p_up = .p_up
  )
}

for_variant$source_name |>
  unique() |>
  purrr::map(
    .f = fn_upset_plot
  ) ->
p_ups

(p_ups[[1]]$p_up | p_ups[[2]]$p_up) +
  plot_annotation(tag_levels = "A") ->
p_ups_together
p_ups_together

ggsave(
  plot = p_ups_together,
  filename = "upset-all.pdf" |> glue::glue(),
  path = outdir,
  width = 25,
  height = 12,
  device = "pdf"
)


# save to xlsx ------------------------------------------------------------

names(p_ups) <- unique(metadata_anno_depth_variant$disease)


p_ups |>
  purrr::map("v") |>
  writexl::write_xlsx(
    path = file.path(
      outdir,
      "/upset-variants.xlsx"
    )
  )


# metadata_anno |> dplyr::glimpse()
# metadata_anno$anno[[1]]

# heteroplasmy ------------------------------------------------------------

# GSM7080019
metadata_anno |>
  dplyr::glimpse()
metadata_anno$srrdir[[14]]
metadata_anno$anno[[14]] |>
  dplyr::mutate(
    v = glue::glue("{Position}{Ref}>{Alt}")
  ) ->
sel_anno

metadata_anno$coverage[[14]] ->
sel_cov
metadata_anno$hetero[[14]] ->
forplot

forplot |>
  dplyr::group_by(variant) |>
  dplyr::summarise(maf = sum(af, na.rm = T)) |>
  dplyr::arrange(-maf) ->
sort_variant

sel_anno |>
  # dplyr::filter(Haplogroup == "T2b") |>
  dplyr::mutate(fill = ifelse(Haplogroup == "T2b", "red", "white")) |>
  dplyr::mutate(color = ifelse(Haplogroup == "T2b", "white", "black")) |>
  dplyr::mutate(
    variant = factor(v, sort_variant$variant)
  ) |>
  dplyr::arrange(variant) ->
t2b_variant

forplot |>
  dplyr::mutate(
    variant = factor(variant, sort_variant$variant |> unique())
  ) |>
  dplyr::arrange(variant) ->
forplot_t2b

library(ggh4x)
forplot_t2b |>
  # dplyr::filter(variant == "2706A>G") |>
  ggplot(aes(x = celltype, y = af)) +
  geom_col(aes(fill = af)) +
  # facet_wrap(~variant, ncol = 10) +
  ggh4x::facet_wrap2(
    ~variant,
    ncol = 10,
    strip = ggh4x::strip_themed(
      # Horizontal strips
      # background_x = elem_list_rect(fill = c("limegreen", "dodgerblue")),
      background_x = elem_list_rect(
        fill = t2b_variant$fill
      ),
      text_x = elem_list_text(
        colour = t2b_variant$color,
        face = c("bold")
      ),
      by_layer_x = FALSE,
      # # Vertical strips
      # background_y = elem_list_rect(
      #   fill = c("gold", "tomato", "deepskyblue")
      # ),
      # text_y = elem_list_text(angle = c(0, 90)),
      # by_layer_y = FALSE
    )
  ) +
  scale_fill_gradient2(
    low = "white",
    mid = "red",
    high = "#3B0049",
    midpoint = 0.5
  ) +
  theme(
    panel.background = element_blank(),
    panel.grid = element_blank(),
    # axis.ticks = element_blank(),
    axis.title = element_blank(),
    axis.text = element_text(
      color = "black",
      # size = 18
    ),
    legend.position = "none ",
    plot.title = element_text(
      size = 16,
      hjust = 0.5
    ),
    # strip.background = element_rect(
    #   # fill = NA,
    #   # fill = t2b_variant$color,
    #   color = "black",
    # ),
    # strip.text = element_text(
    #   # color = "black",
    #   # color = forplot_t2b$t2b,
    #   size = 10,
    #   face = "bold"
    # ),
    axis.line = element_line(
      color = "black"
    ),
    axis.text.x = element_text(
      angle = 45,
      hjust = 1
    )
  ) ->
p
p

ggsave(
  filename = "GSM7080019-cluster-bar.pdf",
  plot = p,
  device = "pdf",
  path = outdir,
  width = 20,
  height = 10
)

unique(forplot$variant)[c(6, 10:15, 17, 18, 19:27, 30:36, 38, 40, 43:48, 50, 52, 54, 55, 58, 59, 61:65, 70)] ->
sel_v

forplot |>
  dplyr::filter(variant %in% sel_v) |>
  dplyr::mutate(
    variant = factor(variant, sel_v)
  ) |>
  # dplyr::filter(variant == "2706A>G") |>
  ggplot(aes(x = celltype, y = af)) +
  geom_col(aes(fill = af)) +
  facet_wrap(~variant, ncol = 9) +
  scale_fill_gradient2(
    low = "white",
    mid = "red",
    high = "#3B0049",
    midpoint = 0.5
  ) +
  theme(
    panel.background = element_blank(),
    panel.grid = element_blank(),
    # axis.ticks = element_blank(),
    axis.title = element_blank(),
    axis.text = element_text(
      color = "black",
      # size = 18
    ),
    legend.position = "none ",
    plot.title = element_text(
      size = 16,
      hjust = 0.5
    ),
    strip.background = element_rect(
      fill = NA,
      color = "black",
    ),
    strip.text = element_text(
      color = "black",
      size = 10,
      face = "bold"
    ),
    axis.line = element_line(
      color = "black"
    ),
    axis.text.x = element_text(
      angle = 45,
    )
  ) ->
p
p


ggsave(
  filename = "GSM7080053-cluster-bar-sel.pdf",
  plot = p,
  device = "pdf",
  path = outdir,
  width = 17,
  height = 8
)


tibble::tibble(
  variant = sel_v
) |>
  dplyr::mutate(pos = gsub(">|A|G|C|T", "", variant)) |>
  dplyr::mutate(pos = as.integer(pos)) |>
  dplyr::left_join(sel_cov, by = "pos") |>
  dplyr::mutate(variant = factor(variant, sel_v)) |>
  ggplot(aes(x = celltype, y = count)) +
  geom_col(aes(fill = count)) +
  facet_wrap(~variant, ncol = 9) +
  scale_fill_gradient(
    low = "white",
    # mid = "red",
    high = "#3B0049",
    # midpoint = 0.5
  ) +
  theme(
    panel.background = element_blank(),
    panel.grid = element_blank(),
    # axis.ticks = element_blank(),
    axis.title = element_blank(),
    axis.text = element_text(
      color = "black",
      # size = 18
    ),
    legend.position = "none ",
    plot.title = element_text(
      size = 16,
      hjust = 0.5
    ),
    strip.background = element_rect(
      fill = NA,
      color = "black",
    ),
    strip.text = element_text(
      color = "black",
      size = 10,
      face = "bold"
    ),
    axis.line = element_line(
      color = "black"
    )
  ) ->
p_c
p_c

ggsave(
  filename = "GSM7080053-cluster-bar-sel-cov.pdf",
  plot = p_c,
  device = "pdf",
  path = outdir,
  width = 17,
  height = 8
)


# filter ------------------------------------------------------------------

metadata_anno |>
  dplyr::mutate(
    hetero_filter = purrr::map(
      .x = hetero,
      .f = \(.d) {
        .d |>
          dplyr::group_by(variant) |>
          dplyr::summarise(
            mean_af = mean(af, na.rm = T)
          ) ->
        .dd


        .dd$mean_af |> hist()
        .dd |>
          dplyr::filter(mean_af < 0.8) |>
          dplyr::filter(mean_af > 0.05) ->
        .ddd

        .d |>
          dplyr::filter(variant %in% .ddd$variant)
      }
    )
  ) ->
metadata_anno_filter



metadata_anno_filter |>
  dplyr::select(
    srrid,
    genotype,
    disease = disease,
    dia = disease,
    age,
    Sex = gender,
    hetero_filter
  ) |>
  dplyr::mutate(
    dia = factor(dia, levels = c("Healthy Control", "Alzheimers Disease"))
  ) |>
  dplyr::mutate(
    Sex = factor(Sex),
    genotype = factor(genotype),
    disease = factor(disease)
  ) |>
  dplyr::mutate(dia = disease) |>
  dplyr::mutate(Age = age) |>
  dplyr::mutate(
    dia = factor(dia, levels = c("Healthy Control", "Alzheimers Disease"))
  ) |>
  tidyr::unnest(cols = hetero_filter) ->
metadata_anno_filter_sel


metadata_anno_filter_sel |>
  dplyr::filter(!is.na(af)) |>
  dplyr::filter(af >= 0.05) |>
  dplyr::mutate(
    variant = factor(variant)
  ) |>
  dplyr::count(variant, celltype) |>
  plotme::count_to_sunburst()



metadata_anno_filter_sel |>
  dplyr::filter(!is.na(af)) |>
  dplyr::filter(af >= 0.05) |>
  dplyr::mutate(
    variant = factor(variant)
  ) |>
  dplyr::select(variant, celltype) |>
  dplyr::distinct() |>
  dplyr::group_by(variant) |>
  dplyr::count() |>
  dplyr::ungroup() |>
  dplyr::filter(n >= 8) |>
  dplyr::select(variant) |>
  dplyr::distinct() |>
  head(20) ->
sv

metadata_anno_filter_sel |>
  dplyr::filter(!is.na(af)) |>
  dplyr::filter(af >= 0.05) |>
  dplyr::filter(variant %in% sv$variant) |>
  ggplot(aes(x = af, fill = celltype)) +
  geom_density() +
  facet_wrap(~variant)


metadata_anno_filter_sel |>
  # dplyr::filter(variant == "8303A>G") |>
  dplyr::filter(variant %in% c(
    # "1602C>A", "1604G>A",
    # "1610A>T", "1670A>G",
    "2442T>C",
    "2517A>T", "2617A>G",
    "3173G>A", "3176A>T", "3178T>A",
    "7526A>G",
    "8303A>G",
    "8362T>G",
    "10413A>G"
  )) |>
  dplyr::mutate(
    variant = factor(variant)
  ) |>
  ggplot(aes(x = age, y = af)) +
  geom_point(aes(color = dia), show.legend = FALSE) +
  geom_smooth(method = "loess", se = FALSE, color = "black", linetype = 21) +
  geom_smooth(aes(color = dia), method = "glm", se = FALSE) +
  ggsci::scale_color_jama(
    name = "Disease type"
  ) +
  scale_y_continuous(
    limits = c(0, 1)
  ) +
  facet_grid(
    rows = vars(variant),
    cols = vars(celltype),
    switch = "y"
  ) +
  theme(
    panel.background = element_blank(),
    panel.grid = element_blank(),
    # axis.ticks = element_blank(),
    axis.title = element_blank(),
    axis.text = element_text(
      color = "black",
      # size = 18
    ),
    legend.position = "bottom",
    legend.background = element_blank(),
    legend.key = element_blank(),
    plot.title = element_text(
      size = 16,
      hjust = 0.5
    ),
    strip.background = element_rect(
      fill = NA,
      color = "black",
    ),
    strip.text = element_text(
      color = "black",
      size = 10,
      face = "bold"
    ),
    axis.line = element_line(
      color = "black"
    )
  ) ->
p_cell_age
p_cell_age
ggsave(
  filename = "Celltype age.pdf",
  plo = p_cell_age,
  device = "pdf",
  width = 15,
  height = 8,
  path = outdir
)
the_v <- c(
  "7526A>G",
  "8303A>G"
)



metadata_anno_filter_sel |>
  dplyr::filter(variant == "7526A>G") |>
  dplyr::filter(celltype == "DC") |>
  dplyr::filter(dia == "Alzheimers Disease") |>
  ggstatsplot::ggscatterstats(
    x = age,
    y = af,
    color = dia
  )

metadata_anno_filter_sel |>
  tidyr::nest(.by = c("variant", "celltype")) |>
  dplyr::mutate(age_cor = purrr::map(
    .x = data,
    .f = \(.d) {
      tryCatch(
        expr = {
          if (nrow(.d) < 6) {
            return(NULL)
          }
          cor.test(~ af + age, data = .d) |>
            broom::tidy() |>
            dplyr::select(
              cor = estimate,
              pval = p.value
            )
        },
        error = \(e) {
          NULL
        }
      )
    }
  )) |>
  dplyr::mutate(
    sex_t = purrr::map(
      .x = data,
      .f = \(.d) {
        tryCatch(
          expr = {
            t.test(af ~ Sex, .d) |>
              broom::tidy() |>
              dplyr::select(pval = p.value)
          },
          error = \(e) {
            NULL
          }
        )
      }
    )
  ) |>
  dplyr::mutate(
    disease_t = purrr::map(
      .x = data,
      .f = \(.d) {
        tryCatch(
          expr = {
            t.test(af ~ disease, .d) |>
              broom::tidy() |>
              dplyr::select(pval = p.value)
          },
          error = \(e) {
            NULL
          }
        )
      }
    )
  ) ->
metadata_anno_filter_sel_test

metadata_anno_filter_sel_test |>
  dplyr::select(variant, celltype, age_cor) |>
  tidyr::unnest(cols = age_cor) |>
  dplyr::filter(pval < 0.05, abs(cor) > 0.5) ->
select_v_by_age


metadata_anno_filter_sel |>
  dplyr::filter(variant == select_v_by_age$variant) |>
  dplyr::mutate(
    variant = factor(variant)
  ) |>
  ggplot(aes(x = age, y = af)) +
  geom_point(aes(color = dia), show.legend = FALSE) +
  geom_smooth(method = "loess", se = FALSE, color = "black", linetype = 21) +
  geom_smooth(aes(color = dia), method = "glm", se = FALSE) +
  ggsci::scale_color_jama(
    name = "Disease type"
  ) +
  scale_y_continuous(
    limits = c(0, 1)
  ) +
  facet_grid(
    rows = vars(variant),
    cols = vars(celltype),
    switch = "y"
  ) +
  theme(
    panel.background = element_blank(),
    panel.grid = element_blank(),
    # axis.ticks = element_blank(),
    axis.title = element_blank(),
    axis.text = element_text(
      color = "black",
      # size = 18
    ),
    legend.position = "bottom",
    legend.background = element_blank(),
    legend.key = element_blank(),
    plot.title = element_text(
      size = 16,
      hjust = 0.5
    ),
    strip.background = element_rect(
      fill = NA,
      color = "black",
    ),
    strip.text = element_text(
      color = "black",
      size = 10,
      face = "bold"
    ),
    axis.line = element_line(
      color = "black"
    )
  ) ->
p_cell_age
p_cell_age
ggsave(
  filename = "Celltype age.pdf",
  plo = p_cell_age,
  device = "pdf",
  width = 15,
  height = 8,
  path = outdir
)

metadata_anno_filter_sel_test |>
  dplyr::select(variant, sex_t) |>
  tidyr::unnest(cols = sex_t) |>
  dplyr::filter(pval < 0.05) ->
select_v_by_sex



metadata_anno_filter_sel |>
  dplyr::filter(variant == select_v_by_sex$variant) |>
  dplyr::mutate(
    variant = factor(variant)
  ) |>
  ggplot(aes(
    x = Sex,
    y = af,
  )) +
  # geom_violin() +
  geom_boxplot(
    aes(color = Sex),
    width = 0.5,
    show.legend = T
  ) +
  # geom_point(position = position_jitter(width = 0.3)) +
  ggsci::scale_color_aaas(
    name = "Sex"
  ) +
  scale_x_discrete(
    limits = c("male", "female"),
    labels = c("Male", "Female")
  ) +
  scale_y_continuous(
    limits = c(0, 1)
  ) +
  facet_grid(
    rows = vars(variant),
    cols = vars(celltype),
    switch = "y"
  ) +
  theme(
    panel.background = element_blank(),
    panel.grid = element_blank(),
    # axis.ticks = element_blank(),
    axis.title = element_blank(),
    axis.text = element_text(
      color = "black",
      # size = 18
    ),
    legend.position = "bottom",
    legend.background = element_blank(),
    legend.key = element_blank(),
    plot.title = element_text(
      size = 16,
      hjust = 0.5
    ),
    strip.background = element_rect(
      fill = NA,
      color = "black",
    ),
    strip.text = element_text(
      color = "black",
      size = 10,
      face = "bold"
    ),
    axis.line = element_line(
      color = "black"
    )
  )



metadata_anno_filter_sel_test |>
  dplyr::select(variant, celltype, disease_t) |>
  tidyr::unnest(cols = disease_t) |>
  dplyr::filter(pval < 0.01)

metadata_anno_filter_sel_test |>
  dplyr::select(variant, disease_t) |>
  tidyr::unnest(cols = disease_t) |>
  dplyr::filter(pval < 0.01) ->
select_v_by_disease



metadata_anno_filter_sel |>
  dplyr::filter(variant == select_v_by_disease$variant) |>
  dplyr::mutate(
    variant = factor(variant)
  ) |>
  ggplot(aes(
    x = dia,
    y = af,
  )) +
  # geom_violin() +
  geom_boxplot(
    aes(color = dia),
    width = 0.5,
    show.legend = T
  ) +
  # geom_point(position = position_jitter(width = 0.3)) +
  scale_y_continuous(
    limits = c(0, 1)
  ) +
  ggsci::scale_color_jama(
    name = "Disease type"
  ) +
  facet_grid(
    rows = vars(variant),
    cols = vars(celltype),
    switch = "y"
  ) +
  theme(
    panel.background = element_blank(),
    panel.grid = element_blank(),
    # axis.ticks = element_blank(),
    axis.title = element_blank(),
    axis.text = element_text(
      color = "black",
      # size = 18
    ),
    legend.position = "bottom",
    legend.background = element_blank(),
    legend.key = element_blank(),
    plot.title = element_text(
      size = 16,
      hjust = 0.5
    ),
    strip.background = element_rect(
      fill = NA,
      color = "black",
    ),
    strip.text = element_text(
      color = "black",
      size = 10,
      face = "bold"
    ),
    axis.line = element_line(
      color = "black"
    )
  ) ->
p_cell_dia
p_cell_dia
ggsave(
  filename = "Celltype disease.pdf",
  plo = p_cell_dia,
  device = "pdf",
  width = 15,
  height = 8,
  path = outdir
)
# footer ------------------------------------------------------------------

# future::plan(future::sequential)

# save image --------------------------------------------------------------
save.image(
  file = file.path(
    outdir,
    "03-merge-meta-variant.rda"
  )
)
load(file.path(
  outdir,
  "03-merge-meta-variant.rda"
))
