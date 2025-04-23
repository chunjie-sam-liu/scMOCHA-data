#!/usr/bin/env Rscript
# Metainfo ----------------------------------------------------------------
# @AUTHOR: Chun-Jie Liu
# @CONTACT: chunjie.sam.liu.at.gmail.com
# @DATE: 2025-04-22 14:48:15
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
  # mt_exons_df <- "/home/liuc9/github/scMOCHA/fasta/mt_exons.df.rds.gz"

  LENGTH <- 16569
  rCRS <- Biostrings::readDNAStringSet("/home/liuc9/github/scMOCHA-data/config/rCRS.MT.fasta")
  gtf_gene_df <- readr::read_rds("/home/liuc9/github/scMOCHA-data/config/mtdna_genes_dloop.rds.gz")


  library(gggenes)
  ggplot(gtf_gene_df, aes(xmin = start, xmax = end, y = seqnames)) +
    # geom_gene_arrow() +
    geom_gene_arrow(
      aes(
        fill = TYPE
      ),
      arrowhead_height = unit(3, "mm"), arrowhead_width = unit(1, "mm")
    ) +
    scale_fill_brewer(
      palette = "Set1",
      name = "Gene type",
      labels = c("D-Loop", "MT rRNA", "MT tRNA", "Protein coding")
    ) +
    ggrepel::geom_text_repel(
      aes(
        x = (start + end) / 2,
        label = gene_name,
        color = TYPE
      ),
      # fill = "white",
      # nudge_x =1,
      # nudge_y = -0.1,
      size = 3,
      show.legend = F,
      max.overlaps = Inf,
    ) +
    scale_color_brewer(palette = "Set1") +
    scale_x_continuous(
      limits = c(0, LENGTH),
      breaks = c(seq(0, LENGTH, 1000), LENGTH),
      labels = c(seq(0, LENGTH, 1000), LENGTH),
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
      axis.ticks.x = element_line(color = "black"),
      axis.line.x = element_line(color = "black"),
      axis.text.x = element_text(
        vjust = -1,
      ),
    )
}
# load data ---------------------------------------------------------------
basedir <- "/home/liuc9/github/scMOCHA-data/data"
foundation_out <- file.path(basedir, "scfoundation/out")
outdir <- "/home/liuc9/github/scMOCHA-data/stats/stats/zzz"


gseids <- c(
  "GSE155673",
  "GSE157344",
  "GSE149689",
  "GSE171555",
  "GSE155223",
  "GSE163668",
  "GSE175524",
  "GSE206283",
  "GSE226598",
  "GSE261140",
  "GSE279945",
  "GSE214865",
  "GSE220189",
  "GSE233844",
  "GSE175499",
  "GSE149313",
  "GSE154386",
  "GSE159117",
  "GSE188632",
  "GSE166992",
  "GSE162117",
  "GSE226602",
  "GSE161354",
  "GSE235050",
  "GSE181279",
  # scfoundation2
  "GSE143353",
  "GSE148215",
  "GSE163314",
  "GSE163633",
  "GSE164690",
  "GSE167825",
  "GSE174125",
  "GSE184703",
  "GSE153421",
  # "GSE168453",
  "GSE147794"
)

# not PBMC
gseids_not_pbmc <- c(
  "GSE168453", # Pool samples together, not individual samples
  "GSE163668", # Some are pooled samples, not individual samples. Whole blood, not PBMC, including others like Red blood cells, Platelets and Plasma
  "GSE163633", # Not all PBMC, some are Mucosa-derived T cells / Myloid cells, Squamous Cell Carcinoma(SCC) -derived T cells / Myloid Cells
  "GSE157344", # It’s not PBMC, but the blood or bronchoalveloar lavage samples
  "GSE163314", # Half of the samples are from Colon, not PBMC
  "GSE164690", # Most of the samples are from tumors (Head and neck squamous cell carcinoma (HNSCC)), not from PBMC
  "GSE148215" # The samples are human embryonic cells, not PBMC cells
)
# enrich cells
gseids_enrich_cells <- c(
  "GSE167825", # CD8T cell enriched
  "GSE175524", # B cell enriched
  "GSE261140" # CD8T cell enriched
)

gseids_tobe_excluded <- c(
  gseids_not_pbmc
)

gse_dataset_metadata_full <- readr::read_rds(
  file.path(foundation_out, "gse_dataset_metadata_full.rds")
) |>
  dplyr::filter(
    !gseid %in% gseids_tobe_excluded
  )

pcc <- readr::read_tsv(file = "https://raw.githubusercontent.com/chunjie-sam-liu/chunjie-sam-liu.life/master/public/data/pcc.tsv") |>
  dplyr::arrange(cancer_types)


# thegseid <- "GSE168453"
# body --------------------------------------------------------------------

tibble::tibble(
  gseid = gseids
) |>
  dplyr::filter(
    !gseid %in% gseids_tobe_excluded
  ) |>
  dplyr::mutate(
    anno = parallel::mclapply(
      X = gseid,
      FUN = function(.gseid) {
        log_info("Loading {.gseid}... ({which(gseids == .gseid)}/{length(gseids)})")
        .anno <- readr::read_rds(
          file.path(basedir, .gseid, "out", glue::glue("{.gseid}.scmocha.out.rds.gz"))
        )
        log_success("Loaded {.gseid}! ({which(gseids == .gseid)}/{length(gseids)})")
        return(.anno)
      },
      mc.cores = 10
    )
  ) ->
gse_data_loaded

gse_data_loaded |>
  tidyr::unnest(cols = anno) ->
gse_data

# body --------------------------------------------------------------------

gse_data |>
  dplyr::select(gseid, srrid, chemistry, anno, hetero, haplo_violin, somatic_variant) ->
for_hetero


gse_data |>
  dplyr::select(gseid, srrid, chemistry, anno, hetero, haplo_variant, haplo_violin, somatic_variant, celltype_ratio) |>
  dplyr::left_join(
    gse_dataset_metadata_full |> dplyr::select(-gseid),
    by = c("srrid" = "srrid")
  ) |>
  dplyr::mutate(
    disease = dplyr::case_when(
      disease %in% c("Alzheimer's Disease", "Healthy", "COVID-19", "Unknown") ~ disease,
      TRUE ~ "Other"
    )
  ) |>
  dplyr::mutate(
    disease = factor(
      disease,
      levels = c(
        "Alzheimer's Disease",
        "COVID-19",
        "Healthy",
        "Unknown",
        "Other"
      )
    )
  ) |>
  dplyr::mutate(
    Chemistry = factor(
      Chemistry,
      levels = c("SC3Pv2", "SC3Pv3", "SC5P-R2", "SC5P-PE")
    )
  ) |>
  dplyr::arrange(disease, Chemistry) ->
gse_data_haplo_variant



# ! all variants --------------------------------------------------------------------

gse_data_haplo_variant$somatic_variant |>
  purrr::map(
    .f = \(.x) {
      .x$somatic
    }
  ) |>
  purrr::reduce(
    union
  ) ->
somatic_variants

gse_data_haplo_variant |>
  dplyr::select(gseid, srrid, chemistry, haplo_variant) |>
  tidyr::unnest(cols = haplo_variant) ->
all_variants

all_variants |>
  dplyr::select(Position, variant, aachange, Disease) |>
  dplyr::distinct() |>
  dplyr::arrange(Position) ->
variant_type


all_variants |>
  dplyr::count(variant) |>
  dplyr::left_join(
    variant_type,
    by = "variant"
  ) |>
  dplyr::mutate(
    issomatic = ifelse(variant %in% somatic_variants, "somatic", "other"),
  ) ->
variant_count



variant_count |>
  # dplyr::filter(color == "black") |>
  dplyr::filter(
    issomatic == "somatic"
  ) |>
  dplyr::arrange(
    desc(n)
  ) |>
  ggplot(aes(
    x = Position,
    y = n
  )) +
  geom_bar(stat = "identity", color = "red") +
  geom_text(
    data = variant_count |>
      # dplyr::filter(color == "black") |>
      dplyr::filter(
        issomatic == "somatic"
      ) |>
      dplyr::arrange(
        desc(n)
      ) |>
      head(5) |>
      dplyr::mutate(
        label = glue::glue("{variant}({aachange})")
      ),
    aes(
      label = label,
    ),
    color = "black",
    size = 3,
    vjust = -0.5,
    hjust = 0.5,
    show.legend = FALSE
  ) +
  scale_x_continuous(
    expand = expansion(mult = c(0.01, 0)),
    limits = c(1, 17000),
    breaks = seq(0, 17000, 2000),
    labels = seq(0, 17000, 2000)
  ) +
  scale_y_continuous(
    expand = c(0.01, 0),
    label = scales::label_number()
  ) +
  # scale_fill_identity(
  #   name = "Sample"
  # ) +
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
    # legend.position = c(0.8, 0.5),
    legend.position = "none",
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
  labs(
    y = "# of Samples",
  ) ->
variant_count_plot
variant_count_plot
p_mtdna <- fn_plot_mtdna()

wrap_plots(
  variant_count_plot,
  p_mtdna,
  ncol = 1,
  heights = c(15, 1)
)

ggsave(
  filename = file.path(outdir, "somatic_variant_distribution_bar_samples.pdf"),
  plot = wrap_plots(
    variant_count_plot,
    p_mtdna,
    ncol = 1,
    heights = c(15, 1)
  ),
  width = 23,
  height = 12,
  dpi = 300
)

# body --------------------------------------------------------------------


# footer ------------------------------------------------------------------

# future: :plan(future: :sequential)

# save image --------------------------------------------------------------
