#!/usr/bin/env Rscript
# Metainfo ----------------------------------------------------------------
# @AUTHOR: Chun-Jie Liu
# @CONTACT: chunjie.sam.liu.at.gmail.com
# @DATE: 2026-02-23 12:30:01
# @DESCRIPTION: this script is used for ...

# Reproducibility ----------------------------------------------------------
set.seed(1)
# Library -----------------------------------------------------------------

suppressMessages({
  load_pkg(jutils)
})

# Args --------------------------------------------------------------------

# s: string, i: integer, f: float, !: boolean, @: array, %: list
VERSION = "v0.0.1"

GetoptLong.options(help_style = "two-column")

# default: default value specified here.

nthread = 8
GetoptLong(
  "nthread=i",
  "Number of threads to use",
  "verbose",
  "Enable verbose logging"
)


# Logger ------------------------------------------------------------------

log_layout(layout_glue_colors)

if (isTRUE(verbose)) {
  log_threshold(TRACE)
  log_info("Verbose mode enabled")
} else {
  log_threshold(INFO)
}


# Source ---------------------------------------------------------------------
source("high-res/00-colors.R")

# Load data ---------------------------------------------------------------
load_pkg(jutils)
dotenv()
suppressMessages({
  conflicted::conflict_prefer("filter", "dplyr")
})

variant_annotation <- import(
  path(Sys.getenv("OUTDIR")) /
    "VARIANT-ANNOTATION-TABLE-APOGEE2.xlsx"
)

cleandatadir <- path(Sys.getenv("CLEANDATADIR"))

cluster_variant <- import(
  path(Sys.getenv("OUTDIR")) / "ALLVARIANT-ALLSAMPLES-CLUSTERAF.xlsx"
)
#

outdir <- path(Sys.getenv("OUTDIR"))
outdirnotuse <- path(
  Sys.getenv("OUTDIRNOTUSE")
)
ALLVARIANTS <- import(
  outdir / "SAMPLE-VARIANT-CLASSIFICATION-CLUSTER-BULK-AF.xlsx"
) |>
  filter(variant_type %in% c("hete", "homo"))
METAFULL <- import(outdir / "SAMPLES-METADATA-FULL.xlsx")

METAFULL |>
  dplyr::filter(
    disease %in% c("Healthy", "Alzheimer's Disease")
  ) |>
  filter(Chemistry == "SC5P-PE") |>
  select(gseid, srrid, Chemistry, disease) |>
  mutate(
    disease = factor(disease, levels = c("Healthy", "Alzheimer's Disease"))
  ) -> admeta

ad_srrid <- admeta$srrid

admeta |>
  left_join(
    ALLVARIANTS,
    by = c("gseid", "srrid")
  ) |>
  select(
    -c(Chemistry, Haplogroup, Verbose_haplogroup)
  ) |>
  mutate(
    disease = factor(disease, levels = c("Healthy", "Alzheimer's Disease"))
  ) -> admeta_af

ad_forttest_ttest <- import(
  outdirnotuse / "AD" / "AD-variant-af-ttest.qs"
)


# Conn ---------------------------------------------------------------
conn <- db_conn(
  Sys.getenv("DUCKDB_PATH"),
  readonly = TRUE
)
tbl_ls(conn)
tbl_allvariants_cell <- dplyr::tbl(
  conn,
  "allvariants_cell"
) |>
  filter(srrid %in% ad_srrid) |>
  filter(variant_type %in% c("colorful", "black"))
# Function ----------------------------------------------------------------

# Main --------------------------------------------------------------------

ad_forttest_ttest |>
  select(-data) |>
  filter(variant %in% admeta_af$variant) |>
  unnest(t) |>
  dplyr::filter(p.value < 0.05) |>
  dplyr::mutate(
    plog10p = -log10(p.value),
    est = abs(estimate),
  ) |>
  dplyr::mutate(
    rank = plog10p * est,
  ) |>
  dplyr::arrange(
    desc(rank)
  ) |>
  dplyr::rename(
    ad = "Alzheimer's Disease",
  ) |>
  dplyr::filter(
    ad >= 5,
    Healthy >= 5
  ) -> top_variants


top_variants

library(ggh4x)
library(ggbeeswarm)
library(ggnewscale)

color_celltype_bulk <- c(
  # "Pseudo-bulk" = "red",
  color_celltype
)


top_variants |>
  dplyr::mutate(
    log10p = -log10(p.value),
  ) |>
  dplyr::mutate(
    celltype = gsub(celltype, pattern = "_", replacement = " "),
    celltype = ifelse(celltype == "bulk", "Pseudo-bulk", celltype),
  ) |>
  dplyr::mutate(
    celltype = factor(
      celltype,
      levels = names(color_celltype_bulk)
    ),
  ) |>
  dplyr::mutate(
    point_color = ifelse(
      p.value < 0.05 & abs(estimate) > 0.03,
      "red",
      "black"
    ),
  ) |>
  dplyr::mutate(
    label = ifelse(
      p.value < 0.05 & abs(estimate) > 0.03,
      variant,
      NA
    ),
  ) -> forplot_test


forplot_test |>
  ggplot(aes(
    x = estimate,
    y = log10p,
  )) +
  geom_point(
    aes(
      color = point_color,
    ),
  ) +
  scale_color_identity() +
  ggrepel::geom_text_repel(
    aes(
      label = label,
    ),
    size = 3,
    show.legend = FALSE,
    # nudge_x = 0.1,
    # nudge_y = 0.1,
    # segment.size = 0.5,
    segment.color = "black",
    max.overlaps = Inf
  ) +
  geom_hline(
    yintercept = -log10(0.05),
    linetype = 20,
    color = "red"
  ) +
  geom_vline(
    xintercept = 0,
    linetype = 20,
    color = "red"
  ) +
  ggh4x::facet_grid2(
    ~celltype,
    strip = ggh4x::strip_themed(
      background_x = ggh4x::elem_list_rect(
        fill = color_celltype_bulk,
        color = NA
      ),
      text_x = ggh4x::elem_list_text(
        colour = "white",
        face = c("bold")
      )
    )
  ) +
  theme(
    plot.margin = margin(t = 0.2, b = 0.1, l = 0.1, r = 0.2, unit = "cm"),
    # panel.background = element_blank(),
    panel.background = element_rect(
      fill = NA,
      color = "black",
      linewidth = 0.5
    ),
    panel.grid = element_blank(),
    axis.line.y.left = element_line(color = "black"),
    axis.line.x.bottom = element_line(color = "black"),
    axis.ticks.x = element_blank(),
    # axis.text.x = element_blank(),
    # axis.line.x = element_blank(),
    # axis.title.x = element_blank(),
    legend.position = "top",
    legend.key = element_blank(),
    axis.title.y = element_text(color = "black", size = 16),
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
    y = "-log10(p-value)",
    x = "Effect Size (AD - Healthy)",
  ) -> p_variant_boxplot_af_sc_ttest

ggsave(
  p_variant_boxplot_af_sc_ttest,
  filename = outdirnotuse / "AD" / "AD-variant-af-sc-ttest.pdf",
  width = 13,
  height = 3.5,
)

topvariants <- c(
  "263A>G",
  "16311T>C",
  "8794C>T",
  "1736A>G",
  "4824A>G",
  "663A>G",
  "8362T>G",
  "4175G>A",
  "5513G>A",
  "11704C>T",
  "1397T>A",
  "1397T>A",
  "1670A>G",
  "3173G>A",
  "3176A>T",
  "3178T>A"
)

ad_forttest_ttest |>
  select(
    variant,
    celltype,
    data
  ) |>
  filter(variant %in% topvariants) |>
  unnest(data) |>
  dplyr::mutate(
    celltype = gsub(celltype, pattern = "_", replacement = " "),
    celltype = ifelse(celltype == "bulk", "Pseudo-bulk", celltype),
  ) |>
  dplyr::mutate(
    celltype = factor(
      celltype,
      levels = names(color_celltype_bulk)
    ),
  ) -> forplot_


forplot_ |>
  ggplot(aes(x = disease)) +
  ggh4x::facet_grid2(
    variant ~ celltype,
    strip = ggh4x::strip_themed(
      background_x = ggh4x::elem_list_rect(
        fill = color_celltype_bulk,
        color = NA
      ),
      text_x = ggh4x::elem_list_text(
        colour = "white",
        face = c("bold")
      ),
      background_y = ggh4x::elem_list_rect(
        fill = "black",
        color = NA
      ),
      text_y = ggh4x::elem_list_text(
        colour = "white",
        face = c("bold")
      )
    )
  ) +
  geom_violin(
    aes(
      y = af,
      fill = disease
    ),
    alpha = 0.7,
    size = 1,
    color = NA
  ) +
  scale_fill_manual(
    values = color_disease,
    name = "Disease",
    # labels = c("AD", "Healthy")
  ) +
  ggbeeswarm::geom_quasirandom(
    aes(
      y = af,
      color = disease
    ),
    size = 1,
    dodge.width = .75,
    alpha = 1,
    show.legend = FALSE
  ) +
  scale_color_manual(
    values = color_disease,
    name = "Disease",
    # labels = c("AD", "Healthy")
  ) +
  ggsignif::geom_signif(
    aes(
      y = af,
    ),
    comparisons = list(
      c("Alzheimer's Disease", "Healthy")
    ),
    y_position = 0.8
  ) +
  theme(
    plot.margin = margin(t = 0.2, b = 0.1, l = 0.1, r = 0.2, unit = "cm"),
    # panel.background = element_blank(),
    panel.background = element_rect(
      fill = NA,
      color = "black",
      linewidth = 0.5
    ),
    panel.grid = element_blank(),
    axis.line.y.left = element_line(color = "black"),
    axis.line.x.bottom = element_line(color = "black"),
    axis.ticks.x = element_blank(),
    axis.text.x = element_blank(),
    # axis.line.x = element_blank(),
    axis.title.x = element_blank(),
    legend.position = "top",
    legend.key = element_blank(),
    axis.title.y = element_text(color = "black", size = 16),
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
    y = "Allele Frequency",
  ) -> p_variant_boxplot_af

ggsave(
  p_variant_boxplot_af,
  filename = outdirnotuse / "AD" / "AD-variant-af-boxplot.pdf",
  width = 16,
  height = 20,
  device = cairo_pdf
)

# Save  --------------------------------------------------------------

# Session info -------------------------------------------------------------
if (isTRUE(verbose)) {
  sessionInfo()
}
