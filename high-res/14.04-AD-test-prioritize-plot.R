#!/usr/bin/env Rscript
# Metainfo ----------------------------------------------------------------
# @AUTHOR: Chun-Jie Liu
# @CONTACT: chunjie.sam.liu.at.gmail.com
# @DATE: 2026-03-02 14:15:16
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

ad_srrid <- admeta_af$srrid |> unique()
ad_variant <- admeta_af$variant |> unique()


ad_variant_ttest_cluster <- import(
  outdirnotuse / "AD" / "AD-variant-af-ttest-cluster.qs"
)
ad_variant_ttest_cell <- import(
  outdirnotuse / "AD" / "AD-variant-af-ttest-cell.qs"
)


color_celltype_bulk <- color_celltype


# Conn ---------------------------------------------------------------

# Function ----------------------------------------------------------------

theme_ad_panel <- function() {
  theme(
    plot.margin = margin(t = 0.2, b = 0.1, l = 0.1, r = 0.2, unit = "cm"),
    panel.background = element_rect(
      fill = NA,
      color = "black",
      linewidth = 0.5
    ),
    panel.grid = element_blank(),
    axis.line.y.left = element_line(color = "black"),
    axis.line.x.bottom = element_line(color = "black"),
    legend.position = "top",
    legend.key = element_blank(),
    axis.title.y = element_text(color = "black", size = 16),
    axis.text.y = element_text(color = "black"),
    legend.text = element_text(size = 14, color = "black"),
    legend.title = element_text(size = 16, colour = "black"),
    strip.background = element_blank(),
    strip.text = element_text(size = 8, color = "black", face = "bold")
  )
}

strip_celltype <- function() {
  ggh4x::strip_themed(
    background_x = ggh4x::elem_list_rect(
      fill = color_celltype_bulk,
      color = NA
    ),
    text_x = ggh4x::elem_list_text(colour = "white", face = "bold")
  )
}


plot_ttest_scatter <- function(ttest_data, label = "cluster") {
  ttest_data |>
    select(-data) |>
    unnest(t) |>
    filter(p.value < 0.05) |>
    mutate(
      plog10p = -log10(p.value),
      est = abs(estimate),
      rank = plog10p * est,
    ) |>
    arrange(desc(rank)) |>
    rename(ad = "Alzheimer's Disease") |>
    filter(ad >= 5, Healthy >= 5) |>
    mutate(
      log10p = -log10(p.value),
      celltype = gsub(celltype, pattern = "_", replacement = " "),
      celltype = ifelse(celltype == "bulk", "Pseudo-bulk", celltype),
      celltype = factor(celltype, levels = names(color_celltype_bulk)),
      point_color = ifelse(
        p.value < 0.05 & abs(estimate) > 0.03,
        "red",
        "black"
      ),
      label = ifelse(p.value < 0.05 & abs(estimate) > 0.03, variant, NA),
    ) -> forplot

  forplot |>
    ggplot(aes(x = estimate, y = log10p)) +
    geom_point(aes(color = point_color)) +
    scale_color_identity() +
    ggrepel::geom_text_repel(
      aes(label = label),
      size = 3,
      show.legend = FALSE,
      segment.color = "black",
      max.overlaps = Inf
    ) +
    geom_hline(yintercept = -log10(0.05), linetype = 20, color = "red") +
    geom_vline(xintercept = 0, linetype = 20, color = "red") +
    ggh4x::facet_grid2(~celltype, strip = strip_celltype()) +
    theme_ad_panel() +
    theme(axis.ticks.x = element_blank()) +
    labs(
      title = paste("AD variant t-test —", label, "level"),
      y = "-log10(p-value)",
      x = "Effect Size (AD - Healthy)",
    )
}


# Main --------------------------------------------------------------------

# Save  --------------------------------------------------------------

# Session info -------------------------------------------------------------
if (isTRUE(verbose)) {
  sessionInfo()
}
