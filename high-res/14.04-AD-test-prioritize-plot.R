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


color_celltype_bulk <- c(color_celltype, "Pseudo-bulk" = "grey30")


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


plot_ttest_dotmatrix <- function(ttest_data, top_n = 20, label = "cluster") {
  # 1. Unnest, filter
  processed <- ttest_data |>
    select(-data) |>
    unnest(t) |>
    filter(!is.na(p.value)) |>
    rename(ad = "Alzheimer's Disease") |>
    filter(ad >= 5, Healthy >= 5) |>
    mutate(
      celltype = gsub(celltype, pattern = "_", replacement = " "),
      celltype = ifelse(celltype == "bulk", "Pseudo-bulk", celltype),
      celltype = factor(celltype, levels = names(color_celltype_bulk))
    )

  # 2. Rank ALL significant variants; identify top N
  all_ranked <- processed |>
    filter(p.value < 0.05) |>
    group_by(variant) |>
    summarise(
      global_rank = max(-log10(p.value), na.rm = TRUE),
      .groups = "drop"
    ) |>
    arrange(desc(global_rank))

  all_variant_order <- all_ranked$variant # best → worst
  top_variant_order <- head(all_variant_order, top_n)
  n_total <- length(all_variant_order)

  # 3. Build y-axis labels: top N get full annotation, rest get bare variant ID
  anno_subset <- variant_annotation |>
    select(variant, Locus, prediction_class, Disease) |>
    distinct()

  make_label <- function(v, is_top) {
    idx <- match(v, anno_subset$variant)
    locus <- if (!is.na(idx)) anno_subset$Locus[idx] else NA_character_
    base <- if (!is.na(locus) & locus != "") paste0(v, " [", locus, "]") else v
    if (!is_top) {
      return(v)
    }
    pred <- if (!is.na(idx)) {
      anno_subset$prediction_class[idx]
    } else {
      NA_character_
    }
    dis <- if (!is.na(idx)) anno_subset$Disease[idx] else NA_character_
    parts <- c(
      if (!is.na(pred) & pred != "") pred else NULL,
      if (!is.na(dis) & dis != "") dis else NULL
    )
    if (length(parts) > 0) {
      paste0(base, "\n", paste(parts, collapse = " | "))
    } else {
      base
    }
  }

  all_labels <- mapply(
    make_label,
    all_variant_order,
    all_variant_order %in% top_variant_order,
    USE.NAMES = FALSE
  )

  # 4. Prepare forplot (all significant variants)
  forplot <- processed |>
    filter(variant %in% all_variant_order) |>
    mutate(
      variant_label = all_labels[match(variant, all_variant_order)],
      variant_label = factor(variant_label, levels = rev(all_labels)),
      dot_size = ifelse(p.value < 0.05, -log10(p.value), NA_real_)
    )

  # 5. Highlight positions for top N rows
  # rev(all_labels) puts best-ranked at top = highest numeric level indices
  top_y_pos <- seq(n_total - top_n + 1, n_total)

  # 6. Draw
  ggplot(forplot, aes(x = celltype, y = variant_label)) +
    annotate(
      "rect",
      xmin = -Inf,
      xmax = Inf,
      ymin = top_y_pos - 0.5,
      ymax = top_y_pos + 0.5,
      fill = "#FFF3CD",
      alpha = 0.8
    ) +
    geom_point(
      data = dplyr::filter(forplot, is.na(dot_size)),
      size = 0.8,
      color = "grey85",
      shape = 16
    ) +
    geom_point(
      data = dplyr::filter(forplot, !is.na(dot_size)),
      aes(size = dot_size, color = estimate),
      shape = 16
    ) +
    scale_size_continuous(
      range = c(1, 8),
      name = "-log10(p-value)",
      guide = guide_legend(override.aes = list(color = "grey40"))
    ) +
    scale_color_gradient2(
      low = "steelblue3",
      mid = "white",
      high = "firebrick3",
      midpoint = 0,
      name = "Effect size\n(AD \u2212 Healthy)"
    ) +
    theme_ad_panel() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
      axis.text.y = element_text(size = 7, color = "grey40"),
      axis.title = element_blank(),
      legend.position = "right"
    ) +
    labs(
      title = paste("AD variants \u2014", label, "level"),
      subtitle = paste0(
        "Top ",
        top_n,
        " highlighted (yellow); all significant shown; ",
        "dot size = -log10(p), color = effect direction"
      )
    )
}


plot_ttest_scatter <- function(ttest_data, label = "cluster") {
  ttest_data |>
    select(-data) |>
    # filter(variant %in% ad_variant) |>
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
        p.value < 0.05 & abs(estimate) > 0.05 | p.value < 10^-2.5,
        "red",
        "black"
      ),
      label = ifelse(
        p.value < 0.05 & abs(estimate) > 0.05 | p.value < 10^-2.5,
        variant,
        NA
      ),
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
    ) -> p
  forplot |> filter(point_color == "red") |> pull(variant) -> topvariants
  list(
    p = p,
    topvariants = topvariants
  )
}


plot_af_violin <- function(
  ttest_data,
  variants,
  label = "cluster",
  anno_title = NULL
) {
  ttest_data |>
    select(variant, celltype, data) |>
    filter(variant %in% variants) |>
    unnest(data) |>
    mutate(
      celltype = gsub(celltype, pattern = "_", replacement = " "),
      celltype = ifelse(celltype == "bulk", "Pseudo-bulk", celltype),
      celltype = factor(celltype, levels = names(color_celltype_bulk)),
    ) -> forplot

  forplot |>
    ggplot(aes(x = disease)) +
    ggh4x::facet_grid2(
      variant ~ celltype,
      strip = ggh4x::strip_themed(
        background_x = ggh4x::elem_list_rect(
          fill = color_celltype_bulk,
          color = NA
        ),
        text_x = ggh4x::elem_list_text(colour = "white", face = "bold"),
        background_y = ggh4x::elem_list_rect(fill = "black", color = NA),
        text_y = ggh4x::elem_list_text(colour = "white", face = "bold")
      )
    ) +
    geom_violin(
      aes(y = af, fill = disease),
      alpha = 0.7,
      color = NA
    ) +
    scale_fill_manual(values = color_disease, name = "Disease") +
    ggbeeswarm::geom_quasirandom(
      aes(y = af, color = disease),
      size = 1,
      dodge.width = 0.75,
      alpha = 1,
      show.legend = FALSE
    ) +
    scale_color_manual(values = color_disease, name = "Disease") +
    ggsignif::geom_signif(
      aes(y = af),
      comparisons = list(c("Alzheimer's Disease", "Healthy")),
      y_position = 0.8
    ) +
    theme_ad_panel() +
    theme(
      axis.ticks.x = element_blank(),
      axis.text.x = element_blank(),
      axis.title.x = element_blank(),
    ) +
    labs(
      title = if (!is.null(anno_title)) {
        paste0("AD variant AF distribution — ", label, "\n", anno_title)
      } else {
        paste("AD variant AF distribution —", label, "level")
      },
      y = "Allele Frequency",
    )
}


# Main --------------------------------------------------------------------

suppressMessages({
  load_pkg(ggh4x)
})
p_scatter_cluster <- plot_ttest_scatter(
  ttest_data = ad_variant_ttest_cluster,
  label = "cluster"
)

ggsave(
  p_scatter_cluster$p,
  filename = outdirnotuse / "AD" / "AD-variant-af-ttest-cluster.pdf",
  width = 13,
  height = 3.5
)

p_violin_cluster <- plot_af_violin(
  ttest_data = ad_variant_ttest_cluster,
  variants = p_scatter_cluster$topvariants,
  label = "cluster"
)

ggsave(
  p_violin_cluster,
  filename = outdirnotuse /
    "AD" /
    "AD-variant-af-violin-cluster.pdf",
  width = 16,
  height = 20,
  device = cairo_pdf
)


p_scatter_cluster$topvariants |>
  map(
    ~ {
      anno <- variant_annotation |>
        filter(variant == .x) |>
        slice(1)
      locus <- if (nrow(anno) > 0 && !is.na(anno$Locus) && anno$Locus != "") {
        anno$Locus
      } else {
        NA_character_
      }
      pred <- if (
        nrow(anno) > 0 &&
          !is.na(anno$prediction_class) &&
          anno$prediction_class != ""
      ) {
        anno$prediction_class
      } else {
        NA_character_
      }
      dis <- if (nrow(anno) > 0 && !is.na(anno$Disease) && anno$Disease != "") {
        anno$Disease
      } else {
        NA_character_
      }

      var_label <- paste0(
        .x,
        if (!is.na(locus)) paste0(" [", locus, "]") else ""
      )
      anno_parts <- na.omit(c(pred, dis))
      anno_title <- if (length(anno_parts) > 0) {
        paste(anno_parts, collapse = " | ")
      } else {
        NULL
      }

      plot_af_violin(
        ttest_data = ad_variant_ttest_cluster,
        variants = .x,
        label = var_label,
        anno_title = anno_title
      )
    }
  ) -> p_list


{
  pdf(
    outdirnotuse / "AD" / "AD-variant-af-violin-cluster-individual.pdf",
    width = 8,
    height = 4
  )
  p_list |> map(print)
  dev.off()
}

# Save  --------------------------------------------------------------

# Session info -------------------------------------------------------------
if (isTRUE(verbose)) {
  sessionInfo()
}
