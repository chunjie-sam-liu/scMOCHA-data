#!/usr/bin/env Rscript
# Metainfo ----------------------------------------------------------------
# @AUTHOR: Chun-Jie Liu
# @CONTACT: chunjie.sam.liu.at.gmail.com
# @DATE: 2026-03-23 00:00:00
# @DESCRIPTION: Prioritized cluster-level COVID-19 variant t-test plots

# Reproducibility ----------------------------------------------------------
set.seed(1)

# Library ------------------------------------------------------------------

suppressMessages({
  load_pkg(jutils)
})

# Args ---------------------------------------------------------------------

VERSION = "v0.0.1"

GetoptLong.options(help_style = "two-column")

nthread = 8
GetoptLong(
  "nthread=i",
  "Number of threads to use",
  "verbose",
  "Enable verbose logging"
)

# Logger -------------------------------------------------------------------

log_layout(layout_glue_colors)

if (isTRUE(verbose)) {
  log_threshold(TRACE)
  log_info("Verbose mode enabled")
} else {
  log_threshold(INFO)
}

# Load data ----------------------------------------------------------------

load_pkg(jutils)
dotenv()
source(path(Sys.getenv("HIGHRESDIR"), "00-colors.R"))

variant_annotation <- import(
  path(Sys.getenv("OUTDIR")) / "VARIANT-ANNOTATION-TABLE-APOGEE2.xlsx"
)

outdirnotuse <- path(Sys.getenv("OUTDIRNOTUSE"))
covid_outdir <- outdirnotuse / "COVID19"

covid_variant_ttest_cluster <- import(
  covid_outdir / "COVID19-variant-af-ttest-cluster.qs2"
)

# Packages -----------------------------------------------------------------

suppressMessages({
  load_pkg(ggh4x, ggbeeswarm, ggnewscale, ggsignif)
})

# Constants ----------------------------------------------------------------

color_celltype_bulk <- c("Pseudo-bulk" = "red", color_celltype)

# Functions ----------------------------------------------------------------

theme_covid_panel <- function() {
  theme(
    plot.margin = margin(t = 0.2, b = 0.1, l = 0.1, r = 0.2, unit = "cm"),
    panel.background = element_rect(fill = NA, color = "black", linewidth = 0.5),
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

normalize_celltype <- function(x) {
  x <- gsub(pattern = "_", replacement = " ", x = x)
  ifelse(tolower(x) == "bulk", "Pseudo-bulk", x)
}

plot_ttest_dotmatrix <- function(ttest_data, top_n = 20, label = "cluster") {
  processed <- ttest_data |>
    select(-data) |>
    unnest(t) |>
    filter(!is.na(p.value)) |>
    filter(Healthy >= 5, `COVID-19` >= 5) |>
    mutate(
      celltype = normalize_celltype(celltype),
      celltype = factor(celltype, levels = names(color_celltype_bulk))
    )

  all_ranked <- processed |>
    filter(p.value < 0.05) |>
    group_by(variant) |>
    summarise(global_rank = max(-log10(pmax(p.value, 1e-300))), .groups = "drop") |>
    arrange(desc(global_rank))

  if (nrow(all_ranked) == 0) {
    log_warn("No significant variants found; falling back to overall ranking")
    all_ranked <- processed |>
      group_by(variant) |>
      summarise(global_rank = max(-log10(pmax(p.value, 1e-300))), .groups = "drop") |>
      arrange(desc(global_rank))
  }

  all_variant_order <- all_ranked$variant
  top_variant_order <- head(all_variant_order, top_n)
  n_total <- length(all_variant_order)

  anno_subset <- variant_annotation |>
    select(variant, Locus, prediction_class, Disease) |>
    distinct()

  make_label <- function(v, is_top) {
    idx <- match(v, anno_subset$variant)
    locus <- if (!is.na(idx)) anno_subset$Locus[idx] else NA_character_
    base <- if (!is.na(locus) && locus != "") paste0(v, " [", locus, "]") else v
    if (!is_top) {
      return(v)
    }
    pred <- if (!is.na(idx)) anno_subset$prediction_class[idx] else NA_character_
    dis <- if (!is.na(idx)) anno_subset$Disease[idx] else NA_character_
    parts <- c(
      if (!is.na(pred) && pred != "") pred else NULL,
      if (!is.na(dis) && dis != "") dis else NULL
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

  forplot <- processed |>
    filter(variant %in% all_variant_order) |>
    mutate(
      variant_label = all_labels[match(variant, all_variant_order)],
      variant_label = factor(variant_label, levels = rev(all_labels)),
      dot_size = ifelse(p.value < 0.05, -log10(pmax(p.value, 1e-300)), NA_real_),
      effect = estimate2 - estimate1
    )

  top_y_pos <- seq(max(1, n_total - length(top_variant_order) + 1), n_total)

  p <- ggplot(forplot, aes(x = celltype, y = variant_label)) +
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
      aes(size = dot_size, color = effect),
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
      name = "Effect size\n(COVID-19 - Healthy)"
    ) +
    theme_covid_panel() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
      axis.text.y = element_text(size = 7, color = "grey40"),
      axis.title = element_blank(),
      legend.position = "right"
    ) +
    labs(
      title = paste("COVID-19 variants -", label, "level"),
      subtitle = paste0(
        "Top ",
        length(top_variant_order),
        " highlighted (yellow); all significant shown; dot size = -log10(p), color = effect direction"
      )
    )

  list(plot = p, topvariants = top_variant_order)
}

plot_af_violin <- function(ttest_data, variants, label = "cluster", anno_title = NULL) {
  ttest_data |>
    select(variant, celltype, data) |>
    filter(variant %in% variants) |>
    unnest(data) |>
    mutate(
      celltype = normalize_celltype(celltype),
      celltype = factor(celltype, levels = names(color_celltype_bulk))
    ) |>
    ggplot(aes(x = disease)) +
    ggh4x::facet_grid2(
      variant ~ celltype,
      strip = ggh4x::strip_themed(
        background_x = ggh4x::elem_list_rect(fill = color_celltype_bulk, color = NA),
        text_x = ggh4x::elem_list_text(colour = "white", face = "bold"),
        background_y = ggh4x::elem_list_rect(fill = "black", color = NA),
        text_y = ggh4x::elem_list_text(colour = "white", face = "bold")
      )
    ) +
    geom_violin(aes(y = af, fill = disease), alpha = 0.7, color = NA) +
    scale_y_continuous(limits = c(0, 1)) +
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
      comparisons = list(c("COVID-19", "Healthy")),
      y_position = 0.8
    ) +
    theme_covid_panel() +
    theme(
      axis.ticks.x = element_blank(),
      axis.text.x = element_blank(),
      axis.title.x = element_blank()
    ) +
    labs(
      title = if (!is.null(anno_title)) {
        paste0("COVID-19 variant AF distribution - ", label, "\n", anno_title)
      } else {
        paste("COVID-19 variant AF distribution -", label, "level")
      },
      y = "Allele Frequency"
    )
}

variant_annotation_title <- function(variant_id) {
  anno <- variant_annotation |>
    filter(variant == variant_id) |>
    slice(1)

  locus <- if (nrow(anno) > 0 && !is.na(anno$Locus) && anno$Locus != "") {
    anno$Locus
  } else {
    NA_character_
  }
  pred <- if (nrow(anno) > 0 && !is.na(anno$prediction_class) && anno$prediction_class != "") {
    anno$prediction_class
  } else {
    NA_character_
  }
  dis <- if (nrow(anno) > 0 && !is.na(anno$Disease) && anno$Disease != "") {
    anno$Disease
  } else {
    NA_character_
  }

  var_label <- paste0(variant_id, if (!is.na(locus)) paste0(" [", locus, "]") else "")
  anno_parts <- na.omit(c(pred, dis))
  anno_title <- if (length(anno_parts) > 0) paste(anno_parts, collapse = " | ") else NULL

  list(label = var_label, anno_title = anno_title)
}

# Main ---------------------------------------------------------------------

dotmatrix <- plot_ttest_dotmatrix(
  ttest_data = covid_variant_ttest_cluster,
  top_n = 20,
  label = "cluster"
)

ggsave(
  filename = covid_outdir / "COVID19-variant-prioritized-cluster.pdf",
  plot = dotmatrix$plot,
  width = 14,
  height = max(8, length(dotmatrix$topvariants) * 0.7),
  device = cairo_pdf
)

p_list <- dotmatrix$topvariants |>
  map(
    \(.x) {
      anno <- variant_annotation_title(.x)
      plot_af_violin(
        ttest_data = covid_variant_ttest_cluster,
        variants = .x,
        label = anno$label,
        anno_title = anno$anno_title
      )
    }
  )

pdf(
  covid_outdir / "COVID19-variant-af-violin-cluster-individual.pdf",
  width = 8,
  height = 4
)
invisible(p_list |> walk(print))
invisible(dev.off())

if (isTRUE(verbose)) {
  sessionInfo()
}
