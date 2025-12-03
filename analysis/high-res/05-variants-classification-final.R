#!/usr/bin/env Rscript
# Metainfo ----------------------------------------------------------------
# @AUTHOR: Chun-Jie Liu
# @CONTACT: chunjie.sam.liu.at.gmail.com
# @DATE: 2025-12-02 16:23:08
# @DESCRIPTION: this script is used for ...

# Library -----------------------------------------------------------------

load_pkg(jutils)

# args --------------------------------------------------------------------

# s: string, i: integer, f: float, !: boolean, @: array, %: list
GetoptLong.options(help_style = "two-column")
VERSION = "v0.0.1"

# default: default value specified here.

verbose = TRUE

GetoptLong("verbose!", "print messages")


logger::log_threshold(logger::TRACE)
logger::log_layout(logger::layout_glue_colors)

# header ------------------------------------------------------------------

# load data ---------------------------------------------------------------
gse_data_variant_classification_clusteraf_bulkaf <- import(
  "/home/liuc9/github/scMOCHA-data/analysis/high-res/MANUSCRIPTFIGURES/SAMPLE-VARIANT-CLASSIFICATION-CLUSTER-BULK-AF.xlsx"
)
# load conn ---------------------------------------------------------------

# src ---------------------------------------------------------------------

# function ----------------------------------------------------------------
#' Very important function
fn_plot_cell_af_depth_forplot <- function(thevariant, thesrrid) {
  source("analysis/00-colors.R")
  conn_all_hetero_af <- DBI::dbConnect(
    duckdb::duckdb(),
    "/home/liuc9/github/scMOCHA-data/analysis/zzz/clean-data/all_hetero_af.cell.duckdb.1.2.1",
    read_only = TRUE
  )

  colorcode <- setNames(names(color_variantcell), color_variantcell)

  dplyr::tbl(
    conn_all_hetero_af,
    "allvariants_cell"
  ) |>
    dplyr::filter(
      srrid == thesrrid,
      variant == thevariant
    ) |>
    dplyr::collect() |>
    dplyr::mutate(
      variant_type = dplyr::case_match(
        variant_type,
        "colorful" ~ "red",
        "black" ~ "darkblue",
        "white" ~ "white",
        "grey" ~ "gray",
        NA ~ "white"
      )
    ) |>
    dplyr::mutate(
      variant_type = factor(
        variant_type,
        levels = color_variantcell
      )
    ) |>
    dplyr::arrange(
      variant_type,
      -af
    ) -> forplot_

  forplot_ |>
    dplyr::mutate(
      barcode = factor(
        barcode,
        levels = forplot_$barcode
      )
    ) |>
    dplyr::mutate(
      celltype = gsub(
        "_",
        " ",
        celltype
      )
    ) |>
    dplyr::mutate(
      celltype = factor(
        celltype,
        names(color_celltype)
      )
    ) |>
    dplyr::mutate(
      af = ifelse(
        af < 0.01,
        NA_real_,
        af
      )
    ) |>
    # dplyr::mutate(
    #   variant_type = as.character(variant_type),
    # ) |>
    dplyr::mutate(
      depth = log2(depth + 1) # log2 transform to reduce skewness
    ) |>
    dplyr::mutate(
      cellvarianttype = colorcode[variant_type]
    ) |>
    dplyr::mutate(
      cellvarianttype = factor(
        cellvarianttype,
        levels = colorcode
      )
    ) -> forplot
  DBI::dbDisconnect(conn_all_hetero_af, shutdown = TRUE)
  forplot
}
fn_plot_cell_af_somatic_variant_af <- function(forplot, thetheme) {
  fn_xy_breaks_limits(forplot$af, step = 0.2) -> .ybl

  forplot |>
    ggplot(aes(
      x = barcode,
      y = af,
      fill = af
    )) +
    geom_col() +
    geom_hline(
      aes(yintercept = 0.05),
      linetype = 20,
      color = "red"
    ) +
    geom_hline(
      aes(yintercept = 0.1),
      linetype = 21,
      color = "black"
    ) +
    scale_fill_gradient2(
      name = "Allele Frequency",
      high = "#FDE725FF",
      mid = "#21908CFF",
      low = "#440154FF"
    ) +
    scale_y_continuous(
      limits = .ybl$limits,
      breaks = c(.ybl$breaks, 0.05, 0.1) |> unique() |> sort(),
      labels = \(b) {
        dplyr::case_when(
          b == 0.1 ~ "gnomAD cutoff 10%",
          b == 0.05 ~ "our cutoff 5%",
          TRUE ~ scales::percent_format(accuracy = 1)(b)
        )
      },
      expand = expansion(mult = c(0, 0)),
    ) +
    theme(
      panel.background = element_blank(),
      panel.grid = element_blank(),
      axis.line = element_line(colour = "black"),
      axis.ticks.x = element_blank(),
      axis.text.x = element_blank(),
      axis.title.x = element_blank(),
      plot.margin = margin(t = 0, b = 0, unit = "cm"),
    ) +
    labs(
      y = "Allele Frequency",
    )
}
fn_plot_cell_af_somatic_variant_celltype <- function(forplot, thetheme) {
  forplot |>
    ggplot(aes(
      x = barcode,
      y = 1,
      fill = celltype
    )) +
    geom_col() +
    scale_fill_manual(
      name = "Cell Type",
      values = color_celltype,
    ) +
    scale_y_continuous(expand = expansion(mult = c(0, 0)), ) +
    thetheme +
    # theme(panel.background = element_rect(color = "red")) +
    labs(
      y = "Cell Type",
    )
}
fn_plot_cell_af_somatic_variant_depth <- function(forplot, thetheme) {
  fn_xy_breaks_limits(forplot$depth, step = 1) -> .ybl_depth
  forplot |>
    ggplot(aes(
      x = barcode,
      y = depth,
      fill = depth
    )) +
    geom_col() +
    geom_hline(
      aes(yintercept = log2(10 + 1)),
      linetype = 21,
      color = "black"
    ) +
    scale_fill_gradient(
      name = "log2(depth + 1)",
      high = "gold",
      low = "white"
    ) +
    scale_y_continuous(
      limits = .ybl_depth$limits,
      breaks = c(.ybl_depth$breaks, log2(10 + 1)) |> unique() |> sort(),
      labels = \(b) {
        dplyr::case_when(
          b == log2(10 + 1) ~ "cutoff 10",
          TRUE ~ scales::label_number()(b)
        )
      },
      expand = expansion(mult = c(0, 0)),
    ) +
    theme(
      panel.background = element_blank(),
      panel.grid = element_blank(),
      axis.line = element_line(colour = "black"),
      axis.ticks.x = element_blank(),
      axis.text.x = element_blank(),
      axis.title.x = element_blank(),
      plot.margin = margin(t = 0, b = 0, unit = "cm"),
    ) +
    labs(
      y = "Log2(Depth + 1)",
    )
}
fn_plot_cell_af_somatic_variant_cell <- function(forplot, thetheme) {
  forplot |>
    dplyr::mutate(
      variant_type = as.character(variant_type),
    ) |>
    ggplot(aes(
      x = barcode,
      y = 1,
      fill = variant_type
    )) +
    geom_col() +
    scale_y_continuous(
      expand = expansion(mult = c(0, 0)),
    ) +
    scale_fill_identity(
      guide = "legend",
      name = "Variant cell",
      breaks = c("red", "darkblue", "gray", "white"),
      labels = c(
        "Heteroplasmy",
        "Sufficient reads",
        "No sufficient reads",
        "No reads"
      )
    ) +
    thetheme +
    labs(
      y = "Variant cells",
    )
}
fn_plot_cell_af_cellvarianttype <- function(forplot) {
  source("analysis/00-colors.R")
  colorcode <- setNames(names(color_variantcell), color_variantcell)

  forplot |>
    dplyr::mutate(
      color = as.character(variant_type),
    ) |>
    dplyr::count(color) |>
    dplyr::mutate(
      varianttype = colorcode[color]
    ) |>
    dplyr::mutate(
      varianttype = factor(
        varianttype,
        levels = colorcode
      )
    ) |>
    dplyr::arrange(varianttype) -> cellvarianttype
}
fn_plot_cell_af_somatic_variant <- function(forplot_) {
  source("analysis/00-colors.R")

  colorcode <- setNames(names(color_variantcell), color_variantcell)
  thetheme <- theme(
    panel.background = element_blank(),
    panel.grid = element_blank(),
    axis.ticks = element_blank(),
    axis.text = element_blank(),
    axis.title.x = element_blank(),
    plot.margin = margin(t = 0, b = 0, unit = "cm"),
  )

  forplot_ |>
    dplyr::mutate(
      barcode = factor(
        barcode,
        levels = forplot_$barcode
      )
    ) |>
    dplyr::mutate(
      celltype = gsub(
        "_",
        " ",
        celltype
      )
    ) |>
    dplyr::mutate(
      celltype = factor(
        celltype,
        names(color_celltype)
      )
    ) |>
    dplyr::mutate(
      af = ifelse(
        af < 0.01,
        NA_real_,
        af
      )
    ) |>
    dplyr::mutate(
      cellvarianttype = colorcode[variant_type]
    ) |>
    dplyr::mutate(
      cellvarianttype = factor(
        cellvarianttype,
        levels = colorcode
      )
    ) -> forplot

  fn_plot_cell_af_somatic_variant_af(forplot, thetheme) -> p_af
  fn_plot_cell_af_somatic_variant_depth(forplot, thetheme) -> p_depth
  fn_plot_cell_af_somatic_variant_cell(forplot, thetheme) -> p_variant_cells
  fn_plot_cell_af_somatic_variant_celltype(forplot, thetheme) -> p_celltype

  .gseid <- unique(forplot$gseid)
  .srrid <- unique(forplot$srrid)
  .variant <- unique(forplot$variant)

  wrap_plots(
    p_af,
    plot_spacer(),
    p_depth,
    plot_spacer(),
    p_variant_cells,
    plot_spacer(),
    p_celltype,
    ncol = 1,
    heights = c(15, -1.05, 15, -1.05, 10, -1.05, 10),
    guides = "collect"
  ) +
    plot_annotation(
      title = glue::glue(
        "Variant {.variant} in {.gseid}-{.srrid}"
      ),
      theme = theme(
        plot.title = element_text(
          hjust = 0.5,
          size = 16,
          face = "bold"
        )
      )
    ) -> p_all

  p_all
}


# body --------------------------------------------------------------------

gse_data_variant_classification_clusteraf_bulkaf |>
  dplyr::filter(variant_type %in% c("homo", "hete")) -> df_variants


#
#
# ? Variant --------------------------------------------------------------------
#
#

gse_data_variant_classification_clusteraf_bulkaf |>
  dplyr::select(variant_type, variant) |>
  dplyr::distinct() -> variants_classification_unique

variants_classification_unique |>
  dplyr::mutate(isvariant = TRUE) |>
  tidyr::pivot_wider(
    names_from = variant_type,
    values_from = isvariant,
    values_fill = FALSE
  ) |>
  dplyr::rename(
    Homoplasmic = homo,
    Heteroplasmic = hete,
    Ethnicity = haplo,
    Somatic = somatic
  ) -> dt_allvariants_euler

color_variant_type <- c(
  "haplo" = "#BEBADAFF",
  "homo" = "#8DD3C7FF",
  "hete" = "#FFFFB3FF",
  "somatic" = "red"
)

library(eulerr)
plot(
  euler(
    dt_allvariants_euler[, c(
      "Ethnicity",
      "Homoplasmic",
      "Heteroplasmic",
      "Somatic"
    )],
    # shape = "ellipse",
    control = list(extraopt = FALSE)
  ),
  quantities = list(type = c("counts"), font = 3),
  labels = list(fontfamily = "serif"),
  edges = list(lty = 3),
  # fills = c("#BEBADAFF", "#8DD3C7FF", "#FFFFB3FF", "red")
  fills = color_variant_type
) -> p_euler


{
  outdir = "/home/liuc9/github/scMOCHA-data/analysis/high-res/MANUSCRIPTFIGURES/notuse"
  pdf(
    path(
      outdir,
      "Variant-type-Euler.pdf"
    ),
    width = 8,
    height = 6
  )
  print(p_euler)
  dev.off()
}
#
#
# ? variant is both hete and homo --------------------------------------------------------------------
#
#

df_variants |>
  tidyr::nest(.by = variant, .key = "sample") -> df_variants_sample


df_variants_sample |>
  dplyr::mutate(
    a = parallel::mclapply(
      sample,
      function(.x) {
        # .x <- df_variants_sample$sample[[1]]
        .x |>
          dplyr::count(variant_type) |>
          tidyr::pivot_wider(
            names_from = variant_type,
            values_from = n,
            values_fill = NA_integer_,
            names_prefix = "n_"
          )
      },
      mc.cores = 10
    )
  ) |>
  dplyr::select(-sample) |>
  tidyr::unnest(a) -> df_variants_sample_count

df_variants_sample_count |>
  dplyr::filter(!is.na(n_hete), !is.na(n_homo)) |>
  dplyr::mutate(total = (n_hete + n_homo))

thevariant <- "2442T>C"
thevariant <- "8362T>G"
df_variants |>
  dplyr::filter(variant == thevariant) |>
  dplyr::arrange(variant_type) |>
  dplyr::count(variant_type) |>
  tibble::deframe() -> n_hete_homo

df_variants |>
  dplyr::filter(variant == thevariant) |>
  dplyr::arrange(variant_type) |>
  dplyr::mutate(
    srrid = forcats::fct_reorder(srrid, Bulk, .desc = TRUE)
  ) |>
  tidyr::pivot_longer(
    cols = c(B, CD4_T, CD8_T, DC, Mono, NK, other, Bulk, other_T),
    names_to = "celltype",
    values_to = "af"
  ) |>
  dplyr::mutate(
    celltype = factor(
      celltype,
      levels = c(
        "Bulk",
        "B",
        "CD4_T",
        "CD8_T",
        "DC",
        "Mono",
        "NK",
        "other",
        "other_T"
      )
    )
  ) |>
  ggplot(aes(
    x = celltype,
    y = srrid,
    fill = af
  )) +
  geom_tile() +
  scale_fill_gradient(
    name = "AF",
    low = "white",
    high = "red"
  ) +
  theme_classic() +
  labs(
    x = "Celltype",
    y = "Sample"
  ) +
  labs(
    title = glue::glue(
      "Variant {thevariant} Allele Frequency\n{n_hete_homo['hete']} samples heteroplasmic and {n_hete_homo['homo']} samples homoplasmic"
    ),
  ) +
  theme(
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
  ) -> p_both_hete_homo

{
  outdir <- "/home/liuc9/github/scMOCHA-data/analysis/high-res/MANUSCRIPTFIGURES/notuse"
  ggsave(
    filename = glue::glue("Example-variant-{thevariant}-Hete-Homo.pdf"),
    plot = p_both_hete_homo,
    path = outdir,
    width = 6,
    height = 11
  )
}


#
#
# ? variant multiple --------------------------------------------------------------------
#
#

df_variants_sample_count |>
  dplyr::mutate(
    coord = parallel::mclapply(
      X = variant,
      FUN = \(.v) {
        # .v <- gse_data_variant_classification_clusteraf_bulkaf$variant[[1]]
        pos <- gsub(
          pattern = "(\\d+)([A-Z])>([A-Z])",
          replacement = "\\1",
          x = .v
        ) |>
          as.integer()
        ref <- gsub(
          pattern = "(\\d+)([A-Z])>([A-Z])",
          replacement = "\\2",
          x = .v
        )
        alt <- gsub(
          pattern = "(\\d+)([A-Z])>([A-Z])",
          replacement = "\\3",
          x = .v
        )
        data.table(
          seqnames = "MT",
          start = pos,
          end = pos,
          ref = ref,
          alt = alt
        )
      },
      mc.cores = 10
    )
  ) |>
  tidyr::unnest(
    cols = coord
  ) -> df_variants_sample_count_coord


df_variants$gseid |>
  unique() |>
  length()
df_variants$srrid |>
  unique() |>
  length()
df_variants_sample_count_coord$variant |>
  unique() |>
  length()
df_variants_sample_count_coord$start |> unique() |> length()

df_variants_sample_count_coord |>
  dplyr::arrange(start) |>
  dplyr::group_by(start) |>
  dplyr::filter(dplyr::n() > 2)


# thevariant <- "2442T>C"
thevariants1 <- c("13761A>G", "13761A>C", "13761A>T")
thevariants2 <- c("7585A>T", "7585A>C", "7585A>G")
thevariants <- c(thevariants1, thevariants2)
df_variants |>
  dplyr::filter(variant %in% thevariants)

df_variants |>
  dplyr::filter(variant %in% thevariants) |>
  dplyr::arrange(variant_type) |>
  dplyr::mutate(
    srrid = forcats::fct_reorder(srrid, Bulk, .desc = TRUE)
  ) |>
  tidyr::pivot_longer(
    cols = c(B, CD4_T, CD8_T, DC, Mono, NK, other, Bulk, other_T),
    names_to = "celltype",
    values_to = "af"
  ) |>
  dplyr::filter(celltype == "Bulk") |>
  ggplot(aes(
    x = variant,
    y = srrid,
    fill = af
  )) +
  geom_tile() +
  scale_fill_gradient(
    name = "AF",
    low = "white",
    high = "red"
  ) +
  theme_classic() +
  labs(
    x = "Variant",
    y = "Sample"
  ) +
  labs(
    title = glue::glue(
      "Variant in Bulk AF"
    ),
  ) +
  theme(
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
  ) -> p_multiple_variants_bulk

{
  outdir <- "/home/liuc9/github/scMOCHA-data/analysis/high-res/MANUSCRIPTFIGURES/notuse"
  ggsave(
    filename = glue::glue("Example-variant-multi.pdf"),
    plot = p_multiple_variants_bulk,
    path = outdir,
    width = 6,
    height = 8
  )
}


#
#
# ? somatic ethnicity --------------------------------------------------------------------
#
#

df_variants |>
  dplyr::select()

dt_allvariants_euler |>
  dplyr::filter(Somatic) |>
  dplyr::filter(Ethnicity) -> dt_allvariants_euler_somatic_ethnicity

df_variants |>
  dplyr::filter(variant %in% dt_allvariants_euler_somatic_ethnicity$variant) |>
  tidyr::nest(.by = variant, .key = "sample") |>
  dplyr::mutate(
    n = parallel::mclapply(
      sample,
      function(.x) {
        # .x <- df_variants_sample$sample[[1]]
        .x |>
          dplyr::count(variant_type) |>
          tidyr::pivot_wider(
            names_from = variant_type,
            values_from = n,
            values_fill = NA_integer_,
            names_prefix = "n_"
          )
      },
      mc.cores = 10
    )
  ) |>
  tidyr::unnest(cols = n) |>
  dplyr::filter(!is.na(n_hete), !is.na(n_homo)) |>
  as.data.table()


thevariant <- "9449C>T"
df_variants |>
  dplyr::filter(variant == thevariant) |>
  dplyr::arrange(variant_type) |>
  dplyr::count(variant_type) |>
  tibble::deframe() -> n_hete_homo


df_variants |>
  dplyr::filter(variant == thevariant) |>
  dplyr::arrange(variant_type) |>
  dplyr::mutate(
    srrid = forcats::fct_reorder(srrid, Bulk, .desc = TRUE)
  ) |>
  tidyr::pivot_longer(
    cols = c(B, CD4_T, CD8_T, DC, Mono, NK, other, Bulk, other_T),
    names_to = "celltype",
    values_to = "af"
  ) |>
  dplyr::mutate(
    celltype = factor(
      celltype,
      levels = c(
        "Bulk",
        "B",
        "CD4_T",
        "CD8_T",
        "DC",
        "Mono",
        "NK",
        "other",
        "other_T"
      )
    )
  ) -> forplot
forplot |>
  ggplot(aes(
    x = celltype,
    y = srrid,
    fill = af
  )) +
  geom_tile() +
  scale_fill_gradient(
    name = "AF",
    low = "white",
    high = "red"
  ) +
  theme_classic() +
  labs(
    x = "Celltype",
    y = "Sample"
  ) +
  theme(
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
    axis.text.x = element_text(face = "bold", size = 12),
    axis.title.x = element_text(face = "bold", size = 12)
  ) -> p_both_hete_homo

gse_data_variant_classification_clusteraf_bulkaf |>
  dplyr::filter(variant == thevariant) |>
  dplyr::select(srrid, variant_type) |>
  dplyr::mutate(v = 1) |>
  tidyr::pivot_wider(
    names_from = variant_type,
    values_from = v,
    values_fill = NA_integer_
  ) |>
  dplyr::mutate(
    srrid = factor(srrid, forplot$srrid |> levels())
  ) |>
  tidyr::pivot_longer(
    cols = -srrid,
    names_to = "variant_type",
    values_to = "value"
  ) |>
  ggplot(
    aes(
      x = variant_type,
      y = srrid,
      fill = value
    )
  ) +
  geom_tile() +
  scale_x_discrete(
    labels = c(
      "homo" = "Homoplasmic",
      "haplo" = "Ethnicity",
      "hete" = "Heteroplasmic",
      "somatic" = "Somatic"
    )
  ) +
  scale_fill_gradient(
    name = "Presence",
    low = "white",
    high = "blue"
  ) +
  theme_classic() +
  labs(
    x = "Variant Type",
    y = "Sample"
  ) +
  theme(
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
    axis.text.x = element_text(face = "bold", size = 12),
    axis.title.x = element_text(face = "bold", size = 12)
  ) -> p_variant_type


{
  wrap_plots(
    p_both_hete_homo,
    p_variant_type,
    ncol = 2,
    widths = c(1, 1),
    guides = "collect"
  ) +
    plot_annotation(
      title = glue::glue(
        "Variant {thevariant} Allele Frequency\n{n_hete_homo['hete']} samples heteroplasmic, {n_hete_homo['homo']} samples homoplasmic and 2 samples somatic"
      ),
      theme = theme(
        plot.title = element_text(size = 14, face = "bold", hjust = 0.5)
      )
    ) -> p_collect
}
{
  outdir <- "/home/liuc9/github/scMOCHA-data/analysis/high-res/MANUSCRIPTFIGURES/notuse"
  ggsave(
    filename = glue::glue("Example-variant-{thevariant}-somatic.pdf"),
    plot = p_collect,
    path = outdir,
    width = 15,
    height = 8
  )
}


fn_plot_cell_af_depth_forplot(
  thevariant = thevariant,
  thesrrid = "GSM7080018"
) |>
  fn_plot_cell_af_somatic_variant() -> p_somatic_variant1
fn_plot_cell_af_depth_forplot(
  thevariant = thevariant,
  thesrrid = "GSM7080039"
) |>
  fn_plot_cell_af_somatic_variant() -> p_somatic_variant2
fn_plot_cell_af_depth_forplot(
  thevariant = thevariant,
  thesrrid = "GSM4697622"
) |>
  fn_plot_cell_af_somatic_variant() -> p_somatic_variant3


{
  plotlist <- list(
    p_somatic_variant1,
    p_somatic_variant2,
    p_somatic_variant3
  )
  outdir <- "/home/liuc9/github/scMOCHA-data/analysis/high-res/MANUSCRIPTFIGURES/notuse"

  pdf(
    path(
      outdir,
      glue::glue("Example-variant-{thevariant}-somatic-GSM-multipage.pdf")
    ),
    width = 13,
    height = 8
  )

  for (i in seq_along(plotlist)) {
    cli_alert_info("Plotting page P{i}...")
    print(plotlist[[i]])
  }

  dev.off()
}

# footer ------------------------------------------------------------------

# save image --------------------------------------------------------------
