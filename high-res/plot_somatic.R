dotenv(".env")
source(path(
  Sys.getenv("HIGHRESDIR"),
  "00-colors.R"
))

fn_plot_cell_af_depth_forplot <- function(thevariant, thesrrid) {
  conn <- conn_db(
    Sys.getenv("DUCKDB_PATH"),
    readonly = TRUE
  )
  colorcode <- setNames(names(color_variantcell), color_variantcell)

  dplyr::tbl(
    conn,
    "allvariants_cell"
  ) |>
    dplyr::filter(
      srrid == thesrrid,
      variant == thevariant
    ) |>
    as.data.table() |>
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
    ) |>
    dplyr::arrange(
      cellvarianttype,
      celltype,
      -af
    ) -> forplot_
  forplot_ |>
    dplyr::mutate(
      barcode = factor(
        barcode,
        levels = forplot_$barcode
      )
    ) -> forplot

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
fn_plot_cell_number_of_cells <- function(forplot) {
  forplot |>
    dplyr::count(variant_type, celltype) |>
    dplyr::mutate(
      variant_type = as.character(variant_type),
    ) |>
    ggplot(aes(
      x = variant_type,
      y = celltype,
      fill = variant_type
    )) +
    geom_tile() +
    # scale_fill_gradient(
    #   name = "Number of cells",
    #   low = "white",
    #   high = "red",
    #   na.value = "grey"
    # ) +
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
    geom_text(
      aes(
        label = n,
        color = variant_type
      ),
      size = 4,
      show.legend = FALSE
    ) +
    scale_color_manual(
      values = c(
        "red" = "white",
        "darkblue" = "white",
        "gray" = "black",
        "white" = "black"
      )
    ) +
    scale_x_discrete(
      limits = c("red", "darkblue", "gray", "white"),
    ) +
    theme_bw() +
    theme(
      panel.grid = element_blank(),
      legend.position = "none",
      axis.text = element_text(face = "bold"),
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      axis.text.y = element_text(size = 12),
      axis.title = element_blank(),
    ) +
    coord_fixed(
      ratio = 1
    )
}

fn_plot_cell_af_somatic_variant <- function(forplot) {
  source(path(
    Sys.getenv("HIGHRESDIR"),
    "00-colors.R"
  ))

  colorcode <- setNames(names(color_variantcell), color_variantcell)
  thetheme <- theme(
    panel.background = element_blank(),
    panel.grid = element_blank(),
    axis.ticks = element_blank(),
    axis.text = element_blank(),
    axis.title.x = element_blank(),
    plot.margin = margin(t = 0, b = 0, unit = "cm"),
  )

  fn_plot_cell_af_somatic_variant_af(forplot, thetheme) -> p_af
  fn_plot_cell_af_somatic_variant_depth(forplot, thetheme) -> p_depth
  fn_plot_cell_af_somatic_variant_cell(forplot, thetheme) -> p_variant_cells
  fn_plot_cell_af_somatic_variant_celltype(forplot, thetheme) -> p_celltype
  fn_plot_cell_number_of_cells(forplot) -> p_heatmap

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
  ) -> p_main

  wrap_plots(
    wrap_elements(p_main),
    p_heatmap,
    ncol = 2,
    widths = c(3, 1)
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

fn_plot_somatic <- function(thevariant, thesrrid) {
  fn_plot_cell_af_depth_forplot(
    thevariant = thevariant,
    thesrrid = thesrrid
  ) -> forplot

  fn_plot_cell_af_somatic_variant(forplot)
}

\() {
  outdir <- "/home/liuc9/github/scMOCHA-data/analysis/high-res/MANUSCRIPTFIGURES-notuse"

  pdf(
    path(
      outdir,
      glue::glue("Example-variant-{thesrrid}-somatic-GSM-multipage.pdf")
    ),
    width = 16,
    height = 8
  )

  # thesrrid = "GSM7493841"
  # thevariant <- "6440C>A"

  print(fn_plot_somatic(
    thevariant = thevariant,
    thesrrid = thesrrid
  ))

  dev.off()
}

\() {
  outdir <- "/home/liuc9/github/scMOCHA-data/analysis/high-res/MANUSCRIPTFIGURES-notuse"

  thesrrid = "GSM7437874"
  thevariant <- "7757G>A"

  pdf(
    path(
      outdir,
      glue::glue(
        "Example-variant-{thesrrid}-{thevariant}-somatic-GSM-multipage.pdf"
      )
    ),
    width = 16,
    height = 8
  )

  print(fn_plot_somatic(
    thevariant = thevariant,
    thesrrid = thesrrid
  ))

  dev.off()
}

\() {
  outdir <- "/home/liuc9/github/scMOCHA-data/analysis/high-res/MANUSCRIPTFIGURES-notuse"

  thesrrid = "GSM7493835"
  thevariant <- "4175G>A"

  pdf(
    path(
      outdir,
      glue::glue(
        "Example-variant-{thesrrid}-{thevariant}-somatic-GSM-multipage.pdf"
      )
    ),
    width = 16,
    height = 8
  )

  print(fn_plot_somatic(
    thevariant = thevariant,
    thesrrid = thesrrid
  ))

  dev.off()
}
