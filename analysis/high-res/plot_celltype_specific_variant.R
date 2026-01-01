conn <- DBI::dbConnect(
  duckdb::duckdb(),
  "/home/liuc9/github/scMOCHA-data/analysis/zzz/clean-data/all_hetero_af.cell.duckdb.1.2.1",
  readonly = TRUE
)
DBI::dbListTables(conn)

tbl_allvariants_cell <- dplyr::tbl(conn, "allvariants_cell")
# dplyr::tbl(conn, "allvariants_af_cell")
# dplyr::tbl(conn, "all_hetero_af_cell")

load_pkg(
  ggdist
)
source("/home/liuc9/github/scMOCHA-data/analysis/high-res/00-colors.R")


# thevariant <- "7428G>A"
# thegseid <- "GSE161354"
# thesrrid <- "GSM4905217"

fn_plot_ggdist <- function(
  thevariant,
  thegseid,
  thesrrid
) {
  tbl_allvariants_cell |>
    dplyr::filter(
      srrid == thesrrid,
      variant == thevariant,
      variant_type %in% c("colorful", "black")
    ) |>
    as.data.table() -> .d

  .variant <- .d$variant[1]
  .gseid <- .d$gseid[1]
  .srrid <- .d$srrid[1]

  .d |>
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
        levels = names(color_celltype) |> rev()
      )
    ) -> forplot_

  forplot_ |>
    ggplot(aes(
      x = af,
      y = celltype,
      fill = celltype
    )) +
    ggdist::stat_halfeye(
      scale = 1
    ) +
    ggdist::stat_interval(
      show.legend = FALSE,
    ) +
    stat_summary(
      geom = "point",
      fun = median,
      show.legend = FALSE
    ) +
    scale_fill_manual(
      values = color_celltype,
      na.value = "grey50"
    ) +
    scale_color_manual(
      values = MetBrewer::met.brewer("VanGogh3")
    ) +
    # scale_color_brewer() +
    guides(col = "none") +
    ggridges::theme_ridges() +
    # ggridges::stat_density_ridges(
    #   quantile_lines = TRUE, quantiles = 2
    # ) +
    theme(
      legend.position = "none",
      plot.title = element_text(
        hjust = 0.5,
        size = 14,
        color = "black",
        face = "bold"
      ),
      axis.title.y = element_blank()
    ) +
    labs(
      title = glue::glue("m.{.variant}\n({.gseid}-{.srrid})"),
      x = "Heteroplasmy Level",
      y = "Cell Type"
    )
}

# fn_plot_ggdist(
#   thevariant = thevariant,
#   thegseid = thegseid,
#   thesrrid = thesrrid
# )

fn_plot_joy <- function(
  thevariant,
  thegseid,
  thesrrid
) {
  tbl_allvariants_cell |>
    dplyr::filter(
      srrid == thesrrid,
      variant == thevariant,
      variant_type %in% c("colorful", "black")
    ) |>
    as.data.table() -> .d

  .variant <- .d$variant[1]
  .gseid <- .d$gseid[1]
  .srrid <- .d$srrid[1]

  .d |>
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
        levels = names(color_celltype) |> rev()
      )
    ) -> forplot_

  forplot_ |>
    dplyr::count(celltype) -> n_celltype

  n_celltype <- n_celltype |>
    dplyr::mutate(
      celltype = factor(celltype, levels = levels(forplot_$celltype))
    ) |>
    dplyr::arrange(celltype)

  forplot_ |>
    ggplot(aes(
      x = af,
      y = celltype,
      fill = celltype
    )) +
    ggridges::geom_density_ridges(
      scale = 2,
      alpha = 0.8,
      rel_min_height = 0.01,
      size = 0.1
    ) +
    scale_fill_manual(
      values = color_celltype,
      na.value = "grey50",
      name = "Cell Type"
    ) +
    scale_y_discrete(
      labels = function(x) {
        n_counts <- n_celltype$n[match(x, n_celltype$celltype)]
        glue("{x} (n={n_counts})")
      }
    ) +
    ggridges::theme_ridges() +
    theme(
      legend.position = "none",
      plot.title = element_text(
        hjust = 0.5,
        size = 14,
        color = "black",
        face = "bold"
      ),
      panel.background = element_blank(),
      panel.grid = element_line(colour = "grey", linetype = "dashed"),
      panel.grid.major = element_line(
        colour = "grey",
        linetype = "dashed",
        size = 0.2
      ),
      axis.title.y = element_blank()
    ) +
    labs(
      title = glue::glue("m.{.variant}\n({.gseid}-{.srrid})"),
      x = "Heteroplasmy Level",
      y = "Cell Type"
    )
}

# fn_plot_joy(
#   thevariant = thevariant,
#   thegseid = thegseid,
#   thesrrid = thesrrid
# )

fn_plot_hist <- function(
  thevariant,
  thegseid,
  thesrrid
) {
  tbl_allvariants_cell |>
    dplyr::filter(
      # gseid == thegseid,
      srrid == thesrrid,
      variant == thevariant,
      variant_type %in% c("colorful", "black")
    ) |>
    dplyr::collect() -> .d
  .variant <- .d$variant[1]
  .gseid <- .d$gseid[1]
  .srrid <- .d$srrid[1]

  .d |>
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
        levels = names(color_celltype)
      )
    ) -> forplot_

  forplot_ |>
    dplyr::count(celltype) |>
    dplyr::mutate(
      label = glue::glue(
        "n={scales::label_comma()(n)}"
      )
    ) -> .forlabel

  fn_xy_breaks_limits(forplot_$af, step = 0.2) -> .xbl
  # fn_xy_breaks_limits(forplot_$celltype, step = 1) -> .ybl

  forplot_ |>
    ggplot(aes(
      x = af,
      # y = celltype,
      fill = celltype
    )) +
    geom_histogram(binwidth = 0.1) +
    scale_fill_manual(
      values = color_celltype,
    ) +
    geom_text(
      data = .forlabel,
      aes(
        x = 0.5,
        y = Inf,
        label = label
      ),
      vjust = 1.5
    ) +
    scale_y_continuous(
      labels = scales::label_number(accuracy = 1),
      expand = expansion(mult = c(0.01, 0.01)),
    ) +
    scale_x_continuous(
      # limits = .xbl$limits,
      breaks = seq(0, 1, 0.2),
      labels = scales::label_number(accuracy = 0.1),
      expand = expansion(add = c(0.01, 0.01)),
    ) +
    theme(
      legend.position = "none",
      # axis.text.y = element_blank(),
      axis.line = element_line(color = "black"),
      strip.background = element_rect(
        fill = "white",
        color = "black"
      ),
      plot.title = element_text(
        hjust = 0.5,
        size = 14,
        color = "black",
        face = "bold"
      ),
      panel.background = element_blank(),
    ) +
    ggh4x::facet_wrap2(
      ~celltype,
      nrow = 1,
      # ncol = 8,
      strip.position = "top",
      strip = ggh4x::strip_themed(
        background_x = ggh4x::elem_list_rect(
          fill = color_celltype
        ),
        text_x = ggh4x::elem_list_text(
          colour = "white",
          face = c("bold")
        ),
        by_layer_y = TRUE,
      ),
      scales = "free_y",
    ) +
    labs(
      title = glue::glue("m.{.variant}\n({.gseid}-{.srrid})"),
      x = glue::glue("Heteroplasmy Level"),
      y = "Cell Count"
    )
}

# fn_plot_hist(
#   thevariant = thevariant,
#   thegseid = thegseid,
#   thesrrid = thesrrid
# )

fn_plot_cumulative_fraction <- function(
  thevariant,
  thegseid,
  thesrrid
) {
  tbl_allvariants_cell |>
    dplyr::filter(
      srrid == thesrrid,
      variant == thevariant,
      variant_type %in% c("colorful", "black")
    ) |>
    dplyr::collect() -> .d

  .variant <- .d$variant[1]
  .gseid <- .d$gseid[1]
  .srrid <- .d$srrid[1]

  .d |>
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
        levels = names(color_celltype)
      )
    ) -> forplot_

  forplot_ |>
    dplyr::count(celltype) |>
    dplyr::mutate(
      label = glue::glue(
        "n={scales::label_comma()(n)}"
      )
    ) -> .forlabel

  # Create labels with cell counts for legend
  celltype_labels <- setNames(
    .forlabel$label |> as.character(),
    .forlabel$celltype |> as.character()
  )

  forplot_ |>
    ggplot(aes(
      x = af, # Convert to percentage
      color = celltype
    )) +
    stat_ecdf(
      geom = "step",
      size = 1.2
    ) +
    scale_color_manual(
      values = color_celltype,
      labels = celltype_labels,
      na.translate = FALSE
    ) +
    scale_x_continuous(
      limits = c(0, 1),
      breaks = seq(0, 1.0, 0.25),
      expand = expansion(add = c(0.01, 0.01)),
    ) +
    scale_y_continuous(
      limits = c(0, 1),
      breaks = seq(0, 1, 0.25),
      labels = scales::label_number(accuracy = 0.01),
      expand = expansion(add = c(0.01, 0.01)),
    ) +
    theme_classic() +
    theme(
      legend.position = "right",
      legend.title = element_blank(),
      legend.text = element_text(size = 10),
      plot.title = element_text(
        hjust = 0.5,
        size = 14,
        color = "black",
        face = "bold"
      ),
      axis.line = element_line(color = "black", size = 0.8),
      axis.ticks = element_line(color = "black", size = 0.8),
      panel.background = element_blank(),
    ) +
    labs(
      title = glue::glue("m.{.variant}\n({.gseid}-{.srrid})"),
      x = "Heteroplasmy Level",
      y = "Cumulative fraction"
    ) +
    coord_fixed(ratio = 1)
}

# fn_plot_cumulative_fraction(
#   thevariant = thevariant,
#   thegseid = thegseid,
#   thesrrid = thesrrid
# )

fn_plot_joy_celltype_level2_level3 <- function(
  thevariant,
  thegseid,
  thesrrid,
  thecelltype,
  thecelltype_prefix,
  thecelltype_level
) {
  # thevariant <- thevariant_celltype_df$thevariant[1]
  # thegseid <- thevariant_celltype_df$thegseid[1]
  # thesrrid <- thevariant_celltype_df$thesrrid[1]
  # thecelltype <- thevariant_celltype_df$thecelltype[1]
  # thecelltype_prefix <- thevariant_celltype_df$thecelltype_prefix[1]
  # thecelltype_level <- thevariant_celltype_df$thecelltype_level[1]
  celltypedetail <- import(
    "/mnt/isilon/u01_project/large-scale/liuc9/raw/{thegseid}/final/{thesrrid}/sc_azimuth_celltype.csv" |>
      glue::glue()
  )

  tbl_allvariants_cell |>
    dplyr::filter(
      gseid == thegseid,
      srrid == thesrrid,
      variant == thevariant,
      variant_type %in% c("colorful", "black")
    ) |>
    dplyr::select(-celltype) |>
    dplyr::collect() |>
    dplyr::left_join(
      celltypedetail,
      by = c("barcode")
    ) |>
    dplyr::rename(
      plotcelltype = "celltype_{thecelltype_level}" |> glue::glue(),
    ) -> thevariant_data

  thevariant_data |>
    dplyr::filter(
      celltype == thecelltype,
    ) |>
    dplyr::filter(grepl(thecelltype_prefix, plotcelltype)) |>
    dplyr::mutate(
      plotcelltype = factor(
        plotcelltype
      )
    ) -> forplot

  levels(forplot$plotcelltype)

  color_celltype_detail <- log(seq(
    1,
    exp(1),
    length.out = length(levels(forplot$plotcelltype))
  )) |>
    purrr::map_chr(
      ~ prismatic::clr_lighten(
        color_celltype[thecelltype],
        .x
      )
    )
  names(color_celltype_detail) <- levels(forplot$plotcelltype)

  n_celltype <- forplot |>
    dplyr::count(plotcelltype) |>
    dplyr::arrange(plotcelltype)

  forplot |>
    ggplot(aes(
      x = af,
      y = plotcelltype,
      fill = plotcelltype
    )) +
    ggridges::geom_density_ridges(
      scale = 2,
      alpha = 0.8,
      rel_min_height = 0.01,
      size = 0.1
    ) +
    scale_fill_manual(
      values = color_celltype_detail,
      na.value = "grey50"
    ) +
    scale_y_discrete(
      labels = function(x) {
        n_counts <- n_celltype$n[match(x, n_celltype$plotcelltype)]
        glue("{x} (n={n_counts})")
      }
    ) +
    ggridges::theme_ridges() +
    theme(
      legend.position = "none",
      plot.title = element_text(
        hjust = 0.5,
        size = 14,
        color = "black",
        face = "bold"
      ),
      panel.background = element_blank(),
      panel.grid = element_line(colour = "grey", linetype = "dashed"),
      panel.grid.major = element_line(
        colour = "grey",
        linetype = "dashed",
        size = 0.2
      ),
      axis.title.y = element_blank()
    ) +
    labs(
      title = "{thecelltype}-{thecelltype_level}-{thevariant}\n({thegseid}-{thesrrid})" |>
        glue::glue(),
      x = "Heteroplasmy Level",
      y = "Cell Type"
    )
}

fn_plot_joy_celltype_detail <- function(
  thevariant,
  thegseid,
  thesrrid
) {
  tibble::tibble(
    thevariant = thevariant,
    thegseid = thegseid,
    thesrrid = thesrrid,
    thecelltype = c(
      c("B", "CD4 T", "CD8 T", "other T") |> rep(times = 2),
      c("NK", "DC", "Mono", "other") |> rep(times = 2)
    ),
    thecelltype_prefix = c(
      c("B", "CD4", "CD8", "") |> rep(times = 2),
      c("NK", "DC", "Mono", "") |> rep(times = 2)
    ),
    thecelltype_level = c("l2", "l3") |> rep(each = 4) |> rep(times = 2)
  ) -> thevariant_celltype_df

  thevariant_celltype_df |>
    dplyr::mutate(
      # p = parallel::mcmapply(
      p = mapply(
        FUN = fn_plot_joy_celltype_level2_level3,
        thevariant = thevariant,
        thegseid = thegseid,
        thesrrid = thesrrid,
        thecelltype = thecelltype,
        thecelltype_prefix = thecelltype_prefix,
        thecelltype_level = thecelltype_level,
        # mc.cores = 5,
        SIMPLIFY = FALSE
      )
    ) -> plot_thevariant_celltype_list

  plot_thevariant_celltype_list |>
    dplyr::pull(p) |>
    wrap_plots(
      ncol = 4
    ) +
    plot_layout(
      guides = "collect",
    )
}

# fn_plot_joy_celltype_detail(
#   thevariant = thevariant,
#   thegseid = thegseid,
#   thesrrid = thesrrid
# )
