# pseudo-bulk
source("/home/liuc9/github/scMOCHA-data/analysis/high-res/00-colors.R")
ALLVARIANTSFORPLOT <- import(
  path("/home/liuc9/github/scMOCHA-data/analysis/high-res/MANUSCRIPTFIGURES") /
    "SAMPLE-VARIANT-CLASSIFICATION-CLUSTER-BULK-AF.xlsx"
) |>
  dplyr::filter(variant_type %in% c("homo", "hete"))

# Connection will be opened inside functions for parallel safety
# thevariant <- "3173G>A"
fn_plot_hetero_pseudo_bulk <- function(thevariant) {
  .d <- ALLVARIANTSFORPLOT |>
    dplyr::filter(variant == thevariant)

  color_celltype_bulk <- c(
    "Pseudo-bulk" = "red",
    color_celltype
  )

  .variant <- thevariant
  .n_srr <- unique(.d$srrid) |> length()
  .n_gse <- unique(.d$gseid) |> length()

  .d |>
    dplyr::select(gseid, srrid, variant, B:Bulk) |>
    tidyr::pivot_longer(
      cols = B:Bulk,
      names_to = "barcode",
      values_to = "af"
    ) -> .all_hetero_af_thevariant

  .all_hetero_af_thevariant |>
    dplyr::mutate(
      barcode = gsub(barcode, pattern = "_", replacement = " "),
      barcode = ifelse(barcode == "Bulk", "Pseudo-bulk", barcode),
    ) |>
    dplyr::mutate(
      barcode = factor(
        barcode,
        levels = names(color_celltype_bulk)
      ),
    ) -> .forplot

  .forplot |>
    dplyr::filter(barcode == "Pseudo-bulk") |>
    dplyr::arrange(-af) -> .rank_pseudo_bulk

  fn_xy_breaks_limits(
    c(0, .forplot$af),
    step = 0.1,
    max = FALSE
  ) -> .ybl

  .forplot |>
    dplyr::mutate(
      srrid = factor(srrid, levels = .rank_pseudo_bulk$srrid),
    ) |>
    ggplot(aes(x = srrid, y = af, fill = barcode)) +
    geom_col() +
    ggh4x::facet_grid2(
      ~barcode,
      strip = ggh4x::strip_themed(
        background_x = ggh4x::elem_list_rect(
          fill = color_celltype_bulk,
          color = NA
        ),
        text_x = ggh4x::elem_list_text(
          colour = "white",
          face = c("bold")
        )
      ),
      switch = "x",
    ) +
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
    scale_fill_manual(
      name = "Cell type",
      values = color_celltype_bulk
    ) +
    scale_y_continuous(
      name = "Heteroplasmy frequency",
      limits = .ybl$limits,
      breaks = c(.ybl$breaks, 0.05, 0.1) |> unique() |> sort(),
      expand = expansion(mult = c(0.005, 0.03)),
      labels = \(b) {
        dplyr::case_when(
          b == 0.1 ~ "gnomAD cutoff 10%",
          b == 0.05 ~ "our cutoff 5%",
          TRUE ~ scales::percent_format(accuracy = 1)(b)
        )
      },
    ) +

    theme(
      plot.margin = margin(t = 0.2, b = 0.1, l = 0.1, r = 0.2, unit = "cm"),
      panel.grid = element_blank(),
      panel.background = element_blank(),
      axis.text.x = element_blank(),
      axis.title = element_text(size = 12, color = "black", face = "bold"),
      axis.ticks.x = element_blank(),
      axis.line = element_line(color = "black"),
      # legend.position = c(0.2, 0.6),
      legend.position = "none",
      strip.placement = "outside",
      plot.title = element_text(hjust = 0.5, size = 16, face = "bold")
    ) +
    labs(
      x = "Individuals",
      title = glue::glue(
        "Heteroplasmy frequency of {.variant} in {.n_gse} projects and {.n_srr} samples"
      ),
    )
}
# fn_plot_hetero_pseudo_bulk(thevariant)

fn_plot_variant_ratio <- function(thevariant) {
  .m <- ALLVARIANTSFORPLOT |>
    dplyr::filter(variant == thevariant)

  .srrids <- unique(.m$srrid)

  conn <- DBI::dbConnect(
    duckdb::duckdb(),
    "/home/liuc9/github/scMOCHA-data/analysis/zzz/clean-data/all_hetero_af.cell.duckdb.1.2.1",
    readonly = TRUE
  )
  tbl_allvariants_cell <- dplyr::tbl(conn, "allvariants_cell")

  tbl_allvariants_cell |>
    dplyr::filter(
      variant == thevariant,
      srrid %in% .srrids
    ) |>
    as.data.table() -> .dt

  DBI::dbDisconnect(conn)
  colorcode <- setNames(names(color_variantcell), color_variantcell)

  color_celltype_bulk <- c(
    "Pseudo-bulk" = "red",
    color_celltype
  )

  .dt |>
    dplyr::count(
      gseid,
      srrid,
      variant_type,
      name = "count"
    ) |>
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
      varianttype = colorcode[variant_type]
    ) |>
    dplyr::mutate(
      variant_type = factor(
        variant_type,
        levels = color_variantcell
      )
    ) -> .d_dt

  .n_gse <- unique(.d_dt$gseid) |> length()
  .n_srr <- unique(.d_dt$srrid) |> length()
  .n_cells <- sum(.d_dt$count)
  .n_cells_hete <- sum(
    .d_dt |>
      dplyr::filter(varianttype == "Heteroplasmy") |>
      dplyr::pull(count),
    na.rm = TRUE
  )
  .variant <- thevariant

  .d_dt |>
    dplyr::mutate(
      varianttype = factor(
        varianttype,
        levels = names(color_variantcell) |>
          rev()
      )
    ) |>
    dplyr::group_by(
      srrid,
    ) |>
    dplyr::mutate(
      ratio = count / sum(count, na.rm = TRUE)
    ) |>
    dplyr::ungroup() -> .dd

  .dd |>
    dplyr::filter(
      varianttype == "Heteroplasmy"
    ) |>
    dplyr::arrange(ratio) |>
    dplyr::pull(srrid) -> rank_srrid

  .dd |>
    dplyr::mutate(
      srrid = factor(
        srrid,
        levels = rank_srrid
      )
    ) -> .d_forplot

  #plot count
  fn_plot_variant_ratio_count(.d_forplot, rank_srrid) -> p_count
  # plot ratio
  fn_plot_variant_ratio_ratio(.d_forplot, rank_srrid) -> p_ratio
  # paf
  fn_plot_variant_ratio_paf(.m, rank_srrid) -> p_haf_list
  p_haf <- p_haf_list$p
  .mean_af <- p_haf_list$mean_af

  wrap_plots(
    p_haf,
    p_ratio,
    p_count,
    ncol = 1,
    heights = c(0.8, 1, 1),
    guides = "collect"
  ) +
    plot_annotation(
      title = "{.variant} in {.n_gse} projects and  {scales::label_comma()(.n_srr)} samples" |>
        glue::glue(),
      subtitle = "{scales::label_percent(accuracy = 0.01)(.n_cells_hete/.n_cells)} ({scales::label_comma()(.n_cells_hete)}/{scales::label_comma()(.n_cells)}) cells with average HAF {scales::label_number(accuracy=0.01)(.mean_af)}" |>
        glue::glue(),
      theme = theme(
        plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
        plot.subtitle = element_text(hjust = 0.5, size = 14, face = "bold"),
      )
    )
}

fn_plot_variant_ratio_count <- function(.d_forplot, rank_srrid) {
  .d_forplot |>
    dplyr::group_by(srrid) |>
    dplyr::summarise(
      count = sum(count),
      .groups = "drop"
    ) -> .m

  fn_xy_breaks_limits(.m$count, step = 2000, max = FALSE) -> .count_ybl
  # count
  .d_forplot |>
    dplyr::mutate(
      srrid = factor(
        srrid,
        levels = rank_srrid
      )
    ) |>
    ggplot(aes(
      x = srrid,
      y = count,
    )) +
    geom_col(
      aes(
        fill = varianttype
      ),
      position = "stack"
    ) +
    scale_fill_manual(
      name = "Variant cell",
      values = color_variantcell
    ) +
    scale_x_discrete(
      limits = rank_srrid,
    ) +
    scale_y_continuous(
      limits = .count_ybl$limits,
      breaks = .count_ybl$breaks,
      expand = expansion(mult = c(0.005, 0.03)),
      labels = scales::label_comma()
    ) +
    theme(
      # panel.background = element_blank(),
      panel.grid = element_blank(),
      # axis.ticks = element_blank(),
      axis.line = element_line(color = "black"),
      axis.text.x = element_text(angle = 45, hjust = 1),
      axis.title.x = element_blank(),
      plot.margin = margin(t = 0, b = 0, unit = "cm"),
    ) +
    labs(y = "# of cell")
}
fn_plot_variant_ratio_ratio <- function(.d_forplot, rank_srrid) {
  .d_forplot |>
    dplyr::mutate(
      srrid = factor(
        srrid,
        levels = rank_srrid
      )
    ) |>
    ggplot(aes(
      x = srrid,
      y = ratio,
    )) +
    geom_col(
      aes(
        fill = varianttype
      ),
      position = "stack"
    ) +
    scale_fill_manual(
      name = "Variant cell",
      values = color_variantcell
    ) +
    scale_y_continuous(
      expand = expansion(add = c(0.005, 0.01)),
      labels = scales::percent_format(accuracy = 1)
    ) +
    scale_x_discrete(
      limits = rank_srrid,
    ) +
    theme(
      panel.grid = element_blank(),
      # axis.ticks = element_blank(),
      # axis.text.x = element_text(angle = 45, hjust = 1),
      axis.text.x = element_blank(),
      axis.title.x = element_blank(),
      plot.margin = margin(t = 0, b = 0, unit = "cm"),
      axis.line = element_line(color = "black"),
    ) +
    labs(y = "Cell ratio")
}
fn_plot_variant_ratio_paf <- function(.m, rank_srrid) {
  .m |>
    dplyr::select(gseid, srrid, af = Bulk) -> .bulk_forplot

  fn_xy_breaks_limits(c(0, .bulk_forplot$af), step = 0.1, max = FALSE) -> .ybl

  .bulk_forplot |>
    dplyr::mutate(
      srrid = factor(
        srrid,
        levels = rank_srrid
      )
    ) |>
    ggplot(aes(x = srrid, y = af, fill = "gold")) +
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
    scale_fill_manual(
      name = "Cell type",
      # values = color_celltype_bulk
      values = "gold"
    ) +
    scale_y_continuous(
      name = "Pseudo-bulk HAF",
      limits = .ybl$limits,
      breaks = c(.ybl$breaks, 0.05, 0.1) |> unique() |> sort(),
      labels = \(b) {
        dplyr::case_when(
          b == 0.1 ~ "gnomAD cutoff 10%",
          b == 0.05 ~ "our cutoff 5%",
          TRUE ~ scales::percent_format(accuracy = 1)(b)
        )
      },
      expand = expansion(mult = c(0.01, 0.01)),
    ) +
    theme(
      panel.grid = element_blank(),
      # axis.ticks = element_blank(),
      # axis.text.x = element_text(angle = 45, hjust = 1),
      axis.text.x = element_blank(),
      axis.title.x = element_blank(),
      plot.margin = margin(t = 0, b = 0, unit = "cm"),
      axis.line = element_line(color = "black"),
      panel.background = element_blank(),
      legend.position = "none",
    ) -> .p
  mean(.bulk_forplot$af) -> .mean_af
  list(
    p = .p,
    mean_af = .mean_af
  )
}

# fn_plot_variant_ratio(thevariant)
