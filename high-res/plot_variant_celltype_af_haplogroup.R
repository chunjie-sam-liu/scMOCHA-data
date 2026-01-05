dotenv(".env")
fn_plot_variant_celltype_af_haplogroup <- function(
  thevariant,
  vtype = NULL
) {
  df <- import(
    path(
      Sys.getenv("HIGHRESDIR"),
      "MANUSCRIPTFIGURES/SAMPLE-VARIANT-CLASSIFICATION-CLUSTER-BULK-AF.xlsx"
    )
  )
  df |>
    dplyr::filter(variant == thevariant) |>
    dplyr::arrange(variant_type) |>
    (\(x) {
      if (!is.null(vtype)) dplyr::filter(x, variant_type %in% vtype) else x
    })() |>
    dplyr::mutate(
      srrid = forcats::fct_reorder(srrid, Bulk, .desc = TRUE)
    ) -> df_thevariant

  df_thevariant |>
    dplyr::count(variant_type) |>
    tibble::deframe() -> n_hete_homo

  df_thevariant |>
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
          "other_T",
          "NK",
          "DC",
          "Mono",
          "other"
        )
      )
    ) -> forplot

  forplot |>
    dplyr::select(Haplogroup, Verbose_haplogroup) |>
    dplyr::distinct() |>
    dplyr::mutate(
      Haplogroup_s = purrr::map_chr(
        .x = Haplogroup,
        .f = \(.x) {
          # if (stringr::str_starts(.x, "L")) {
          #   gsub("L", "L0", .x)
          # }
          gsub("\\d+.*", "", .x)
        }
      )
    ) |>
    dplyr::mutate(
      color_haplogroup = color(color_haplogroup[Haplogroup_s])
    ) |>
    dplyr::mutate(
      color_verbose_haplogroup = ifelse(
        Haplogroup == Verbose_haplogroup,
        color_haplogroup,
        prismatic::clr_lighten(color_haplogroup, 0.5)
      )
    ) |>
    dplyr::select(-Haplogroup_s) -> haplo_colors

  c(
    haplo_colors$color_haplogroup,
    haplo_colors$color_verbose_haplogroup
  ) -> haplo_color_vector
  names(haplo_color_vector) <- c(
    haplo_colors$Haplogroup,
    haplo_colors$Verbose_haplogroup
  )

  forplot |>
    dplyr::select(srrid, Haplogroup, Verbose_haplogroup) |>
    tidyr::pivot_longer(
      cols = c(Haplogroup, Verbose_haplogroup),
      names_to = "type",
      values_to = "Haplogroup"
    ) -> forplot_haplogroup

  forplot_haplogroup |>
    ggplot(aes(
      x = type,
      y = srrid,
      fill = Haplogroup
    )) +
    geom_tile(
      show.legend = FALSE
    ) +
    scale_fill_manual(
      values = haplo_color_vector
    ) +
    geom_text(
      data = forplot_haplogroup |> dplyr::distinct(),
      aes(
        label = Haplogroup,
      ),
      color = "black",
      fontface = "bold"
    ) +
    scale_x_discrete(
      expand = c(0, 0)
    ) +
    theme_classic() +
    theme(
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      axis.title.x = element_blank(),
      axis.line = element_blank(),
      axis.title.y = element_blank(),
      axis.text.y = element_text(
        face = "bold",
        size = 12
      ),
      axis.ticks.y = element_blank(),
    ) -> p_haplogroup

  forplot |>
    ggplot(aes(
      x = celltype,
      y = srrid,
      fill = af
    )) +
    geom_tile() +
    geom_text(
      aes(
        label = ifelse(af >= 0.01, sprintf("%.2f", af), "")
      ),
      color = "black",
      fontface = "bold"
    ) +
    scale_fill_gradient(
      name = "AF",
      low = "white",
      high = "red"
    ) +
    scale_x_discrete(
      expand = c(0, 0)
    ) +
    theme_classic() +
    labs(
      x = "Celltype",
      y = "Sample"
    ) +
    theme(
      plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
      axis.text.x = element_text(
        face = "bold",
        size = 12,
        angle = 15,
        hjust = 1
      ),
      # axis.title.x = element_text(face = "bold", size = 12),
      axis.title = element_blank(),
      axis.line.y = element_blank(),
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
    ) -> p_af

  df_thevariant |>
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
      expand = c(0, 0),
      limits = c("homo", "haplo", "hete", "somatic"),
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
      axis.text.x = element_text(
        face = "bold",
        size = 12,
        angle = 15,
        hjust = 1
      ),
      # axis.title.x = element_text(face = "bold", size = 12),
      axis.title = element_blank(),
      axis.line.y = element_blank(),
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank()
    ) -> p_variant_type

  wrap_plots(
    p_haplogroup,
    p_af,
    p_variant_type,
    ncol = 3,
    widths = c(0.2, 1, 0.8),
    guides = "collect"
  ) +
    plot_annotation(
      title = glue::glue(
        "Variant {thevariant} Allele Frequency"
      ),
      theme = theme(
        plot.title = element_text(size = 14, face = "bold", hjust = 0.5)
      )
    ) -> p_collect_all

  wrap_plots(
    p_haplogroup,
    p_af,
    ncol = 2,
    widths = c(0.2, 1),
    guides = "collect"
  ) +
    plot_annotation(
      title = glue::glue(
        "Variant {thevariant} Allele Frequency"
      ),
      theme = theme(
        plot.title = element_text(size = 14, face = "bold", hjust = 0.5)
      )
    ) -> p_collect_2
  list(
    p_collect_2 = p_collect_2,
    p_collect_all = p_collect_all
  )
}
