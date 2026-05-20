load_pkg(jutils)
dotenv(".env")
source(path(
  Sys.getenv("HIGHRESDIR"),
  "00-colors.R"
))
outdir <- path(Sys.getenv("OUTDIR"))

METAFULL <- import(outdir / "SAMPLES-METADATA-FULL.xlsx")
ALLVARIANTS <- import(
  outdir / "SAMPLE-VARIANT-CLASSIFICATION-CLUSTER-BULK-AF.xlsx"
) |>
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
  )

SOMATIC_VARIANTS <- ALLVARIANTS |>
  dplyr::filter(variant_type == "somatic")


fn_plot_somatic_dist <- function() {
  SOMATIC_VARIANTS |>
    dplyr::mutate(srridvariant = glue("{srrid}_{variant}")) |>
    dplyr::select(
      srridvariant,
      c(B, CD4_T, CD8_T, DC, Mono, NK, other, Bulk, other_T)
    ) |>
    tibble::column_to_rownames("srridvariant") |>
    as.matrix() -> mat_somatic_af

  # Identify dominant cell type for each variant (excluding Bulk)
  mat_specific <- mat_somatic_af[, c(
    "B",
    "CD4_T",
    "CD8_T",
    "other_T",
    "NK",
    "DC",
    "Mono",
    "other"
  )]

  # For each variant, find the cell type with maximum AF
  dominant_celltype <- apply(mat_specific, 1, function(x) {
    if (max(x, na.rm = TRUE) < 0.01) {
      return("Bulk")
    }
    colnames(mat_specific)[which.max(x)]
  })

  # Create ordering: group by dominant cell type, then by AF within each group
  celltype_order <- c(
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

  variant_df <- data.frame(
    srridvariant = rownames(mat_somatic_af),
    dominant_celltype = dominant_celltype,
    max_af = apply(mat_specific, 1, max, na.rm = TRUE),
    stringsAsFactors = FALSE
  )

  variant_df$dominant_celltype <- factor(
    variant_df$dominant_celltype,
    levels = celltype_order
  )
  variant_ordered <- variant_df[
    order(variant_df$dominant_celltype, -variant_df$max_af),
  ]
  clustered_order <- variant_ordered$srridvariant

  SOMATIC_VARIANTS |>
    dplyr::mutate(srridvariant = glue("{srrid}_{variant}")) |>
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
    ) |>
    dplyr::mutate(
      srridvariant = factor(
        srridvariant,
        levels = clustered_order
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
    dplyr::filter(!is.na(Haplogroup)) |>
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
    dplyr::select(srridvariant, srrid, Haplogroup, Verbose_haplogroup) |>
    dplyr::mutate(
      srridvariant = factor(
        srridvariant,
        levels = clustered_order
      )
    ) |>
    tidyr::pivot_longer(
      cols = c(Haplogroup, Verbose_haplogroup),
      names_to = "type",
      values_to = "Haplogroup"
    ) -> forplot_haplogroup

  forplot_haplogroup |>
    ggplot(aes(
      x = type,
      y = srridvariant,
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
    scale_y_discrete(
      labels = function(x) sub("_.*", "", x)
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

  gtf_gene_df <- import(
    path(Sys.getenv("REPODIR"), "config/mtdna_genes_dloop.qs")
  )

  forplot |>
    dplyr::select(srridvariant, variant) |>
    dplyr::distinct() |>
    dplyr::mutate(
      pos = as.integer(gsub("(\\d+)[A-Z]>[A-Z]", "\\1", variant))
    ) |>
    dplyr::mutate(
      gene_df = purrr::map(
        .x = pos,
        .f = \(.x) {
          gtf_gene_df |>
            dplyr::filter(start <= .x & end >= .x) |>
            dplyr::select(gene_name, TYPE, COLOR)
        }
      )
    ) |>
    tidyr::unnest(cols = gene_df) -> forplot_variant_gene

  # Create color vector: COLOR is for TYPE, lighten for different gene_name
  type_color_df <- forplot_variant_gene |>
    dplyr::select(TYPE, COLOR) |>
    dplyr::distinct()

  # For each TYPE, create lightened colors for each gene_name
  gene_colors_list <- forplot_variant_gene |>
    dplyr::select(TYPE, gene_name, COLOR) |>
    dplyr::distinct() |>
    dplyr::arrange(TYPE, gene_name) |>
    dplyr::group_by(TYPE) |>
    dplyr::mutate(
      gene_idx = dplyr::row_number(),
      n_genes = dplyr::n()
    ) |>
    dplyr::ungroup() |>
    dplyr::mutate(
      gene_color = purrr::pmap_chr(
        .l = list(
          .color = COLOR,
          .idx = gene_idx
        ),
        .f = \(.color, .idx) {
          if (.idx == 1) {
            .color
          } else {
            as.character(prismatic::clr_lighten(.color, shift = 0.05 * .idx))
          }
        }
      )
    )

  # Create named vectors
  type_color_vector <- type_color_df$COLOR
  names(type_color_vector) <- type_color_df$TYPE

  gene_color_vector <- gene_colors_list$gene_color
  names(gene_color_vector) <- gene_colors_list$gene_name

  all_color_vector <- c(type_color_vector, gene_color_vector)

  forplot_variant_gene |>
    dplyr::select(srridvariant, gene_name, TYPE, COLOR) |>
    dplyr::mutate(
      srridvariant = factor(srridvariant, levels = clustered_order)
    ) |>
    tidyr::pivot_longer(
      cols = c(TYPE, gene_name),
      names_to = "type",
      values_to = "label"
    ) |>
    dplyr::mutate(
      type = factor(type, levels = c("TYPE", "gene_name"))
    ) -> forplot_gene

  forplot_gene |>
    ggplot(aes(
      x = type,
      y = srridvariant,
      fill = label
    )) +
    geom_tile(
      show.legend = FALSE
    ) +
    scale_fill_manual(
      values = all_color_vector
    ) +
    geom_text(
      data = forplot_gene |> dplyr::distinct(),
      aes(label = label),
      color = "black",
      fontface = "bold",
      size = 3
    ) +
    scale_x_discrete(
      expand = c(0, 0)
    ) +
    scale_y_discrete(
      labels = function(x) sub(".*_", "", x)
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
        size = 10
      ),
      axis.ticks.y = element_blank(),
    ) -> p_gene

  forplot |>
    ggplot(aes(
      x = celltype,
      y = srridvariant,
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
    scale_y_discrete(
      labels = function(x) sub(".*_", "", x)
    ) +
    theme_classic() +
    labs(
      x = "Celltype",
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

  SOMATIC_VARIANTS$variant |> unique() |> length() -> n_variants
  SOMATIC_VARIANTS$srrid |> unique() |> length() -> n_samples

  wrap_plots(
    p_haplogroup,
    p_gene,
    p_af,
    ncol = 3,
    widths = c(0.2, 0.3, 1),
    guides = "collect"
  ) +
    plot_annotation(
      title = glue::glue(
        "Somatic variants distribution across cell types and haplogroups\n({n_variants} variants across {n_samples} samples)"
      ),
      theme = theme(
        plot.title = element_text(size = 14, face = "bold", hjust = 0.5)
      )
    ) -> p_collect

  pdf(
    file = path(
      Sys.getenv("OUTDIRNOTUSE")
    ) /
      "HOTSPOTS-SOMATIC-CLUSTER.pdf",
    width = 15,
    height = 35
  )
  print(p_collect)
  dev.off()

  png(
    file = path(
      Sys.getenv("OUTDIRNOTUSE")
    ) /
      "HOTSPOTS-SOMATIC-CLUSTER.png",
    width = 15,
    height = 35,
    units = "in",
    res = 300
  )
  print(p_collect)
  dev.off()
}

# fn_plot_somatic_dist()

fn_save_somatic_table <- function() {
  load_pkg(openxlsx, gt)

  # Prepare data same as in fn_plot_somatic_dist
  SOMATIC_VARIANTS |>
    dplyr::mutate(srridvariant = glue("{srrid}_{variant}")) |>
    dplyr::select(
      srridvariant,
      c(B, CD4_T, CD8_T, DC, Mono, NK, other, Bulk, other_T)
    ) |>
    tibble::column_to_rownames("srridvariant") |>
    as.matrix() -> mat_somatic_af

  mat_specific <- mat_somatic_af[, c(
    "B",
    "CD4_T",
    "CD8_T",
    "other_T",
    "NK",
    "DC",
    "Mono",
    "other"
  )]

  dominant_celltype <- apply(mat_specific, 1, function(x) {
    if (max(x, na.rm = TRUE) < 0.01) {
      return("Bulk")
    }
    colnames(mat_specific)[which.max(x)]
  })

  celltype_order <- c(
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

  variant_df <- data.frame(
    srridvariant = rownames(mat_somatic_af),
    dominant_celltype = dominant_celltype,
    max_af = apply(mat_specific, 1, max, na.rm = TRUE),
    stringsAsFactors = FALSE
  )
  variant_df$dominant_celltype <- factor(
    variant_df$dominant_celltype,
    levels = celltype_order
  )
  variant_ordered <- variant_df[
    order(variant_df$dominant_celltype, -variant_df$max_af),
  ]
  clustered_order <- variant_ordered$srridvariant

  # Get gene annotation

  gtf_gene_df <- import(
    path(Sys.getenv("REPODIR"), "config/mtdna_genes_dloop.qs")
  )

  # Prepare haplogroup colors
  SOMATIC_VARIANTS |>
    dplyr::select(Haplogroup, Verbose_haplogroup) |>
    dplyr::distinct() |>
    dplyr::filter(!is.na(Haplogroup)) |>
    dplyr::mutate(
      Haplogroup_s = gsub("\\d+.*", "", Haplogroup)
    ) |>
    dplyr::mutate(
      color_haplogroup = color_haplogroup[Haplogroup_s]
    ) |>
    dplyr::mutate(
      color_verbose_haplogroup = ifelse(
        Haplogroup == Verbose_haplogroup,
        color_haplogroup,
        as.character(prismatic::clr_lighten(color_haplogroup, 0.5))
      )
    ) |>
    dplyr::select(-Haplogroup_s) -> haplo_colors_df

  # Create final table data
  SOMATIC_VARIANTS |>
    dplyr::mutate(srridvariant = glue("{srrid}_{variant}")) |>
    dplyr::mutate(
      pos = as.integer(gsub("(\\d+)[A-Z]>[A-Z]", "\\1", variant))
    ) |>
    dplyr::mutate(
      gene_info = purrr::map(
        .x = pos,
        .f = \(.x) {
          gtf_gene_df |>
            dplyr::filter(start <= .x & end >= .x) |>
            dplyr::select(gene_name, TYPE, COLOR) |>
            dplyr::slice(1)
        }
      )
    ) |>
    tidyr::unnest(cols = gene_info) |>
    dplyr::left_join(
      haplo_colors_df,
      by = c("Haplogroup", "Verbose_haplogroup")
    ) |>
    dplyr::mutate(
      srridvariant = factor(srridvariant, levels = clustered_order)
    ) |>
    dplyr::arrange(srridvariant) |>
    dplyr::select(
      srrid,
      variant,
      Haplogroup,
      Verbose_haplogroup,
      color_haplogroup,
      color_verbose_haplogroup,
      TYPE,
      gene_name,
      COLOR,
      Bulk,
      B,
      CD4_T,
      CD8_T,
      other_T,
      NK,
      DC,
      Mono,
      other
    ) -> table_data_

  table_data_ |>
    dplyr::select(
      color_haplogroup,
      color_verbose_haplogroup,
      COLOR
    ) -> table_data_colors

  table_data_ |>
    dplyr::select(
      Sample = srrid,
      Haplogroup,
      Haplogroup_V = Verbose_haplogroup,
      Variant = variant,
      `Gene type` = TYPE,
      `Gene name` = gene_name,
      Bulk:other
    ) -> table_data
  table_data_colnames <- colnames(table_data)

  # ========== Save to Excel with formatting ==========
  wb <- openxlsx::createWorkbook()
  openxlsx::addWorksheet(wb, "Somatic Variants")

  # Write data
  openxlsx::writeData(wb, "Somatic Variants", table_data)

  # Header style
  header_style <- openxlsx::createStyle(
    textDecoration = "bold",
    halign = "center",
    valign = "center",
    border = "Bottom",
    borderColour = "black"
  )

  openxlsx::addStyle(
    wb,
    "Somatic Variants",
    header_style,
    rows = 1,
    cols = 1:ncol(table_data),
    gridExpand = TRUE
  )

  # Apply cell background colors for haplogroups
  which(table_data_colnames == "Haplogroup") -> haplogroup_col_idx
  for (i in seq_len(nrow(table_data))) {
    # Haplogroup column
    hg_color <- table_data_colors$color_haplogroup[i]
    if (!is.na(hg_color) && nchar(hg_color) > 0 && grepl("^#", hg_color)) {
      hg_style <- openxlsx::createStyle(
        fgFill = gsub("FF$", "", hg_color),
        halign = "center",
        fontColour = "black"
      )
      openxlsx::addStyle(
        wb,
        "Somatic Variants",
        hg_style,
        rows = i + 1,
        cols = haplogroup_col_idx
      )
    }

    # Haplogroup_V column
    vhg_color <- table_data_colors$color_verbose_haplogroup[i]
    which(table_data_colnames == "Haplogroup_V") -> verbose_haplogroup_col_idx
    if (!is.na(vhg_color) && nchar(vhg_color) > 0 && grepl("^#", vhg_color)) {
      vhg_style <- openxlsx::createStyle(
        fgFill = gsub("FF$", "", vhg_color),
        halign = "center",
        fontColour = "black"
      )
      openxlsx::addStyle(
        wb,
        "Somatic Variants",
        vhg_style,
        rows = i + 1,
        cols = verbose_haplogroup_col_idx
      )
    }

    # Gene type and Gene name columns
    type_color <- table_data_colors$COLOR[i]
    which(table_data_colnames == "Gene type") -> type_col_idx
    which(table_data_colnames == "Gene name") -> gene_col_idx
    if (
      !is.na(type_color) && nchar(type_color) > 0 && grepl("^#", type_color)
    ) {
      type_style <- openxlsx::createStyle(
        fgFill = gsub("FF$", "", type_color),
        halign = "center",
        fontColour = "black"
      )
      openxlsx::addStyle(
        wb,
        "Somatic Variants",
        type_style,
        rows = i + 1,
        cols = type_col_idx
      )
      openxlsx::addStyle(
        wb,
        "Somatic Variants",
        type_style,
        rows = i + 1,
        cols = gene_col_idx
      )
    }

    # AF columns (columns celltypes) - gradient from white to red
    af_cols <- c(
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

    for (celltype in af_cols) {
      which(table_data_colnames == celltype) -> celltype_col_idx
      af_val <- table_data[[celltype]][i]
      if (!is.na(af_val) && af_val > 0) {
        # Create gradient: white (0) to red (max ~0.5)
        intensity <- min(af_val / 0.3, 1)
        r <- 255
        g <- round(255 * (1 - intensity))
        b <- round(255 * (1 - intensity))
        hex_color <- sprintf("#%02X%02X%02X", r, g, b)
        af_style <- openxlsx::createStyle(
          fgFill = hex_color,
          halign = "center",
          numFmt = "0.00"
        )
        openxlsx::addStyle(
          wb,
          "Somatic Variants",
          af_style,
          rows = i + 1,
          cols = celltype_col_idx
        )
      }
    }
  }

  # Set column widths
  openxlsx::setColWidths(
    wb,
    "Somatic Variants",
    cols = 1:ncol(table_data),
    widths = "auto"
  )

  # Save Excel
  excel_path <- path(
    Sys.getenv("OUTDIRNOTUSE")
  ) /
    "HOTSPOTS-SOMATIC-CLUSTER.xlsx"
  openxlsx::saveWorkbook(wb, excel_path, overwrite = TRUE)
  log_info("Saved Excel to: {excel_path}")

  # ========== Create gt table ==========
  table_data |>
    gt::gt() |>
    gt::tab_header(
      title = gt::md("**Somatic Variants Distribution**"),
      subtitle = gt::md(
        "*mtDNA somatic variants across cell types and haplogroups*"
      )
    ) |>
    gt::tab_spanner(
      label = gt::md("**Haplogroup**"),
      columns = c(Haplogroup, Haplogroup_V)
    ) |>
    gt::tab_spanner(
      label = gt::md("**Gene Annotation**"),
      columns = c(`Gene type`, `Gene name`)
    ) |>
    gt::tab_spanner(
      label = gt::md("**Allele Frequency**"),
      columns = c(Bulk, B, CD4_T, CD8_T, other_T, NK, DC, Mono, other)
    ) |>
    gt::fmt_number(
      columns = c(Bulk, B, CD4_T, CD8_T, other_T, NK, DC, Mono, other),
      decimals = 3
    ) |>
    gt::cols_align(align = "center", columns = gt::everything()) |>
    gt::cols_label(
      Sample = "Sample",
      Variant = "Variant",
      Haplogroup = "Haplogroup",
      Haplogroup_V = "Verbose",
      `Gene type` = "Type",
      `Gene name` = "Gene"
    ) -> gt_table

  # Apply haplogroup colors row by row
  for (i in seq_len(nrow(table_data))) {
    hg <- table_data$Haplogroup[i]
    hg_color <- table_data_colors$color_haplogroup[i]
    vhg_color <- table_data_colors$color_verbose_haplogroup[i]
    type_color <- table_data_colors$COLOR[i]

    if (!is.na(hg_color)) {
      gt_table <- gt_table |>
        gt::tab_style(
          style = gt::cell_fill(color = hg_color),
          locations = gt::cells_body(columns = Haplogroup, rows = i)
        )
    }
    if (!is.na(vhg_color)) {
      gt_table <- gt_table |>
        gt::tab_style(
          style = gt::cell_fill(color = vhg_color),
          locations = gt::cells_body(columns = c(Haplogroup_V), rows = i)
        )
    }
    if (!is.na(type_color)) {
      gt_table <- gt_table |>
        gt::tab_style(
          style = gt::cell_fill(color = type_color),
          locations = gt::cells_body(
            columns = c(`Gene type`, `Gene name`),
            rows = i
          )
        )
    }
  }

  # Apply AF gradient colors
  af_cols <- c(
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
  for (col in af_cols) {
    gt_table <- gt_table |>
      gt::data_color(
        columns = !!rlang::sym(col),
        palette = c("white", "#FFCCCC", "#FF6666", "#FF0000"),
        domain = c(0, 0.5),
        na_color = "white"
      )
  }

  # Style the table
  gt_table <- gt_table |>
    gt::tab_options(
      table.font.size = gt::px(12),
      heading.title.font.size = gt::px(18),
      heading.subtitle.font.size = gt::px(14),
      column_labels.font.weight = "bold",
      table.border.top.width = gt::px(2),
      table.border.bottom.width = gt::px(2),
      heading.border.bottom.width = gt::px(2)
    ) |>
    gt::opt_row_striping() |>
    gt::opt_table_outline()

  # Save gt table as HTML
  html_path <- path(
    Sys.getenv("OUTDIRNOTUSE")
  ) /
    "HOTSPOTS-SOMATIC-CLUSTER.html"
  gt::gtsave(gt_table, html_path)
  log_info("Saved gt table HTML to: {html_path}")

  return(list(excel = excel_path, html = html_path, data = table_data))
}

fn_save_somatic_table()
