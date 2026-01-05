dotenv(".env")

fn_umap_coord <- function(.x) {
  .col_names <- c("UMAP_1", "UMAP_2")

  if ("ref.umap" %in% names(.x@reductions)) {
    .umap <- .x@reductions$ref.umap@cell.embeddings |>
      data.table::as.data.table()
    colnames(.umap) <- .col_names
    .tsne <- NULL
  } else {
    .umap <- .x@reductions$umap@cell.embeddings |> data.table::as.data.table()
    colnames(.umap) <- .col_names
    .tsne <- .x@reductions$tsne@cell.embeddings |> data.table::as.data.table()
    colnames(.tsne) <- .col_names
  }

  # .umap
  .x@meta.data |>
    dplyr::select(
      celltype
    ) |>
    data.table::as.data.table() -> .xx

  .xxx <- dplyr::bind_cols(.umap, .xx) |>
    dplyr::mutate(barcode = rownames(.x@meta.data))

  .xxx
}

fn_load_hetero <- function(.filename) {
  # .filename <- file.path(
  #   sc_dir,
  #   "mgatk_out/final",
  #   "sc.cell_heteroplasmic_df.tsv.gz"
  # )

  data.table::fread(input = .filename) -> .d

  data.table::setnames(.d, "V1", "barcode")
  .d <- data.table::melt(
    .d,
    id.vars = "barcode",
    variable.name = "variant",
    value.name = "af"
  )
  .d[, pos := gsub(pattern = ">|[AGCT]", replacement = "", x = variant)]
  .d[, pos := as.integer(pos)]
  .d
}

fn_load_coverage <- function(.filename) {
  data.table::fread(
    input = .filename,
    sep = ",",
    col.names = c("pos", "barcode", "depth")
  ) -> .d
  .d
}

fn_plot_vaf_featureplot <- function(.thevariant, sc, .cell_annotation = NULL) {
  pcc <- readr::read_tsv(
    file = "https://raw.githubusercontent.com/chunjie-sam-liu/chunjie-sam-liu.life/master/public/data/pcc.tsv"
  ) |>
    dplyr::arrange(cancer_types)
  source(path(
    Sys.getenv("HIGHRESDIR"),
    "00-colors.R"
  ))
  sc$cell_hetero_coverage |>
    dplyr::filter(variant == .thevariant) -> vhc

  sc$umap_coord |>
    dplyr::left_join(vhc, by = "barcode") -> vhc_umap

  if (is.null(.cell_annotation)) {
    vhc_umap |>
      dplyr::arrange(af) |>
      ggplot(aes(x = UMAP_1, y = UMAP_2)) +
      geom_point(aes(color = af)) +
      scale_color_gradient2(
        name = "AF",
        low = "grey",
        mid = "grey50",
        high = "#F02415",
        midpoint = 0.01,
      ) +
      theme_bw() +
      labs(
        title = .thevariant
      ) +
      theme(
        plot.title = element_text(
          color = "black",
          face = "bold",
          hjust = 0.5
        )
      ) -> p_feature
    return(p_feature)
  } else {
    vhc_umap |>
      dplyr::arrange(af) |>
      ggplot(aes(x = UMAP_1, y = UMAP_2)) +
      geom_point(aes(color = celltype)) +
      scale_color_manual(
        name = "Cell Type",
        # values = pcc$color
        values = color_celltype
      ) +
      theme_bw() +
      labs() +
      theme(
        plot.title = element_blank()
      ) -> p_annotation
    return(p_annotation)
  }
}


fn_load_by_path <- function(thepath) {
  library(Seurat)
  azimuth_file <- file.path(thepath, "sc_azimuth.rds.gz")
  cell_hetero_file <- file.path(
    thepath,
    "cell.cell_heteroplasmic_df_raw.tsv.gz"
  )
  cell_coverage_file <- file.path(thepath, "cell.coverage.txt.gz")

  sc <- readr::read_rds(azimuth_file)
  sc$umap_coord <- fn_umap_coord(.x = sc$sc_azimuth)
  sc$cell_hetero <- fn_load_hetero(cell_hetero_file)
  sc$cell_coverage <- fn_load_coverage(cell_coverage_file)
  sc$cell_hetero_coverage <- sc$cell_hetero |>
    dplyr::left_join(sc$cell_coverage, by = c("barcode", "pos"))

  sc
}

fn_plot_vaf_featureplot_multi <- function(.thevariants, sc) {
  purrr::map(
    .thevariants,
    fn_plot_vaf_featureplot,
    sc = sc
  ) -> variant_plots

  fn_plot_vaf_featureplot(
    .thevariants[[1]],
    sc,
    .cell_annotation = TRUE
  ) -> p_annotation

  c(variant_plots, list(p_annotation)) |>
    wrap_plots(ncol = 4) +
    guide_area() +
    plot_layout(guides = "collect") &
    theme(
      legend.direction = "vertical",
      legend.box = "horizontal",
      # plot.background = element_rect(color = "black"),
    )
}

fn_load_count <- function(thepath, type = c("cluster", "cell")) {
  type <- match.arg(type)

  pattern <- if (type == "cluster") {
    "*cluster.*.txt.gz*"
  } else {
    "*cell.*.txt.gz*"
  }

  tibble::tibble(
    path = list.files(
      thepath,
      pattern,
      full.names = T
    )
  ) |>
    dplyr::filter(!grepl("coverage", x = path)) |>
    dplyr::mutate(d = purrr::map(path, data.table::fread)) |>
    dplyr::mutate(n = basename(path)) |>
    dplyr::mutate(n = gsub(paste0(type, ".|.txt.gz"), "", n)) |>
    dplyr::select(n, d) |>
    tidyr::unnest(cols = d) |>
    as.data.table() |>
    dplyr::mutate(nv = V3 + V4) |>
    dplyr::select(
      gt = n,
      pos = V1,
      group = V2,
      fw = V3,
      rv = V4,
      nv
    ) -> cluster_n

  fasta <- Biostrings::readDNAStringSet(
    "/home/liuc9/github/scMOCHA/fasta/rCRS.chrM.fasta"
  )

  fasta$chrM |>
    as.data.table() |>
    tibble::rownames_to_column(var = "pos") |>
    dplyr::rename(ref = x) |>
    dplyr::mutate(posref = glue::glue("{pos}{ref}")) |>
    dplyr::mutate(pos = as.integer(pos)) -> fasta_df
  # data.table::fwrite(
  #   fasta_df,
  #   file = "/home/liuc9/github/scMOCHA-data/config/rCRS.MT.fasta.csv",
  #   sep = ","
  # )

  cluster_n |>
    dtplyr::lazy_dt() |>
    dplyr::left_join(fasta_df, by = "pos") |>
    dplyr::mutate(gt = factor(gt, levels = c("A", "G", "C", "T"))) |>
    as.data.table() -> cluster_n_temp

  cluster_n_temp[, ratio := nv / sum(nv), by = .(group, pos)]

  cluster_n_temp |>
    dplyr::mutate(
      label = glue::glue(
        "total = {nv} \n forward = {fw} \n reverse = {rv} \n ratio = {round(ratio, 3) * 100}%"
      )
    ) -> cluster_n_forplot

  cluster_n_forplot
}

fn_plot_count <- function(cluster_n_forplot, thepos, group_sel = NA) {
  pcc <- readr::read_tsv(
    file = "https://raw.githubusercontent.com/chunjie-sam-liu/chunjie-sam-liu.life/master/public/data/pcc.tsv"
  ) |>
    dplyr::arrange(cancer_types)

  if (!all(is.na(group_sel))) {
    cluster_n_forplot |>
      dplyr::filter(group %in% group_sel) -> cluster_n_forplot
  }

  if (length(unique(cluster_n_forplot$group)) > 10) {
    # stop("The number of unique groups exceeds 50.")
    cluster_n_forplot |>
      dplyr::filter(pos == thepos) |>
      dplyr::arrange(-nv) |>
      dplyr::slice(1:10) -> cluster_n_forplot
  }

  cluster_n_forplot |>
    dplyr::filter(pos %in% thepos) |>
    dplyr::mutate(pos = as.character(pos)) -> cluster_n_forplot_

  gt <- factor(c("A", "G", "C", "T"), levels = c("A", "G", "C", "T"))
  posref = cluster_n_forplot_$posref |> unique()
  group = cluster_n_forplot_$group |> unique()

  posref_df <- data.table::data.table(
    gt = rep(gt, each = length(posref)),
    posref = rep(posref, length(gt)),
    group = rep(group, each = length(gt) * length(posref))
  )

  cluster_n_forplot_ |>
    dplyr::full_join(posref_df, by = c("gt", "posref", "group")) |>
    ggplot(aes(x = posref, y = gt)) +
    geom_tile(aes(fill = ratio)) +
    geom_text(aes(label = label)) +
    scale_fill_gradient(
      low = "white",
      high = "red",
      na.value = "white"
    ) +
    ggh4x::facet_wrap2(
      ~group,
      ncol = 8,
      strip.position = "top",
      strip = ggh4x::strip_themed(
        background_x = ggh4x::elem_list_rect(
          fill = pcc$color
        ),
        text_x = ggh4x::elem_list_text(
          colour = "white",
          face = c("bold")
        ),
        by_layer_y = FALSE,
      )
    ) +
    scale_x_discrete(
      expand = expansion(mult = c(0.01, 0.01))
    ) +
    scale_y_discrete(
      expand = expansion(mult = c(0, 0))
    ) +
    theme(
      panel.background = element_blank(),
      panel.grid = element_blank(),
      # axis.ticks = element_blank(),
      axis.title = element_blank(),
      axis.text = element_text(
        color = "black",
        size = 18
      ),
      legend.position = "none ",
      plot.title = element_text(
        size = 16,
        hjust = 0.5
      ),
      strip.background = element_rect(
        fill = NA,
        color = "black",
      ),
      strip.text = element_text(
        color = "black",
        size = 14,
        face = "bold"
      ),
      axis.line = element_line(
        color = "black"
      )
    ) -> p_tile
  p_tile
}

fn_plot_count_multi <- function(cluster_n_forplot, theposes, group_sel = NA) {
  theposes |>
    purrr::map(
      ~ fn_plot_count(cluster_n_forplot, thepos = ., group_sel = NA)
    ) |>
    wrap_plots(ncol = 1) +
    plot_layout(guides = "collect") &
    theme(
      legend.justification = "left",
      legend.position = "none",
    ) -> p_count
  p_count
}

fn_plot_mtdna <- function() {
  mt_exons_df <- "/home/liuc9/github/scMOCHA/fasta/mt_exons.df.rds.gz"

  gtf_gene_df <-
    readr::read_rds(
      file = mt_exons_df
    )
  library(gggenes)
  ggplot(gtf_gene_df, aes(xmin = start, xmax = end, y = seqnames)) +
    # geom_gene_arrow() +
    geom_gene_arrow(
      aes(
        fill = gene_biotype
      ),
      arrowhead_height = unit(3, "mm"),
      arrowhead_width = unit(1, "mm")
    ) +
    scale_fill_brewer(
      palette = "Set1",
      name = "Gene type",
      labels = c("MT rRNA", "MT tRNA", "Protein coding")
    ) +
    ggrepel::geom_text_repel(
      aes(x = (start + end) / 2, label = gene_name, color = gene_biotype),
      # fill = "white",
      # nudge_x =1,
      # nudge_y = -0.1,
      size = 3,
      show.legend = F,
      max.overlaps = Inf,
    ) +
    scale_color_brewer(palette = "Set1") +
    scale_x_continuous(
      limits = c(0, 17000),
      breaks = seq(0, 17000, 1000),
      labels = seq(0, 17000, 1000),
      expand = expansion(mult = c(0, 0.01)),
    ) +
    scale_y_discrete(
      expand = expansion(mult = c(0, 0), add = c(0, 0))
    ) +
    # theme_genes() +
    theme(
      legend.position = "bottom",
      axis.title = element_blank(),
      axis.text.y = element_blank(),
      # axis.text.x = element_text(size = 14),
      # legend.text = element_text(size = 14),
      panel.background = element_blank(),
      panel.grid = element_blank(),
      axis.ticks.y = element_blank(),
      axis.line.x = element_line(color = "black"),
      axis.text.x = element_text(
        vjust = -1,
      ),
    ) +
    coord_cartesian(xlim = c(0, 17000)) -> pg
  pg
}

fn_plot_mtdna <- function() {
  # mt_exons_df <- "/home/liuc9/github/scMOCHA/fasta/mt_exons.df.rds.gz"

  LENGTH <- 16569
  rCRS <- Biostrings::readDNAStringSet(
    "/home/liuc9/github/scMOCHA-data/config/rCRS.MT.fasta"
  )
  gtf_gene_df <- readr::read_rds(
    "/home/liuc9/github/scMOCHA-data/config/mtdna_genes_dloop.rds.gz"
  )

  library(gggenes)
  ggplot(gtf_gene_df, aes(xmin = start, xmax = end, y = seqnames)) +
    # geom_gene_arrow() +
    geom_gene_arrow(
      aes(
        fill = TYPE
      ),
      arrowhead_height = unit(3, "mm"),
      arrowhead_width = unit(1, "mm")
    ) +
    scale_fill_brewer(
      palette = "Set1",
      name = "Gene type",
      labels = c("D-Loop", "MT rRNA", "MT tRNA", "Protein coding")
    ) +
    ggrepel::geom_text_repel(
      aes(
        x = (start + end) / 2,
        label = gene_name,
      ),
      color = "black",
      # fill = "white",
      # nudge_x =1,
      # nudge_y =0.001,
      size = 3,
      show.legend = F,
      max.overlaps = Inf,
    ) +
    scale_color_brewer(palette = "Set1") +
    scale_x_continuous(
      limits = c(0, LENGTH),
      breaks = c(seq(0, LENGTH, 1000), LENGTH),
      labels = c(seq(0, LENGTH, 1000), LENGTH),
      expand = expansion(mult = c(0, 0.01)),
    ) +
    scale_y_discrete(
      expand = expansion(mult = c(0, 0), add = c(0, 0))
    ) +
    # theme_genes() +
    theme(
      legend.position = "bottom",
      axis.title = element_blank(),
      axis.text.y = element_blank(),
      # axis.text.x = element_text(size = 14),
      # legend.text = element_text(size = 14),
      panel.background = element_blank(),
      panel.grid = element_blank(),
      axis.ticks.y = element_blank(),
      axis.ticks.x = element_line(color = "black"),
      axis.line.x = element_line(color = "black"),
      axis.text.x = element_text(
        vjust = -1,
      ),
    )
}

fn_plot_coverage <- function(thepath, theposes = NULL) {
  .cluster_coverage <- fn_load_coverage(file.path(
    thepath,
    "cluster.coverage.txt.gz"
  ))
  pcc <- readr::read_tsv(
    file = "https://raw.githubusercontent.com/chunjie-sam-liu/chunjie-sam-liu.life/master/public/data/pcc.tsv"
  ) |>
    dplyr::arrange(cancer_types)

  if (!is.null(theposes)) {
    .plot_vline <- geom_vline(
      xintercept = theposes,
      linetype = "dashed",
      color = "red",
      linetype = 21
    )
  } else {
    .plot_vline <- NULL
  }

  .cluster_coverage |>
    ggplot(aes(x = pos, y = depth, fill = barcode)) +
    geom_bar(stat = "identity", show.legend = FALSE) +
    scale_x_continuous(
      limits = c(0, 17000),
      breaks = seq(0, 17000, 1000),
      labels = seq(0, 17000, 1000),
      expand = expansion(mult = c(0, 0.01)),
    ) +
    .plot_vline +
    scale_y_continuous(
      expand = c(0.01, 0),
      # limits = c(0, 520000),
      label = scales::label_number()
    ) +
    scale_fill_manual(
      name = "Cell type",
      values = pcc$color
    ) +
    ggh4x::facet_wrap2(
      ~barcode,
      ncol = 1,
      strip.position = "right",
      strip = ggh4x::strip_themed(
        background_y = ggh4x::elem_list_rect(
          fill = pcc$color
        ),
        text_y = ggh4x::elem_list_text(
          colour = "white",
          face = c("bold")
        ),
        by_layer_y = FALSE,
      )
    ) +
    theme(
      plot.margin = margin(t = 0, b = 0, unit = "cm"),
      panel.background = element_blank(),
      panel.grid = element_blank(),
      axis.line.y.left = element_line(color = "black"),
      # axis.line.x.bottom = element_line(color = "black"),
      axis.ticks.x = element_blank(),
      axis.text.x = element_blank(),
      axis.line.x = element_blank(),
      axis.title.x = element_blank(),
      legend.position = c(0.8, 0.5),
      legend.key = element_blank(),
      axis.title.y = element_text(color = "black"),
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
    coord_cartesian(xlim = c(0, 17000)) +
    labs(y = "Depth") -> p_mt_depth_celltype

  .cluster_coverage |>
    dplyr::group_by(pos) |>
    dplyr::summarise(depth = sum(depth, na.rm = T)) |>
    ggplot(aes(x = pos, y = depth)) +
    geom_bar(fill = "grey", stat = "identity", show.legend = FALSE) +
    .plot_vline +
    scale_x_continuous(
      limits = c(0, 17000),
      breaks = seq(0, 17000, 1000),
      labels = seq(0, 17000, 1000),
      expand = expansion(mult = c(0, 0.01)),
    ) +
    scale_y_continuous(
      expand = c(0.01, 0),
      # limits = c(0, 520000),
      label = scales::label_number()
    ) +
    theme(
      plot.margin = margin(t = 0, b = 0, unit = "cm"),
      panel.background = element_blank(),
      panel.grid = element_blank(),
      axis.line.y.left = element_line(color = "black"),
      # axis.line.x.bottom = element_line(color = "black"),
      axis.ticks.x = element_blank(),
      axis.text.x = element_blank(),
      axis.line.x = element_blank(),
      axis.title.x = element_blank(),
      legend.position = c(0.8, 0.5),
      legend.key = element_blank(),
      axis.title.y = element_text(color = "black"),
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
    coord_cartesian(xlim = c(0, 17000)) +
    labs(y = "Depth") -> p_mt_depth_allcell

  list(
    p_mt_depth_celltype = p_mt_depth_celltype,
    p_mt_depth_allcell = p_mt_depth_allcell
  )
}

fn_load_meta <- function(.filename) {
  data.table::fread(
    input = .filename,
    sep = "\t"
  ) |>
    dplyr::rename(
      barcode = cellbarcode
    ) |>
    dplyr::select(-orig.ident)
}

fn_load_cluster <- function(.filename) {
  data.table::fread(
    input = .filename,
    sep = "\t",
    col.names = c("barcode", "tag", "celltype")
  ) |>
    dplyr::arrange(celltype) |>
    dplyr::mutate(celltype = factor(celltype)) |>
    dplyr::select(-tag)
}

fn_forheatmap_plot <- function(.af, .coverage, .meta) {
  # print(.meta)
  .af |>
    dplyr::select(barcode, cluster, dplyr::contains(">")) |>
    tidyr::pivot_longer(
      cols = -c(barcode, cluster),
      names_to = "variant",
      values_to = "af"
    ) |>
    dplyr::group_by(barcode, cluster) |>
    dplyr::summarise(s_af = sum(af, na.rm = T)) |>
    dplyr::ungroup() |>
    dplyr::arrange(cluster, -s_af) -> .rank

  .af |>
    dplyr::select(barcode, dplyr::contains(">")) |>
    tidyr::pivot_longer(
      cols = -barcode,
      names_to = "variant",
      values_to = "af"
    ) |>
    dplyr::mutate(
      pos = gsub(pattern = "([[:digit:]]*).*", "\\1", variant) |>
        as.numeric()
    ) |>
    dplyr::left_join(
      .coverage,
      by = c("barcode", "pos")
    ) |>
    tidyr::replace_na(
      replace = list(
        af = 0
      )
    ) |>
    dplyr::mutate(af = ifelse(is.na(depth), NA, af)) |>
    # dplyr::mutate(af = ifelse(depth < log2(10), -0.1, af)) |>
    dplyr::mutate(af = ifelse(depth < 10, -0.1, af)) |>
    dplyr::arrange(pos) -> .forplot

  .coverage |>
    dplyr::group_by(barcode) |>
    dplyr::summarise(sum_depth = sum(depth, na.rm = TRUE)) -> .coverage_cell

  list(
    rank = .rank,
    forplot = .forplot,
    meta = .meta,
    coverage_cell = .coverage_cell
  )
}

fn_plot_hotspots <- function(thepath, thevariants = NULL) {
  pcc <- readr::read_tsv(
    file = "https://raw.githubusercontent.com/chunjie-sam-liu/chunjie-sam-liu.life/master/public/data/pcc.tsv"
  ) |>
    dplyr::arrange(cancer_types)

  .cell_anno <- data.table::fread(
    ifelse(
      file.exists(file.path(thepath, "variant_annotation.tsv")),
      file.path(thepath, "variant_annotation.tsv"),
      file.path(thepath, "cell_variant_annotation.tsv")
    )
  ) |>
    dplyr::mutate(variant = glue::glue("{Position}{Ref}>{Alt}"))

  cluster_umap <- fn_load_cluster(
    .filename = file.path(thepath, "barcode_cluster.tsv")
  )

  cell_hetero_raw <- fn_load_hetero(
    .filename = file.path(thepath, "cell.cell_heteroplasmic_df_raw.tsv.gz")
  ) |>
    dplyr::filter(af > 0.05) |>
    dplyr::filter(variant %in% .cell_anno$variant)

  cell_raw_cluster_af <- cluster_umap |>
    dplyr::left_join(cell_hetero_raw, by = "barcode") |>
    dplyr::rename(cluster = celltype) |>
    tidyr::pivot_wider(
      names_from = variant,
      values_from = af
    )

  cell_coverage <- fn_load_coverage(
    .filename = file.path(thepath, "cell.coverage.txt.gz")
  )

  metadata <- fn_load_meta(
    .filename = file.path(thepath, "cell_meta_data.tsv")
  )

  .forplot <- fn_forheatmap_plot(
    .af = cell_raw_cluster_af,
    .coverage = cell_coverage,
    .meta = metadata
  )

  if (!is.null(thevariants)) {
    .cell_anno <- .cell_anno |>
      dplyr::filter(variant %in% thevariants)
  }

  .forplot$forplot |>
    dplyr::left_join(
      .forplot$rank |> dplyr::select(-s_af),
      by = "barcode"
    ) |>
    dplyr::filter(variant %in% .cell_anno$variant) |>
    dplyr::mutate(
      af = ifelse(af == 0, NA_real_, af)
    ) |>
    dplyr::mutate(
      af = ifelse(depth < 10, NA_real_, af)
    ) |>
    dplyr::mutate(
      depth_log2 = log2(depth + 1)
    ) -> .forplot_cluster_cell_variant

  .forplot_cluster_cell_variant |>
    dplyr::filter(af > 0) -> .theforplot

  .theforplot |>
    dplyr::select(variant, pos) |>
    dplyr::distinct() |>
    dplyr::arrange(pos) -> .sort_variant

  .forplot_cluster_cell_variant |>
    dplyr::group_by(cluster, variant) |>
    dplyr::summarise(
      mean_cluster_variant_af = mean(af, na.rm = T),
      sum_cluster_variant_depth = sum(depth, na.rm = T),
      max_cluster_variant_depth = max(depth, na.rm = T)
    ) |>
    dplyr::ungroup() |>
    dplyr::mutate(
      sum_cluster_variant_depth_log2 = log2(sum_cluster_variant_depth + 1)
    ) |>
    dplyr::filter(variant %in% .sort_variant$variant) -> .cluster_variant_af

  .c <- unique(.cluster_variant_af$cluster)
  .v <- unique(.cluster_variant_af$variant)

  tibble::tibble(
    cluster = rep(.c, each = length(.v)),
    variant = rep(.v, length(.c))
  ) |>
    dplyr::left_join(
      .cluster_variant_af,
      by = c("cluster", "variant")
    ) |>
    dplyr::mutate(
      rect_type = ifelse(max_cluster_variant_depth >= 10, "white", "grey")
    ) |>
    dplyr::mutate(
      variant = factor(variant, .sort_variant$variant)
    ) -> .no_depth

  .cell_anno |>
    dplyr::filter(variant %in% .sort_variant$variant) |>
    dplyr::mutate(fill = ifelse(!is.na(Haplogroup), "#3B0049", "white")) |>
    dplyr::mutate(color = ifelse(!is.na(Haplogroup), "white", "black")) |>
    dplyr::mutate(
      variant = factor(variant, .sort_variant$variant)
    ) |>
    dplyr::arrange(variant) -> .haplo_variant

  .theforplot |>
    dplyr::inner_join(
      .cluster_variant_af,
      by = c("cluster", "variant")
    ) |>
    dplyr::mutate(
      variant = factor(variant, .sort_variant$variant |> unique())
    ) |>
    dplyr::arrange(variant) -> .haplo_forplot

  library(ggh4x)
  library(ggbeeswarm)
  library(ggnewscale)

  .cl <- as.character(unique(.haplo_forplot$cluster)[[1]])
  .haplo_forplot |>
    dplyr::select(variant, pos, cluster, mean_cluster_variant_af) |>
    dplyr::distinct() |>
    dplyr::filter(cluster == .cl) -> .forlabel

  ggplot() +
    ggrepel::geom_text_repel(
      data = .forlabel,
      aes(
        x = pos,
        y = mean_cluster_variant_af,
        label = variant,
      ),
      size = 3,
      # direction = "y",
      # nudge_x = -0.2,
      direction = "y",
      # hjust = "right",
      segment.size = 0,
    ) +
    scale_x_continuous(
      limits = c(0, 17000),
      breaks = seq(0, 17000, 1000),
      labels = seq(0, 17000, 1000),
      expand = expansion(mult = c(0, 0.01)),
    ) +
    theme(
      plot.margin = margin(t = 0, b = 0, unit = "cm"),
      panel.background = element_blank(),
      panel.grid = element_blank(),
      # axis.line.y.left = element_line(color = "black"),
      axis.ticks.x = element_blank(),
      axis.text.x = element_blank(),
      axis.line.x = element_blank(),
      axis.title.x = element_blank(),
      legend.position = "right",
      legend.key = element_blank(),
      axis.title.y = element_blank(),
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      ,
      axis.line.y = element_blank(),
    ) +
    coord_cartesian(xlim = c(0, 17000)) -> p_label

  .haplo_forplot |>
    ggplot() +
    ggh4x::facet_wrap2(
      ~cluster,
      ncol = 1,
      strip.position = "right",
      strip = ggh4x::strip_themed(
        background_y = ggh4x::elem_list_rect(
          fill = pcc$color
        ),
        text_y = ggh4x::elem_list_text(
          colour = "white",
          face = c("bold")
        ),
        by_layer_y = FALSE,
      )
    ) +
    ggbeeswarm::geom_quasirandom(
      aes(
        x = pos,
        y = af,
        color = af
      ),
      size = 1,
      dodge.width = .75,
      alpha = .5,
      varwidth = TRUE
    ) +
    scale_color_gradient2(
      name = "AF",
      low = "white",
      mid = "red",
      high = "#3B0049",
      midpoint = 0.5,
    ) +
    scale_x_continuous(
      limits = c(0, 17000),
      breaks = seq(0, 17000, 1000),
      labels = seq(0, 17000, 1000),
      expand = expansion(mult = c(0, 0.01)),
    ) +
    scale_y_continuous(
      expand = c(0.01, 0),
      limits = c(0, 1),
    ) +
    theme(
      plot.margin = margin(t = 0, b = 0, unit = "cm"),
      panel.background = element_blank(),
      panel.grid = element_blank(),
      axis.line.y.left = element_line(color = "black"),
      # axis.line.x.bottom = element_line(color = "black"),
      axis.ticks.x = element_blank(),
      axis.text.x = element_blank(),
      axis.line.x = element_blank(),
      axis.title.x = element_blank(),
      legend.position = "right",
      legend.key = element_blank(),
      axis.title.y = element_text(color = "black"),
      # axis.title.y = element_blank(),
      axis.text.y = element_text(color = "black"),
      # legend.text = element_text(
      #   size = 14,
      #   color = "black"
      # ),
      # legend.title = element_text(
      #   size = 16,
      #   colour = "black"
      # ),
      # strip.background = element_blank(),
      # strip.text = element_text(
      #   # size = 8,
      #   color = "black",
      #   face = "bold"
      # )
    ) +
    coord_cartesian(xlim = c(0, 17000)) +
    labs(y = "AF") -> p_af_cell

  wrap_plots(
    p_label,
    p_af_cell,
    ncol = 1,
    heights = c(1, 8)
  )

  # p_af_cell
}


fn_plot_all <- function(thepath, thevariants = thevariants) {
  log_info("Start to plot ", thepath)
  # ! parse --------------------------------------------------------------------

  gsmid <- basename(thepath)
  gseid <- basename(dirname(dirname(thepath)))

  theposes <- thevariants |>
    purrr::map(~ gsub(pattern = "[>|AGCT]", "", x = .)) |>
    purrr::map_int(as.integer)

  # load data ---------------------------------------------------------------
  # load sc
  sc <- fn_load_by_path(thepath)
  # load count

  # vaf cell umap -------------------------------------------------------------

  fn_plot_vaf_featureplot_multi(
    .thevariants = thevariants,
    sc = sc
  ) -> p_vaf_feature
}

\() {
  rawdir <- path("/mnt/isilon/u01_project/large-scale/liuc9/raw")
  thgseid <- "GSE233844"
  thesrrid <- "GSM7437874"
  thevariant <- "7757G>A"
  fn_plot_all(
    thepath = rawdir / thgseid / "final" / thesrrid,
    thevariants = thevariant
  ) -> p_vaf_feature
  outdir <- path(
    "/home/liuc9/github/scMOCHA-data/analysis/high-res/MANUSCRIPTFIGURES-notuse"
  )
  ggsave(
    filename = "{thevariant}-{thgseid}-{thesrrid}.pdf" |>
      glue::glue(),
    path = outdir,
    plot = p_vaf_feature,
    width = 12,
    height = 4,
  )
}


\() {
  rawdir <- path("/mnt/isilon/u01_project/large-scale/liuc9/raw")
  thgseid <- "GSE149689"
  thesrrid <- "GSM4509015"
  thevariant <- "9033A>G"
  fn_plot_all(
    thepath = rawdir / thgseid / "final" / thesrrid,
    thevariants = thevariant
  ) -> p_vaf_feature
  outdir <- path(
    "/home/liuc9/github/scMOCHA-data/analysis/high-res/MANUSCRIPTFIGURES-notuse"
  )
  ggsave(
    filename = "{thevariant}-{thgseid}-{thesrrid}.pdf" |>
      glue::glue(),
    path = outdir,
    plot = p_vaf_feature,
    width = 12,
    height = 4,
  )
}
