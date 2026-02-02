#!/usr/bin/env Rscript
# Metainfo ----------------------------------------------------------------
# @AUTHOR: Chun-Jie Liu
# @CONTACT: chunjie.sam.liu.at.gmail.com
# @DATE: 2026-01-01 17:17:24
# @DESCRIPTION: this script is used for ...

# Library -----------------------------------------------------------------

load_pkg(jutils)
load_pkg(magrittr)


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
load_pkg(jutils)

conflicted::conflict_prefer("filter", "dplyr")

dotenv(".env")
outdir <- path(Sys.getenv("OUTDIR"))
outdirnotuse <- path(
  Sys.getenv("OUTDIRNOTUSE")
)
cleandatadir <- path(Sys.getenv("CLEANDATADIR"))
scmergedir <- outdirnotuse / "scmerge"
scintegrateddir <- outdirnotuse / "scintegrated"

# Ensure output directories exist
fs::dir_create(scmergedir)
fs::dir_create(scintegrateddir)

# load conn ---------------------------------------------------------------

# src ---------------------------------------------------------------------

# function ----------------------------------------------------------------

fn_de_plot <- function(
  markers,
  .cutoff_pval = 0.05,
  .cutoff_log2fc = 0.25,
  .pct = 0.05
) {
  markers |>
    tibble::rownames_to_column("gene") |>
    dplyr::mutate(
      fdr = -log10(p_val_adj)
    ) |>
    dplyr::mutate(
      fdr = ifelse(
        fdr > -log10(1e-300),
        -log10(1e-300),
        fdr
      )
    ) |>
    dplyr::mutate(
      avg_log2FC = ifelse(
        abs(avg_log2FC) > 100,
        sign(avg_log2FC) * 100,
        avg_log2FC
      )
    ) |>
    dplyr::mutate(
      color = dplyr::case_when(
        p_val_adj < .cutoff_pval &
          (pct.1 >= .pct | pct.2 >= .pct) &
          avg_log2FC > .cutoff_log2fc ~
          "red",
        p_val_adj < .cutoff_pval &
          (pct.1 >= .pct | pct.2 >= .pct) &
          avg_log2FC < -.cutoff_log2fc ~
          "blue",
        TRUE ~ "grey"
      )
    ) -> forplot

  forplot |>
    dplyr::count(color) |>
    print()

  forplot |>
    dplyr::count(color) |>
    tibble::deframe() -> n_color

  forplot |>
    ggplot(aes(
      x = avg_log2FC,
      y = fdr,
      color = color
    )) +
    geom_point(aes()) +
    ggrepel::geom_text_repel(
      data = forplot |>
        dplyr::filter(color != "grey") |>
        dplyr::group_by(color) |>
        dplyr::slice_head(n = 20) |>
        dplyr::ungroup(),
      aes(label = gene),
      size = 3,
      max.overlaps = 20
    ) +
    scale_color_identity() +
    geom_vline(
      xintercept = c(
        -.cutoff_log2fc,
        .cutoff_log2fc
      ),
      linetype = "dashed"
    ) +
    geom_hline(
      yintercept = -log10(.cutoff_pval),
      linetype = "dashed"
    ) +
    theme_classic() +
    labs(
      x = "Fold change(Heteroplasmy/Sufficient reads)",
      y = "FDR",
      subtitle = glue::glue(
        "Up:{
          scales::label_comma()(
          ifelse(is.na(n_color['red']), 0, n_color['red'])
        )
        }, down: {
          scales::label_comma()(
          ifelse(is.na(n_color['blue']), 0, n_color['blue'])
        )
        }; (Cutoff: FDR<{.cutoff_pval}, log2FC>{.cutoff_log2fc}, Pct>{.pct})"
      )
    ) +
    theme(
      plot.title = element_text(
        hjust = 0.5,
        face = "bold",
        color = "black",
        size = 16
      ),
      plot.subtitle = element_text(
        hjust = 0.5,
        # face = "bold",
        color = "black",
        size = 12
      ),
    ) -> p
  list(
    p = p,
    markers = forplot
  )
}
fn_go_enrich <- function(cancer_sgene, ont = c("BP", "CC", "MF")) {
  .ont <- match.arg(ont)

  gobp <- clusterProfiler::enrichGO(
    gene = cancer_sgene |> unique(),
    OrgDb = "org.Hs.eg.db",
    keyType = "SYMBOL",
    # keyType = "ENSEMBL",
    ont = .ont
  )

  gobp
}

fn_plot_go <- function(.go, .topn = Inf, .ont = c("BP", "CC", "MF")) {
  if (is.null(.go) || nrow(.go) == 0) {
    p <- ggplot() +
      theme_void() +
      labs(
        title = "No significant GO terms found"
      ) +
      theme(
        plot.title = element_text(
          hjust = 0.5,
          face = "bold",
          color = "black",
          size = 16
        )
      )
    return(p)
  }

  .ont <- match.arg(.ont)

  base_fill <- c("BP" = "#AE1700", "CC" = "#DF8F44FF", "MF" = "#00A1D5FF")
  ont_fullname <- c(
    "BP" = "Biological Process",
    "CC" = "Cellular Component",
    "MF" = "Molecular Function"
  )

  .ont_fill <- base_fill[.ont]
  x_label <- ont_fullname[.ont]

  .go |>
    tibble::as_tibble() |>
    dplyr::mutate(
      Description = stringr::str_wrap(
        stringr::str_to_sentence(string = Description),
        width = 60
      )
    ) |>
    dplyr::mutate(adjp = -log10(p.adjust)) |>
    dplyr::select(ID, Description, adjp, Count, geneID) |>
    dplyr::arrange(adjp, Count) |>
    dplyr::mutate(
      Description = factor(Description, levels = Description)
    ) -> .go_bp_for_plot

  if (!is.infinite(.topn)) {
    .go_bp_for_plot |>
      tail(.topn) -> .go_bp_for_plot
  }

  .go_bp_for_plot |>
    ggplot(aes(x = Description, y = adjp)) +
    geom_col(fill = .ont_fill, color = NA, width = 0.7) +
    geom_text(aes(label = Count), hjust = 1, color = "white", size = 5) +
    labs(y = "-log10(Adj. P value)", x = x_label) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
    coord_flip() +
    theme(
      panel.background = element_rect(fill = NA),
      panel.grid = element_blank(),
      axis.line.x = element_line(color = "black"),
      axis.line.y = element_line(color = "black"),
      # axis.title.y = element_blank(),
      axis.text.y = element_text(color = "black", size = 13, hjust = 1),
      axis.ticks.length.y = unit(3, units = "mm"),
      axis.text.x = element_text(color = "black", size = 12),
      axis.title = element_text(colour = "black", size = 16, face = "bold")
    )
}

fn_variant_go <- function(markers, .variant) {
  # .variant <- variants[[1]]
  .corr_pos <- markers |>
    dplyr::filter(color == "red") |>
    dplyr::mutate(
      genename = gene
    )
  .corr_neg <- markers |>
    dplyr::filter(color == "blue") |>
    dplyr::mutate(
      genename = gene
    )

  .pos_bp <- fn_go_enrich(cancer_sgene = unique(.corr_pos$genename), "BP")
  .pos_cc <- fn_go_enrich(cancer_sgene = unique(.corr_pos$genename), "CC")
  .pos_mf <- fn_go_enrich(cancer_sgene = unique(.corr_pos$genename), "MF")

  .neg_bp <- fn_go_enrich(cancer_sgene = unique(.corr_neg$genename), "BP")
  .neg_cc <- fn_go_enrich(cancer_sgene = unique(.corr_neg$genename), "CC")
  .neg_mf <- fn_go_enrich(cancer_sgene = unique(.corr_neg$genename), "MF")

  tibble::tibble(
    # variant = .variant,
    pos_bp = list(.pos_bp),
    pos_cc = list(.pos_cc),
    pos_mf = list(.pos_mf),
    neg_bp = list(.neg_bp),
    neg_cc = list(.neg_cc),
    neg_mf = list(.neg_mf),
    pos_bp_plot = list(
      fn_plot_go(.pos_bp, 20, "BP") +
        labs(title = "{.variant} POS BP" |> glue::glue()) +
        theme(
          plot.title = element_text(size = 20)
        )
    ),
    pos_cc_plot = list(
      fn_plot_go(.pos_cc, 20, "CC") +
        labs(title = "{.variant} POS CC" |> glue::glue()) +
        theme(
          plot.title = element_text(size = 20)
        )
    ),
    pos_mf_plot = list(
      fn_plot_go(.pos_mf, 20, "MF") +
        labs(title = "{.variant} POS MF" |> glue::glue()) +
        theme(
          plot.title = element_text(size = 20)
        )
    ),
    neg_bp_plot = list(
      fn_plot_go(.neg_bp, 20, "BP") +
        labs(title = "{.variant} NEG BP" |> glue::glue()) +
        theme(
          plot.title = element_text(size = 20)
        )
    ),
    neg_cc_plot = list(
      fn_plot_go(.neg_cc, 20, "CC") +
        labs(title = "{.variant} NEG CC" |> glue::glue()) +
        theme(
          plot.title = element_text(size = 20)
        )
    ),
    neg_mf_plot = list(
      fn_plot_go(.neg_mf, 20, "MF") +
        labs(title = "{.variant} NEG MF" |> glue::glue()) +
        theme(
          plot.title = element_text(size = 20)
        )
    )
  )
}


fn_de_ <- function(
  thevariant,
  sc,
  .ident.1,
  .ident.2,
  .group.by,
  .prefix,
  .labs,
  .celltype = NULL,
  .tmpdir = "deg_merge_new"
) {
  # .ident.1 <- glue::glue("{hetero_label} high")
  # .ident.2 <- glue::glue("{hetero_label} low")
  # .group.by <- "cellvarianttype2"
  # .prefix <- "hetero_high_vs_low"
  # .labs <- labs(
  #   x = "Fold change {hetero_label} High vs Low" |> glue::glue(),
  #   y = "FDR",
  #   # title = "m.{thevariant}" |> glue::glue()
  #   title = "Markers: {hetero_label} High vs Low (m.{thevariant})" |>
  #     glue::glue()
  # )
  meta.data <- sc@meta.data

  .outdir <- outdirnotuse / "deg" / thevariant / .tmpdir

  if (!is.null(.celltype)) {
    .outdir <- outdirnotuse / "deg" / thevariant / .tmpdir / .celltype
  }
  dir_create(.outdir)
  DefaultAssay(sc) <- "RNA"

  markers_hetero_high_vs_low <- Seurat::FindMarkers(
    object = sc,
    ident.1 = .ident.1,
    ident.2 = .ident.2,
    assay = "RNA",
    slot = "data",
    test.use = "wilcox",
    group.by = .group.by,
    latent.vars = "srrid",
    features = Seurat::VariableFeatures(sc)
  )

  export(
    markers_hetero_high_vs_low,
    file.path(
      .outdir,
      "markers.{.prefix}.{thevariant}.qs" |>
        glue::glue()
    )
  )

  fn_de_plot(
    markers_hetero_high_vs_low,
    .cutoff_pval = 0.05,
    .cutoff_log2fc = 0.25,
    .pct = 0.05
  ) -> p_hetero_high_vs_low
  p_hetero_high_vs_low$p + .labs -> p_hetero_high_vs_low$p

  ggsave(
    filename = "markers.{.prefix}.{thevariant}.pdf" |>
      glue::glue() |>
      fs::path_sanitize(),
    plot = p_hetero_high_vs_low$p,
    path = .outdir,
    device = "pdf",
    width = 10,
    height = 6
  )
  p_hetero_high_vs_low
}
fn_go_ <- function(
  thevariant,
  p_hetero_high_vs_low,
  .prefix,
  .celltype = NULL,
  .tmpdir = "go_merge_new"
) {
  # .prefix <- "hetero_high_vs_low"

  fn_variant_go(
    p_hetero_high_vs_low$markers,
    thevariant
  ) -> p_go_hetero_high_vs_low

  .outdir <- path(
    outdirnotuse / "deg" / thevariant / .tmpdir
  )
  if (!is.null(.celltype)) {
    .outdir <- path(
      outdirnotuse / "deg" / thevariant / .tmpdir / .celltype
    )
  }
  dir_create(.outdir)

  export(
    p_go_hetero_high_vs_low,
    file.path(
      .outdir,
      "markers.{.prefix}.{thevariant}.go.qs" |>
        glue::glue()
    )
  )

  tibble::tibble(
    pn = c("pos", "neg") |> rep(each = 3),
    t = c("bp", "cc", "mf") |> rep(times = 2)
  ) |>
    dplyr::mutate(
      saveimage = purrr::map2(
        .x = pn,
        .y = t,
        .f = \(.x, .y) {
          .p <- p_go_hetero_high_vs_low[[glue::glue("{.x}_{.y}_plot")]][[1]]
          .filename <- "markers.{.prefix}.{thevariant}.go.{.x}_{.y}_plot.pdf" |>
            glue::glue() |>
            fs::path_sanitize()
          ggsave(
            filename = .filename,
            plot = .p +
              labs(
                title = glue::glue(
                  "m.{thevariant} {ifelse(is.null(.celltype), '', .celltype)}"
                )
              ),
            path = .outdir,
            device = "pdf",
            width = 10,
            height = 6
          )
        }
      )
    )
}

fn_load_sc <- function(thevariant) {
  library(Seurat)
  sc <- import(
    outdirnotuse / "scintegrated" / glue::glue("sc.{thevariant}.integrated.qs")
  )
  sc
}

fn_variant_ <- function(
  thevariant,
  sc,
  .vs = c(0.5, 0.5),
  .celltype = NULL
) {
  # thevariant <- "3727T>C"
  log_info(
    "Start analysis for {thevariant} with {.vs[1]} and {.vs[2]}"
  )

  sc@meta.data |>
    dplyr::count(cellvarianttype) |>
    dplyr::arrange(cellvarianttype) |>
    tibble::deframe() -> n_cellvarianttype

  if (all(.vs == c("Heteroplasmy", "Sufficient reads"))) {
    log_info("Special case for {thevariant} Heteroplasmy vs Sufficient reads")
    log_info(
      "Start DE for  {thevariant} Heteroplasmy vs Sufficient reads"
    )
    p_hetero_vs_sufficient <- fn_de_(
      thevariant = thevariant,
      sc = sc,
      .ident.1 = "Heteroplasmy",
      .ident.2 = "Sufficient reads",
      .group.by = "cellvarianttype",
      .prefix = "hetero_vs_sufficient",
      .labs <- labs(
        x = "Fold change Heteroplasmy (n={
          scales::label_comma()(n_cellvarianttype['Heteroplasmy'])
        }) vs Sufficient reads (n={
          scales::label_comma()(n_cellvarianttype['Sufficient reads'])
        })" |>
          glue::glue(),
        y = "FDR",
        title = "Markers: Heteroplasmy vs Sufficient Reads (m.{thevariant}) {ifelse(is.na(.celltype), '', .celltype)}" |>
          glue::glue()
      ),
      .celltype = .celltype
    )

    if (!is.null(p_hetero_vs_sufficient)) {
      log_info(
        "GO enrichment for {thevariant} Heteroplasmy vs Sufficient reads"
      )
      fn_go_(
        thevariant = thevariant,
        p_hetero_high_vs_low = p_hetero_vs_sufficient,
        .prefix = "hetero_vs_sufficient",
        .celltype = .celltype
      )
    }
    log_success(
      "Finished analysis for {thevariant} with {.vs[1]} and {.vs[2]}"
    )
    return(invisible(NULL))
  }

  sc@meta.data |>
    as.data.table() |>
    dplyr::filter(cellvarianttype == "Heteroplasmy") |>
    dplyr::pull(af) |>
    quantile(probs = seq(0, 1, 0.05), na.rm = FALSE) -> .quant

  .high <- .quant[glue::glue("{.vs[1] * 100}%")]
  .low <- .quant[glue::glue("{.vs[2] * 100}%")]

  # scales::label_number(accuracy = 0.01)(median_af)
  .label_high <- glue::glue(
    "High={scales::label_number(accuracy = 1)(.vs[1] * 100)}% AF={scales::label_number(accuracy = 0.01)(.high)}"
  )
  .label_low <- glue::glue(
    "Low={scales::label_number(accuracy = 1)(.vs[2] * 100)}% AF={scales::label_number(accuracy = 0.01)(.low)}"
  )
  hetero_label <- glue::glue(
    "Heteroplasmy ({.label_high}) vs ({.label_low})"
  )

  sc@meta.data |>
    dplyr::mutate(
      cellvarianttype2 = dplyr::case_when(
        cellvarianttype == "Heteroplasmy" &
          af >= .high ~
          glue::glue("{.label_high}"),
        cellvarianttype == "Heteroplasmy" &
          af < .low ~
          glue::glue("{.label_low}"),
        TRUE ~ as.character(cellvarianttype)
      )
    ) -> sc@meta.data

  DefaultAssay(sc) <- "SCT"

  sc@meta.data |>
    dplyr::count(cellvarianttype2) |>
    dplyr::arrange(cellvarianttype2) |>
    tibble::deframe() -> n_cellvarianttype2

  log_info(
    "Start DE for Heteroplasmy {thevariant}  high {.label_high} vs low {.label_low}"
  )

  p_hetero_high_vs_low <- fn_de_(
    thevariant = thevariant,
    sc = sc,
    .ident.1 = glue::glue("{.label_high}"),
    .ident.2 = glue::glue("{.label_low}"),
    .group.by = "cellvarianttype2",
    .prefix = hetero_label,
    .labs = labs(
      x = glue::glue(
        "Fold change Heteroplasmy ({.label_high} n={
          scales::label_comma()(n_cellvarianttype2[.label_high])
        }) vs Low ({.label_low} n={
          scales::label_comma()(n_cellvarianttype2[.label_low])
        })"
      ),
      y = "FDR",
      # title = "m.{thevariant}" |> glue::glue()
      title = "Markers: {hetero_label} (m.{thevariant}) {ifelse(is.null(.celltype), '', .celltype)}" |>
        glue::glue()
    ),
    .celltype = .celltype
  )

  if (!is.null(p_hetero_high_vs_low)) {
    log_info(
      "GO enrichment for  {thevariant} Heteroplasmy high {.label_high} vs low {.label_low}"
    )
    fn_go_(
      thevariant = thevariant,
      p_hetero_high_vs_low = p_hetero_high_vs_low,
      .prefix = hetero_label,
      .celltype = .celltype
    )
  }
  log_success(
    "Finished analysis for {thevariant} with {.vs[1]} and {.vs[2]}"
  )
}

fn_variant_vaf_ <- function(
  thevariant,
  sc,
  .vs = 0.4,
  .celltype = NULL
) {
  # thevariant <- "3727T>C"
  log_info(
    "Start analysis for {thevariant} with {.vs}"
  )

  sc@meta.data |>
    dplyr::count(cellvarianttype) |>
    dplyr::arrange(cellvarianttype) |>
    tibble::deframe() -> n_cellvarianttype

  sc@meta.data |>
    as.data.table() |>
    dplyr::filter(cellvarianttype == "Heteroplasmy") |>
    dplyr::pull(af) |>
    quantile(probs = seq(0, 1, 0.05), na.rm = FALSE) -> .quant

  .high <- .vs
  .low <- .vs

  # scales::label_number(accuracy = 0.01)(median_af)
  .label_high <- glue::glue(
    "High (AF>={scales::label_number(accuracy = 0.01)(.high)})"
  )
  .label_low <- glue::glue(
    "Low (AF<{scales::label_number(accuracy = 0.01)(.low)}) and wildtype"
  )
  hetero_label <- glue::glue(
    "{.label_high} vs {.label_low}"
  )

  sc@meta.data |>
    dplyr::mutate(
      cellvarianttype2 = dplyr::case_when(
        cellvarianttype == "Heteroplasmy" &
          af >= .high ~
          glue::glue("{.label_high}"),
        cellvarianttype == "Sufficient reads" ~
          glue::glue("{.label_low}"),
        cellvarianttype == "Heteroplasmy" &
          af < .low ~
          glue::glue("{.label_low}"),
        TRUE ~ as.character(cellvarianttype)
      )
    ) -> sc@meta.data

  # DefaultAssay(sc) <- "SCT"

  sc@meta.data |>
    dplyr::count(cellvarianttype2) |>
    dplyr::arrange(cellvarianttype2) |>
    tibble::deframe() -> n_cellvarianttype2

  log_info(
    "Start DE for Heteroplasmy {thevariant}  high {.label_high} (n={scales::label_comma()(n_cellvarianttype2[.label_high])}) vs low {.label_low} (n={scales::label_comma()(n_cellvarianttype2[.label_low])})"
  )

  p_hetero_high_vs_low <- fn_de_(
    thevariant = thevariant,
    sc = sc,
    .ident.1 = glue::glue("{.label_high}"),
    .ident.2 = glue::glue("{.label_low}"),
    .group.by = "cellvarianttype2",
    .prefix = hetero_label,
    .labs = labs(
      x = glue::glue(
        "Fold change ({.label_high} n={
          scales::label_comma()(n_cellvarianttype2[.label_high])
        }) vs ({.label_low} n={
          scales::label_comma()(n_cellvarianttype2[.label_low])
        })"
      ),
      y = "FDR",
      # title = "m.{thevariant}" |> glue::glue()
      title = "Markers: {hetero_label} (m.{thevariant}) {ifelse(is.null(.celltype), '', .celltype)}" |>
        glue::glue()
    ),
    .celltype = .celltype,
    .tmpdir = "deg_merge_vaf"
  )

  if (!is.null(p_hetero_high_vs_low)) {
    log_info(
      "GO enrichment for  {thevariant} Heteroplasmy high {.label_high} vs low {.label_low}"
    )
    fn_go_(
      thevariant = thevariant,
      p_hetero_high_vs_low = p_hetero_high_vs_low,
      .prefix = hetero_label,
      .celltype = .celltype,
      .tmpdir = "go_merge_vaf"
    )
  }
  log_success(
    "Finished analysis for {thevariant} with {hetero_label}"
  )
}

fn_variant_cell_ <- function(thevariant, sc, .vs = c(0.5, 0.5)) {
  library(Seurat)
  sc$predicted.celltype.l1 |> unique() -> celltypes

  log_info(
    "Start celltype-specific analysis for {thevariant} with {.vs[1]} and {.vs[2]}"
  )

  lapply(
    celltypes,
    function(.celltype) {
      sc_sub <- subset(
        sc,
        subset = predicted.celltype.l1 == .celltype
      )
      log_info(
        "Start celltype-specific analysis for {thevariant} with {.vs[1]} and {.vs[2]} in {.celltype}"
      )
      fn_variant_(
        thevariant = thevariant,
        sc = sc_sub,
        .vs = .vs,
        .celltype = .celltype
      )
      log_success(
        "Finished celltype-specific analysis for {thevariant} with {.vs[1]} and {.vs[2]} in {.celltype}"
      )
    }
  )
  log_success(
    "Finished celltype-specific analysis for {thevariant} with {.vs[1]} and {.vs[2]}"
  )
}

fn_variant_cell_vaf_ <- function(thevariant, sc, .vs = 0.4) {
  library(Seurat)
  sc$predicted.celltype.l1 |> unique() -> celltypes

  log_info(
    "Start celltype-specific analysis for {thevariant} with {.vs}"
  )

  lapply(
    celltypes,
    function(.celltype) {
      sc_sub <- subset(
        sc,
        subset = predicted.celltype.l1 == .celltype
      )
      log_info(
        "Start celltype-specific analysis for {thevariant} with {.vs} in {.celltype}"
      )
      fn_variant_vaf_(
        thevariant = thevariant,
        sc = sc_sub,
        .vs = .vs,
        .celltype = .celltype
      )
      log_success(
        "Finished celltype-specific analysis for {thevariant} with {.vs} in {.celltype}"
      )
    }
  )
  log_success(
    "Finished celltype-specific analysis for {thevariant} with {.vs}"
  )
}


fn_main <- function(thevariant) {
  vaf_cutoff <- c(0.4, 0.5, 0.6, 0.7, 0.8)
  cli_alert_info("Processing variant {thevariant}")

  sc <- fn_load_sc(thevariant)
  cli_alert_success("Loaded sc for variant {thevariant}")

  cli_alert_info("Running variant analysis for {thevariant} with vss")

  pbmclapply(
    X = vaf_cutoff,
    FUN = \(.vs) {
      cli_alert_info(
        "Processing variant {thevariant} with vaf cutoff {.strong {(.vs)}}"
      )
      tryCatch(
        expr = {
          fn_variant_vaf_(
            thevariant = thevariant,
            sc = sc,
            .vs = .vs,
            .celltype = NULL
          )
        },
        error = function(e) {
          cli_alert_danger(
            "Error in variant {thevariant} with vaf cutoff {.strong {(.vs)}}: {e$message}"
          )
        }
      )
    },
    mc.cores = 2
  )

  cli_alert_success("Finished variant analysis for {thevariant} with vss")
  pbmclapply(
    X = vaf_cutoff,
    FUN = \(.vs) {
      tryCatch(
        expr = {
          fn_variant_cell_vaf_(
            thevariant = thevariant,
            sc = sc,
            .vs = .vs
          )
        },
        error = function(e) {
          cli_alert_danger(
            "Error in celltype-specific variant {thevariant} with vaf cutoff {.strong {(.vs)}}: {e$message}"
          )
        }
      )
    },
    mc.cores = 2
  )
}
# body --------------------------------------------------------------------

thevariants <- c(
  "3173G>A",
  "3176A>T",
  "3178T>A",
  "3727T>C",
  "3728C>T",
  "13271T>C",
  "14063T>C",
  "14831G>A",
  "1643A>G",
  "3667T>G",
  "4175G>A",
  "5513G>A",
  "7065G>A",
  "9025G>A",
  "9237G>A",
  "10398A>G"
)

thevariants <- c()


thevariant <- "4175G>A"

# fn_main(thevariant)

thevariants <- c(
  # "4175G>A",
  c(
    "14082C>G",
    "15169A>G",
    "3240C>G",
    "7757G>A",
    "3173G>A",
    "3176A>T",
    "3178T>A",
    "9025G>A",
    "9237G>A",
    "10398A>G"
  )
)

# thats for all variants, don't run or run once
thevariants |>
  purrr::map(
    .f = \(thevariant) {
      fn_main(thevariant)
    }
  ) -> res_all_variants
# footer ------------------------------------------------------------------

# save image --------------------------------------------------------------
