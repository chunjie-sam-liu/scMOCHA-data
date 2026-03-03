#!/usr/bin/env Rscript
# Metainfo ----------------------------------------------------------------
# @AUTHOR: Chun-Jie Liu
# @CONTACT: chunjie.sam.liu.at.gmail.com
# @DATE: 2026-03-03 13:40:58
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

outdir <- path(Sys.getenv("OUTDIR"))
outdirnotuse <- path(Sys.getenv("OUTDIRNOTUSE"))
ALLVARIANTS <- import(
  outdir / "SAMPLE-VARIANT-CLASSIFICATION-CLUSTER-BULK-AF.xlsx"
) |>
  filter(variant_type %in% c("hete", "homo"))
METAFULL <- import(outdir / "SAMPLES-METADATA-FULL.xlsx")
METAFULL |>
  dplyr::filter(
    disease %in% c("Healthy", "Alzheimer's Disease"),
    Chemistry == "SC5P-PE"
  ) |>
  dplyr::select(gseid, srrid, Chemistry, disease) -> admeta

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

cluster_variant <- import(
  path(Sys.getenv("OUTDIR")) / "ALLVARIANT-ALLSAMPLES-CLUSTERAF.qs"
)

expr <- import(
  "/mnt/isilon/u01_project/large-scale/liuc9/raw/zzz/db/EXPR/gse_srrid_celltype_gene_expr.fst"
) |>
  dplyr::inner_join(
    admeta |> select(srrid, disease),
    by = c("srrid")
  ) |>
  tidyr::pivot_longer(
    cols = -c(gseid, srrid, genename, disease),
    names_to = "celltype",
    values_to = "expr"
  ) |>
  tidyr::nest(
    .by = c(genename, celltype),
    .key = "expr"
  )


color_celltype_bulk <- c("Pseudo-bulk" = "red", color_celltype)
variants <- c(
  "13592C>T",
  "5031G>T",
  "8362T>G"
)

# Conn ---------------------------------------------------------------

# Function ----------------------------------------------------------------

fn_go_enrich <- function(cancer_sgene, ont = c("BP", "CC", "MF")) {
  .ont <- match.arg(ont)

  .genes <- unique(cancer_sgene)
  if (length(.genes) == 0) {
    return(NULL)
  }

  tryCatch(
    clusterProfiler::enrichGO(
      gene = .genes,
      OrgDb = "org.Hs.eg.db",
      keyType = "SYMBOL",
      ont = .ont
    ),
    error = function(e) NULL
  )
}

fn_plot_go <- function(.go, .topn = Inf, .ont = c("BP", "CC", "MF")) {
  .ont <- match.arg(.ont)

  if (is.null(.go)) {
    return(NULL)
  }
  if (nrow(as.data.frame(.go)) == 0) {
    return(NULL)
  }

  base_fill <- c("BP" = "#AE1700", "CC" = "#DF8F44FF", "MF" = "#00A1D5FF")
  ont_fullname <- c(
    "BP" = "Biological Process",
    "CC" = "Cellular Component",
    "MF" = "Molecular Function"
  )

  .ont_fill <- base_fill[.ont]
  x_label <- ont_fullname[.ont]

  .go %>%
    tibble::as_tibble() %>%
    dplyr::mutate(
      Description = stringr::str_wrap(
        stringr::str_to_sentence(string = Description),
        width = 60
      )
    ) %>%
    dplyr::mutate(adjp = -log10(p.adjust)) %>%
    dplyr::select(ID, Description, adjp, Count, geneID) %>%
    dplyr::arrange(adjp, Count) %>%
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

fn_load_corr <- function(.variant) {
  import(
    outdirnotuse /
      "AD" /
      "corr" /
      "ad-celltype-variant-af-{.variant}-corr.csv" |> glue(),
    lazy = FALSE
  ) |>
    dplyr::filter(pval < 0.05) |>
    dplyr::arrange(desc(corr)) |>
    dplyr::filter(abs(corr) > 0.3)
}

fn_variant_go <- function(.variant) {
  # .variant <- variants[[1]]
  .corr_all <- fn_load_corr(.variant)

  # Full cartesian product: celltype × direction × ont — all run in parallel
  .grid <- tidyr::expand_grid(
    celltype = names(color_celltype_bulk),
    direction = c("pos", "neg"),
    ont = c("BP", "CC", "MF")
  )
  nthread <- 27
  pbmclapply(
    X = seq_len(nrow(.grid)),
    FUN = function(.i) {
      .ct <- .grid$celltype[.i]
      .dir <- .grid$direction[.i]
      .ont <- .grid$ont[.i]

      tryCatch(
        {
          .corr <- .corr_all |> dplyr::filter(celltype == .ct)
          .genes <- if (.dir == "pos") {
            .corr |>
              dplyr::filter(corr > 0.3) |>
              dplyr::pull(genename) |>
              unique()
          } else {
            .corr |>
              dplyr::filter(corr < -0.3) |>
              dplyr::pull(genename) |>
              unique()
          }

          .go <- fn_go_enrich(.genes, .ont)
          .title <- glue::glue("{.variant} | {.ct}")
          .plot <- fn_plot_go(.go, 20, .ont)
          if (!is.null(.plot)) {
            .plot <- .plot +
              labs(title = .title) +
              theme(plot.title = element_text(size = 20))
          }

          tibble::tibble(
            celltype = .ct,
            direction = .dir,
            ont = tolower(.ont),
            go = list(.go),
            plot = list(.plot)
          )
        },
        error = function(e) {
          message(glue::glue(
            "WARN [{.ct} {.dir} {.ont}]: {conditionMessage(e)}"
          ))
          tibble::tibble(
            celltype = .ct,
            direction = .dir,
            ont = tolower(.ont),
            go = list(NULL),
            plot = list(NULL)
          )
        }
      )
    },
    mc.cores = nthread
  ) |>
    dplyr::bind_rows() |>
    (\(.flat) {
      .flat |>
        dplyr::mutate(
          col_go = paste0(direction, "_", ont),
          col_plot = paste0(direction, "_", ont, "_plot")
        ) -> .f

      .go_wide <- .f |>
        dplyr::select(celltype, col_go, go) |>
        tidyr::pivot_wider(names_from = col_go, values_from = go)

      .plot_wide <- .f |>
        dplyr::select(celltype, col_plot, plot) |>
        tidyr::pivot_wider(names_from = col_plot, values_from = plot)

      dplyr::left_join(.go_wide, .plot_wide, by = "celltype")
    })() -> .go_plot_wide

  export(
    .go_plot_wide,
    outdirnotuse /
      "AD" /
      "corr" /
      "go" /
      "ad-variant-{.variant}-go-enrichment.qs" |> glue()
  )
}

variants[-1] |>
  purrr::walk(
    fn_variant_go
  )

# Main --------------------------------------------------------------------

# Save  --------------------------------------------------------------

# Session info -------------------------------------------------------------
if (isTRUE(verbose)) {
  sessionInfo()
}
