#!/usr/bin/env Rscript
# Metainfo ----------------------------------------------------------------
# @AUTHOR: Chun-Jie Liu
# @CONTACT: chunjie.sam.liu.at.gmail.com
# @DATE: 2026-03-23 00:00:00
# @DESCRIPTION: COVID-19 GO enrichment from variant-expression correlations

# Reproducibility ----------------------------------------------------------
set.seed(1)

# Library -----------------------------------------------------------------

suppressMessages({
  load_pkg(jutils)
})

# Args --------------------------------------------------------------------

VERSION = "v0.0.1"

GetoptLong.options(help_style = "two-column")

nthread = 30
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

# Load data ----------------------------------------------------------------

dotenv()
source(path(Sys.getenv("HIGHRESDIR"), "00-colors.R"))
suppressMessages({
  conflicted::conflict_prefer("filter", "dplyr")
})

outdirnotuse <- path(Sys.getenv("OUTDIRNOTUSE"))
covid_outdir <- outdirnotuse / "COVID19"
go_outdir <- covid_outdir / "corr" / "go"
dir_create(go_outdir)

variants <- import(covid_outdir / "COVID19-variant-top-ttest-cluster-variants.qs2")
color_celltype_bulk <- c("Pseudo-bulk" = "red", color_celltype)

# Functions ----------------------------------------------------------------

normalize_qs2_output <- function(.path_qs2) {
  .path_qs2 <- as.character(.path_qs2)
  .path_qs <- sub("\\.qs2$", ".qs", .path_qs2)
  if (file.exists(.path_qs)) {
    file.copy(.path_qs, .path_qs2, overwrite = TRUE)
    invisible(file.remove(.path_qs))
  }
}

fn_go_enrich <- function(genes, universe, ont = c("BP", "CC", "MF")) {
  ont <- match.arg(ont)
  genes <- unique(genes)
  universe <- unique(universe)
  if (length(genes) == 0) {
    return(NULL)
  }

  tryCatch(
    clusterProfiler::enrichGO(
      gene = genes,
      universe = universe,
      OrgDb = "org.Hs.eg.db",
      keyType = "SYMBOL",
      ont = ont
    ),
    error = function(e) NULL
  )
}

fn_plot_go <- function(go_result, top_n = 20, ont = c("BP", "CC", "MF")) {
  ont <- match.arg(ont)
  if (is.null(go_result)) {
    return(NULL)
  }
  if (nrow(as.data.frame(go_result)) == 0) {
    return(NULL)
  }

  base_fill <- c("BP" = "#AE1700", "CC" = "#DF8F44FF", "MF" = "#00A1D5FF")
  ont_fullname <- c(
    "BP" = "Biological Process",
    "CC" = "Cellular Component",
    "MF" = "Molecular Function"
  )

  go_result |>
    tibble::as_tibble() |>
    mutate(
      Description = stringr::str_wrap(
        stringr::str_to_sentence(Description),
        width = 60
      ),
      adjp = -log10(p.adjust)
    ) |>
    select(ID, Description, adjp, Count, geneID) |>
    arrange(adjp, Count) |>
    mutate(Description = factor(Description, levels = Description)) |>
    tail(top_n) |>
    ggplot(aes(x = Description, y = adjp)) +
    geom_col(fill = base_fill[[ont]], color = NA, width = 0.7) +
    geom_text(aes(label = Count), hjust = 1, color = "white", size = 5) +
    labs(y = "-log10(Adj. P value)", x = ont_fullname[[ont]]) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
    coord_flip() +
    theme(
      panel.background = element_rect(fill = NA),
      panel.grid = element_blank(),
      axis.line.x = element_line(color = "black"),
      axis.line.y = element_line(color = "black"),
      axis.text.y = element_text(color = "black", size = 13, hjust = 1),
      axis.ticks.length.y = unit(3, units = "mm"),
      axis.text.x = element_text(color = "black", size = 12),
      axis.title = element_text(colour = "black", size = 16, face = "bold")
    )
}

fn_load_corr <- function(.variant) {
  corr_path <- covid_outdir / "corr" / glue("covid-celltype-variant-af-{.variant}-corr.csv")
  if (!file_exists(corr_path)) {
    log_warn("Missing correlation file for {.variant}: {corr_path}")
    return(tibble())
  }

  import(corr_path, lazy = FALSE) |>
    filter(adj_pval < 0.05) |>
    arrange(desc(corr)) |>
    filter(abs(corr) > 0.3)
}

fn_variant_go <- function(.variant) {
  out_path <- go_outdir / glue("covid-variant-{.variant}-go-enrichment.qs2")
  if (file_exists(out_path)) {
    log_info("Skipping {.variant}; GO enrichment output already exists")
    return(invisible(NULL))
  }

  log_info("Processing GO enrichment for {.variant}")
  corr_all <- fn_load_corr(.variant)
  if (nrow(corr_all) == 0) {
    log_warn("No significant adjusted correlations for {.variant}")
    return(invisible(NULL))
  }

  grid <- expand_grid(
    celltype = names(color_celltype_bulk),
    direction = c("pos", "neg"),
    ont = c("BP", "CC", "MF")
  )

  pbmclapply(
    X = seq_len(nrow(grid)),
    FUN = function(i) {
      ct <- grid$celltype[[i]]
      direction <- grid$direction[[i]]
      ont <- grid$ont[[i]]

      tryCatch(
        {
          corr_ct <- corr_all |> filter(celltype == ct)
          universe <- corr_ct |> pull(genename) |> unique()
          genes <- if (direction == "pos") {
            corr_ct |> filter(corr > 0.3) |> pull(genename) |> unique()
          } else {
            corr_ct |> filter(corr < -0.3) |> pull(genename) |> unique()
          }

          go <- fn_go_enrich(genes, universe, ont)
          plot <- fn_plot_go(go, 20, ont)
          if (!is.null(plot)) {
            plot <- plot +
              labs(title = glue("{.variant} | {ct}")) +
              theme(plot.title = element_text(size = 20))
          }

          tibble(
            celltype = ct,
            direction = direction,
            ont = tolower(ont),
            go = list(go),
            plot = list(plot)
          )
        },
        error = function(e) {
          log_warn("GO failed for {.variant} {ct} {direction} {ont}: {conditionMessage(e)}")
          tibble(
            celltype = ct,
            direction = direction,
            ont = tolower(ont),
            go = list(NULL),
            plot = list(NULL)
          )
        }
      )
    },
    mc.cores = nthread
  ) |>
    bind_rows() |>
    (
      function(flat) {
        flat |>
          mutate(
            col_go = paste0(direction, "_", ont),
            col_plot = paste0(direction, "_", ont, "_plot")
          ) -> flat2

        go_wide <- flat2 |>
          select(celltype, col_go, go) |>
          pivot_wider(names_from = col_go, values_from = go)

        plot_wide <- flat2 |>
          select(celltype, col_plot, plot) |>
          pivot_wider(names_from = col_plot, values_from = plot)

        left_join(go_wide, plot_wide, by = "celltype")
      }
    )() -> go_plot_wide

  export(go_plot_wide, out_path)
  normalize_qs2_output(out_path)
}

# Main ---------------------------------------------------------------------

variants |>
  walk(fn_variant_go)

if (isTRUE(verbose)) {
  sessionInfo()
}
