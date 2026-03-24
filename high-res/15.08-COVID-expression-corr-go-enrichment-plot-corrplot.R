#!/usr/bin/env Rscript
# Metainfo ----------------------------------------------------------------
# @AUTHOR: Chun-Jie Liu
# @CONTACT: chunjie.sam.liu.at.gmail.com
# @DATE: 2026-03-23 00:00:00
# @DESCRIPTION: COVID-19 variant-expression corrplot workflow

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

load_pkg(ggrepel)
dotenv()
source(path(Sys.getenv("HIGHRESDIR"), "00-colors.R"))
suppressMessages({
  conflicted::conflict_prefer("filter", "dplyr")
})

outdir <- path(Sys.getenv("OUTDIR"))
outdirnotuse <- path(Sys.getenv("OUTDIRNOTUSE"))
covid_outdir <- outdirnotuse / "COVID19"
expr_path <- Sys.getenv(
  "SCMOCHA_EXPR_PATH",
  unset = "/mnt/isilon/u01_project/large-scale/liuc9/raw/zzz/db/EXPR/gse_srrid_celltype_gene_expr.fst"
)
rootdir <- path_dir(path(Sys.getenv("HIGHRESDIR")))

allvariants <- import(
  outdir / "SAMPLE-VARIANT-CLASSIFICATION-CLUSTER-BULK-AF.xlsx"
) |>
  filter(variant_type %in% c("hete", "homo"))
metafull <- import(outdir / "SAMPLES-METADATA-FULL.xlsx")

covid_meta <- metafull |>
  filter(disease %in% c("Healthy", "COVID-19")) |>
  select(gseid, srrid, disease) |>
  mutate(disease = factor(disease, levels = c("Healthy", "COVID-19")))

cluster_variant <- import(outdir / "ALLVARIANT-ALLSAMPLES-CLUSTERAF.qs")

expr <- import(expr_path) |>
  inner_join(covid_meta |> select(srrid, disease), by = "srrid") |>
  pivot_longer(
    cols = -c(gseid, srrid, genename, disease),
    names_to = "celltype",
    values_to = "expr"
  ) |>
  nest(.by = c(genename, celltype), .key = "expr")

genebiotype <- import(rootdir / "config" / "Homo_sapiens.GRCh38.107.gtf.id_name_length_genetype.fst") |>
  distinct(gene_name, .keep_all = TRUE)

color_celltype_bulk <- c("Pseudo-bulk" = "red", color_celltype)
variants <- import(covid_outdir / "COVID19-variant-top-ttest-cluster-variants.qs2")

# Functions ----------------------------------------------------------------

theme_cor <- function() {
  theme(
    plot.title = element_text(size = rel(1.3), vjust = 2, hjust = 0.5, lineheight = 0.8),
    axis.title.x = element_text(face = "bold", size = 16),
    axis.title.y = element_text(face = "bold", size = 16, angle = 90),
    axis.text = element_text(size = rel(1.1)),
    axis.text.x = element_text(hjust = 0.5, vjust = 0, size = 14),
    axis.text.y = element_text(vjust = 0.5, hjust = 0, size = 14),
    axis.line = element_line(colour = "black"),
    axis.ticks = element_line(colour = "black"),
    legend.title = element_text(size = rel(1.1), face = "bold"),
    legend.text = element_text(size = rel(1.1), face = "bold"),
    legend.background = element_blank(),
    legend.key = element_rect(fill = NA, colour = NA),
    strip.text = element_text(size = rel(1.3)),
    panel.background = element_blank(),
    complete = TRUE
  )
}

normalize_qs2_output <- function(.path_qs2) {
  .path_qs2 <- as.character(.path_qs2)
  .path_qs <- sub("\\.qs2$", ".qs", .path_qs2)
  if (file.exists(.path_qs)) {
    file.copy(.path_qs, .path_qs2, overwrite = TRUE)
    invisible(file.remove(.path_qs))
  }
}

fn_corr_path <- function(.variant) {
  covid_outdir / "corr" / glue("covid-celltype-variant-af-{.variant}-corr.csv")
}

fn_load_corr <- function(.variant) {
  import(fn_corr_path(.variant), lazy = FALSE) |>
    filter(adj_pval < 0.05) |>
    arrange(desc(corr)) |>
    filter(abs(corr) > 0.3)
}

fn_load_corr_all <- function(.variant) {
  import(fn_corr_path(.variant), lazy = FALSE)
}

# Main ---------------------------------------------------------------------

variants_available <- variants[
  purrr::map_lgl(variants, \(.v) file_exists(fn_corr_path(.v)))
]

if (length(variants_available) < length(variants)) {
  missing_variants <- setdiff(variants, variants_available)
  log_warn(
    "Skipping {length(missing_variants)} variant(s) with no corr file: {paste(missing_variants, collapse = ', ')}"
  )
}

variants_available |>
  walk(
    \(.variant) {
      log_info("Processing corrplot for {.variant}")
      safe_variant <- gsub(">", "_", .variant, fixed = TRUE)
      plotdir <- path(covid_outdir, "corr", "corrplot", safe_variant)
      dir_create(plotdir)

      variant_af <- cluster_variant |>
        select(celltype, gseid, srrid, af = all_of(.variant)) |>
        filter(srrid %in% covid_meta$srrid) |>
        nest(.by = celltype, .key = "af")

      corr_sig <- fn_load_corr(.variant)
      corr_all <- fn_load_corr_all(.variant)

      celltypes_available <- intersect(unique(corr_sig$celltype), unique(variant_af$celltype))

      celltypes_available |>
        walk(
          \(.ct) {
            safe_ct <- gsub(" ", "_", .ct)
            variant_af_ct <- variant_af |> filter(celltype == .ct)

            corr_all_ct <- corr_all |>
              filter(celltype == .ct, adj_pval < 0.05) |>
              arrange(desc(corr))

            if (nrow(corr_all_ct) == 0 || nrow(variant_af_ct) == 0) {
              log_warn("No data for {.variant} {.ct}; skipping")
              return(invisible(NULL))
            }

            corr_all_ct |>
              filter(!grepl("ENSG0", genename)) |>
              left_join(genebiotype, by = c("genename" = "gene_name")) |>
              filter(`COVID-19` >= 20, Healthy >= 20) |>
              tibble::rowid_to_column() -> corr_filtered

            if (nrow(corr_filtered) == 0) {
              return(invisible(NULL))
            }

            corr_filtered |>
              mutate(
                label = ifelse(abs(corr) > 0.5 & Gene_type == "protein_coding", genename, NA_character_)
              ) -> corr_labeled

            p_waterfall <- corr_labeled |>
              ggplot(aes(x = rowid, y = corr)) +
              geom_col() +
              ggrepel::geom_text_repel(aes(label = label)) +
              scale_x_continuous(expand = expansion(mult = c(0.01, 0.02))) +
              scale_y_continuous(expand = expansion(mult = c(0.02, 0.01))) +
              theme_cor() +
              labs(
                x = "Gene rank",
                y = glue("Corr. between gene expr and {.variant} AF in {.ct}")
              )

            ggsave(
              filename = glue("corr_gene_{safe_variant}_in_{safe_ct}.pdf"),
              path = plotdir,
              plot = p_waterfall,
              width = 10,
              height = 6,
              dpi = 300
            )

            selected_genes <- corr_labeled |>
              filter(!is.na(label)) |>
              pull(genename)

            if (length(selected_genes) == 0) {
              log_warn("No top genes (|r|>0.5) for scatter in {.variant} {.ct}")
              return(invisible(NULL))
            }

            corr_selected <- corr_all_ct |>
              filter(genename %in% selected_genes)

            expr |>
              filter(celltype == .ct) |>
              left_join(variant_af_ct, by = "celltype") |>
              inner_join(corr_selected, by = c("genename", "celltype")) |>
              arrange(corr) -> expr_variant

            if (nrow(expr_variant) == 0) {
              return(invisible(NULL))
            }

            scatterdir <- path(plotdir, "scatter", safe_ct)
            dir_create(scatterdir)

            expr_variant |>
              mutate(
                corr_plot = parallel::mcmapply(
                  FUN = function(.expr, .af, .genename, .corr, .pval) {
                    tryCatch(
                      {
                        .expr |>
                          left_join(.af, by = c("gseid", "srrid")) -> expr_af

                        expr_af |>
                          ggplot(aes(x = af, y = expr)) +
                          geom_point(aes(color = disease)) +
                          scale_color_manual(values = color_disease, name = "Disease") +
                          geom_smooth(color = "red", method = "lm") +
                          scale_x_continuous(expand = expansion(mult = c(0.01, 0.01))) +
                          scale_y_continuous(expand = expansion(mult = c(0.01, 0))) +
                          theme_cor() +
                          theme(legend.position = "bottom") +
                          labs(
                            title = human_read_latex_pval(
                              x = human_read(.pval),
                              s = glue("R={round(.corr, 3)}")
                            ),
                            x = glue("{.variant} Allele frequency"),
                            y = glue("{.genename} normalized gene expression")
                          )
                      },
                      error = function(e) NULL
                    )
                  },
                  .expr = expr,
                  .af = af,
                  .genename = genename,
                  .corr = corr,
                  .pval = pval,
                  mc.cores = nthread,
                  SIMPLIFY = FALSE
                )
              ) -> expr_variant_plot

            expr_variant_plot |>
              filter(!purrr::map_lgl(corr_plot, is.null)) |>
              mutate(
                saved = purrr::pmap(
                  list(genename, corr_plot),
                  function(genename, corr_plot) {
                    ggsave(
                      filename = glue("{genename}.pdf"),
                      path = scatterdir,
                      plot = corr_plot,
                      width = 7,
                      height = 6,
                      dpi = 300
                    )
                  }
                )
              )

            export(
              expr_variant_plot |> select(-corr_plot),
              path(plotdir, glue("{safe_variant}_{safe_ct}_corr.qs2"))
            )
            normalize_qs2_output(path(plotdir, glue("{safe_variant}_{safe_ct}_corr.qs2")))
          }
        )

      topcorrgenes_all <- corr_all |>
        filter(pval < 0.05) |>
        filter(adj_pval < 0.05) |>
        filter(!grepl("ENSG0", genename)) |>
        left_join(genebiotype, by = c("genename" = "gene_name")) |>
        filter(`COVID-19` >= 20, Healthy >= 20, abs(corr) > 0.5, Gene_type == "protein_coding") |>
        pull(genename) |>
        unique()

      if (length(topcorrgenes_all) > 0) {
        corr_all |>
          filter(genename %in% topcorrgenes_all) |>
          mutate(
            celltype = gsub("_", " ", celltype),
            celltype = factor(celltype, levels = rev(names(color_celltype_bulk))),
            genename = factor(genename, levels = topcorrgenes_all),
            mark = case_when(
              adj_pval < 0.001 ~ "***",
              adj_pval < 0.01 ~ "**",
              adj_pval < 0.05 ~ "*",
              TRUE ~ ""
            )
          ) -> tile_data

        p_tile <- tile_data |>
          ggplot(aes(x = genename, y = celltype)) +
          geom_tile(aes(fill = corr)) +
          geom_text(aes(label = mark), size = 5, color = "black") +
          scale_fill_gradient2(low = "#00fefe", mid = "white", high = "#fe0000") +
          theme_cor() +
          theme(
            axis.text.y = element_text(hjust = 1, size = 14, face = "bold"),
            axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, face = "bold"),
            axis.line = element_blank(),
            axis.ticks = element_blank(),
            legend.position = "right",
            axis.title = element_blank()
          ) +
          guides(
            fill = guide_legend(
              title = "Pearson's correlation (r)",
              title.position = "right",
              title.theme = element_text(angle = -90, size = 14, family = "Times"),
              title.vjust = -0.5,
              title.hjust = 0.5,
              label = TRUE,
              label.position = "left",
              label.theme = element_text(size = 14, family = "Times"),
              label.hjust = 0.5,
              label.vjust = 0.5,
              keywidth = 1,
              keyheight = 1.8,
              reverse = TRUE
            )
          ) +
          coord_fixed(ratio = 1) +
          labs(x = "Gene", y = "Cell type")

        ggsave(
          filename = glue("variant_topcorrgenes_{safe_variant}_tileplot.pdf"),
          path = plotdir,
          plot = p_tile,
          width = 12,
          height = 5,
          dpi = 300
        )
      } else {
        log_warn("No top correlated genes for {.variant} heatmap")
      }
    }
  )

if (isTRUE(verbose)) {
  sessionInfo()
}
