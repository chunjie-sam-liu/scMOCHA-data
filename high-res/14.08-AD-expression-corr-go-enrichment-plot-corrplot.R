#!/usr/bin/env Rscript
# Metainfo ----------------------------------------------------------------
# @AUTHOR: Chun-Jie Liu
# @CONTACT: chunjie.sam.liu.at.gmail.com
# @DATE: 2026-03-03 14:06:42
# @DESCRIPTION: Plot scatter correlation, waterfall rank, and heatmap tile
#   for AD variant-expression correlations (adapted from 13-corrplot)

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


# Source -------------------------------------------------------------------
source("high-res/00-colors.R")

# Load data ---------------------------------------------------------------
load_pkg(jutils, ggrepel)
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
  dplyr::left_join(ALLVARIANTS, by = c("gseid", "srrid")) |>
  dplyr::select(-c(Chemistry, Haplogroup, Verbose_haplogroup)) |>
  dplyr::mutate(
    disease = factor(disease, levels = c("Healthy", "Alzheimer's Disease"))
  ) -> admeta_af

cluster_variant <- import(
  path(Sys.getenv("OUTDIR")) / "ALLVARIANT-ALLSAMPLES-CLUSTERAF.qs"
)

expr <- import(
  "/mnt/isilon/u01_project/large-scale/liuc9/raw/zzz/db/EXPR/gse_srrid_celltype_gene_expr.fst"
) |>
  dplyr::inner_join(
    admeta |> dplyr::select(srrid, disease),
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

genebiotype <- import(
  "config/Homo_sapiens.GRCh38.107.gtf.id_name_length_genetype.fst"
)

color_celltype_bulk <- c("Pseudo-bulk" = "red", color_celltype)

variants <- c(
  "13592C>T",
  "5031G>T",
  "8362T>G"
)

variants <- import(
  outdirnotuse / "AD" / "AD-variant-top-ttest-cluster-variants.qs"
)

# Function ----------------------------------------------------------------

theme_cor <- function() {
  theme(
    plot.title = element_text(
      size = rel(1.3),
      vjust = 2,
      hjust = 0.5,
      lineheight = 0.8
    ),
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

fn_load_corr <- function(.variant) {
  import(
    outdirnotuse /
      "AD" /
      "corr" /
      "ad-celltype-variant-af-{.variant}-corr.csv" |> glue::glue(),
    lazy = FALSE
  ) |>
    dplyr::filter(pval < 0.05) |>
    dplyr::arrange(dplyr::desc(corr)) |>
    dplyr::filter(abs(corr) > 0.3)
}

fn_load_corr_all <- function(.variant) {
  import(
    outdirnotuse /
      "AD" /
      "corr" /
      "ad-celltype-variant-af-{.variant}-corr.csv" |> glue::glue(),
    lazy = FALSE
  )
}

# Main --------------------------------------------------------------------

variants |>
  purrr::walk(\(.variant) {
    log_info("Processing variant: {.variant}")
    .safe_variant <- gsub(">", "_", .variant)
    .plotdir <- path(outdirnotuse, "AD", "corr", "corrplot", .safe_variant)
    fs::dir_create(.plotdir)

    # -- 1. Load variant AF per celltype --
    .v <- cluster_variant |>
      dplyr::select(celltype, gseid, srrid, af = dplyr::all_of(.variant)) |>
      dplyr::filter(srrid %in% admeta$srrid) |>
      tidyr::nest(.by = celltype, .key = "af")

    # -- 2. Load correlation results --
    .corr_sig <- fn_load_corr(.variant)
    .corr_all <- fn_load_corr_all(.variant)

    # -- 3. Per-celltype: scatter plots, waterfall barplot --
    celltypes_available <- intersect(
      unique(.corr_sig$celltype),
      unique(.v$celltype)
    )

    celltypes_available |>
      purrr::walk(\(.ct) {
        .safe_ct <- gsub(" ", "_", .ct)
        log_info("  Celltype: {.ct}")

        .v_ct <- .v |> dplyr::filter(celltype == .ct)

        # --- 3a. Waterfall barplot: gene rank vs correlation ---
        .corr_all_ct <- .corr_all |>
          dplyr::filter(celltype == .ct, pval < 0.05) |>
          dplyr::arrange(dplyr::desc(corr))

        if (nrow(.corr_all_ct) == 0 || nrow(.v_ct) == 0) {
          log_warn("  No data for {.ct}, skipping")
          return(invisible(NULL))
        }

        .corr_all_ct |>
          dplyr::filter(!grepl("ENSG0", genename)) |>
          dplyr::left_join(
            genebiotype,
            by = c("genename" = "gene_name")
          ) |>
          dplyr::filter(
            `Alzheimer's Disease` >= 20,
            Healthy >= 20
          ) |>
          tibble::rowid_to_column() -> .corr_filtered

        if (nrow(.corr_filtered) == 0) {
          return(invisible(NULL))
        }

        # Identify top genes: abs(corr) > 0.5, protein_coding
        .corr_filtered |>
          dplyr::mutate(
            label = ifelse(
              abs(corr) > 0.5 & Gene_type == "protein_coding",
              genename,
              NA_character_
            )
          ) -> .corr_labeled

        # Waterfall plot
        .corr_labeled |>
          ggplot(aes(x = rowid, y = corr)) +
          geom_col() +
          ggrepel::geom_text_repel(aes(label = label)) +
          scale_x_continuous(
            expand = expansion(mult = c(0.01, 0.02))
          ) +
          scale_y_continuous(
            expand = expansion(mult = c(0.02, 0.01))
          ) +
          theme_cor() +
          labs(
            x = "Gene rank",
            y = glue::glue(
              "Corr. between gene expr and {.variant} AF in {.ct}"
            )
          ) -> .p_waterfall

        ggsave(
          filename = glue::glue(
            "corr_gene_{.safe_variant}_in_{.safe_ct}.pdf"
          ),
          path = .plotdir,
          plot = .p_waterfall,
          width = 10,
          height = 6,
          dpi = 300
        )

        # --- 3b. Scatter plots for selected top genes (abs(corr) > 0.5) ---
        .selected_genes <- .corr_labeled |>
          dplyr::filter(!is.na(label)) |>
          dplyr::pull(genename)

        if (length(.selected_genes) == 0) {
          log_warn("  No top genes (|r|>0.5) for scatter in {.ct}")
          return(invisible(NULL))
        }

        log_info("  {length(.selected_genes)} top genes for scatter in {.ct}")

        # Get corr info for selected genes
        .corr_selected <- .corr_all_ct |>
          dplyr::filter(genename %in% .selected_genes)

        expr |>
          dplyr::filter(celltype == .ct) |>
          dplyr::left_join(.v_ct, by = "celltype") |>
          dplyr::inner_join(
            .corr_selected,
            by = c("genename", "celltype")
          ) |>
          dplyr::arrange(corr) -> .expr_v

        if (nrow(.expr_v) == 0) {
          return(invisible(NULL))
        }

        .scatterdir <- path(.plotdir, "scatter", .safe_ct)
        fs::dir_create(.scatterdir)

        .expr_v |>
          dplyr::mutate(
            corr_plot = parallel::mcmapply(
              FUN = \(.expr, .af, .genename, .corr, .pval) {
                tryCatch(
                  {
                    # .expr <- .expr_v$expr[[1]]
                    # .af <- .expr_v$af[[1]]
                    # .genename <- .expr_v$genename[[1]]
                    # .corr <- .expr_v$corr[[1]]
                    # .pval <- .expr_v$pval[[1]]
                    # .expr already has disease column from expr creation
                    .expr |>
                      dplyr::left_join(
                        .af,
                        by = c("gseid", "srrid")
                      ) -> .expr_af

                    .expr_af |>
                      ggplot(aes(x = af, y = expr)) +
                      geom_point(
                        aes(color = disease),
                        # position = position_jitter(
                        #   width = 0.1,
                        #   height = 0.1
                        # )
                      ) +
                      scale_color_manual(
                        values = color_disease,
                        name = "Disease"
                      ) +
                      geom_smooth(color = "red", method = "lm") +
                      scale_x_continuous(
                        expand = expansion(mult = c(0.01, 0.01))
                      ) +
                      scale_y_continuous(
                        expand = expansion(mult = c(0.01, 0))
                      ) +
                      theme_cor() +
                      theme(legend.position = "bottom") +
                      labs(
                        title = human_read_latex_pval(
                          x = human_read(.pval),
                          s = glue::glue("R={round(.corr, 3)}")
                        ),
                        x = glue::glue("{.variant} Allele frequency"),
                        y = glue::glue(
                          "{.genename} normalized gene expression"
                        )
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
          ) -> .expr_v_plot

        # Save scatter plots
        .expr_v_plot |>
          dplyr::filter(!purrr::map_lgl(corr_plot, is.null)) |>
          dplyr::mutate(
            a = purrr::pmap(
              list(genename, corr_plot),
              \(.g, .p) {
                ggsave(
                  filename = glue::glue("{.g}.pdf"),
                  path = .scatterdir,
                  plot = .p,
                  width = 7,
                  height = 6,
                  dpi = 300
                )
              }
            )
          )

        export(
          .expr_v_plot |> dplyr::select(-corr_plot),
          path(
            .plotdir,
            glue::glue("{.safe_variant}_{.safe_ct}_corr.qs")
          )
        )
      }) # end celltype walk

    # -- 4. Heatmap tile: top correlated genes across ALL cell types --
    .corr_all |>
      dplyr::filter(pval < 0.05) |>
      dplyr::filter(!grepl("ENSG0", genename)) |>
      dplyr::left_join(
        genebiotype,
        by = c("genename" = "gene_name")
      ) |>
      dplyr::filter(
        `Alzheimer's Disease` >= 20,
        Healthy >= 20,
        abs(corr) > 0.5,
        Gene_type == "protein_coding"
      ) |>
      dplyr::pull(genename) |>
      unique() -> .topcorrgenes_all

    if (length(.topcorrgenes_all) > 0) {
      .corr_all |>
        dplyr::filter(genename %in% .topcorrgenes_all) |>
        dplyr::mutate(
          celltype = gsub("_", " ", celltype),
          celltype = factor(
            celltype,
            levels = rev(names(color_celltype_bulk))
          ),
          genename = factor(
            genename,
            levels = .topcorrgenes_all
          ),
          mark = dplyr::case_when(
            pval < 0.001 ~ "***",
            pval < 0.01 ~ "**",
            pval < 0.05 ~ "*",
            TRUE ~ ""
          )
        ) -> .tile_data

      .tile_data |>
        ggplot(aes(x = genename, y = celltype)) +
        geom_tile(aes(fill = corr)) +
        geom_text(aes(label = mark), size = 5, color = "black") +
        scale_fill_gradient2(
          low = "#00fefe",
          mid = "white",
          high = "#fe0000"
        ) +
        theme_cor() +
        theme(
          axis.text.y = element_text(
            hjust = 1,
            size = 14,
            face = "bold"
          ),
          axis.text.x = element_text(
            angle = 45,
            hjust = 1,
            vjust = 1,
            face = "bold"
          ),
          axis.line = element_blank(),
          axis.ticks = element_blank(),
          legend.position = "right",
          axis.title = element_blank()
        ) +
        guides(
          fill = guide_legend(
            title = "Pearson's correlation (r)",
            title.position = "right",
            title.theme = element_text(
              angle = -90,
              size = 14,
              family = "Times"
            ),
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
        labs(x = "Gene", y = "Cell type") -> .p_tile

      ggsave(
        filename = glue::glue(
          "variant_topcorrgenes_{.safe_variant}_tileplot.pdf"
        ),
        path = .plotdir,
        plot = .p_tile,
        width = 12,
        height = 5,
        dpi = 300
      )
    } else {
      log_warn("No top correlated genes for {.variant} heatmap")
    }

    log_info("Done: {.variant}")
  })

# Save  --------------------------------------------------------------

# Session info -------------------------------------------------------------
if (isTRUE(verbose)) {
  sessionInfo()
}
