#!/usr/bin/env Rscript
# Metainfo ----------------------------------------------------------------
# @AUTHOR: Chun-Jie Liu
# @CONTACT: chunjie.sam.liu.at.gmail.com
# @DATE: 2026-01-01 12:00:39
# @DESCRIPTION: this script is used for ...

# Library -----------------------------------------------------------------

load_pkg(jutils)
dotenv(".env")

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
dotenv(".env")
outdir <- path(Sys.getenv("OUTDIR"))
outdirnotuse <- path(
  Sys.getenv("OUTDIRNOTUSE")
)
cleandatadir <- path(Sys.getenv("CLEANDATADIR"))

METAFULL <- import(outdir / "SAMPLES-METADATA-FULL.xlsx")
ALLVARIANTS <- import(
  outdir / "SAMPLE-VARIANT-CLASSIFICATION-CLUSTER-BULK-AF.xlsx"
)

HOMO_HETE_VARIANTS <- ALLVARIANTS |>
  dplyr::filter(variant_type %in% c("homo", "hete"))


ALLVARIANTS_TEST <- import(
  outdir / "VARIANT-KRUSKAL-WALLIS-TEST.xlsx"
)

# load conn ---------------------------------------------------------------

# src ---------------------------------------------------------------------

# function ----------------------------------------------------------------

# body --------------------------------------------------------------------
\() {
  ALLVARIANTS_TEST |>
    ggplot(aes(
      x = statistic
    )) +
    geom_histogram(
      aes(y = after_stat(density)),
      bins = 100,
      fill = "grey50",
      color = "black",
      alpha = 0.5
    ) +
    geom_vline(
      xintercept = 55,
      linetype = "dashed",
      color = "red"
    ) +
    geom_text(
      data = data.frame(
        x = 55,
        y = 0.02,
        label = "KS statistic = 55",
        vjust = -1
      ),
      aes(
        x = 55,
        y = 0.02,
        label = "KS statistic = 55",
        vjust = -1
      ),
      color = "red",
      size = 4
    ) +
    theme_bw() +
    labs(
      x = "Kruskal-Wallis H statistic",
      y = "Density",
      title = "Distribution of Kruskal-Wallis statistic for all variants"
    ) -> plot_ks_statistic
  ggsave(
    file.path(
      outdirnotuse,
      "KRUSKAL-WALLIS-STATISTIC-DISTRIBUTION.pdf"
    ),
    plot = plot_ks_statistic,
    width = 8,
    height = 6
  )
}

ALLVARIANTS_TEST |>
  dplyr::mutate(
    log10p = -log10(p.value),
    log10p = ifelse(is.infinite(log10p), 300, log10p)
  ) -> ALLVARIANTS_TEST_logp

fn_xy_breaks_limits(
  ALLVARIANTS_TEST_logp$statistic,
  max = FALSE
) -> xbl
fn_xy_breaks_limits(
  ALLVARIANTS_TEST_logp$log10p,
  max = FALSE
) -> ybl

ggplot(
  ALLVARIANTS_TEST_logp,
  aes(
    x = statistic,
    y = log10p
  )
) +
  geom_point(
    # aes(size = parameter),
    alpha = 0.5,
    shape = 16,
  ) +
  ggsci::scale_color_aaas() +
  geom_hline(
    yintercept = -log10(0.05),
    linetype = "dashed",
    color = "red"
  ) +
  # geom_abline(
  #   slope = 1,
  #   intercept = 0,
  #   linetype = "dashed",
  #   color = "red",
  # ) +
  ggrepel::geom_text_repel(
    data = ALLVARIANTS_TEST_logp |>
      dplyr::filter(
        log10p > 60,
        statistic > 300
        # variant %in% c("4175G>A", "3727T>C", "3728G>A", "3664G>A", "3243A>G")
      ),
    aes(label = variant),
    size = 3,
    max.overlaps = 20,
    show.legend = FALSE,
    # nudge_x = 2,
    nudge_y = 1
  ) +
  scale_x_continuous(
    limits = xbl$limits,
    labels = scales::label_number(
      accuracy = 1
    ),
    breaks = xbl$breaks,
    expand = expansion(mult = c(0.01, 0.01))
  ) +
  scale_y_continuous(
    limits = ybl$limits,
    labels = scales::label_number(
      accuracy = 1
    ),
    breaks = ybl$breaks,
    expand = expansion(mult = c(0.01, 0.01))
  ) +
  theme_bw() +
  theme(
    panel.background = element_blank(),
    panel.grid = element_line(colour = "grey", linetype = "dashed"),
    panel.grid.major = element_line(
      colour = "grey",
      linetype = "dashed",
      size = 0.2
    ),
    plot.title = element_text(
      hjust = 0.5,
      size = 16,
      face = "bold",
      color = "black"
    ),
    plot.subtitle = element_text(hjust = 0.5),
    legend.position = c(0.25, 0.75),
    legend.background = element_rect(
      fill = NA,
      color = NA
    ),
    legend.box.background = element_rect(
      fill = NA,
      color = NA
    ),
    axis.title = element_text(size = 12, face = "bold", color = "black"),
  ) +
  labs(
    x = "Kruskal-Wallis H Statistic",
    y = "-log10(P-value)",
    title = "Cell type specific variant",
  ) -> plot_ks_statistic_vs_logp

ggsave(
  file.path(
    outdir,
    "VARIANT-KRUSKAL-WALLIS-STATISTIC-vs-LOG10P-VALUE.pdf"
  ),
  plot = plot_ks_statistic_vs_logp,
  width = 8,
  height = 6
)


#
#
# ? replot --------------------------------------------------------------------
#
#

ALLVARIANTS_TEST_logp |>
  dplyr::filter(variant == "4175G>A") |>
  arrange(desc(log10p)) |>
  select(gseid, srrid, variant, statistic, log10p)

ALLVARIANTS_TEST_logp |>
  arrange(desc(log10p))

ALLVARIANTS_TEST_logp |>
  select(gseid, srrid, variant, statistic, log10p) |>
  dplyr::filter(!is.na(log10p)) |>
  dplyr::filter(
    log10p > -log10(0.05),
    statistic > 20
  ) |>
  group_by() |>
  tidyr::nest(.by = variant, .key = "test") |>
  mutate(
    n = purrr::map_int(test, nrow),
    mean_statistic = purrr::map_dbl(
      test,
      ~ mean(.x$statistic)
    ),
    mean_log10p = purrr::map_dbl(
      test,
      ~ mean(.x$log10p)
    )
  ) |>
  select(-test) |>
  arrange(-mean_log10p) -> ALLVARIANTS_TEST_logp_replot


fn_xy_breaks_limits(
  ALLVARIANTS_TEST_logp_replot$mean_statistic,
  max = FALSE
) -> xbl
fn_xy_breaks_limits(
  ALLVARIANTS_TEST_logp_replot$mean_log10p,
  max = FALSE
) -> ybl


ALLVARIANTS_TEST_logp_replot |>
  arrange(desc(mean_log10p)) |>
  dplyr::filter(n > 10)

ggplot(
  ALLVARIANTS_TEST_logp_replot,
  aes(
    x = mean_statistic,
    y = mean_log10p,
    size = n
  )
) +
  geom_point(
    # aes(size = parameter),
    alpha = 0.5,
    shape = 16,
  ) +
  ggsci::scale_color_aaas() +
  geom_hline(
    yintercept = -log10(0.05),
    linetype = "dashed",
    color = "red"
  ) +
  ggrepel::geom_text_repel(
    data = ALLVARIANTS_TEST_logp_replot |>
      dplyr::filter(
        mean_log10p > 60,
        mean_statistic > 300
        # variant %in% c("4175G>A", "3727T>C", "3728G>A", "3664G>A", "3243A>G")
      ),
    aes(label = variant),
    size = 3,
    max.overlaps = 20,
    show.legend = FALSE,
    # nudge_x = 2,
    nudge_y = 1
  ) +
  scale_x_continuous(
    limits = xbl$limits,
    labels = scales::label_number(
      accuracy = 1
    ),
    breaks = xbl$breaks,
    expand = expansion(mult = c(0.01, 0.01))
  ) +
  scale_y_continuous(
    limits = ybl$limits,
    labels = scales::label_number(
      accuracy = 1
    ),
    breaks = ybl$breaks,
    expand = expansion(mult = c(0.01, 0.01))
  ) +
  theme_bw() +
  theme(
    panel.background = element_blank(),
    panel.grid = element_line(colour = "grey", linetype = "dashed"),
    panel.grid.major = element_line(
      colour = "grey",
      linetype = "dashed",
      size = 0.2
    ),
    plot.title = element_text(
      hjust = 0.5,
      size = 16,
      face = "bold",
      color = "black"
    ),
    plot.subtitle = element_text(hjust = 0.5),
    legend.position = c(0.25, 0.75),
    legend.background = element_rect(
      fill = NA,
      color = NA
    ),
    legend.box.background = element_rect(
      fill = NA,
      color = NA
    ),
    axis.title = element_text(size = 12, face = "bold", color = "black"),
  ) +
  labs(
    x = "Kruskal-Wallis H Statistic",
    y = "-log10(P-value)",
    title = "Cell type specific variant",
  ) -> plot_ks_statistic_vs_logp_replot

# ggsave(
#   file.path(
#     outdir,
#     "VARIANT-KRUSKAL-WALLIS-STATISTIC-vs-LOG10P-VALUE.pdf"
#   ),
#   plot = plot_ks_statistic_vs_logp,
#   width = 8,
#   height = 6
# )

plot_folder <- "https://chop365-my.sharepoint.com/:f:/r/personal/liuc9_chop_edu/Documents/01-CHOP/01-ProjectOn/Collaboration/01-Mito/Sharing-for-manuscript/Draft-figures/celltype-specific-each?csf=1&web=1&e=sUkCKq"

ALLVARIANTS_TEST_logp |>
  select(
    gseid,
    srrid,
    Haplogroup,
    Verbose_haplogroup,
    variant_type,
    variant,
    statistic,
    p.value,
    log10p
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
          pos = pos,
          ref = ref,
          alt = alt
        )
      },
      mc.cores = 10
    )
  ) |>
  tidyr::unnest(
    cols = coord
  ) |>
  mutate(
    url = glue(
      '=HYPERLINK("{plot_folder}", "DETAIL-{pos}{ref}>{alt}-{gseid}-{srrid}.pdf")'
    )
  ) -> dd
load_pkg(openxlsx2)

wb <- wb_workbook() |>
  wb_add_worksheet(
    sheet = "KWT"
  ) |>
  wb_add_data(
    sheet = "KWT",
    x = dd
  )

wb |>
  wb_add_plot(
    sheet = "KWT",
    dims = anchor_cell,
    width = 5,
    file_type = "png",
    units = "in",
    res = 300,
  ) |>
  wb_save(
    file = "test.xlsx"
  )
# footer ------------------------------------------------------------------

# save image --------------------------------------------------------------
