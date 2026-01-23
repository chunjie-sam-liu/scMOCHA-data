#!/usr/bin/env Rscript
# Metainfo ----------------------------------------------------------------
# @AUTHOR: Chun-Jie Liu
# @CONTACT: chunjie.sam.liu.at.gmail.com
# @DATE: 2026-01-21 22:09:20
# @DESCRIPTION: this script is used for ...

# Library -----------------------------------------------------------------

load_pkg(jutils)

# args --------------------------------------------------------------------

# s: string, i: integer, f: float, !: boolean, @: array, %: list
GetoptLong.options(help_style = "two-column")
VERSION = "v0.0.1"

# default: default value specified here.

verbose = TRUE

GetoptLong("verbose!", "print messages")


logger::log_threshold(logger::TRACE)
logger::log_layout(logger::layout_glue_colors)


# src ---------------------------------------------------------------------

# header ------------------------------------------------------------------

# load data ---------------------------------------------------------------
load_pkg(jutils)
dotenv()

outdir <- path(Sys.getenv("OUTDIR"))

color_celltype <- c(
  "homo hetero" = "blue",
  "Bimodal" = "red",
  "Mono" = "#A6D854FF",
  "B" = "#66C2A5FF",
  "T cell" = "#E5C494FF",
  "CD8 T" = "#8DA0CBFF",
  "NK" = "#FFD92FFF",
  "other" = "black"
)

variant_annotation <- import(
  outdir / "VARIANT-ANNOTATION-TABLE.xlsx"
)

thevariants <- tibble(
  variant = c(
    "709G>A", # homo hetero
    "14905G>A", # Bimodal
    "8138A>G", # Bimodal
    "2011G>A", # Bimodal
    "7751T>C", # Bimodal
    "7609T>C", # Bimodal
    "4813T>C", # Mono
    "7159T>C", # B cell
    "7833T>C", # T cell
    "10500G>A", # T cell
    "10097A>G", # T cell
    "8005T>C", # CD8 T
    "7850G>A", # CD8 T
    "9033A>G", # NK
    "7757G>A", # NK
    "9390A>G", # NK
    "6374T>C", # NK
    "10236A>G", # NK
    "1474G>A", # NK
    "9609T>C", # NK
    "2636G>A", # NK
    "15612G>A", # NK
    "2343G>A", # NK
    "7837T>C", # NK
    "6928T>C", # NK
    "2666T>C" # NK
  ),
  type = c(
    "homo hetero", # homo hetero
    "Bimodal",
    "Bimodal",
    "Bimodal",
    "Bimodal",
    "Bimodal", # Bimodal
    "Mono", # Mono
    "B", # B cell
    "T cell",
    "T cell",
    "T cell", # T cell
    "CD8 T",
    "CD8 T", # CD8 T
    "NK",
    "NK",
    "NK",
    "NK",
    "NK",
    "NK",
    "NK",
    "NK",
    "NK",
    "NK",
    "NK",
    "NK",
    "NK" # NK
  )
)

ALLVARIANTS_TEST <- import(
  outdir / "VARIANT-KRUSKAL-WALLIS-TEST.xlsx"
)

ALLVARIANTS_TEST_SIG <- ALLVARIANTS_TEST |>
  dplyr::filter(p.value < 0.05, statistic > 20)


thevariants |>
  left_join(variant_annotation, by = "variant") |>
  mutate(Locus = gsub(",Humanin", "", Locus)) -> thevariants_annotation

thevariants_annotation |> print(n = Inf)

ALLVARIANTS_TEST |>
  dplyr::mutate(
    log10p = -log10(p.value),
    log10p = ifelse(is.infinite(log10p), 300, log10p)
  ) -> ALLVARIANTS_TEST_logp


ALLVARIANTS_TEST_logp |>
  select(gseid, srrid, variant, statistic, log10p) |>
  dplyr::filter(!is.na(log10p)) |>
  # dplyr::filter(
  #   log10p > -log10(0.05),
  #   statistic > 20
  # ) |>
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

ggsci::pal_npg()(10) |> color()
thevariants_annotation$Locus |> unique() |> sort()


color_locus <- c(
  "12S" = "#E64B35FF",
  "16S" = "#4DBBD5FF",
  "ATPase6" = "#00A087FF",
  "COI" = "#3C5488FF",
  "COII" = "#F39B7FFF",
  "COIII" = "#8491B4FF",
  "Cytb" = "#91D1C2FF",
  "ND2" = "#B09C85FF",
  "ND3" = "#7E6148FF",
  "ND4L" = "#DC0000FF"
)

ALLVARIANTS_TEST_logp_replot |>
  left_join(thevariants_annotation, by = "variant") |>
  mutate(
    type = ifelse(is.na(type), "other", type)
  ) |>
  mutate(type = factor(type, names(color_celltype))) |>
  ggplot(
    aes(
      x = mean_statistic,
      y = mean_log10p,
      size = n
    )
  ) +
  geom_point(
    data = ALLVARIANTS_TEST_logp_replot |>
      # left_join(thevariants, by = "variant") |>
      left_join(thevariants_annotation, by = "variant") |>
      dplyr::filter(!is.na(type)) |>
      mutate(
        label = glue("{variant} - {aachange}")
      ),
    aes(color = Locus),
    alpha = 0.5,
    shape = 16,
  ) +
  geom_hline(
    yintercept = -log10(0.05),
    linetype = "dashed",
    color = "red"
  ) +
  ggrepel::geom_text_repel(
    data = ALLVARIANTS_TEST_logp_replot |>
      # left_join(thevariants, by = "variant") |>
      left_join(thevariants_annotation, by = "variant") |>
      dplyr::filter(!is.na(type)) |>
      mutate(
        label = glue("{variant} - {aachange}")
      ),
    aes(label = label, color = Locus),
    size = 3,
    # max.overlaps = 20,
    show.legend = FALSE,
    # nudge_x = 2,
    # nudge_y = 1
  ) +
  scale_size(name = "# of samples") +
  scale_color_manual(
    name = "Cell type specific",
    values = color_locus
  ) +
  # ggsci::scale_color_npg(
  #   name = "Locus"
  # ) +
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
    legend.position = "right",
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

ggsave(
  file.path(
    outdir,
    "VARIANT-KRUSKAL-WALLIS-STATISTIC-vs-LOG10P-VALUE-CANDIDATE-LOCUS.pdf"
  ),
  plot = plot_ks_statistic_vs_logp_replot,
  width = 11,
  height = 6
)

# load conn ---------------------------------------------------------------

# function ----------------------------------------------------------------

# body --------------------------------------------------------------------

# footer ------------------------------------------------------------------

# save image --------------------------------------------------------------
