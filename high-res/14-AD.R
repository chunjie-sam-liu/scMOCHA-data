#!/usr/bin/env Rscript
# Metainfo ----------------------------------------------------------------
# @AUTHOR: Chun-Jie Liu
# @CONTACT: chunjie.sam.liu.at.gmail.com
# @DATE: 2026-02-09 19:46:38
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
dotenv(".env")
suppressMessages({
  conflicted::conflict_prefer("filter", "dplyr")
})
cleandatadir <- path(Sys.getenv("CLEANDATADIR"))


outdir <- path(Sys.getenv("OUTDIR"))
outdirnotuse <- path(
  Sys.getenv("OUTDIRNOTUSE")
)
ALLVARIANTS <- import(
  outdir / "SAMPLE-VARIANT-CLASSIFICATION-CLUSTER-BULK-AF.xlsx"
) |>
  filter(variant_type %in% c("hete", "homo"))
METAFULL <- import(outdir / "SAMPLES-METADATA-FULL.xlsx")


variant_annotation <- import(
  path(Sys.getenv("OUTDIR")) /
    "VARIANT-ANNOTATION-TABLE-APOGEE2.xlsx"
)

# Conn ---------------------------------------------------------------

# Function ----------------------------------------------------------------

# Main --------------------------------------------------------------------
METAFULL |>
  dplyr::filter(
    disease %in% c("Healthy", "Alzheimer's Disease")
  ) |>
  filter(Chemistry == "SC5P-PE") |>
  select(gseid, srrid, Chemistry, disease) -> admeta


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


admeta_af |>
  nest(
    .by = c(variant, disease)
  ) |>
  mutate(
    nsamples = map_int(data, nrow),
  ) |>
  select(-data) -> variant_nsamples


variant_nsamples |>
  pivot_wider(
    names_from = disease,
    values_from = nsamples,
    values_fill = 0
  ) |>
  mutate(
    total = `Healthy` + `Alzheimer's Disease`
  ) |>
  arrange(-total) -> variant_nsamples_wide

variant_nsamples_wide |>
  export(outdirnotuse / "AD" / "AD-variant-sample-counts.xlsx")
admeta_af |>
  export(
    outdirnotuse / "AD" / "AD-variant-af.xlsx"
  )
admeta_af |>
  export(
    outdirnotuse / "AD" / "AD-variant-af.qs"
  )


variant_nsamples_wide |>
  mutate(
    abs = abs(`Alzheimer's Disease` - `Healthy`),
    ratio = (`Alzheimer's Disease` + 1) / (`Healthy` + 1)
  ) |>
  arrange(
    -total,
    -ratio - abs,
  ) -> variant_nsamples_wide_sorted

variant_nsamples_wide_sorted |> filter(variant == "11536C>T")
variant_nsamples_wide_sorted |> group_by(variant) |> count() |> filter(n > 1)

variant_nsamples_wide_sorted |>
  filter(Healthy > 10 & `Alzheimer's Disease` > 10)

load_pkg(stringr)
variant_annotation |>
  filter(
    !is.na(Disease)
  ) |>
  filter(variant %in% admeta_af$variant) |>
  mutate(
    label = glue("{variant}\n{Disease}"),
    disease_category = case_when(
      str_detect(
        Disease,
        regex("\\bAD\\b|\\bPD\\b|Alzheimer|Parkinson", ignore_case = TRUE)
      ) ~ "AD / PD",
      str_detect(Disease, regex("LHON", ignore_case = TRUE)) ~ "LHON",
      str_detect(
        Disease,
        regex("T2D|diabetes|metabolic", ignore_case = TRUE)
      ) ~ "Metabolic / Diabetes",
      str_detect(
        Disease,
        regex("altitude|VO2|exercise|EXIT|cyclic vomiting", ignore_case = TRUE)
      ) ~ "Exercise / Altitude",
      str_detect(
        Disease,
        regex(
          "DEAF|SNHL|hearing|dystonia|encephalomyopathy|Mitochondrial|neuropathy|Respiratory Chain",
          ignore_case = TRUE
        )
      ) ~ "Neurological",
      TRUE ~ "Other"
    )
  ) -> ad_variant_annotation

color_disease_category <- c(
  "AD / PD" = "#ff1e05",
  "LHON" = "#b300ff",
  "Metabolic / Diabetes" = "#fe7702",
  "Exercise / Altitude" = "#00d0ff",
  "Neurological" = "#0099ff",
  "Other" = "grey55"
)


fn_xy_breaks_limits(variant_nsamples_wide_sorted$total, step = 10) -> xyb

variant_nsamples |>
  mutate(
    variant = factor(
      variant,
      levels = variant_nsamples_wide_sorted$variant |> unique()
    )
  ) |>
  ggplot(aes(
    x = variant,
    y = nsamples
  )) +
  geom_col(
    aes(fill = disease)
  ) +
  scale_fill_manual(
    name = "Disease",
    values = color_disease
  ) +
  geom_hline(
    yintercept = c(2, xyb$breaks),
    color = "grey80",
    linetype = "dashed"
  ) +
  ggrepel::geom_label_repel(
    data = variant_nsamples_wide_sorted |>
      inner_join(
        ad_variant_annotation |> select(variant, label, disease_category),
        by = "variant"
      ) |>
      mutate(
        variant = factor(
          variant,
          levels = variant_nsamples_wide_sorted$variant |> unique()
        )
      ),
    aes(
      x = variant,
      y = total,
      label = label,
      color = disease_category,
      segment.colour = after_scale(colour)
    ),
    inherit.aes = FALSE,
    size = 2.5,
    lineheight = 0.9,
    label.size = 0.2,
    label.padding = unit(0.2, "lines"),
    label.r = unit(0.1, "lines"),
    box.padding = 0.4,
    point.padding = 0.3,
    min.segment.length = 0,
    segment.size = 0.3,
    segment.curvature = -0.2,
    direction = "both",
    nudge_y = max(xyb$breaks) * 0.15,
    max.overlaps = Inf,
    fill = alpha("white", 0.85)
  ) +
  scale_color_manual(
    name = "Variant annotation",
    values = color_disease_category
  ) +
  scale_x_discrete(
    expand = expansion(add = c(3, 1)),
    name = "Variant"
  ) +
  scale_y_continuous(
    breaks = xyb$breaks,
    limits = xyb$limits,
    labels = scales::comma_format(),
    expand = expansion(mult = c(0.01, 0.05)),
    name = "Number of samples"
  ) +
  guides(
    fill = guide_legend(ncol = 1, order = 1),
    color = guide_legend(ncol = 1, order = 2)
  ) +
  theme(
    panel.grid = element_blank(),
    panel.background = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.line = element_line(color = "black"),
    legend.position = c(0.8, 0.8),
    legend.box = "horizontal"
  ) -> p_ad_variant_dis


{
  ggsave(
    p_ad_variant_dis,
    filename = outdirnotuse / "AD" / "AD-variant-distribution.pdf",
    width = 17,
    height = 8
  )
}


#
#
# ? test --------------------------------------------------------------------
#
#

\() {
  olddata <- import(
    "/home/liuc9/github/scMOCHA-data/analysis/zzz/plot-ad/admeta_sc5p_variant_type_af_ttest.qs"
  )

  olddata |>
    filter(variant == "3173G>A") |>
    tidyr::unnest(data) |>
    rename(old_af = af) |>
    rename(celltype = barcode) |>
    mutate(celltype = ifelse(celltype == "bulk", "Bulk", celltype)) |>
    mutate(disease = as.character(disease)) -> olddata_3173

  admeta_af |>
    pivot_longer(
      cols = c(B:Bulk),
      names_to = "celltype",
      values_to = "af"
    ) |>
    mutate(disease = as.character(disease)) |>
    filter(variant == "3173G>A") -> admeta_af_long

  admeta_af_long |>
    left_join(
      olddata_3173 |> select(gseid, srrid, celltype, old_af),
      by = c("gseid", "srrid", "celltype")
    ) -> admeta_af_long_joined

  admeta_af_long_joined |> mutate(diff = abs(af - old_af)) |> arrange(-diff)

  admeta_af_long_joined |>
    ggplot(aes(
      x = old_af,
      y = af
    )) +
    geom_point(aes(color = celltype)) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed") -> p

  ggsave(
    p,
    filename = outdirnotuse / "AD" / "AD-variant-af-comparison.pdf",
    width = 6,
    height = 6
  )
}


admeta_af |>
  pivot_longer(
    cols = c(B:Bulk),
    names_to = "celltype",
    values_to = "af"
  ) |>
  nest(
    .by = c(variant, variant_type, celltype)
  ) -> admeta_af_celltype

admeta_af_celltype


admeta_af_celltype |>
  # filter(variant == "3173G>A")
  # head(2) |>
  dplyr::mutate(
    t = map(
      .x = data,
      .f = \(.x) {
        # .x <- a$data[[1]]
        .x |>
          dplyr::count(disease) |>
          tidyr::pivot_wider(
            names_from = disease,
            values_from = n
          ) -> .xx
        tryCatch(
          expr = {
            t.test(
              af ~ disease,
              data = .x,
              var.equal = TRUE
            ) |>
              broom::tidy() |>
              dplyr::select(
                estimate,
                estimate1,
                estimate2,
                p.value,
                conf.low,
                conf.high
              ) |>
              dplyr::bind_cols(
                .xx
              )
          },
          error = function(e) {
            message("Error: ", conditionMessage(e))
            return(NULL)
          }
        )
      }
    )
  ) -> admeta_af_ttest


admeta_af_ttest |>
  select(-data) |>
  unnest(t) |>
  dplyr::filter(p.value < 0.05) |>
  dplyr::mutate(
    plog10p = -log10(p.value),
    est = abs(estimate),
  ) |>
  dplyr::mutate(
    rank = plog10p * est,
  ) |>
  dplyr::arrange(
    desc(rank)
  ) |>
  dplyr::rename(
    ad = "Alzheimer's Disease",
  ) |>
  dplyr::filter(
    ad >= 5,
    Healthy >= 5
  )


# Save  --------------------------------------------------------------

# Session info -------------------------------------------------------------
if (isTRUE(verbose)) {
  sessionInfo()
}
