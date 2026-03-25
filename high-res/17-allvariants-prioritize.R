#!/usr/bin/env Rscript
# Metainfo ----------------------------------------------------------------
# @AUTHOR: Chun-Jie Liu
# @CONTACT: chunjie.sam.liu.at.gmail.com
# @DATE: 2026-03-25 00:01:05
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
  dplyr::filter(variant_type %in% c("hete", "homo"))
METAFULL <- import(outdir / "SAMPLES-METADATA-FULL.xlsx")
variant_annotation <- import(outdir / "VARIANT-ANNOTATION-TABLE-APOGEE2.xlsx")

# Source ---------------------------------------------------------------------

# Conn ---------------------------------------------------------------

# Function ----------------------------------------------------------------
fn_count_prediction_class <- function(df, dataname = "NA") {
  df |>
    mutate(
      prediction_class = case_when(
        prediction_class %in%
          c("Likely pathogenic", "Pathogenic") ~ "Pathogenic",
        prediction_class %in%
          c("Benign", "Likely benign") ~ "Benign",
        prediction_class %in% c("VUS", "VUS-", "VUS+") ~ "VUS",
        TRUE ~ as.character(prediction_class)
      )
    ) |>
    mutate(
      prediction_class = if_else(
        !is.na(Disease),
        "MITOMAP reported disease",
        as.character(prediction_class)
      )
    ) |>
    # count(prediction_class)
    mutate(
      prediction_class = factor(
        prediction_class,
        levels = c(
          "MITOMAP reported disease",
          "Pathogenic",
          "VUS",
          "Benign",
          "No prediction"
        )
      )
    ) |>
    count(prediction_class) |>
    pivot_wider(
      names_from = prediction_class,
      values_from = n,
      values_fill = 0
    ) |>
    add_column(
      dataname = dataname,
      .before = 1
    )
}

# Main --------------------------------------------------------------------
ALLVARIANTS |>
  left_join(
    METAFULL,
    by = c("gseid", "srrid")
  ) |>
  nest(
    .by = c(variant),
    .key = "samples"
  )

variant_annotation |>
  fn_count_prediction_class(
    dataname = "All annotated variants"
  ) |>
  export(
    outdirnotuse /
      "allvariants-prioritize" /
      "allvariants-prioritize-prediction-class-count.xlsx"
  )

variant_annotation |>
  left_join(
    ALLVARIANTS |>
      left_join(
        METAFULL |> select(-Haplogroup),
        by = c("gseid", "srrid")
      ) |>
      nest(
        .by = c(variant),
        .key = "samples"
      ),
    by = "variant"
  ) |>
  mutate(
    prediction_class_new = case_when(
      prediction_class %in%
        c("Likely pathogenic", "Pathogenic") ~ "Pathogenic",
      prediction_class %in%
        c("Benign", "Likely benign") ~ "Benign",
      prediction_class %in% c("VUS", "VUS-", "VUS+") ~ "VUS",
      TRUE ~ as.character(prediction_class)
    )
  ) |>
  mutate(
    prediction_class_new = if_else(
      !is.na(Disease),
      "MITOMAP reported disease",
      as.character(prediction_class_new)
    )
  ) -> variant_annotation_with_samples

variant_annotation_with_samples |>
  export(
    outdirnotuse /
      "allvariants-prioritize" /
      "allvariants-prioritize-prediction-class-count-with-samples.qs"
  )


variant_annotation_with_samples |>
  # count(prediction_class)
  mutate(
    prediction_class_new = factor(
      prediction_class_new,
      levels = c(
        "MITOMAP reported disease",
        "Pathogenic",
        "VUS",
        "Benign",
        "No prediction"
      )
    )
  ) |>
  mutate(
    n_homo_hete = pbmclapply(
      X = samples,
      FUN = \(.x) {
        .x |>
          count(variant_type) |>
          deframe() -> .xvt

        .x |> count(disease) |> deframe() -> .xd

        tibble(
          n_hete = ifelse(is.na(.xvt["hete"]), 0, .xvt["hete"]),
          n_homo = ifelse(is.na(.xvt["homo"]), 0, .xvt["homo"]),
          n_ad = ifelse(
            is.na(.xd["Alzheimer's Disease"]),
            0,
            .xd["Alzheimer's Disease"]
          ),
          n_covid = ifelse(
            is.na(.xd["COVID-19"]),
            0,
            .xd["COVID-19"]
          ),
          n_other = ifelse(
            is.na(.xd["Other"]),
            0,
            .xd["Other"]
          ),
          n_unknown = ifelse(
            is.na(.xd["Unknown"]),
            0,
            .xd["Unknown"]
          ),
          n_health = ifelse(
            is.na(.xd["Healthy"]),
            0,
            .xd["Healthy"]
          )
        )
      },
      mc.cores = 10
    )
  ) |>
  unnest(
    n_homo_hete
  ) -> variant_annotation_with_samples_n


variant_annotation_with_samples_n |>
  select(
    -contains("APOGEE"),
    -prediction_category,
    -prediction_class,
    -ntchange,
    -`Gnomad Frequency`,
    -pathogenicity_assessment
  ) -> variant_annotation_with_samples_n_clean


export(
  variant_annotation_with_samples_n_clean,
  outdirnotuse /
    "allvariants-prioritize" /
    "allvariants-prioritize-variant-annotation-with-samples-n-clean.qs"
)

variant_annotation_with_samples_n_clean |>
  select(
    variant,
    aachange,
    aachange_group,
    Locus,
    Locus_group,
    Disease,
    Status,
    prediction_class_new,
    contains("n_")
  ) |>
  group_by(
    prediction_class_new
  ) -> variant_annotation_with_samples_n_clean_group


variant_annotation_with_samples_n_clean_group |>
  group_split(
    .keep = TRUE
  ) |>
  set_names(
    variant_annotation_with_samples_n_clean_group |> group_keys() |> pull()
  ) |>
  as.list() -> variant_annotation_with_samples_n_clean_group_list

variant_annotation_with_samples_n_clean_group_list |>
  export(
    outdirnotuse /
      "allvariants-prioritize" /
      "allvariants-prioritize-variant-annotation-with-samples-n-clean-group.xlsx"
  )


variant_annotation_with_samples_n_clean_group_list |>
  export(
    outdirnotuse /
      "allvariants-prioritize" /
      "allvariants-prioritize-variant-annotation-with-samples-n-clean-group.qs"
  )


# Save  --------------------------------------------------------------

# Session info -------------------------------------------------------------
if (isTRUE(verbose)) {
  sessionInfo()
}
