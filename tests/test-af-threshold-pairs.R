#!/usr/bin/env Rscript

load_script_env <- function(path) {
  env <- new.env(parent = globalenv())
  source(path, local = env)
  env
}


assert_identical <- function(actual, expected, label) {
  if (!identical(actual, expected)) {
    stop(
      paste0(
        label,
        "\nExpected: ", paste(deparse(expected), collapse = " "),
        "\nActual: ", paste(deparse(actual), collapse = " ")
      ),
      call. = FALSE
    )
  }
}


assert_true <- function(value, label) {
  if (!isTRUE(value)) {
    stop(label, call. = FALSE)
  }
}


run_threshold_tests <- function(env, script_label) {
  assert_true(
    exists("fn_parse_cutoff_specs", envir = env, inherits = FALSE),
    paste(script_label, "must define fn_parse_cutoff_specs")
  )

  specs <- env$fn_parse_cutoff_specs(
    cutoffs_raw = "0.6",
    cutoff_pairs_raw = "0.6:0.4,0.8:0.2"
  )

  assert_identical(length(specs), 3L, paste(script_label, "spec count"))

  symmetric_spec <- specs[[1]]
  pair_spec <- specs[[2]]

  assert_identical(
    env$fn_assign_af_group(
      af_values = c(0.59, 0.6, 0.61),
      high_cutoff = symmetric_spec$high_cutoff,
      low_cutoff = symmetric_spec$low_cutoff,
      high_inclusive = symmetric_spec$high_inclusive,
      low_inclusive = symmetric_spec$low_inclusive
    ),
    c("low_af", NA_character_, "high_af"),
    paste(script_label, "symmetric grouping")
  )

  assert_identical(
    env$fn_assign_af_group(
      af_values = c(0.39, 0.4, 0.5, 0.6, 0.61),
      high_cutoff = pair_spec$high_cutoff,
      low_cutoff = pair_spec$low_cutoff,
      high_inclusive = pair_spec$high_inclusive,
      low_inclusive = pair_spec$low_inclusive
    ),
    c("low_af", "low_af", NA_character_, "high_af", "high_af"),
    paste(script_label, "asymmetric grouping")
  )

  assert_identical(
    env$fn_cutoff_dir(pair_spec),
    "cutoff-high-0.6-low-0.4",
    paste(script_label, "asymmetric directory")
  )

  assert_identical(
    env$fn_cutoff_comparison_label(pair_spec),
    "AF >= 0.6 vs AF <= 0.4",
    paste(script_label, "asymmetric comparison label")
  )
}


script_paths <- c(
  "high-res/17.05-af-cutoff-comparison-by-disease.R",
  "high-res/17.07-af-cutoff-comparison-by-disease-gse166992.R"
)

for (script_path in script_paths) {
  env <- load_script_env(script_path)
  run_threshold_tests(env, basename(script_path))
}

message("threshold-pair tests passed")
