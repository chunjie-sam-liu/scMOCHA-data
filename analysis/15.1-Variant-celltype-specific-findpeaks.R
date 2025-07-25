#!/usr/bin/env Rscript
# Metainfo ----------------------------------------------------------------
# @AUTHOR: Chun-Jie Liu
# @CONTACT: chunjie.sam.liu.at.gmail.com
# @DATE: 2025-06-20 19:17:51
# @DESCRIPTION: filename
# @VERSION: v0.0.1

# Library -----------------------------------------------------------------

suppressPackageStartupMessages(library(magrittr))
library(ggplot2)
library(patchwork)
library(prismatic)
library(paletteer)
library(data.table)
# library(rlang)
library(glue)
library(parallel)
library(GetoptLong)
library(logger)

# args --------------------------------------------------------------------

# s: string, i: integer, f: float, !: boolean
# @: array
# %: hash
# default: default value specified here.
verbose <- FALSE
spec <- "
Usage: Rscript foorbar.R [options]
Options:

<verbose!> Print messages
"

GetoptLong.options(help_style = "two-column")
GetoptLong(spec, template_control = list(opt_width = 21))

# header ------------------------------------------------------------------

# future: :plan(future: :multisession, workers = 10)

# load data ---------------------------------------------------------------
cleandatadir <- "/home/liuc9/github/scMOCHA-data/data/zzz/clean-data"
dbdir <- "/home/liuc9/github/scMOCHA-data/analysis/zzz/db"
ks_test_dir <- file.path(dbdir, "all_hetero_af.cell.ks_test")
plotdir <- "/home/liuc9/github/scMOCHA-data/analysis/zzz/plot-celltype-specific-variant"

# gseid_srrid_ks_load <- import(
#   file.path(
#     ks_test_dir,
#     "a_gseid_srrid_ks_load.nocellaf.qs"
#   )
# )

# ALLVARIANTS <- import(file.path(
#   "/home/liuc9/github/scMOCHA-data/analysis/zzz/clean-data/", "all_variant.qs"
# ))

# META <- import("/home/liuc9/github/scMOCHA-data/analysis/zzz/clean-data/gse_dataset_metadata_full.sex_pred.qs") |>
#   dplyr::select(
#     gseid, srrid, Age_new, Age_group,
#     Haplogroup,
#     disease, Chemistry, sex_pred
#   ) |>
#   dplyr::mutate(
#     Haplogroup = purrr::map_chr(
#       .x = Haplogroup,
#       .f = \(.x) {
#         # if (stringr::str_starts(.x, "L")) {
#         #   gsub("L", "L0", .x)
#         # }
#         gsub("\\d+.*", "", .x)
#       }
#     )
#   )

# load conn --------------------------------------------------------------------

conn <- DBI::dbConnect(
  duckdb::duckdb(),
  "/home/liuc9/github/scMOCHA-data/analysis/zzz/clean-data/all_hetero_af.cell.duckdb.1.2.1"
)
DBI::dbListTables(conn)


tbl_all_hetero_af_cell <- dplyr::tbl(
  conn,
  "all_hetero_af_cell"
)
tbl_all_hetero_altdepth_cell <- dplyr::tbl(
  conn,
  "all_hetero_altdepth_cell"
)
tbl_all_hetero_sumdepth_cell <- dplyr::tbl(
  conn,
  "all_hetero_sumdepth_cell"
)
tbl_barcode <- dplyr::tbl(
  conn,
  "barcode"
)
tbl_allvariants <- dplyr::tbl(
  conn,
  "allvariants"
)
tbl_meta <- dplyr::tbl(conn, "meta") |>
  dplyr::select(
    gseid,
    srrid,
    Age_new,
    Age_group,
    Haplogroup,
    disease,
    Chemistry,
    sex_pred
  )

# DBI::dbDisconnect(conn, shutdown = TRUE)

tbl_allvariants_celltype_peaks <- dplyr::tbl(
  conn,
  "allvariants_celltype_peaks"
)

# src ---------------------------------------------------------------------
source("./analysis/00-colors.R")

# function ----------------------------------------------------------------

fn_plot_ggdist <- function(
  thevariant,
  thegseid,
  thesrrid
) {
  library(ggdist)
  tbl_all_hetero_af_cell |>
    dplyr::filter(
      # gseid == thegseid,
      srrid == thesrrid,
      variant == thevariant,
      af > 0
    ) |>
    dplyr::collect() -> .d
  # thevariant <- "7833T>C"
  .variant <- .d$variant[1]
  .gseid <- .d$gseid[1]
  .srrid <- .d$srrid[1]

  .d |>
    dplyr::filter(af > 0) |>
    dplyr::mutate(
      celltype = gsub(
        "_",
        " ",
        celltype
      )
    ) |>
    dplyr::mutate(
      celltype = factor(
        celltype,
        levels = names(color_celltype) |> rev()
      )
    ) -> forplot_

  forplot_ |>
    ggplot(aes(
      x = af,
      y = celltype,
      fill = celltype
    )) +
    ggdist::stat_halfeye(
      scale = 1
    ) +
    ggdist::stat_interval(
      show.legend = FALSE,
    ) +
    stat_summary(geom = "point", fun = median, show.legend = FALSE) +
    scale_fill_manual(
      values = color_celltype,
      na.value = "grey50"
    ) +
    scale_color_manual(values = MetBrewer::met.brewer("VanGogh3")) +
    # scale_color_brewer() +
    guides(col = "none") +
    ggridges::theme_ridges() +
    # ggridges::stat_density_ridges(
    #   quantile_lines = TRUE, quantiles = 2
    # ) +
    theme(
      legend.position = "none",
      plot.title = element_text(
        hjust = 0.5,
        size = 16
      ),
    ) +
    labs(
      title = paste0(.variant, "\n(", .gseid, "-", .srrid, ")"),
      x = "Allele Frequency",
      y = "Cell Type"
    )
}

fn_plot_joy <- function(
  thevariant,
  thegseid,
  thesrrid
) {
  tbl_all_hetero_af_cell |>
    dplyr::filter(
      # gseid == thegseid,
      srrid == thesrrid,
      variant == thevariant,
      af > 0
    ) |>
    dplyr::collect() -> .d
  # thevariant <- "7833T>C"
  .variant <- .d$variant[1]
  .gseid <- .d$gseid[1]
  .srrid <- .d$srrid[1]

  .d |>
    dplyr::filter(af > 0) |>
    dplyr::mutate(
      celltype = gsub(
        "_",
        " ",
        celltype
      )
    ) |>
    dplyr::mutate(
      celltype = factor(
        celltype,
        levels = names(color_celltype) |> rev()
      )
    ) -> forplot_

  forplot_ |>
    ggplot(aes(
      x = af,
      y = celltype,
      fill = celltype
    )) +
    ggridges::geom_density_ridges(
      # scale = 3,
      # alpha = 0.8,
      rel_min_height = 0.01,
      size = 0.1
    ) +
    scale_fill_manual(
      values = color_celltype,
      na.value = "grey50"
    ) +
    ggridges::theme_ridges() +
    # ggridges::stat_density_ridges(
    #   quantile_lines = TRUE, quantiles = 2
    # ) +
    theme(
      legend.position = "none",
      plot.title = element_text(
        hjust = 0.5,
        size = 16
      ),
    ) +
    labs(
      title = paste0(.variant, "\n(", .gseid, "-", .srrid, ")"),
      x = "Allele Frequency",
      y = "Cell Type"
    )
}

fn_plot_joy_celltype_level2_level3 <- function(
  thevariant,
  thegseid,
  thesrrid,
  thecelltype,
  thecelltype_prefix,
  thecelltype_level
) {
  tbl_barcode |>
    dplyr::filter(
      srrid == thesrrid
    ) |>
    dplyr::collect() -> celltypedetail

  tbl_all_hetero_af_cell |>
    dplyr::filter(
      gseid == thegseid,
      srrid == thesrrid,
      variant == thevariant,
      af > 0
    ) |>
    dplyr::select(-celltype) |>
    dplyr::collect() |>
    dplyr::left_join(
      celltypedetail,
      by = c("barcode")
    ) |>
    dplyr::rename(
      plotcelltype = "celltype_{thecelltype_level}" |> glue::glue(),
    ) -> thevariant_data

  thevariant_data |>
    dplyr::filter(
      celltype == thecelltype,
      af > 0
    ) |>
    dplyr::filter(grepl(thecelltype_prefix, plotcelltype)) |>
    dplyr::mutate(
      plotcelltype = factor(
        plotcelltype
      )
    ) -> forplot

  levels(forplot$plotcelltype)

  color_celltype_detail <- log(seq(
    1,
    exp(1),
    length.out = length(levels(forplot$plotcelltype))
  )) |>
    purrr::map_chr(
      ~ prismatic::clr_lighten(
        color_celltype[thecelltype],
        .x
      )
    )
  names(color_celltype_detail) <- levels(forplot$plotcelltype)

  forplot |>
    ggplot(aes(
      x = af,
      y = plotcelltype,
      fill = plotcelltype
    )) +
    ggridges::geom_density_ridges(
      # scale = 3,
      # alpha = 0.8,
      rel_min_height = 0.01,
      size = 0.1
    ) +
    scale_fill_manual(
      values = color_celltype_detail,
      na.value = "grey50"
    ) +
    ggridges::theme_ridges() +
    theme(
      legend.position = "none",
      plot.title = element_text(
        hjust = 0.5,
        # size = 16
      ),
    ) +
    labs(
      title = "{thecelltype}-{thecelltype_level}-{thevariant}\n({thegseid}-{thesrrid})" |>
        glue::glue(),
      x = "Allele Frequency",
      y = "Cell Type"
    )
}

fn_plot_joy_celltype_detail <- function(
  thevariant,
  thegseid,
  thesrrid
) {
  tibble::tibble(
    thevariant = thevariant,
    thegseid = thegseid,
    thesrrid = thesrrid,
    thecelltype = c(
      c("B", "CD4 T", "CD8 T", "other T") |> rep(times = 2),
      c("NK", "DC", "Mono", "other") |> rep(times = 2)
    ),
    thecelltype_prefix = c(
      c("B", "CD4", "CD8", "") |> rep(times = 2),
      c("NK", "DC", "Mono", "") |> rep(times = 2)
    ),
    thecelltype_level = c("l2", "l3") |> rep(each = 4) |> rep(times = 2)
  ) -> thevariant_celltype_df

  thevariant_celltype_df |>
    dplyr::mutate(
      p = parallel::mcmapply(
        thevariant = thevariant,
        thegseid = thegseid,
        thesrrid = thesrrid,
        thecelltype = thecelltype,
        thecelltype_prefix = thecelltype_prefix,
        thecelltype_level = thecelltype_level,
        FUN = fn_plot_joy_celltype_level2_level3,
        mc.cores = 5,
        SIMPLIFY = FALSE
      )
    ) -> plot_thevariant_celltype_list
}

find_major_peaks <- function(
  x,
  adjust = 1,
  threshold_factor = 1.5,
  minpeakdistance_frac = 0.05,
  plot = FALSE
) {
  d <- density(x, adjust = adjust)
  x_vals <- d$x
  y_vals <- d$y

  threshold_auto <- mean(y_vals) + threshold_factor * sd(y_vals)
  minpeakdistance_auto <- round(length(x_vals) * minpeakdistance_frac)

  peaks <- pracma::findpeaks(
    y_vals,
    threshold = threshold_auto,
    minpeakdistance = minpeakdistance_auto
  )

  if (is.null(peaks)) {
    warning("No major peaks detected.")
    return(data.frame(Peak_X = numeric(0), Peak_Y = numeric(0)))
  }

  peak_x <- x_vals[peaks[, 2]]
  peak_y <- y_vals[peaks[, 2]]

  if (plot) {
    plot(d, main = "Density with Major Peaks")
    points(peak_x, peak_y, col = "red", pch = 19)
    text(peak_x, peak_y, labels = round(peak_x, 2), pos = 3, col = "red")
  }

  return(data.frame(Peak_X = peak_x, Peak_Y = peak_y))
}

find_density_peaks <- function(
  x,
  bw = "nrd0",
  min_height = 0.01,
  min_prominence = 0.1,
  min_distance = 0.05,
  ...
) {
  # Calculate density
  dx <- density(x, bw = bw, ...)

  # Find local maxima
  n <- length(dx$y)
  peaks <- c()

  # Check each point (excluding endpoints)
  for (i in 2:(n - 1)) {
    if (
      dx$y[i] > dx$y[i - 1] && dx$y[i] > dx$y[i + 1] && dx$y[i] >= min_height
    ) {
      peaks <- c(peaks, i)
    }
  }

  # Filter peaks by prominence and distance
  if (length(peaks) > 0) {
    peaks_info <- data.frame(
      idx = peaks,
      x = dx$x[peaks],
      y = dx$y[peaks]
    )

    # Calculate prominence for each peak
    peaks_info$prominence <- 0
    for (i in 1:nrow(peaks_info)) {
      peak_idx <- peaks_info$idx[i]
      peak_height <- peaks_info$y[i]

      # Find minimum heights on both sides
      left_min <- min(dx$y[1:peak_idx])
      right_min <- min(dx$y[peak_idx:n])
      baseline <- max(left_min, right_min)

      peaks_info$prominence[i] <- peak_height - baseline
    }

    # Filter by prominence
    peaks_info <- peaks_info[peaks_info$prominence >= min_prominence, ]

    # Filter by minimum distance between peaks
    if (nrow(peaks_info) > 1) {
      # Sort by height (keep strongest peaks when too close)
      peaks_info <- peaks_info[order(peaks_info$y, decreasing = TRUE), ]

      keep <- rep(TRUE, nrow(peaks_info))
      for (i in 1:(nrow(peaks_info) - 1)) {
        if (keep[i]) {
          for (j in (i + 1):nrow(peaks_info)) {
            if (abs(peaks_info$x[i] - peaks_info$x[j]) < min_distance) {
              keep[j] <- FALSE
            }
          }
        }
      }
      peaks_info <- peaks_info[keep, ]
    }

    # Sort by density value (highest first) and remove helper columns
    peaks_info <- peaks_info[
      order(peaks_info$y, decreasing = TRUE),
      c("x", "y", "prominence")
    ]
  } else {
    peaks_info <- data.frame(
      x = numeric(0),
      y = numeric(0),
      prominence = numeric(0)
    )
  }

  return(list(
    peaks = peaks_info,
    density_obj = dx
  ))
}

find_cell_density_peak <- function(x, pheight = 0.1, pprominence = 0.2) {
  dx <- density(x)
  # Find multiple peaks with adaptive filtering for bimodal distributions
  max_density <- max(dx$y)
  adaptive_min_height <- max_density * pheight # 10% of maximum density
  adaptive_min_prominence <- max_density * pprominence # 20% of maximum density

  peaks_result <- find_density_peaks(
    x,
    min_height = adaptive_min_height,
    min_prominence = adaptive_min_prominence,
    min_distance = 0.2
  )

  data.table(
    npeaks = nrow(peaks_result$peaks),
    peakmin = min(peaks_result$peaks$x, na.rm = TRUE),
    peaks = peaks_result$peaks |> jsonlite::toJSON() |> as.character()
  )
}

find_variant_celltype_peak <- function(thevariant) {
  tbl_all_hetero_af_cell |>
    dplyr::filter(
      variant == thevariant,
      af > 0,
    ) |>
    as.data.table() -> thevariant_data

  thevariant_data |>
    tidyr::nest(
      .by = c(variant, gseid, srrid, celltype),
      .key = "data"
    ) |>
    dplyr::mutate(
      peaks = parallel::mcmapply(
        .x = data,
        FUN = \(.x) {
          if (nrow(.x) < 30) {
            return(
              NULL
            )
          }
          tryCatch(
            find_cell_density_peak(.x$af),
            error = function(e) {
              return(
                NULL
              )
            }
          )
        },
        SIMPLIFY = FALSE,
        mc.cores = 5
      )
    ) |>
    dplyr::select(-data, -variant) |>
    tidyr::unnest(cols = peaks)
}

extract_haplogroup <- function(haplo) {
  library(stringr)
  level4 <- haplo
  level3 <- str_extract(haplo, "^[A-Z]+\\d*")
  level2 <- str_extract(haplo, "^[A-Z]+")

  level1 <- dplyr::case_when(
    str_detect(level2, "^L") ~ "L",
    str_detect(level2, "^M") ~ "M",
    str_detect(level2, "^(C|D|Z|G|E|Q|M7|M8|M9)") ~ "M",
    str_detect(level2, "^N") ~ "N",
    str_detect(level2, "^(R|H|J|T|U|K|B|F|V|A|Y)") ~ "R",
    str_detect(level2, "^(W|X)") ~ "Other",
    TRUE ~ "Unknown"
  )

  data.table(HG1 = level1, HG2 = level2, HG3 = level3, HG4 = level4)
}

# body --------------------------------------------------------------------

tbl_allvariants |>
  dplyr::filter(issomatic == "heteroplasmic") |>
  dplyr::select(variant) |>
  as.data.table() -> allvariants
# ! don't run, it only run once
# allvariants |>
#   dplyr::mutate(
#     peaks = parallel::mcmapply(
#       thevariant = variant,
#       FUN = find_variant_celltype_peak,
#       mc.cores = 10,
#       SIMPLIFY = FALSE
#     )
#   ) |>
#   tidyr::unnest(cols = peaks) ->
# allvariants_celltype_peaks

# export(
#   allvariants_celltype_peaks,
#   file.path(
#     "/home/liuc9/github/scMOCHA-data/analysis/zzz/clean-data",
#     "celltype_specific_allvariants_celltype_peaks.qs"
#   ),
#   format = "both"
# )
# DBI::dbWriteTable(
#   conn,
#   "allvariants_celltype_peaks",
#   allvariants_celltype_peaks,
#   overwrite = TRUE,
#   temporary = FALSE
# )

tbl_allvariants_celltype_peaks |>
  dplyr::select(-peaks) |>
  as.data.table() |>
  dplyr::left_join(
    tbl_meta |>
      as.data.table() |>
      dplyr::mutate(
        HG = parallel::mcmapply(
          haplo = Haplogroup,
          FUN = extract_haplogroup,
          mc.cores = 10,
          SIMPLIFY = FALSE
        )
      ) |>
      tidyr::unnest(cols = HG),
    by = c("gseid", "srrid")
  ) -> allvariants_celltype_peaks_meta


# ? peak and age --------------------------------------------------------------------

allvariants_celltype_peaks_meta |>
  dplyr::filter(
    !is.na(Age_new),
    !is.na(peakmin)
  ) |>
  tidyr::nest(
    .by = c(variant, celltype)
  ) |>
  # head(10) |>
  dplyr::mutate(
    cor_age = parallel::mclapply(
      X = data,
      FUN = \(.x) {
        tryCatch(
          {
            .x |>
              dplyr::filter(!is.na(Age_new)) |>
              dplyr::filter(!is.na(peakmin)) -> .x
            if (nrow(.x) > 30) {
              cor.test(
                .x$peakmin,
                .x$Age_new,
                method = "spearman",
                use = "pairwise.complete.obs"
              ) |>
                broom::tidy() |>
                dplyr::select(
                  cor = estimate,
                  pval = p.value
                )
            } else {
              tibble::tibble(cor = NA_real_, pval = NA_real_)
            }
          },
          error = function(e) {
            tibble::tibble(cor = NA_real_, pval = NA_real_)
          }
        )
      },
      mc.cores = 10
    )
  ) |>
  tidyr::unnest(cols = cor_age) -> allvariants_celltype_peaks_meta_cor_age

# export(
#   allvariants_celltype_peaks_meta_cor_age,
#   file.path(
#     "/home/liuc9/github/scMOCHA-data/analysis/zzz/plot-celltype-specific-variant",
#     "000-celltype_specific_allvariants_celltype_peaks_meta_cor_age.qs"
#   )
# )
allvariants_celltype_peaks_meta_cor_age |>
  dplyr::filter(variant == "3727T>C") |>
  dplyr::filter(!is.na(cor)) |>
  tidyr::unnest(cols = data) |>
  ggplot(aes(
    x = Age_new,
    y = peakmin,
  )) +
  geom_point(aes(color = disease)) +
  geom_smooth(method = "lm", se = TRUE, color = "black", linetype = "dashed") +
  ggpubr::stat_cor(
    method = "spearman",
    label.x.npc = "left",
    label.y.npc = "top"
  ) +
  theme_bw() +
  labs(
    title = "Correlation between Age and 3727T>C Peak Minimum",
    x = "Age (years)",
    y = "Peak Minimum"
  ) +
  facet_wrap(
    ~celltype,
    nrow = 1
  )


allvariants_celltype_peaks_meta_cor_age |>
  dplyr::filter(pval < 0.05, abs(cor) > 0.3) |>
  dplyr::mutate(
    n = purrr::map_int(
      .x = data,
      .f = nrow
    ),
  ) |>
  dplyr::filter(n > 100) |>
  dplyr::arrange(-cor) |>
  dplyr::slice(1) |>
  tidyr::unnest(cols = data) |>
  ggplot(aes(
    x = Age_new,
    y = peakmin,
  )) +
  geom_point(aes(color = disease))

allvariants_celltype_peaks_meta_cor_age |>
  dplyr::group_by(variant) |>
  dplyr::summarise(
    ncelltype = dplyr::n(),
    ncelltype_significant = sum(cor > 0.3, na.rm = TRUE),
    ncelltype_significant_neg = sum(cor < -0.3, na.rm = TRUE),
    ncelltype_significant_pos = sum(cor > 0.3, na.rm = TRUE),
    .groups = "drop"
  )


# ? peak and haplogroup --------------------------------------------------------------------
allvariants_celltype_peaks_meta |>
  dplyr::filter(
    !is.na(HG2)
  ) |>
  tidyr::nest(
    .by = c(variant, celltype)
  ) |>
  dplyr::mutate(
    cor_haplo = parallel::mclapply(
      X = data,
      FUN = \(.x) {
        tryCatch(
          {
            if (nrow(.x) > 30) {
              kruskal.test(
                peakmin ~ HG2,
                data = .x
              ) |>
                broom::tidy() |>
                dplyr::select(
                  statistic = statistic,
                  pval = p.value
                )
            } else {
              tibble::tibble(statistic = NA_real_, pval = NA_real_)
            }
          },
          error = function(e) {
            tibble::tibble(statistic = NA_real_, pval = NA_real_)
          }
        )
      },
      mc.cores = 10
    )
  ) |>
  tidyr::unnest(cols = cor_haplo) -> allvariants_celltype_peaks_meta_cor_haplo


# export(
#   allvariants_celltype_peaks_meta_cor_haplo,
#   file.path(
#     "/home/liuc9/github/scMOCHA-data/analysis/zzz/plot-celltype-specific-variant",
#     "000-celltype_specific_allvariants_celltype_peaks_meta_cor_haplo.qs"
#   )
# )

allvariants_celltype_peaks_meta_cor_haplo |>
  dplyr::filter(pval < 0.05) |>
  dplyr::arrange(pval)
allvariants_celltype_peaks_meta_cor_haplo |>
  dplyr::filter(
    variant == "13762T>G",
    celltype == "B"
  ) |>
  tidyr::unnest(cols = data) |>
  ggplot(aes(
    x = HG1,
    y = peakmin,
  )) +
  geom_boxplot() +
  geom_point()

# ? sex --------------------------------------------------------------------

allvariants_celltype_peaks_meta |>
  tidyr::nest(
    .by = c(variant, celltype)
  ) |>
  dplyr::mutate(
    cor_sex = parallel::mclapply(
      X = data,
      FUN = \(.x) {
        tryCatch(
          {
            if (nrow(.x) > 50) {
              t.test(
                peakmin ~ sex_pred,
                .x
              ) |>
                broom::tidy() |>
                dplyr::select(
                  statistic = statistic,
                  pval = p.value
                )
            } else {
              tibble::tibble(statistic = NA_real_, pval = NA_real_)
            }
          },
          error = function(e) {
            tibble::tibble(statistic = NA_real_, pval = NA_real_)
          }
        )
      },
      mc.cores = 10
    )
  ) |>
  tidyr::unnest(cols = cor_sex) -> allvariants_celltype_peaks_meta_cor_sex


# export(
#   allvariants_celltype_peaks_meta_cor_sex,
#   file.path(
#     "/home/liuc9/github/scMOCHA-data/analysis/zzz/plot-celltype-specific-variant",
#     "000-celltype_specific_allvariants_celltype_peaks_meta_cor_sex.qs"
#   )
# )

allvariants_celltype_peaks_meta_cor_sex |>
  dplyr::filter(pval < 0.05) |>
  dplyr::arrange(statistic)

allvariants_celltype_peaks_meta_cor_sex |>
  # dplyr::filter(pval < 0.05) |>
  # dplyr::arrange(-statistic) |>
  # dplyr::slice(1) |>
  dplyr::filter(variant == "6191C>T") |>
  dplyr::mutate(
    celltype = gsub(
      "_",
      " ",
      celltype
    )
  ) |>
  dplyr::mutate(
    celltype = factor(
      celltype,
      levels = names(color_celltype)
    )
  ) |>
  tidyr::unnest(cols = data) |>
  ggplot(aes(
    x = sex_pred,
    y = peakmin,
  )) +
  geom_boxplot() +
  geom_point(aes(color = celltype)) +
  scale_color_manual(
    values = color_celltype,
    na.value = "grey50"
  ) +
  ggsignif::geom_signif(
    aes(
      y = peakmin,
    ),
    comparisons = list(
      c("Female", "Male")
    ),
    y_position = 0.8
  ) +
  facet_wrap(
    ~celltype,
    ncol = 8
  )

# ? don't run below --------------------------------------------------------------------

thevariant <- "3727T>C"
tbl_all_hetero_af_cell |>
  dplyr::filter(
    variant == thevariant,
    af > 0,
  ) |>
  dplyr::left_join(
    tbl_all_hetero_altdepth_cell,
    by = c("gseid", "srrid", "variant", "barcode")
  ) |>
  dplyr::left_join(
    tbl_all_hetero_sumdepth_cell,
    by = c("gseid", "srrid", "variant", "barcode")
  ) |>
  dplyr::left_join(
    tbl_barcode,
    by = c("gseid", "srrid", "barcode", "celltype")
  ) |>
  dplyr::select(
    variant,
    srrid,
    barcode,
    celltype,
    af,
    altdepth,
    sumdepth,
    celltype_l2,
    celltype_l3
  ) |>
  dplyr::collect() -> thevariant_data


thevariant_data |>
  tidyr::nest(
    .by = c(variant, srrid, celltype),
    .key = "data"
  ) |>
  dplyr::mutate(
    peaks = parallel::mclapply(
      X = data,
      FUN = \(.x) {
        tryCatch(
          find_cell_density_peak(.x$af),
          error = function(e) {
            NULL
          }
        )
      },
      mc.cores = 5
    )
  ) -> thevariant_data_peaks

thevariant_data_peaks |>
  dplyr::select(-data) |>
  tidyr::unnest(cols = peaks) |>
  dplyr::left_join(
    tbl_meta |> dplyr::collect(),
    by = c("srrid")
  ) -> a

a |>
  dplyr::select(-peaks) |>
  tidyr::nest(
    .by = c(variant, celltype)
  ) |>
  dplyr::mutate(
    cor_age = purrr::map(
      .x = data,
      .f = \(.x) {
        tryCatch(
          {
            .x |>
              dplyr::filter(!is.na(Age_new)) |>
              dplyr::filter(!is.na(peakmin)) -> .x
            if (nrow(.x) > 30) {
              cor.test(
                .x$peakmin,
                .x$Age_new,
                method = "spearman",
                use = "pairwise.complete.obs"
              ) |>
                broom::tidy() |>
                dplyr::select(
                  cor = estimate,
                  pval = p.value
                )
            } else {
              tibble::tibble(cor = NA_real_, pval = NA_real_)
            }
          },
          error = function(e) {
            tibble::tibble(cor = NA_real_, pval = NA_real_)
          }
        )
      }
    )
  ) |>
  tidyr::unnest(cols = cor_age) |>
  dplyr::filter(pval < 0.05) |>
  dplyr::arrange(-cor) -> aa


aa |>
  dplyr::slice(1) |>
  dplyr::select(data) |>
  tidyr::unnest(cols = data) |>
  dplyr::arrange(-peakmin) -> ddd


thevariant_data_peaks |>
  tidyr::unnest(cols = peaks) |>
  dplyr::select(srrid, celltype, npeaks) |>
  tidyr::pivot_wider(
    names_from = celltype,
    values_from = npeaks,
    values_fill = 0
  ) |>
  dplyr::mutate(
    m = purrr::map2_chr(
      .x = Mono,
      .y = DC,
      .f = \(.x, .y) {
        # 1. .x == 2 and .y == 2, M two peaks
        # 2. .x == 1 or .y == 1, M one peak
        # 3. no peaks, M no peaks
        if (.x == 2 && .y == 2) {
          "M two peak"
        } else if (.x == 1 || .y == 1) {
          "M one peak"
        } else {
          "M no peaks"
        }
      }
    )
  ) |>
  dplyr::select(srrid, m) -> thevariant_data_peaks_rank_m

thevariant_data_peaks |>
  tidyr::unnest(cols = peaks) |>
  dplyr::select(srrid, celltype, npeaks) |>
  tidyr::pivot_wider(
    names_from = celltype,
    values_from = npeaks,
    values_fill = 0
  ) |>
  dplyr::filter(
    (CD4_T == 2 | CD8_T == 2 | NK == 2 | B == 2) &
      (Mono == 1)
  ) -> thevariant_data_peaks_rank_m2

thevariant_data_peaks |>
  tidyr::unnest(cols = peaks) |>
  dplyr::select(
    srrid,
    celltype,
    npeaks
  ) |>
  dplyr::group_by(
    srrid,
  ) |>
  dplyr::summarise(
    sum_peaks = sum(npeaks, na.rm = TRUE),
  ) |>
  dplyr::arrange(sum_peaks) |>
  dplyr::filter(sum_peaks > 0) |>
  dplyr::filter(
    srrid %in% thevariant_data_peaks_rank_m2$srrid
  ) -> thevariant_data_peaks_rank


thevariant_data_peaks |>
  dplyr::filter(celltype == "CD4_T") |>
  tidyr::unnest(cols = peaks) |>
  dplyr::arrange(-peakmin) -> thevariant_data_peaks_rank_b

thevariant_data |>
  # dplyr::filter(
  #   srrid %in% thevariant_data_peaks_rank$srrid,
  # ) |>
  dplyr::filter(
    srrid %in% ddd$srrid,
  ) |>
  # dplyr::filter(celltype == "CD4_T") |>
  dplyr::mutate(
    srrid = factor(
      srrid,
      ddd$srrid
    )
  ) |>
  dplyr::mutate(
    celltype = gsub(
      "_",
      " ",
      celltype
    )
  ) |>
  dplyr::mutate(
    celltype = factor(
      celltype,
      levels = names(color_celltype)
    )
  ) |>
  ggplot(aes(
    x = af,
    y = srrid,
  )) +
  ggridges::geom_density_ridges(
    aes(
      fill = celltype,
      color = celltype
    ),
    rel_min_height = 0.01,
    size = 0.1,
    # from = 0,
    # to = 1,
  ) +
  scale_color_manual(
    values = color_celltype,
    na.value = "grey50"
  ) +
  scale_fill_manual(
    values = color_celltype,
    na.value = "grey50"
  ) +
  scale_y_discrete(expand = c(0, 0)) +
  scale_x_continuous(expand = c(0, 0)) +
  ggridges::theme_ridges(
    grid = FALSE
  ) +
  coord_cartesian(clip = "off") +
  theme(
    legend.position = "none",
    axis.text.y = element_blank(),
    axis.line = element_line(color = "black"),
    strip.background = element_rect(
      fill = "white",
      color = "black"
    ),
  ) +
  facet_wrap(
    ~celltype,
    ncol = 8
  ) +
  labs(
    title = paste0(thevariant, ", pisitive correlated with age in DC cells"),
    x = "Allele Frequency",
    y = "Sample"
  ) -> p
ggsave(
  plot = p,
  filename = file.path(
    plotdir,
    paste0(thevariant, ".variant_celltype_peaks.pdf")
  ),
  width = 12,
  height = 8
)

fn_plot_ggdist(
  thevariant = thevariant,
  thegseid = "GSE235050",
  thesrrid = "GSM7493832"
) -> p_ggdist
ggsave(
  plot = p_ggdist,
  filename = file.path(
    plotdir,
    paste0(thevariant, ".variant_celltype_peaks.ggdist.pdf")
  ),
  width = 12,
  height = 8
)

# footer ------------------------------------------------------------------

# future: :plan(future: :sequential)

# save image --------------------------------------------------------------
DBI::dbDisconnect(conn, shutdown = TRUE)
