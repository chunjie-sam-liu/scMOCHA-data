#!/usr/bin/env Rscript
# Metainfo ----------------------------------------------------------------
# @AUTHOR: Chun-Jie Liu
# @CONTACT: chunjie.sam.liu.at.gmail.com
# @DATE: 2025-06-06 10:58:26
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

gseid_srrid_ks_load <- import(
  file.path(
    ks_test_dir,
    "a_gseid_srrid_ks_load.nocellaf.qs"
  )
)


gseid_srrid_ks_load |>
  dplyr::filter(
    p.value < 0.05,
    statistic > 25
  ) |>
  dplyr::count(variant) |>
  dplyr::arrange(-n) |>
  ggplot(aes(
    x = n
  )) +
  geom_histogram(
    aes(y = after_stat(density)),
    bins = 100,
    fill = "grey50",
    color = "black",
    alpha = 0.5
  )


ALLVARIANTS <- import(file.path(
  "/home/liuc9/github/scMOCHA-data/analysis/zzz/clean-data/", "all_variant.qs"
)) |>
  dplyr::filter(
    issomatic == "heteroplasmic"
  )

gseid_srrid_ks_load |>
  dplyr::filter(
    p.value < 0.05,
    statistic > 25
  ) |>
  dplyr::filter(
    variant %in% ALLVARIANTS$variant
  ) ->
gseid_srrid_ks_load_variant

export(
  gseid_srrid_ks_load_variant,
  file.path(
    "/home/liuc9/github/scMOCHA-data/analysis/zzz/clean-data",
    "celltype_specific_variant.qs"
  )
)
export(
  gseid_srrid_ks_load_variant,
  file.path(
    "/home/liuc9/github/scMOCHA-data/analysis/zzz/clean-data",
    "celltype_specific_variant.csv"
  ),
  format = "both"
)


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

# src ---------------------------------------------------------------------
source("./analysis/00-colors.R")

# function ----------------------------------------------------------------

fn_plot_ggdist <- function(
    thevariant,
    thegseid,
    thesrrid) {
  library(ggdist)
  tbl_all_hetero_af_cell |>
    dplyr::filter(
      # gseid == thegseid,
      srrid == thesrrid,
      variant == thevariant,
      af > 0
    ) |>
    dplyr::collect() ->
  .d
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
    ) ->
  forplot_

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
    thesrrid) {
  tbl_all_hetero_af_cell |>
    dplyr::filter(
      # gseid == thegseid,
      srrid == thesrrid,
      variant == thevariant,
      af > 0
    ) |>
    dplyr::collect() ->
  .d
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
    ) ->
  forplot_

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
    thecelltype_level) {
  celltypedetail <- import("/mnt/isilon/u01_project/large-scale/liuc9/raw/{thegseid}/final/{thesrrid}/sc_azimuth_celltype.csv" |> glue::glue())

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
    ) ->
  thevariant_data


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
    ) ->
  forplot

  levels(forplot$plotcelltype)


  color_celltype_detail <- log(seq(1, exp(1), length.out = length(levels(forplot$plotcelltype)))) |>
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
      title = "{thecelltype}-{thecelltype_level}-{thevariant}\n({thegseid}-{thesrrid})" |> glue::glue(),
      x = "Allele Frequency",
      y = "Cell Type"
    )
}

fn_plot_joy_celltype_detail <- function(
    thevariant,
    thegseid,
    thesrrid) {
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
  ) ->
  thevariant_celltype_df

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
    ) ->
  plot_thevariant_celltype_list
}


# body --------------------------------------------------------------------



# ? plot ks statistic--------------------------------------------------------------------


gseid_srrid_ks_load |>
  dplyr::filter(p.value < 0.05) |>
  ggplot(
    aes(x = statistic)
  ) +
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
    x = "KS statistic",
    y = "Density",
    title = "Distribution of KS statistic for all variants"
  ) ->
plot_ks_statistic
ggsave(
  file.path(
    plotdir,
    "ks_statistic_distribution.pdf"
  ),
  plot = plot_ks_statistic,
  width = 8,
  height = 6
)



# ? find examples --------------------------------------------------------------------

gseid_srrid_ks_load |>
  dplyr::filter(
    p.value < 0.05,
    statistic > 20
  ) |>
  dplyr::group_by(
    variant
  ) |>
  dplyr::summarise(
    n = dplyr::n(),
    mean_statistic = mean(statistic, na.rm = TRUE),
    mean_p_value = mean(p.value, na.rm = TRUE),
  ) |>
  dplyr::mutate(
    mean_log10p = -log10(mean_p_value),
  ) ->
variant_count_statistic

variant_count_statistic$mean_log10p |> summary()
variant_count_statistic$mean_statistic |> summary()

variant_count_statistic |>
  ggplot(aes(
    x = mean_statistic,
    y = mean_log10p,
    size = n
  )) +
  ggpointdensity::geom_pointdensity(
    adjust = 0.01,
    # show.legend = FALSE,
  ) +
  # geom_point() +
  viridis::scale_color_viridis() +
  scale_x_continuous(
    limits = c(20, 170),
    labels = scales::number,
    breaks = seq(20, 170, by = 10),
    expand = expansion(add = c(0.01, 0))
  ) +
  scale_y_continuous(
    limits = c(2, 15),
    labels = scales::number,
    breaks = seq(2, 15, by = 1),
    expand = expansion(add = c(0.01, 0))
  ) +
  ggrepel::geom_text_repel(
    data = variant_count_statistic |>
      dplyr::filter(
        variant %in% c(
          "7833T>C",
          # sort by mean_statistic
          "3030A>G",
          "7430A>C",
          "4175G>A",
          "3727T>C",
          "7418C>A",
          "6409T>C",
          "6669C>G",
          "7583T>G",
          "7582C>G",
          "929A>C",
          # sort by mean_log10p
          "4886C>T",
          "1082A>G",
          # "15213T>C",
          "3173G>A",
          "6669C>G",
          "13956A>G",
          "14063T>C",
          "8849T>C",
          "8285C>A",
          "4794G>A",
          "11502T>C",
          # sort by n
          "4175G>A",
          "2289G>T",
          "10645T>G",
          "3584A>C",
          "3030A>G",
          "3520A>C",
          "2193T>A",
          "9076A>C",
          "8072T>G",
          "3577A>C"
        )
      ),
    aes(label = variant),
    size = 3,
    max.overlaps = 20,
    show.legend = FALSE,
    nudge_x = 5,
    nudge_y = 2
  ) +
  labs(
    x = "Mean KS Statistic",
    y = "-log10(Mean P-value)",
    title = "Variant Count vs. Mean KS Statistic and P-value"
  ) +
  theme_bw() +
  guides(
    size = guide_legend(title = "Variant Count"),
    color = guide_colorbar(title = "Density")
  ) +
  theme(
    plot.title = element_text(hjust = 0.5),
  ) ->
plot_variant_count_statistic
ggsave(
  file.path(
    plotdir,
    "ks_variant_count_statistic.pdf"
  ),
  plot = plot_variant_count_statistic,
  width = 11,
  height = 6
)


variant_count_statistic |> dplyr::arrange(-n)

# ? example --------------------------------------------------------------------
variant_list <- c(
  "7833T>C",
  # sort by mean_statistic
  "3030A>G",
  "7430A>C",
  "4175G>A",
  "3727T>C",
  "7418C>A",
  "6409T>C",
  "6669C>G",
  "7583T>G",
  "7582C>G",
  "929A>C",
  # sort by mean_log10p
  # "4886C>T",
  "1082A>G",
  # "15213T>C",
  "3173G>A",
  "6669C>G",
  "13956A>G",
  "14063T>C",
  "8849T>C",
  "8285C>A",
  "4794G>A",
  "11502T>C",
  # sort by n
  "4175G>A",
  "2289G>T",
  "10645T>G",
  "3584A>C",
  "3030A>G",
  "3520A>C",
  "2193T>A",
  "9076A>C",
  "8072T>G",
  "3577A>C"
)
thevariant <- "7833T>C"
thevariant <- "3727T>C"

variant_count_statistic |>
  dplyr::arrange(-mean_statistic) |>
  dplyr::slice(1:400) |>
  dplyr::pull(variant) ->
variant_list


length(variant_list)
variant_list |>
  parallel::mclapply(
    function(thevariant) {
      gseid_srrid_ks_load |>
        dplyr::filter(variant == thevariant) |>
        dplyr::filter(
          p.value < 0.05,
        ) |>
        dplyr::slice(1) ->
      thevariant_data
      .gseid <- thevariant_data$gseid[1]
      .srrid <- thevariant_data$srrid[1]
      .variant <- thevariant_data$variant[1]
      fn_plot_joy(
        thevariant = .variant,
        thegseid = .gseid,
        thesrrid = .srrid
      )
    },
    mc.cores = 10
  ) ->
plot_variants_list

plot_variants_list |>
  wrap_plots(ncol = 20) +
  plot_layout(
    guides = "collect",
  ) ->
plot_variants_joy
ggsave(
  plot = plot_variants_joy,
  filename = file.path(
    plotdir,
    "joy_variants400.pdf"
  ),
  width = 100,
  height = 100,
  limitsize = FALSE
)



# ? single variant joy --------------------------------------------------------------------

thevariant <- "3173G>A"
thevariant <- "7833T>C"
thevariant <- "3727T>C"
thevariants <- c(
  "3173G>A",
  "7833T>C",
  "3727T>C",
  "7159T>C",
  "2666T>C",
  "2193T>A",
  "1474G>A",
  "10500G>A",
  "2666T>C",
  "7833T>C",
  "2763T>C",
  "8005T>C",
  "2871T>C",
  "7757G>A"
)

thevariants |>
  purrr::map(
    function(thevariant) {
      gseid_srrid_ks_load |>
        dplyr::filter(variant == thevariant) |>
        dplyr::filter(
          p.value < 0.05,
          # statistic > 100
        ) |>
        dplyr::mutate(
          p = parallel::mcmapply(
            .variant = variant,
            .gseid = gseid,
            .srrid = srrid,
            FUN = function(.variant, .gseid, .srrid) {
              fn_plot_joy(
                thevariant = .variant,
                thegseid = .gseid,
                thesrrid = .srrid
              )
            },
            mc.cores = 10,
            SIMPLIFY = FALSE
          )
        ) ->
      plot_thevariant_sample_list


      plot_thevariant_sample_list |>
        dplyr::slice(1:100) |>
        dplyr::pull(p) |>
        wrap_plots(ncol = 10) +
        plot_layout(
          guides = "collect",
        ) ->
      plot_thevariant_sample_joy

      ggsave(
        plot = plot_thevariant_sample_joy,
        filename = file.path(
          plotdir,
          "{thevariant}_joy_sample.pdf" |> glue::glue()
        ),
        width = 50,
        height = 50,
        limitsize = FALSE
      )
    }
  )




# ? GSE235050_GSM7493832 azimuth.rda --------------------------------------------------------------------


thevariants <- c(
  "3173G>A",
  "7833T>C",
  "3727T>C",
  "7159T>C",
  "2666T>C",
  "2193T>A",
  "1474G>A",
  "10500G>A",
  "2666T>C",
  "7833T>C",
  "2763T>C",
  "8005T>C",
  "2871T>C"
)

thevariants |>
  purrr::map(
    function(thevariant) {
      gseid_srrid_ks_load |>
        dplyr::filter(variant == thevariant) |>
        dplyr::filter(
          p.value < 0.05,
          # statistic > 100
        ) ->
      thevariant_data
      .gseid <- thevariant_data$gseid[1]
      .srrid <- thevariant_data$srrid[1]
      .variant <- thevariant_data$variant[1]

      fn_plot_joy_celltype_detail(.variant, .gseid, .srrid) ->
      plot_thevariant_celltype_list

      plot_thevariant_celltype_list |>
        dplyr::pull(p) |>
        wrap_plots(
          ncol = 4
        ) +
        plot_layout(
          guides = "collect",
        ) ->
      plot_thevariant_celltype_joy
      ggsave(
        plot = plot_thevariant_celltype_joy,
        filename = file.path(
          plotdir,
          "{thevariant}_joy_celltype.pdf" |> glue::glue()
        ),
        width = 23,
        height = 12,
        limitsize = FALSE
      )
    }
  )






# ? dont run below --------------------------------------------------------------------

# ? META --------------------------------------------------------------------

META <- import("/home/liuc9/github/scMOCHA-data/analysis/zzz/clean-data/gse_dataset_metadata_full.sex_pred.qs") |>
  dplyr::select(
    gseid, srrid, Age_new, Age_group,
    Haplogroup,
    disease, Chemistry, sex_pred
  ) |>
  dplyr::mutate(
    Haplogroup = purrr::map_chr(
      .x = Haplogroup,
      .f = \(.x) {
        # if (stringr::str_starts(.x, "L")) {
        #   gsub("L", "L0", .x)
        # }
        gsub("\\d+.*", "", .x)
      }
    )
  )


thevariant <- "3727T>C"
thegseid <- "GSE235050"
thesrrid <- "GSM7493832"

gseid_srrid_ks_load |>
  dplyr::filter(variant == thevariant) |>
  dplyr::left_join(
    META,
    by = c("gseid", "srrid")
  ) ->
thevariant_meta

library(plotly)
thevariant_meta |>
  dplyr::filter(!is.na(statistic)) |>
  dplyr::filter(statistic > 0) |>
  # dplyr::filter(
  #   Chemistry == "SC5P-PE"
  # ) |>
  # dplyr::filter(
  #   disease %in% c("Healthy", "Alzheimer's Disease")
  # ) |>
  ggplot(aes(
    x = disease,
    y = statistic,
    label = gseid_srrid,
  )) +
  geom_point()

fn_plot_joy(
  thevariant = thevariant,
  # thegseid = thegseid,
  thesrrid = "GSM7493832"
)



tbl_all_hetero_af_cell |>
  dplyr::filter(
    srrid == "GSM7493832",
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
  ) ->
tbl_thevariant_data

# footer ------------------------------------------------------------------

# future: :plan(future: :sequential)

# save image --------------------------------------------------------------
DBI::dbDisconnect(conn, shutdown = TRUE)
