#!/usr/bin/env Rscript
# Metainfo ----------------------------------------------------------------
# @AUTHOR: Chun-Jie Liu
# @CONTACT: chunjie.sam.liu.at.gmail.com
# @DATE: 2025-11-29 23:01:04
# @DESCRIPTION: this script is used for ...

# Library -----------------------------------------------------------------

load_pkg(jutils)

# args --------------------------------------------------------------------

# s: string, i: integer, f: float, !: boolean, @: array, %: list
GetoptLong.options(help_style = "two-column")
VERSION = "v0.0.1"

# default: default value specified here.

verbose = TRUE
# gseid = "GSE235050" # no default gseid, current gseid is for testing
basedir = "/mnt/isilon/u01_project/large-scale/liuc9/raw"

GetoptLong(
  "gseid=s",
  "GSE ID",
  "basedir=s",
  "GSE base directory",
  "verbose!",
  "print messages"
)


logger::log_threshold(logger::TRACE)
logger::log_layout(logger::layout_glue_colors)
# header ------------------------------------------------------------------

# load data ---------------------------------------------------------------
datadir <- path(basedir, gseid)
finaldir <- path(datadir, "final")
outdir <- path(datadir, "out")
dir_create(outdir)

srrid_list <- readr::read_lines(
  file = path(
    datadir,
    "{gseid}.srrid.list" |> glue::glue()
  )
)

# Technical positions -------------------------------------------------------
POS_RNA_EDITING = c(
  585,
  1610,
  3238,
  4271,
  5520,
  7526,
  8303, # tRNA p9
  9999,
  10413,
  12146,
  12274,
  14734,
  15896, # tRNA p9
  295,
  2617,
  13710 # RNA editing
)


POS_MISSALIGNMENT_ERROR = c(
  66:71,
  300:316,
  513:525,
  3106:3107,
  12418:12425,
  16182:16194
)

CUTOFF_NOTRELIABLE = 10 # Cells less than this
CUTOFF_MIN_CELLS = 10
CUTOFF_MIN_READS = 10

# Biological excluding
ETHNICITY = TRUE
HOMOPLASMIC_N_CELLTYPE = 7
CUTOFF_HOMOPLASMIC = 0.95
CUTOFF_HETEROPLASMIC = 0.05
CUTOFF_SOMATIC_IGNORE_N_CELLS = 2
CUTOFF_SOMATIC_MIN_N_CELLS = 10
# load conn ---------------------------------------------------------------

# src ---------------------------------------------------------------------

# function ----------------------------------------------------------------
# fn_plot_clusteraf(.srrdir, af_list, .somatic)
fn_plot_clusteraf <- function(
  .srrdir,
  af_list,
  .somatic
) {
  .somatic |>
    tibble::enframe(name = "variant_type", value = "variant") |>
    tidyr::unnest(cols = variant) -> .somatic_variant_type

  .somatic_variant_type |>
    dplyr::count(variant_type) -> n_hete_homo

  af_list$clusteraf |>
    dplyr::select(variant, celltype, clusteraf) |>
    tidyr::pivot_wider(
      names_from = celltype,
      values_from = clusteraf
    ) |>
    dplyr::left_join(
      af_list$bulkaf |> dplyr::select(variant, Bulk = bulkaf),
      by = "variant"
    ) |>
    dplyr::filter(
      variant %in% .somatic_variant_type$variant
    ) |>
    dplyr::mutate(
      variant = forcats::fct_reorder(variant, Bulk, .desc = TRUE)
    ) |>
    tidyr::pivot_longer(
      cols = -variant,
      names_to = "celltype",
      values_to = "clusteraf"
    ) -> .clusteraf_remain

  af_list$clusteraf |>
    dplyr::select(variant, celltype, cluster_total_reads = total_reads) |>
    tidyr::pivot_wider(
      names_from = celltype,
      values_from = cluster_total_reads
    ) |>
    dplyr::left_join(
      af_list$bulkaf |> dplyr::select(variant, Bulk = total_reads),
      by = "variant"
    ) |>
    dplyr::filter(
      variant %in% .somatic_variant_type$variant
    ) |>
    dplyr::mutate(
      variant = factor(variant, levels = levels(.clusteraf_remain$variant))
    ) |>
    tidyr::pivot_longer(
      cols = -variant,
      names_to = "celltype",
      values_to = "cluster_total_reads"
    ) -> .clusteraf_remain_total_reads

  variant_type_labels <- c(
    "homo" = "Homoplasmic",
    "haplo" = "Ethnicity",
    "hete" = "Heteroplasmic",
    "somatic" = "Somatic",
    "n_cells" = "Low Cells",
    "editing" = "RNA Editing",
    "excluding_pos" = "Mis-alignment"
  )

  # -> forplot

  .clusteraf_remain_total_reads |>
    dplyr::mutate(
      `Depth(log2)` = log2(cluster_total_reads + 1)
    ) |>
    ggplot(aes(
      x = celltype,
      y = variant,
      fill = `Depth(log2)`
    )) +
    geom_tile() +
    scale_fill_gradient(low = "white", high = "gold", na.value = "grey90") +
    theme_classic() +
    labs(
      x = "Depth(log2)",
      y = "Sample"
    ) +
    theme(
      plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
      axis.text.x = element_text(face = "bold", size = 12),
      axis.title.x = element_text(face = "bold", size = 12)
    ) -> p1

  .clusteraf_remain |>
    ggplot(aes(
      x = celltype,
      y = variant,
      fill = clusteraf
    )) +
    geom_tile() +
    scale_fill_gradient(low = "white", high = "red", na.value = "grey90") +
    theme_classic() +
    labs(
      x = "Allele Frequency",
      y = "Sample"
    ) +
    theme(
      plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
      axis.text.x = element_text(face = "bold", size = 12),
      axis.title.x = element_text(face = "bold", size = 12)
    ) -> p2

  .somatic_variant_type |>
    dplyr::mutate(v = 1) |>
    tidyr::pivot_wider(
      names_from = variant_type,
      values_from = v,
      values_fill = NA_integer_
    ) |>
    dplyr::mutate(
      variant = factor(
        variant,
        levels = levels(.clusteraf_remain_total_reads$variant)
      )
    ) |>
    tidyr::pivot_longer(
      cols = -variant,
      names_to = "variant_type",
      values_to = "value"
    ) |>
    ggplot(
      aes(
        x = variant_type,
        y = variant,
        fill = value
      )
    ) +
    geom_tile() +
    scale_x_discrete(
      limits = names(variant_type_labels),
      labels = variant_type_labels
    ) +
    scale_fill_gradient(
      name = "Presence",
      low = "white",
      high = "blue"
    ) +
    theme_classic() +
    labs(
      x = "Variant Type",
      y = "Sample"
    ) +
    theme(
      plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
      axis.text.x = element_text(face = "bold", size = 10),
      axis.title.x = element_text(face = "bold", size = 12)
    ) -> p3

  wrap_plots(
    p1,
    p2,
    p3,
    ncol = 3,
    widths = c(1.2, 1.2, 1.5),
    guides = "collect"
  ) +
    plot_annotation(
      title = glue::glue(
        "{length(unique(.clusteraf_remain$variant))} Variant Allele Frequency and Type\n{.srrdir}"
      ),
      theme = theme(
        plot.title = element_text(size = 14, face = "bold", hjust = 0.5)
      )
    ) -> p_collect
  {
    # outdir <- "/home/liuc9/github/scMOCHA-data/analysis/high-res/MANUSCRIPTFIGURES/notuse"
    ggsave(
      filename = glue::glue("Variant-Type.pdf"),
      plot = p_collect,
      path = .srrdir,
      width = 20,
      height = 15
    )
  }
  p_collect
}

fn_cell_cluster_bulk_af <- function(.srrdir) {
  .cellaf <- import(
    path(.srrdir, "variant_info_from_heatmap.qs")
  ) |>
    dplyr::mutate(
      ref = gsub("\\d+|>.*", "", variant),
      alt = gsub("\\d+.*>", "", variant)
    )

  .cellaf |>
    dplyr::select(
      -c(af, variant_type, variant_in_cell_cluster)
    ) -> .tt_

  .tt_ |>
    tidyr::nest(
      .key = "reads",
      .by = c(variant, celltype, ref, alt)
    ) |>
    dplyr::mutate(
      clusteraf = purrr::pmap(
        .l = list(reads, ref, alt),
        .f = function(reads, ref, alt) {
          .total_reads <- sum(reads$depth, na.rm = TRUE)
          .alt_reads <- sum(
            c(
              reads[[paste0(alt, "F")]],
              reads[[paste0(alt, "R")]]
            ),
            na.rm = TRUE
          )
          .ref_reads <- sum(
            c(
              reads[[paste0(ref, "F")]],
              reads[[paste0(ref, "R")]]
            ),
            na.rm = TRUE
          )
          af <- sum(
            c(
              reads[[paste0(alt, "F")]]
            ),
            na.rm = TRUE
          )
          ar <- sum(
            c(
              reads[[paste0(alt, "R")]]
            ),
            na.rm = TRUE
          )
          rf <- sum(
            c(
              reads[[paste0(ref, "F")]]
            ),
            na.rm = TRUE
          )
          rr <- sum(
            c(
              reads[[paste0(ref, "R")]]
            ),
            na.rm = TRUE
          )
          tibble::tibble(
            clusteraf = .alt_reads / .total_reads,
            reff = rf,
            refr = rr,
            altf = af,
            altr = ar,
            total_reads = .total_reads,
          ) -> .res

          if (.total_reads == 0) {
            return(
              tibble::tibble(
                clusteraf = NA_real_,
                reff = NA_integer_,
                refr = NA_integer_,
                altf = NA_integer_,
                altr = NA_integer_,
                total_reads = NA_integer_
              )
            )
          } else {
            return(.res)
          }
        }
      )
    ) |>
    dplyr::select(variant, celltype, clusteraf) |>
    tidyr::unnest(cols = clusteraf) -> .clusteraf

  .tt_ |>
    tidyr::nest(
      .key = "reads",
      .by = c(variant, ref, alt)
    ) |>
    dplyr::mutate(
      bulkaf = purrr::pmap(
        .l = list(reads, ref, alt),
        .f = function(reads, ref, alt) {
          # reads <- a$reads[[1]]
          # ref <- a$ref[[1]]
          # alt <- a$alt[[1]]
          .total_reads <- sum(reads$depth, na.rm = TRUE)
          .alt_reads <- sum(
            c(
              reads[[paste0(alt, "F")]],
              reads[[paste0(alt, "R")]]
            ),
            na.rm = TRUE
          )
          .ref_reads <- sum(
            c(
              reads[[paste0(ref, "F")]],
              reads[[paste0(ref, "R")]]
            ),
            na.rm = TRUE
          )
          af <- sum(
            c(
              reads[[paste0(alt, "F")]]
            ),
            na.rm = TRUE
          )
          ar <- sum(
            c(
              reads[[paste0(alt, "R")]]
            ),
            na.rm = TRUE
          )
          rf <- sum(
            c(
              reads[[paste0(ref, "F")]]
            ),
            na.rm = TRUE
          )
          rr <- sum(
            c(
              reads[[paste0(ref, "R")]]
            ),
            na.rm = TRUE
          )

          tibble::tibble(
            bulkaf = .alt_reads / .total_reads,
            reff = rf,
            refr = rr,
            altf = af,
            altr = ar,
            total_reads = .total_reads,
          ) -> .res

          if (.total_reads == 0) {
            return(
              tibble::tibble(
                bulkaf = NA_real_,
                reff = NA_integer_,
                refr = NA_integer_,
                altf = NA_integer_,
                altr = NA_integer_,
                total_reads = NA_integer_
              )
            )
          } else {
            return(.res)
          }
        }
      )
    ) |>
    dplyr::select(variant, bulkaf) |>
    tidyr::unnest(cols = bulkaf) -> .bulkaf

  list(
    cellaf = .cellaf,
    clusteraf = .clusteraf,
    bulkaf = .bulkaf
  )
}

fn_depth_all <- function(.srrdir) {
  .depth_read <- data.table::fread(
    path(.srrdir, "possorted_genome_bam.MT.depth"),
    col.names = c("chrom", "pos", "depth")
  )
  .depth_cell <- data.table::fread(
    path(.srrdir, "cell.coverage.txt.gz"),
    col.names = c("pos", "barcode", "depth")
  )
  .barcode_celltype <- data.table::fread(
    path(.srrdir, "barcode_cluster.tsv"),
    col.names = c("barcode", "tag", "celltype")
  )
  .depth_cluster <- .depth_cell |>
    dplyr::left_join(
      .barcode_celltype |> dplyr::select(-tag),
      by = "barcode"
    ) |>
    dplyr::group_by(pos, celltype) |>
    dplyr::summarise(
      depth = sum(depth, na.rm = T)
    ) |>
    dplyr::ungroup() |>
    dplyr::filter(!is.na(celltype)) |>
    dplyr::arrange(celltype, pos)
  .depth <- .depth_cluster |>
    dplyr::group_by(pos) |>
    dplyr::summarise(
      depth = sum(depth, na.rm = T)
    )
  list(
    depth_read = .depth_read,
    depth_cluster = .depth_cluster,
    depth = .depth
  )
}

fn_n_cell <- function(.cellaf) {
  # .cellaf <- af_list$cellaf
  .cellaf |>
    dplyr::filter(af >= CUTOFF_HETEROPLASMIC) |>
    dplyr::filter(depth >= CUTOFF_MIN_READS) |>
    dplyr::count(variant) |>
    dplyr::filter(n < CUTOFF_NOTRELIABLE) |>
    dplyr::pull(variant)
}

fn_homo_hete <- function(.haplo_variant_remain, .clusteraf) {
  # .haplo_variant_remain
  # .clusteraf <- af_list$clusteraf

  .clusteraf |>
    dplyr::filter(
      variant %in% .haplo_variant_remain$variant
    ) |>
    dplyr::filter(total_reads >= CUTOFF_MIN_READS) |>
    dplyr::select(variant, celltype, clusteraf) -> .clusteraf_remain

  .clusteraf_remain |>
    dplyr::count(
      variant,
      name = "n_celltypes_have_reads"
    ) -> .clusteraf_remain_reads

  # .clusteraf_remain |> dplyr::filter(variant == "7584T>G")

  .clusteraf_remain |>
    dplyr::group_by(variant) |>
    dplyr::summarise(
      mean_clusteraf = mean(clusteraf, na.rm = TRUE)
    ) -> .clusteraf_remain_mean

  .clusteraf_remain |>
    dplyr::select(variant, celltype, clusteraf) |>
    dplyr::mutate(
      ishighaf = clusteraf > CUTOFF_HOMOPLASMIC
    ) |>
    dplyr::filter(ishighaf) |>
    dplyr::count(variant, ishighaf, name = "n_celltypes_have_highaf") |>
    dplyr::left_join(
      .clusteraf_remain_reads,
      by = "variant"
    ) |>
    dplyr::left_join(
      .clusteraf_remain_mean,
      by = "variant"
    ) -> .clusteraf_remain_homo_test

  # .clusteraf_remain_homo_test |>
  #   dplyr::filter(
  #     n_celltypes_have_highaf < n_celltypes_have_reads
  #   )

  .clusteraf_remain_homo_test |>
    dplyr::filter(
      n_celltypes_have_highaf >= n_celltypes_have_reads |
        mean_clusteraf >= CUTOFF_HOMOPLASMIC
    ) |>
    dplyr::pull(variant) -> .homo_variant

  .hete_variant <- setdiff(unique(.clusteraf_remain$variant), .homo_variant)

  list(
    homo_variant = .homo_variant,
    hete_variant = .hete_variant
  )
}

fn_somatic <- function(hete, .cellaf) {
  # .cellaf <- af_list$cellaf
  .cellaf |>
    dplyr::filter(
      variant %in% hete
    ) |>
    dplyr::filter(depth >= CUTOFF_MIN_READS) |>
    dplyr::filter(
      variant_type %in% c("colorful", "black")
    ) |>
    dplyr::select(
      variant,
      barcode,
      af,
      depth,
      variant_type,
      celltype
    ) -> .hete

  # .hete |>
  #   dplyr::filter(variant == "3548T>G") |>
  #   dplyr::count(variant, celltype, variant_type) |>
  #   tidyr::pivot_wider(
  #     names_from = variant_type,
  #     values_from = n
  #   )

  .hete |>
    dplyr::count(variant, celltype, variant_type) -> .hete_n_cells

  .hete_n_cells |>
    tidyr::pivot_wider(
      names_from = variant_type,
      values_from = n
    ) -> .hete_n_cells_wide

  .n_colnames <- length(intersect(
    c("black", "colorful"),
    colnames(.hete_n_cells_wide)
  ))
  if (.n_colnames != 2) {
    return(character())
  }

  .hete_n_cells_wide |>
    tidyr::nest(
      .key = "data",
      .by = variant
    ) |>
    dplyr::mutate(
      n = purrr::map(
        .x = data,
        .f = \(.x) {
          data.table(
            n_celltypes_have_reads = sum(
              .x$black >= CUTOFF_MIN_CELLS,
              na.rm = TRUE
            ),
            n_celltypes_have_variant = sum(
              .x$colorful >= CUTOFF_SOMATIC_IGNORE_N_CELLS,
              na.rm = TRUE
            ),
            n_sufficient_cells = sum(
              .x$colorful >= CUTOFF_SOMATIC_MIN_N_CELLS,
              na.rm = TRUE
            )
          )
        }
      )
    ) |>
    tidyr::unnest(cols = n) -> .hete_d

  .hete_d |>
    dplyr::filter(
      n_celltypes_have_variant > 0,
      n_celltypes_have_reads > 0,
      n_sufficient_cells > 0
    ) |>
    dplyr::filter(
      n_celltypes_have_reads - n_celltypes_have_variant > 1
    ) -> .dd
  # .dd$data[[4]]
  .dd$variant
}

fn_variant_classification <- function(.srrdir, .haplo_variant, af_list) {
  # don't use raw variant_somatic.rds file, recompute here
  # .somatic <- import(
  #   path(.srrdir, "variant_somatic.rds")
  # )
  excluding_pos <- .haplo_variant |>
    dplyr::filter(
      Position %in%
        c(
          POS_MISSALIGNMENT_ERROR
        )
    ) |>
    dplyr::pull(variant)

  editing <- .haplo_variant |>
    dplyr::filter(
      Position %in% POS_RNA_EDITING
    ) |>
    dplyr::pull(variant)

  n_cells <- fn_n_cell(af_list$cellaf)

  .haplo_variant_remain <- .haplo_variant |>
    dplyr::filter(
      !variant %in% c(excluding_pos, editing, n_cells)
    )

  haplo <- .haplo_variant_remain |>
    dplyr::filter(Haplogroup != "") |>
    dplyr::pull(variant)

  homo_hete <- fn_homo_hete(.haplo_variant_remain, af_list$clusteraf)

  homo <- homo_hete$homo_variant
  hete <- homo_hete$hete_variant
  somatic <- fn_somatic(hete, af_list$cellaf)

  list(
    homo = homo,
    hete = hete,
    somatic = somatic,
    haplo = haplo,
    n_cells = n_cells,
    editing = editing,
    excluding_pos = excluding_pos
  )
}

# body --------------------------------------------------------------------

#
#
# scmocha.out --------------------------------------------------------------------
#
#

tibble::tibble(
  srrid = srrid_list
) |>
  dplyr::mutate(
    srrdir = path(finaldir, srrid)
  ) |>
  dplyr::mutate(
    dir_exists = file_exists(srrdir)
  ) -> srr_out


srr_out |>
  dplyr::mutate(
    cell_stats = parallel::mclapply(
      X = srrdir,
      FUN = purrr::safely(\(.srrdir) {
        log_info(
          "Start processing {gseid} - {srrid}",
          srrid = basename(.srrdir),
          gseid = basename(dirname(dirname(.srrdir)))
        )
        {
          # gseid <- "GSE235050"
          # srrid <- "GSM7493841"
          # .srrdir <- path(
          #   glue(
          #     "/mnt/isilon/u01_project/large-scale/liuc9/raw/{gseid}/final/{srrid}"
          #   )
          # )
        }
        if (!file_exists(.srrdir)) {
          return(NULL)
        }

        .chemistry <- import(
          path(.srrdir, "chemistry.csv")
        ) |>
          dplyr::pull(name)

        .metrics <- import(
          path(.srrdir, "metrics_summary.csv")
        ) |>
          purrr::map_dfr(~ as.numeric(gsub("[,%]", "", .x)))

        .cs <- import(
          path(.srrdir, "qc_cell_stats.xlsx")
        )

        .celltype_ratio <- import(
          path(.srrdir, "celltype_ratio.tsv")
        )

        depth_list <- fn_depth_all(.srrdir)

        .cva <- import(
          ifelse(
            file_exists(path(.srrdir, "variant_annotation.tsv")),
            path(.srrdir, "variant_annotation.tsv"),
            path(.srrdir, "cell_variant_annotation.tsv")
          )
        ) |>
          dplyr::mutate(
            variant = glue::glue("{Position}{Ref}>{Alt}")
          )

        .hetero <- import(
          path(.srrdir, "cluster.cell_heteroplasmic_df.tsv.gz")
        ) |>
          dplyr::rename(celltype = V1) |>
          tidyr::gather(-celltype, key = variant, value = af) |>
          dplyr::filter(variant %in% .cva$variant)

        .cov <- depth_list$depth_cluster |>
          dplyr::filter(pos %in% .cva$Position)

        .haplo_variant <- import(
          path(.srrdir, "violin_haplo_variant.csv")
        )

        .haplo_violin <- import(
          path(.srrdir, "violin_haplo_forplot.csv")
        )

        af_list <- fn_cell_cluster_bulk_af(.srrdir)

        .somatic <- fn_variant_classification(
          .srrdir,
          .haplo_variant,
          af_list
        )

        log_success(
          "Finished processing {gseid} - {srrid}",
          srrid = basename(.srrdir),
          gseid = basename(dirname(dirname(.srrdir)))
        )

        data.table::data.table(
          chemistry = .chemistry,
          metrics = list(.metrics),
          cell_stats = list(.cs),
          depth_read = list(depth_list$depth_read),
          depth_cluster = list(depth_list$depth_cluster),
          depth = list(depth_list$depth),
          celltype_ratio = list(.celltype_ratio),
          anno = list(.cva),
          hetero = list(.hetero),
          coverage = list(.cov),
          haplo_variant = list(.haplo_variant),
          haplo_violin = list(.haplo_violin),
          somatic_variant = list(.somatic),
          cellaf = list(af_list$cellaf),
          clusteraf = list(af_list$clusteraf),
          bulkaf = list(af_list$bulkaf)
        )
      }),
      mc.cores = 20
    )
  ) |>
  dplyr::mutate(
    cell_stats = purrr::map(cell_stats, "result")
  ) |>
  tidyr::unnest(cols = cell_stats) |>
  as.data.table() -> srr_out_cell_stats

if (nrow(srr_out_cell_stats) == 0) {
  log_warn("{gseid} has no valid srrid processed, exiting now." |> glue::glue())
  quit(save = "no", status = 0)
}

log_success("{gseid} save to {outdir}/{gseid}.scmocha.out.qs" |> glue::glue())

export(
  srr_out_cell_stats,
  path(
    outdir,
    "{gseid}.scmocha.out.qs" |> glue::glue()
  )
)

#
#
# variant --------------------------------------------------------------------
#
#

# srr_out_cell_stats -> variant
srr_out_cell_stats |>
  dplyr::mutate(
    total_reads = purrr::map_dbl(
      .x = metrics,
      .f = \(.x) {
        if (is.null(.x)) {
          return(NA_real_)
        }
        .x$`Number of Reads`
      }
    ),
    depth_read_mean = purrr::map_dbl(
      .x = depth_read,
      .f = \(.x) {
        if (is.null(.x)) {
          return(NA_real_)
        }
        mean(.x$depth, na.rm = T)
      }
    ),
    depth_mean = purrr::map_dbl(
      .x = depth,
      .f = \(.x) {
        if (is.null(.x)) {
          return(NA_real_)
        }
        mean(.x$depth, na.rm = T)
      }
    )
  ) |>
  dplyr::mutate(
    nmut = purrr::map_int(
      .x = somatic_variant,
      .f = \(.x) {
        # if (is.null(.x)) {
        #   return(NA_integer_)
        # }
        # nrow(.x)
        length(c(.x$homo, .x$hete))
      }
    ),
    nmut_variant = purrr::map(
      .x = somatic_variant,
      .f = \(.x) {
        .x |>
          purrr::map_int(length) |>
          tibble::enframe() |>
          tidyr::spread(key = name, value = value) -> .xx
        names(.xx) <- glue::glue("nmut_{names(.xx)}")
        .xx
      }
    )
  ) |>
  # tidyr::unnest(cols = nmut_variant)
  dplyr::mutate(
    haplogroup = purrr::map2(
      .x = anno,
      .y = srrid,
      .f = \(.x, .y) {
        log_info(gseid, " ", .y)
        if (is.null(.x)) {
          return(
            tibble::tibble(
              Haplogroup = NA_character_,
              Verbose_haplogroup = NA_character_
            )
          )
        }
        .x |>
          dplyr::select(Haplogroup, Verbose_haplogroup) |>
          dplyr::filter(!is.na(Haplogroup)) |>
          dplyr::filter(Haplogroup != "") |>
          dplyr::distinct() |>
          dplyr::mutate_all(.funs = as.character) -> .xx

        if (nrow(.xx) == 0) {
          tibble::tibble(
            Haplogroup = NA_character_,
            Verbose_haplogroup = NA_character_
          )
        } else {
          .xx
        }
      }
    )
  ) |>
  tidyr::unnest(cols = haplogroup) |>
  # dplyr::inner_join(
  #   pheno, by = "srrid"
  # ) |>
  tidyr::unnest(
    cols = cell_stats
  ) |>
  dplyr::mutate(
    ratio = round(
      `number of cells after filtering` / `estimated number of cells`,
      2
    )
  ) -> metadata_anno


metadata_anno |>
  tidyr::unnest(cols = nmut_variant) |>
  dplyr::select(
    srrid,
    Chemistry = chemistry,
    Haplogroup = Haplogroup,
    `# of variants` = nmut,
    `# of homoplasmic variants` = nmut_homo,
    `# of heteroplasmic variants` = nmut_hete,
    `# of somatic variants` = nmut_somatic,
    `Median UMI/cell` = `median UMI counts per cell`,
    `Median genes/cell` = `median genes per cell`,
    `# of cells` = `estimated number of cells`,
    `# cells after filter` = `number of cells after filtering`,
    # `Cell ratio` = ratio,
    `Total reads` = total_reads,
    `Depth read mean` = depth_read_mean,
    `Depth mean` = depth_mean
  ) -> metadata_clean

metadata_clean |>
  writexl::write_xlsx(
    path = path(
      outdir,
      "{gseid}.cell_ratio_and_variant_clean.xlsx" |> glue::glue()
    )
  )

log_success(
  "save metadata to {outdir}/{gseid}.cell_ratio_and_variant_clean.xlsx"
)

export(
  x = metadata_clean,
  file = path(
    outdir,
    "{gseid}.cell_ratio_and_variant_clean.csv" |> glue::glue()
  )
)

log_success(
  "save metadata to {outdir}/{gseid}.cell_ratio_and_variant_clean.csv"
)

log_success("{gseid} save to {outdir}/{gseid}.scmocha.out.qs" |> glue::glue())

# footer ------------------------------------------------------------------

# save image --------------------------------------------------------------
