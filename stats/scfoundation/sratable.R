#!/usr/bin/env Rscript
# Metainfo ----------------------------------------------------------------
# @AUTHOR: Chun-Jie Liu
# @CONTACT: chunjie.sam.liu.at.gmail.com
# @DATE: 2025-01-02 12:49:40
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

# src ---------------------------------------------------------------------

# header ------------------------------------------------------------------
log_threshold(TRACE)
log_layout(layout_glue_colors)

# future: :plan(future: :multisession, workers = 10)

# function ----------------------------------------------------------------
fn_celltype_pie_plot <- function(.xxx_celltype) {
  # pcc <- readr::read_tsv(file = "https://raw.githubusercontent.com/chunjie-sam-liu/chunjie-sam-liu.life/master/public/data/pcc.tsv")

  .xxx_celltype |>
    dplyr::ungroup() %>%
    dplyr::mutate(csum = rev(cumsum(rev(n)))) %>%
    dplyr::mutate(pos = n / 2 + dplyr::lead(csum, 1)) %>%
    dplyr::mutate(pos = dplyr::if_else(is.na(pos), n / 2, pos)) %>%
    dplyr::mutate(percentage = n / sum(n)) %>%
    dplyr::arrange(-n) |>
    dplyr::mutate(
      celltype = as.factor(celltype),
      levels = celltype
    ) |>
    ggplot(aes(x = "", y = n, fill = celltype)) +
    geom_bar(
      stat = "identity",
      width = 1,
      color = "white",
      show.legend = FALSE
    )
  # scale_fill_brewer(palette = "Dark2", name = NULL) +
  # scale_fill_manual(
  #   name = NULL,
  #   values = pcc$color,
  # ) +
  # scale_color_manual(
  #   name = NULL,
  #   values = pcc$color
  # ) +
  ggrepel::geom_label_repel(
    aes(
      y = pos,
      # label = glue::glue("{celltype}\n{n} ({scales::percent(percentage)})"),
      fill = celltype,
      # color = celltype
    ),
    size = 6,
    # fill = "white",
    nudge_x = 1,
    show.legend = FALSE,
  ) +
    coord_polar(theta = "y", start = 0) +
    theme_void() +
    theme(
      plot.title = element_text(
        # vjust = -2,
        hjust = 0.5,
        size = 22,
      ),
      legend.position = "none"
    )
}


# load data ---------------------------------------------------------------

basedir <- "/home/liuc9/github/scMOCHA-data/data/scfoundation"

readr::read_lines("/home/liuc9/github/scMOCHA-data/data/scfoundation/cmds.sh") |> gsub(
  pattern = "Rscript /home/liuc9/github/scMOCHA-data/src/01-sra-metadata.R -g | -b /home/liuc9/github/scMOCHA-data/data/scfoundation",
  replacement = "",
  x = _
) -> gseids

# body --------------------------------------------------------------------
sqlite_file <- "/mnt/isilon/xing_lab/liuc9/refdata/sradb/SRAmetadb.sqlite"

sra_con <- DBI::dbConnect(
  RSQLite::SQLite(),
  sqlite_file
)
DBI::dbListTables(sra_con)

sra_table <- dplyr::tbl(
  sra_con, "sra"
)
study_table <- dplyr::tbl(
  sra_con, "study"
)


sra_table |>
  dplyr::filter(study_name %in% gseids) |>
  dplyr::inner_join(study_table, by = "study_accession") |>
  as.data.table() ->
sra_df

sra_df |>
  dplyr::pull(sample_accession) |>
  unique() ->
sample_accessions

sample_table <- dplyr::tbl(
  sra_con, "sample"
)

sample_df <- sample_table |>
  dplyr::filter(
    sample_accession %in% sample_accessions
  ) |>
  as.data.table()

cleaned_sample_df <- sample_df %>%
  dplyr::select(-which(apply(is.na(.), 2, all))) |>
  dplyr::select(sample_accession, sample_attribute)

cleaned_sample_df$sample_attribute[[1]] |>
  stringr::str_split(pattern = "\\|\\|", simplify = T) |>
  stringr::str_split(pattern = ": ") |>
  purrr::map(.f = \(.x) {
    .x[[1]]
  }) |>
  purrr::reduce(.f = c) ->
new_columns

cleaned_sample_df |>
  tidyr::separate(
    col = sample_attribute,
    into = new_columns,
    sep = " \\|\\|"
  ) %>%
  dplyr::mutate_at(
    new_columns,
    ~ stringr::str_remove(., ".*:")
  ) ->
cleaned_sample_df
DBI::dbDisconnect(sra_con)

sra_df |>
  dplyr::inner_join(
    cleaned_sample_df,
    by = "sample_accession"
  ) |>
  dplyr::relocate(
    run_accession,
    .before = 1
  ) |>
  dplyr::mutate(
    source_name = stringr::str_trim(source_name)
  ) ->
cleaned_sample_df_sra


data.table::fwrite(
  x = cleaned_sample_df_sra,
  file = file.path(
    basedir,
    "1.all.metadata.csv" |> glue::glue()
  )
)

cleaned_sample_df_sra |>
  dplyr::filter(library_strategy == "RNA-Seq") ->
cleaned_sample_df_sra_rnaseq

pbmc <- c("PBMCs", "PBMC", "Cryopreserved PBMC", "pbmc", "engineered PBMCs")

cleaned_sample_df_sra_rnaseq |>
  dplyr::filter(!is.na(source_name)) |>
  dplyr::select(experiment_name, source_name) |>
  dplyr::distinct() |>
  # dplyr::filter(grepl("pbmc", source_name, ignore.case = T)) |>
  dplyr::mutate(
    source_name = ifelse(
      source_name %in% pbmc,
      "PBMC",
      source_name
    )
  ) |>
  dplyr::count(source_name) |>
  dplyr::mutate(
    source_name = ifelse(n < 40, "others", source_name)
  ) |>
  dplyr::group_by(source_name) |>
  dplyr::summarise(n = sum(n)) |>
  dplyr::rename(celltype = source_name) |>
  dplyr::arrange(-n) ->
celltype_df
celltype_df |>
  data.table::fwrite(
    file = file.path(
      basedir,
      "1.all.metadata.celltypes.csv" |> glue::glue()
    )
  )


cleaned_sample_df_sra_rnaseq |>
  dplyr::filter(!is.na(source_name)) |>
  dplyr::filter(source_name %in% pbmc) |>
  dplyr::glimpse()

cleaned_sample_df_sra_rnaseq |>
  dplyr::filter(!is.na(source_name)) |>
  dplyr::filter(source_name %in% pbmc) |>
  dplyr::mutate(
    Run = run_accession,
    `Sample Name` = experiment_name,
    Experiment = experiment_accession,
    study_name = study_name
  ) |>
  dplyr::relocate(
    Run,
    `Sample Name`,
    Experiment,
    study_name,
    source_name,
    .before = 1
  ) ->
sratable_pbmc

already_processed <- c(
  "GSE149689", "GSE155223", "GSE155673", "GSE157344", "GSE163668", "GSE166992", "GSE171555", "GSE181279", "GSE226602", "GSE161354",
  "GSE175524",
  "GSE206283",
  "GSE226598",
  "GSE235050",
  "GSE261140",
  "GSE279945"
)

sratable_pbmc |>
  data.table::fwrite(
    file = file.path(
      basedir,
      "1.all.metadata.celltypes.pbmc.csv" |> glue::glue()
    )
  )


sratable_pbmc |>
  dplyr::filter(!study_name %in% already_processed) |>
  dplyr::group_by(study_name) |>
  dplyr::pull(study_name) |>
  unique() ->
new_gseids_pbmc


new_gseids_pbmc <- c("GSE140881", "GSE142595", "GSE149313", "GSE154386", "GSE159117", "GSE162117", "GSE167825", "GSE179566", "GSE188632", "GSE192391")

readr::write_lines(
  x = new_gseids_pbmc,
  file = file.path(
    basedir,
    "1.new_gseids_pbmc.txt" |> glue::glue()
  )
)

sratable_pbmc |>
  dplyr::filter(study_name %in% new_gseids_pbmc) |>
  dplyr::mutate(gseid = study_name) |>
  dplyr::group_by(gseid) |>
  tidyr::nest() |>
  dplyr::ungroup() |>
  dplyr::mutate(
    a = purrr::map2(
      .x = data,
      .y = gseid,
      .f = \(.x, .y) {
        data.table::fwrite(
          x = .x,
          file = file.path(
            basedir,
            .y,
            glue::glue("{.y}.SraRunTable")
          )
        )
      }
    )
  )





# run_accession,experiment_name,experiment_accession

#
# footer ------------------------------------------------------------------

# future: :plan(future: :sequential)

# save image --------------------------------------------------------------
