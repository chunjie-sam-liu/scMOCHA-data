#!/usr/bin/env Rscript
# Metainfo ----------------------------------------------------------------
# @AUTHOR: Chun-Jie Liu
# @CONTACT: chunjie.sam.liu.at.gmail.com
# @DATE: 2025-11-30 00:40:39
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

# header ------------------------------------------------------------------

# load data ---------------------------------------------------------------

basedir <- path("/home/liuc9/github/scMOCHA-data/data")
foundation_out <- path(basedir, "scfoundation/out")
outdir <- path("/home/liuc9/github/scMOCHA-data/analysis/zzz")


#
#
# Publications --------------------------------------------------------------------
#
#

# Create a tibble with gseid and Publication columns

PUBLICATIONS = tibble::tribble(
  ~gseid      , ~Publication                                             ,
  "GSE181279" , NA                                                       ,
  "GSE235050" , "Nat Immunol 2024"                                       ,
  "GSE161354" , "Front Med (Lausanne) 2021"                              ,
  "GSE166992" , NA                                                       ,
  "GSE226602" , NA                                                       ,
  "GSE147794" , "PLoS Pathog 2020"                                       ,
  "GSE167825" , "Proc Natl Acad Sci U S A 2022"                          ,
  "GSE163668" , NA                                                       ,
  "GSE154386" , "PLoS Pathog 2021"                                       ,
  "GSE155223" , NA                                                       ,
  "GSE233844" , "Am J Respir Crit Care Med 2024"                         ,
  "GSE175524" , "Aging (Albany NY) 2023"                                 ,
  "GSE171555" , NA                                                       ,
  "GSE206283" , "Immun Ageing 2025"                                      ,
  "GSE184703" , "Front Immunol 2021"                                     ,
  "GSE188632" , "Front Cell Dev Biol 2021"                               ,
  "GSE153421" , "Protein Cell 2022"                                      ,
  "GSE175499" , "Cell Rep 2021"                                          ,
  "GSE261140" , "Sci Immunol 2024"                                       ,
  "GSE214865" , "Eur Respir J 2023"                                      ,
  "GSE149689" , NA                                                       ,
  "GSE220189" , "Nat Comput Sci 2023"                                    ,
  "GSE157344" , NA                                                       ,
  "GSE226598" , "Cell Rep Med 2023"                                      ,
  "GSE155673" , NA                                                       ,
  "GSE279945" , "Advances in Neural Information Processing Systems 2024" ,
  "GSE159117" , "J Exp Med 2022"                                         ,
  "GSE174125" , "Cell Rep 2022"                                          ,
  "GSE162117" , "Blood cancer discovery 2021"                            ,
  "GSE163314" , "Arthritis Res Ther 2021"                                ,
  "GSE149313" , "mBio 2021"
)

SEXPRED = import(
  "/home/liuc9/github/scMOCHA-data/analysis/zzz/clean-data/gse_srrid_srrdir_sex.qs"
) |>
  dplyr::select(
    srrid,
    SEXPRED = sex
  )

# load conn ---------------------------------------------------------------

# src ---------------------------------------------------------------------

# function ----------------------------------------------------------------

# body --------------------------------------------------------------------

#
#
# GSE samples excluding --------------------------------------------------------------------
#
#

gseids <- c(
  "GSE155673",
  "GSE157344",
  "GSE149689",
  "GSE171555",
  "GSE155223",
  "GSE163668",
  "GSE175524",
  "GSE206283",
  "GSE226598",
  "GSE261140",
  "GSE279945",
  "GSE214865",
  "GSE220189",
  "GSE233844",
  "GSE175499",
  "GSE149313",
  "GSE154386",
  "GSE159117",
  "GSE188632",
  "GSE166992",
  "GSE162117",
  "GSE226602",
  "GSE161354",
  "GSE235050",
  "GSE181279",
  # scfoundation2
  "GSE143353",
  "GSE148215",
  "GSE163314",
  "GSE163633",
  "GSE164690",
  "GSE167825",
  "GSE174125",
  "GSE184703",
  "GSE153421",
  # "GSE168453",
  "GSE147794"
)

# not PBMC
gseids_not_pbmc <- c(
  "GSE168453", # Pool samples together, not individual samples
  # "GSE163668", # Some are pooled samples, not individual samples. Whole blood, not PBMC, including others like Red blood cells, Platelets and Plasma
  "GSE163633", # Not all PBMC, some are Mucosa-derived T cells / Myloid cells, Squamous Cell Carcinoma(SCC) -derived T cells / Myloid Cells
  # "GSE157344", # It’s not PBMC, but the blood or bronchoalveloar lavage samples
  # "GSE163314", # Half of the samples are from Colon, not PBMC
  "GSE164690", # Most of the samples are from tumors (Head and neck squamous cell carcinoma (HNSCC)), not from PBMC
  "GSE148215" # The samples are human embryonic cells, not PBMC cells
)

# cell line
gseids_cellline <- c(
  "GSE143353"
)

# enrich cells
gseids_enrich_cells <- c(
  "GSE167825", # CD8T cell enriched
  "GSE175524", # B cell enriched
  "GSE261140" # CD8T cell enriched
)
# GSE163314 only keep the PBMC samples
GSE143353_blood <- c(
  "GSM4976993", # Patient 2 Blood
  "GSM4976995", # Patient 3 Blood
  "GSM4976997", # Patient 5 Blood
  "GSM4976999", # Patient 7 Blood
  "GSM4977001", # Patient 21 Blood
  "GSM4977003", # Patient 23 Blood
  "GSM4977005", # Patient 27 Blood
  "GSM4977007" # Patient 33 Blood
)
GSE143353_colon <- c(
  "GSM4976992", # Patient 2 Colon
  "GSM4976994", # Patient 3 Colon
  "GSM4976996", # Patient 5 Colon
  "GSM4976998", # Patient 7 Colon
  "GSM4977000", # Patient 21 Colon
  "GSM4977002", # Patient 23 Colon
  "GSM4977004", # Patient 27 Colon
  "GSM4977006" # Patient 33 Colon
)
GSE163668_blood <- c(
  "GSM4995431",
  "GSM4995432",
  "GSM4995433",
  "GSM4995434",
  "GSM4995435",
  "GSM4995436",
  "GSM4995437",
  "GSM4995438",
  "GSM4995439",
  "GSM4995440",
  "GSM4995441",
  "GSM4995442",
  "GSM4995443",
  "GSM4995444",
  "GSM4995445",
  "GSM4995446",
  "GSM4995447",
  "GSM4995448",
  "GSM4995449",
  "GSM4995450",
  "GSM4995451",
  "GSM4995452",
  "GSM4995453",
  "GSM4995454",
  "GSM4995455",
  "GSM4995456",
  "GSM4995457",
  "GSM4995458",
  "GSM4995459",
  "GSM4995460",
  "GSM4995461",
  "GSM4995462"
)
GSE163668_pooled <- c(
  "GSM4995425", # Pooled 10X GEX libraries for Patients 1, 2, and 3
  "GSM4995426", # Pooled 10X GEX libraries for Patients 5 and 6
  "GSM4995427", # Pooled 10X GEX libraries for Patients 7, 8, and 9
  "GSM4995428", # Pooled 10X GEX libraries for Patient 10
  "GSM4995429", # Pooled 10X GEX libraries for Patients 17, 20, and 21
  "GSM4995430" # Pooled 10X GEX libraries for Patients 50 and 51
)
# update new pooled samples to be excluded
GSE162117_pooled <- c(
  "GSM4933450",
  "GSM4933442",
  "GSM4933449",
  "GSM4933446",
  "GSM4933445",
  "GSM4933448"
)

gseids_tobe_excluded <- c(
  gseids_not_pbmc,
  gseids_cellline,
  "GSE162117" # Pooled samples
)
gsmids_tobe_excluded <- c(
  GSE143353_colon,
  GSE163668_pooled
)


#
#
# metadata --------------------------------------------------------------------
#
#

tibble::tibble(gseid = gseids) |>
  dplyr::filter(
    !gseid %in% gseids_tobe_excluded
  ) -> clean_gseids_df

clean_gseids_df$gseid |> unique() -> clean_gseids


clean_gseids_df |>
  dplyr::mutate(
    anno = purrr::map(
      .x = gseid,
      .f = function(.gseid) {
        out_file <- path(
          basedir,
          .gseid,
          "out",
          glue::glue("{.gseid}.cell_ratio_and_variant_clean.csv")
        )
        cli_alert_info(
          "{out_file}"
        )
        .anno <- import(out_file)

        return(.anno)
      }
    )
  ) |>
  tidyr::unnest(cols = anno) |>
  dplyr::filter(
    !srrid %in% gsmids_tobe_excluded
  ) -> cell_ratio_and_variant_clean


gse_dataset_metadata_full <- import(
  path(foundation_out, "gse_dataset_metadata_full.rds")
) |>
  dplyr::filter(
    !gseid %in% gseids_tobe_excluded
  ) |>
  dplyr::filter(
    !srrid %in% gsmids_tobe_excluded
  ) |>
  dplyr::mutate(
    disease = dplyr::case_when(
      disease %in%
        c("Alzheimer's Disease", "Healthy", "COVID-19", "Unknown") ~ disease,
      TRUE ~ "Other"
    )
  ) |>
  dplyr::mutate(
    disease = factor(
      disease,
      levels = c(
        "Alzheimer's Disease",
        "COVID-19",
        "Healthy",
        "Unknown",
        "Other"
      )
    )
  ) |>
  dplyr::mutate(
    Chemistry = factor(
      Chemistry,
      levels = c("SC3Pv2", "SC3Pv3", "SC5P-R2", "SC5P-PE") |> rev()
    )
  ) |>
  dplyr::select(-Publication) |>
  dplyr::left_join(
    PUBLICATIONS,
    by = "gseid"
  ) |>
  dplyr::left_join(
    SEXPRED,
    by = "srrid"
  ) |>
  as.data.table()

gse_dataset_metadata_full |>
  dplyr::select(
    -dplyr::all_of(
      setdiff(
        intersect(
          colnames(gse_dataset_metadata_full),
          colnames(cell_ratio_and_variant_clean)
        ),
        c("gseid", "srrid")
      )
    )
  ) |>
  dplyr::left_join(
    cell_ratio_and_variant_clean,
    by = c("gseid" = "gseid", "srrid" = "srrid")
  ) -> gse_dataset_metadata_full


# save gse_dataset_metadata_full ----------------------------------------------------------
{
  outdir_clean <- "/home/liuc9/github/scMOCHA-data/data/zzz/clean-data"
  export(
    gse_dataset_metadata_full,
    path(outdir_clean, "gse_dataset_metadata_full.csv"),
    format = "both"
  )
  export(
    gse_dataset_metadata_full,
    path(outdir_clean, "gse_dataset_metadata_full.qs"),
  )

  export(
    gse_dataset_metadata_full,
    path(outdir_clean, "gse_dataset_metadata_full.rds"),
  )
}


#
#
# GSE anno data --------------------------------------------------------------------
#
#

gse_dataset_metadata_full |>
  dplyr::select(gseid) |>
  dplyr::distinct() -> clean_gseids_df

clean_gseids_df$gseid |> unique() -> clean_gseids

clean_gseids_df |>
  dplyr::mutate(
    anno = parallel::mclapply(
      X = gseid,
      FUN = function(.gseid) {
        cli_alert_info(
          "Loading {.gseid}... ({which(clean_gseids == .gseid)}/{length(clean_gseids)})"
        )
        .anno <- import(
          path(
            basedir,
            .gseid,
            "out",
            glue::glue("{.gseid}.scmocha.out.qs")
          )
        )
        cli_alert_success(
          "Loaded {.gseid}! ({which(clean_gseids == .gseid)}/{length(clean_gseids)})"
        )
        return(.anno)
      },
      mc.cores = 10
    )
  ) -> gse_data_loaded


gse_data_loaded |>
  tidyr::unnest(cols = anno) |>
  dplyr::filter(
    !srrid %in% gsmids_tobe_excluded
  ) |>
  dplyr::mutate(
    chemistry = factor(
      chemistry,
      levels = c("SC3Pv2", "SC3Pv3", "SC5P-R2", "SC5P-PE") |> rev()
    )
  ) -> gse_data


# save gse_data
{
  outdir <- "/home/liuc9/github/scMOCHA-data/data/zzz/clean-data"
  # data.table::fwrite(
  #   gse_data,
  #   file.path(outdir, "gse_data.csv"),
  #   sep = ",",
  # )
  cli_alert_info("Exporting gse_data.rds very slow")
  cli_alert_info("Only export once.")
  # export(
  #   gse_data,
  #   path(outdir, "gse_data.rds")
  # )

  cli_h1("Save gse_data.qs")
  cli_h2("Save raw gse_data.raw.qs which include all columns")
  export(
    gse_data,
    path(outdir, "gse_data.raw.qs"),
    preset = "fast"
  )
  export(
    gse_data |> dplyr::select(-c(cellaf)),
    path(outdir, "gse_data.qs"),
    preset = "fast"
  )
  # export(
  #   gse_data |> dplyr::select(gseid, srrid, cellaf, clusteraf, bulkaf),
  #   path(outdir, "gse_data.af.qs"),
  #   preset = "fast"
  # )
  export(
    gse_data |> dplyr::select(gseid, srrid, cellaf),
    path(outdir, "gse_data.cellaf.qs"),
    preset = "fast"
  )
  export(
    gse_data |> dplyr::select(gseid, srrid, clusteraf),
    path(outdir, "gse_data.clusteraf.qs"),
    preset = "fast"
  )
  export(
    gse_data |> dplyr::select(gseid, srrid, bulkaf),
    path(outdir, "gse_data.bulkaf.qs"),
    preset = "fast"
  )

  cli_rule("Save gse_srrid_srrdir.csv and .rds")

  gse_data |>
    dplyr::select(1, 2, 3) |>
    export(
      path(outdir, "gse_srrid_srrdir.csv"),
      format = "both",
    )
  gse_data |>
    dplyr::select(1, 2, 3) |>
    export(
      path(outdir, "gse_srrid_srrdir.rds")
    )
}

# footer ------------------------------------------------------------------

# save image --------------------------------------------------------------
