# ting <- readxl::read_excel("/mnt/isilon/u01_project/large-scale/liuc9/raw/zzz/PT_samples/all_samples_PT_used.xlsx")
# cj <- data.table::fread(
#   file.path(
#     "/mnt/isilon/u01_project/large-scale/liuc9/raw/zzz/used_samples",
#     "gse_srrid_srrdir.cj.csv"
#   )
# )

# colnames(ting)
# ggvenn::ggvenn(
#   data = list(
#     "PT" = unique(ting$Samples),
#     "CJ" = unique(cj$srrid)
#   ),
#   fill_color = c("#FF9999", "#66B3FF"),
#   stroke_color = "white"
# )

# ggvenn::ggvenn(
#   data = list(
#     "PT" = unique(ting$Datasets),
#     "CJ" = unique(cj$gseid)
#   ),
#   fill_color = c("#FF9999", "#66B3FF"),
#   stroke_color = "white"
# )

# setdiff(
#   unique(ting$Datasets),
#   unique(cj$gseid)
# ) ->
# ting_only
# setdiff(
#   unique(cj$gseid),
#   unique(ting$Datasets)
# ) ->
# cj_only



# ! # clean data --------------------------------------------------------------------


basedir <- "/home/liuc9/github/scMOCHA-data/data"
foundation_out <- file.path(basedir, "scfoundation/out")
outdir <- "/home/liuc9/github/scMOCHA-data/stats/stats/zzz"


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


gseids_tobe_excluded <- c(
  gseids_not_pbmc,
  gseids_cellline
)
gsmids_tobe_excluded <- c(
  GSE143353_colon,
  GSE163668_pooled
)


gse_dataset_metadata_full <- readr::read_rds(
  file.path(foundation_out, "gse_dataset_metadata_full.rds")
) |>
  dplyr::filter(
    !gseid %in% gseids_tobe_excluded
  ) |>
  dplyr::filter(
    !srrid %in% gsmids_tobe_excluded
  ) |>
  dplyr::mutate(
    disease = dplyr::case_when(
      disease %in% c("Alzheimer's Disease", "Healthy", "COVID-19", "Unknown") ~ disease,
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
      levels = c("SC3Pv2", "SC3Pv3", "SC5P-R2", "SC5P-PE")
    )
  )

# save gse_dataset_metadata_full
{
  outdir <- "/home/liuc9/github/scMOCHA-data/data/zzz/clean-data"
  data.table::fwrite(
    gse_dataset_metadata_full,
    file.path(outdir, "gse_dataset_metadata_full.csv"),
    sep = ",",
  )
  readr::write_rds(
    gse_dataset_metadata_full,
    file.path(outdir, "gse_dataset_metadata_full.rds")
  )
}



tibble::tibble(
  gseid = gseids
) |>
  dplyr::filter(
    !gseid %in% gseids_tobe_excluded
  ) |>
  dplyr::mutate(
    anno = parallel::mclapply(
      X = gseid,
      FUN = function(.gseid) {
        log_info("Loading {.gseid}... ({which(gseids == .gseid)}/{length(gseids)})")
        .anno <- readr::read_rds(
          file.path(basedir, .gseid, "out", glue::glue("{.gseid}.scmocha.out.rds.gz"))
        )
        log_success("Loaded {.gseid}! ({which(gseids == .gseid)}/{length(gseids)})")
        return(.anno)
      },
      mc.cores = 10
    )
  ) ->
gse_data_loaded

gse_data_loaded |>
  tidyr::unnest(cols = anno) |>
  dplyr::filter(
    !srrid %in% gsmids_tobe_excluded
  ) ->
gse_data

# save gse_data
{
  outdir <- "/home/liuc9/github/scMOCHA-data/data/zzz/clean-data"
  # data.table::fwrite(
  #   gse_data,
  #   file.path(outdir, "gse_data.csv"),
  #   sep = ",",
  # )
  readr::write_rds(
    gse_data,
    file.path(outdir, "gse_data.rds")
  )
  gse_data |>
    dplyr::select(1, 2, 3) |>
    data.table::fwrite(
      file.path(outdir, "gse_srrid_srrdir.csv"),
      sep = ",",
    )
  gse_data |>
    dplyr::select(1, 2, 3) |>
    readr::write_rds(
      file.path(outdir, "gse_srrid_srrdir.rds")
    )
}


# python / home / liuc9 / github / scMOCHA - data / stats / stats / barcode_celltype.py
