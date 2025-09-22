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
outdir <- "/home/liuc9/github/scMOCHA-data/analysis/zzz"


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

conn <- DBI::dbConnect(
  duckdb::duckdb(),
  dbdir = "/mnt/isilon/u01_project/large-scale/liuc9/raw/zzz/clean-data/all_hetero_af.cell.duckdb.1.2.1" |>
    glue::glue(),
  read_only = TRUE
)
tbl_allvariants_cell_fishertest <- dplyr::tbl(
  conn,
  "allvariants_cell_fishertest"
)


tibble::tibble(
  gseid = gseids
) |>
  dplyr::filter(
    !gseid %in% gseids_tobe_excluded
  ) |>
  dplyr::mutate(
    anno = purrr::map(
      .x = gseid,
      .f = function(.gseid) {
        log_info(
          "Loading {.gseid}... ({which(gseids == .gseid)}/{length(gseids)})"
        )
        # .gseid <- "GSE226602"
        .anno <- readr::read_rds(
          file.path(
            basedir,
            .gseid,
            "out",
            glue::glue("{.gseid}.scmocha.out.rds.gz")
          )
        )

        .anno |>
          dplyr::mutate(
            haplo_violin2 = parallel::mcmapply(
              .srrid = srrid,
              .haplo_violin = haplo_violin,
              .somatic_variant = somatic_variant,
              FUN = function(.srrid, .haplo_violin, .somatic_variant) {
                # .srrid <- .anno$srrid[[1]]
                # .haplo_violin <- .anno$haplo_violin[[1]]
                # .somatic_variant <- .anno$somatic_variant[[1]]

                tbl_allvariants_cell_fishertest |>
                  dplyr::filter(srrid == .srrid) |>
                  dplyr::select(
                    barcode,
                    variant,
                    celltype,
                    depth,
                    AFO,
                    ARE,
                    CFO,
                    CRE,
                    GFO,
                    GRE,
                    TFO,
                    TRE,
                    variant_type_fisher_test
                  ) |>
                  dplyr::collect() -> .tt

                .haplo_violin |>
                  dplyr::left_join(
                    .tt |>
                      dplyr::select(barcode, variant, variant_type_fisher_test),
                    by = c("barcode", "variant")
                  ) -> .hv

                .tt |>
                  dplyr::select(
                    -c(
                      barcode,
                      variant_type_fisher_test
                    )
                  ) |>
                  dplyr::mutate(
                    ref = gsub("\\d+|>.*", "", variant),
                    alt = gsub("\\d+.*>", "", variant)
                  ) -> .tt_

                .tt_ |>
                  tidyr::nest(
                    .key = "reads",
                    .by = c(variant, celltype, ref, alt)
                  ) |>
                  dplyr::mutate(
                    clusteraf = purrr::pmap_dbl(
                      .l = list(reads, ref, alt),
                      .f = function(reads, ref, alt) {
                        # reads <- .tt_$reads[[1]]
                        # ref <- .tt_$ref[[1]]
                        # alt <- .tt_$alt[[1]]

                        .total_reads <- sum(reads$depth, na.rm = TRUE)
                        .alt_reads <- sum(
                          c(
                            reads[[paste0(alt, "FO")]],
                            reads[[paste0(alt, "RE")]]
                          ),
                          na.rm = TRUE
                        )
                        .ref_reads <- sum(
                          c(
                            reads[[paste0(ref, "FO")]],
                            reads[[paste0(ref, "RE")]]
                          ),
                          na.rm = TRUE
                        )
                        if (.total_reads == 0) {
                          return(NA_real_)
                        } else {
                          return(.alt_reads / .total_reads)
                        }
                      }
                    )
                  ) |>
                  dplyr::select(variant, celltype, clusteraf) -> .hetero

                .tt_ |>
                  dplyr::select(-celltype) |>
                  tidyr::nest(
                    .key = "reads",
                    .by = c(variant, ref, alt)
                  ) |>
                  dplyr::mutate(
                    bulkaf = purrr::pmap_dbl(
                      .l = list(reads, ref, alt),
                      .f = function(reads, ref, alt) {
                        # reads <- a$reads[[1]]
                        # ref <- a$ref[[1]]
                        # alt <- a$alt[[1]]
                        .total_reads <- sum(reads$depth, na.rm = TRUE)
                        .alt_reads <- sum(
                          c(
                            reads[[paste0(alt, "FO")]],
                            reads[[paste0(alt, "RE")]]
                          ),
                          na.rm = TRUE
                        )
                        .ref_reads <- sum(
                          c(
                            reads[[paste0(ref, "FO")]],
                            reads[[paste0(ref, "RE")]]
                          ),
                          na.rm = TRUE
                        )
                        if (.total_reads == 0) {
                          return(NA_real_)
                        } else {
                          return(.alt_reads / .total_reads)
                        }
                      }
                    )
                  ) |>
                  dplyr::select(variant, bulkaf) -> .bulkaf

                .hv |>
                  dplyr::filter(variant_type_fisher_test == "colorful") |>
                  dplyr::count(variant) |>
                  dplyr::filter(n >= 3) |>
                  dplyr::filter(
                    variant %in%
                      unique(c(
                        .somatic_variant$somatic,
                        .somatic_variant$haplo
                      ))
                  ) |>
                  dplyr::pull(variant) -> .real_variants

                .somatic_variant$strand_bias <- setdiff(
                  unique(c(
                    .somatic_variant$somatic,
                    .somatic_variant$haplo
                  )),
                  .real_variants
                )
                .somatic_variant$somatic <- intersect(
                  .somatic_variant$somatic,
                  .real_variants
                )
                .somatic_variant$haplo <- intersect(
                  .somatic_variant$haplo,
                  .real_variants
                )

                tibble::tibble(
                  haplo_violin_fisher = list(.hv),
                  clusteraf = list(.hetero),
                  bulkaf = list(.bulkaf),
                  somatic_variant_fisher = list(.somatic_variant)
                )
              },
              mc.cores = 30,
              SIMPLIFY = FALSE
            )
          ) |>
          tidyr::unnest(cols = haplo_violin2) -> .anno_new

        log_success(
          "Loaded {.gseid}! ({which(gseids == .gseid)}/{length(gseids)})"
        )
        return(.anno_new)
      }
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
  export(
    gse_data,
    file.path(outdir, "gse_data_fisher.rds")
  )
  export(
    gse_data,
    file.path(outdir, "gse_data_fisher.qs")
  )
  # gse_data |>
  #   dplyr::select(1, 2, 3) |>
  #   export(
  #     file.path(outdir, "gse_srrid_srrdir.csv"),
  #     format = "both",
  #     sep = ",",
  #   )
  # gse_data |>
  #   dplyr::select(1, 2, 3) |>
  #   export(
  #     file.path(outdir, "gse_srrid_srrdir.rds")
  #   )
}
DBI::dbDisconnect(conn, shutdown = TRUE)
# python / home / liuc9 / github / scMOCHA - data / stats / stats / barcode_celltype.py
