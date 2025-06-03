# ! read --------------------------------------------------------------------


tictoc::tic()
d <- readr::read_rds("/scr1/users/liuc9/tmp/gse_data.rds")
tictoc::toc()

tictoc::tic()
d <- qs::qread("/scr1/users/liuc9/tmp/gse_data.qs")
tictoc::toc()

bench::mark(
  read_rds = readr::read_rds("/scr1/users/liuc9/tmp/gse_data.rds"),
  qread = qs::qread("/scr1/users/liuc9/tmp/gse_data.qs")
)


# ! write --------------------------------------------------------------------



# tictoc::tic()
# readr::write_rds(d, "/scr1/users/liuc9/tmp/gse_data.rds")
# tictoc::toc()

# tictoc::tic()
# qs::qsave(d, "/scr1/users/liuc9/tmp/gse_data.qs")
# tictoc::toc()

tictoc::tic()
d <- readr::read_rds("/scr1/users/liuc9/tmp/gse_data.rds")
tictoc::toc()

tictoc::tic()
m <- qs::qread("/scr1/users/liuc9/tmp/gse_data.qs")
tictoc::toc()

bench::mark(
  read_rds = readr::write_rds(d, "/scr1/users/liuc9/tmp/gse_data.rds"),
  qsave = qs::qsave(m, "/scr1/users/liuc9/tmp/gse_data.qs")
)

export(m, "/scr1/users/liuc9/tmp/gse_data.rds", format = "csv")

a <- import("/scr1/users/liuc9/tmp/gse_data.rds")

tictoc::tic()
a <- import("/scr1/users/liuc9/tmp/gse_data.qs")
tictoc::toc()

export(a, "/scr1/users/liuc9/tmp/gse_data.rds")
import("/home/liuc9/github/dotfiles/Renv/renv.yaml")


convert(
  "/home/liuc9/github/scMOCHA-data/stats/stats/zzz/clean-data/gse_srrid_srrdir.rds",
  "/home/liuc9/github/scMOCHA-data/stats/stats/zzz/clean-data/gse_srrid_srrdir.qs"
)

clean_data_dir <- "/home/liuc9/github/scMOCHA-data/stats/stats/zzz/clean-data"
rdss <- list.files(
  clean_data_dir,
  pattern = ".rds$",
)


purrr::map(
  rdss,
  ~ convert(
    file.path(clean_data_dir, .x),
    file.path(clean_data_dir, gsub(".rds$", ".qs", .x))
  )
)

clean_data_dir <- "/home/liuc9/github/scMOCHA-data/stats/stats/zzz/clean-data"
csvs <- list.files(
  clean_data_dir,
  pattern = ".csv$",
)


purrr::map(
  csvs,
  ~ convert(
    file.path(clean_data_dir, .x),
    file.path(clean_data_dir, gsub(".csv$", ".fst", .x))
  )
)


convert(
  "/home/liuc9/github/scMOCHA-data/stats/stats/zzz/clean-data/barcode_celltype.fst",
  "/home/liuc9/github/scMOCHA-data/stats/stats/zzz/clean-data/barcode_celltype.feather"
)

{
  convert(
    "/mnt/isilon/u01_project/large-scale/liuc9/raw/zzz/db/EXPR/gse_srrid_celltype_gene_expr.csv",
    "/mnt/isilon/u01_project/large-scale/liuc9/raw/zzz/db/EXPR/gse_srrid_celltype_gene_expr.fst"
  )


  d <- import("/mnt/isilon/u01_project/large-scale/liuc9/raw/zzz/db/EXPR/gse_srrid_celltype_gene_expr.fst")


  d |>
    tidyr::pivot_longer(
      cols = -c(gseid, srrid, genename),
      names_to = "celltype",
      values_to = "expr"
    ) ->
  dd

  dd |>
    dplyr::group_by(genename, celltype) |>
    tidyr::nest(.key = "expr") |>
    dplyr::ungroup() ->
  ddd

  export(
    ddd,
    "/mnt/isilon/u01_project/large-scale/liuc9/raw/zzz/db/EXPR/gse_srrid_celltype_gene_expr.qs"
  )
}



convert(
  "/home/liuc9/github/scMOCHA-data/config/mtdna_genes_dloop.csv.gz",
  "/home/liuc9/github/scMOCHA-data/config/mtdna_genes_dloop.fst"
)

# feather --------------------------------------------------------------------
a <- import("/mnt/isilon/u01_project/large-scale/liuc9/raw/zzz/clean-data/barcode_celltype.fst")
a$barcode |>
  unique() |>
  length()
nrow(a)




# ! test --------------------------------------------------------------------

old <- import("/scr1/users/liuc9/tmp/scanpy/GSE155673_GSM4712885_celltype_gene_expr.csv")

new <- import("/home/liuc9/github/scMOCHA-data/stats/stats/zzz/db/EXPR/GSE155673_GSM4712885_celltype_gene_expr.csv")

old |>
  dplyr::select(genename, bold = B) |>
  dplyr::left_join(
    new |>
      dplyr::select(genename, newb = B),
    by = "genename"
  )



convert(
  "/mnt/isilon/xing_lab/liuc9/refdata/ensembl/Homo_sapiens.GRCh38.107.gtf.id_name_length_genetype.rds",
  "/mnt/isilon/xing_lab/liuc9/refdata/ensembl/Homo_sapiens.GRCh38.107.gtf.id_name_length_genetype.fst"
)


# ? /home/liuc9/github/scMOCHA-data/config/Homo_sapiens.GRCh38.107.gtf.id_name_length_genetype.fst --------------------------------------------------------------------

d <- import("/home/liuc9/github/scMOCHA-data/config/mtdna_genes_dloop.fst")

d
