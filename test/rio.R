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
  "/home/liuc9/github/scMOCHA-data/analysis/zzz/clean-data/gse_srrid_srrdir.rds",
  "/home/liuc9/github/scMOCHA-data/analysis/zzz/clean-data/gse_srrid_srrdir.qs"
)

clean_data_dir <- "/home/liuc9/github/scMOCHA-data/analysis/zzz/clean-data"
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

clean_data_dir <- "/home/liuc9/github/scMOCHA-data/analysis/zzz/clean-data"
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
  "/home/liuc9/github/scMOCHA-data/analysis/zzz/clean-data/barcode_celltype.fst",
  "/home/liuc9/github/scMOCHA-data/analysis/zzz/clean-data/barcode_celltype.feather"
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

new <- import("/home/liuc9/github/scMOCHA-data/analysis/zzz/db/EXPR/GSE155673_GSM4712885_celltype_gene_expr.csv")

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

phastCons100way <- data.table::fread(
  "/home/liuc9/github/scMOCHA-data/config/chrM.phastCons100way.wigFix"
) |>
  tibble::rowid_to_column() |>
  tibble::add_column(
    seqnames = "MT",
    .before = 1
  ) |>
  dplyr::mutate(
    start = rowid,
    end = rowid
  ) |>
  dplyr::select(
    seqnames,
    start = rowid,
    end = rowid,
    phastCons100wayScore = "fixedStep chrom=chrM start=1 step=1"
  )

export(
  phastCons100way,
  "/home/liuc9/github/scMOCHA-data/config/chrM.phastCons100way.wigFix.qs"
)

import(
  "https://raw.githubusercontent.com/chunjie-sam-liu/scMOCHA-data/main/config/chrM.phastCons100way.wigFix.qs"
)
data.table::fread("https://raw.githubusercontent.com/chunjie-sam-liu/scMOCHA-data/main/config/chrM.phastCons100way.wigFix")
import("https://raw.githubusercontent.com/chunjie-sam-liu/scMOCHA-data/main/config/chrM.phastCons100way.wigFix.qs")


d <- import("/home/liuc9/github/scMOCHA-data/config/mtdna_genes_dloop.rds.gz")

export(
  d,
  "/home/liuc9/github/scMOCHA-data/config/mtdna_genes_dloop.qs"
)

export(
  gnomad,
  "/home/liuc9/github/scMOCHA-data/analysis/zzz/db/gnomad.qs"
)
convert(
  "/home/liuc9/github/scMOCHA-data/config/rCRS.MT.fasta.df.csv",
  "/home/liuc9/github/scMOCHA-data/config/rCRS.MT.fasta.df.fst"
)
convert(
  "/home/liuc9/github/scMOCHA-data/config/Mito-Genome-Loci-MitoMAP-Foswiki.csv",
  "/home/liuc9/github/scMOCHA-data/config/Mito-Genome-Loci-MitoMAP-Foswiki.fst"
)


d <- import("/home/liuc9/github/scMOCHA-data/analysis/zzz/used_samples/Fig_2c.csv")

d |>
  dplyr::count(tissue_site_detail)



d |>
  tidyr::nest(.by = tissue_site_detail) |>
  dplyr::mutate(
    a = purrr::map(
      .x = data,
      .f = \(.x) {
        lm(HFN ~ AGE, data = .x) |>
          broom::tidy() |>
          dplyr::slice(2)
      }
    )
  ) |>
  tidyr::unnest(cols = a) |>
  dplyr::arrange(-estimate)

d$AGE |> hist()


d |>
  dplyr::filter(
    tissue_site_detail == "Adrenal Gland"
  ) |>
  dplyr::arrange(AGE)

conn <- DBI::dbConnect(
  duckdb::duckdb(),
  dbdir = "/home/liuc9/github/scMOCHA-data/analysis/zzz/db/all_hetero_af.cell.ks_test/b_gseid_srrid_ks_load_p0.05_s25.duckdb.1.2.1"
)
DBI::dbListTables(conn)
dplyr::tbl(conn, "gseid_srrid_ks_load_p0.05_s25_unnest")

convert(
  "/home/liuc9/github/scMOCHA-data/config/rCRS.MT.fasta.df.fst",
  "/home/liuc9/github/scMOCHA-data/config/rCRS.MT.fasta.df.csv"
)


a <- import("/home/liuc9/github/scMOCHA-data/analysis/zzz/clean-data/gse_data.qs")

a


a |>
  dplyr::mutate(
    dplyr::across(
      dplyr::where(is.list),
      ~ jsonlite::toJSON(.x, auto_unbox = TRUE, null = "null")
    )
  ) ->
aa
