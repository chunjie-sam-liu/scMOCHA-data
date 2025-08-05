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

  d <- import(
    "/mnt/isilon/u01_project/large-scale/liuc9/raw/zzz/db/EXPR/gse_srrid_celltype_gene_expr.fst"
  )

  d |>
    tidyr::pivot_longer(
      cols = -c(gseid, srrid, genename),
      names_to = "celltype",
      values_to = "expr"
    ) -> dd

  dd |>
    dplyr::group_by(genename, celltype) |>
    tidyr::nest(.key = "expr") |>
    dplyr::ungroup() -> ddd

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
a <- import(
  "/mnt/isilon/u01_project/large-scale/liuc9/raw/zzz/clean-data/barcode_celltype.fst"
)
a$barcode |>
  unique() |>
  length()
nrow(a)


# ! test --------------------------------------------------------------------

old <- import(
  "/scr1/users/liuc9/tmp/scanpy/GSE155673_GSM4712885_celltype_gene_expr.csv"
)

new <- import(
  "/home/liuc9/github/scMOCHA-data/analysis/zzz/db/EXPR/GSE155673_GSM4712885_celltype_gene_expr.csv"
)

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
data.table::fread(
  "https://raw.githubusercontent.com/chunjie-sam-liu/scMOCHA-data/main/config/chrM.phastCons100way.wigFix"
)
import(
  "https://raw.githubusercontent.com/chunjie-sam-liu/scMOCHA-data/main/config/chrM.phastCons100way.wigFix.qs"
)


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


d <- import(
  "/home/liuc9/github/scMOCHA-data/analysis/zzz/used_samples/Fig_2c.csv"
)

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


a <- import(
  "/home/liuc9/github/scMOCHA-data/analysis/zzz/clean-data/gse_data.qs"
)

a


a |>
  dplyr::mutate(
    dplyr::across(
      dplyr::where(is.list),
      ~ jsonlite::toJSON(.x, auto_unbox = TRUE, null = "null")
    )
  ) -> aa


conn <- DBI::dbConnect(
  duckdb::duckdb(),
  dbdir = "/home/liuc9/github/scMOCHA-data/analysis/zzz/db/DUCKDB/cov.duckdb"
)
path <- "/mnt/isilon/u01_project/large-scale/liuc9/raw/zzz/db/PARQUET/covall/**/*.parquet"
sql <- glue(
  "CREATE OR REPLACE VIEW covall AS SELECT * FROM read_parquet('{path}')"
)
DBI::dbExecute(conn, sql)
DBI::dbListTables(conn)
a <- dplyr::tbl(conn, "covall")
dplyr::count(a)
a |>
  dplyr::select(
    gseid,
    srrid,
    barcode,
    base,
    `3173`
  ) |>
  dplyr::filter(
    base == "A"
  )

DBI::dbDisconnect(conn = conn, shutdown = TRUE)

d <- import(
  "/home/liuc9/github/scMOCHA-data/analysis/zzz/clean-data/all_variant.qs"
)

d |>
  dplyr::count(variant, issomatic) |>
  dplyr::filter(n > 1)

gse_data <- import(
  "analysis/zzz/clean-data/gse_data.qs"
)


regions_missalignment_error <- c(
  66:71,
  300:316,
  513:525,
  3106:3107,
  12418:12425,
  16182:16194
)
regions_rare_heteroplasmic_variants <- c(499, 538, 545, 10953, 12684)
variants_tobe_excluded <- c(
  regions_missalignment_error,
  regions_rare_heteroplasmic_variants
)

gse_data_variant_heteroplasmic <- import(
  "/home/liuc9/github/scMOCHA-data/analysis/zzz/clean-data/gse_data_variant_heteroplasmic.qs"
)

DBI::dbListTables(conn)

gse_data_variant_heteroplasmic |>
  dplyr::select(gseid, srrid, heteroplasmic) |>
  dplyr::mutate(
    a = purrr::map_chr(
      heteroplasmic,
      \(.x) {
        # tibble::tibble(
        #   hetero = .x$heteroplasmic_variant |>
        #     jsonlite::toJSON(
        #       auto_unbox = TRUE,
        #       null = "null"
        #     ),
        #   homo = .x$homoplasmic_variant |>
        #     jsonlite::toJSON(
        #       auto_unbox = TRUE,
        #       null = "null"
        #     )
        # )
        .x |> jsonlite::toJSON()
      }
    )
  ) |>
  dplyr::select(
    gseid,
    srrid,
    variant_alltype = a
  ) -> gseid_srrid_variant

DBI::dbWriteTable(
  conn,
  "gseid_srrid_variant",
  gseid_srrid_variant,
  overwrite = TRUE,
  temporary = FALSE
)
DBI::dbListTables(conn)


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


ALLVARIANTS <- import(file.path(
  "/home/liuc9/github/scMOCHA-data/analysis/zzz/clean-data/",
  "all_variant.qs"
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
  ) -> gseid_srrid_ks_load_variant

gseid_srrid_variant_celltype_ks_test <- import(
  file.path(
    ks_test_dir,
    "a_gseid_srrid_ks_load.nocellaf.qs"
  )
)


DBI::dbListTables(conn)
DBI::dbWriteTable(
  conn,
  "gseid_srrid_variant_celltype_ks_test",
  gseid_srrid_variant_celltype_ks_test,
  overwrite = TRUE,
  temporary = FALSE
)


d <- import(
  "/home/liuc9/github/scMOCHA-data/analysis/zzz/clean-data/gse_data.qs"
)


d |> dplyr::glimpse()
d$gseid[[1]] |>
  jsonlite::toJSON(auto_unbox = TRUE, null = "null") |>
  jsonlite::fromJSON()

d$haplo_violin[[1]] |> jsonlite::toJSON(auto_unbox = TRUE, null = "null") -> m

jsonlite::fromJSON(m)

d |>
  # dplyr::select(1:7) |>
  dplyr::mutate_if(
    dplyr::where(is.list),
    ~ {
      parallel::mclapply(
        X = .x,
        FUN = \(.xx) {
          jsonlite::toJSON(.xx, auto_unbox = TRUE, null = "null")
        },
        mc.cores = 20
      )
    }
  ) -> dd

DBI::dbListTables(conn_all_hetero_af)

DBI::dbWriteTable(
  conn_all_hetero_af,
  "gse_data",
  dd,
  overwrite = TRUE,
  temporary = FALSE
)

dplyr::tbl(conn_all_hetero_af, "gse_data") |>
  dplyr::collect() |>
  dplyr::mutate(
    metrics = purrr::map(
      .x = metrics,
      ~ jsonlite::fromJSON(.x)
    )
  ) -> ddd

dplyr::tbl(conn_all_hetero_af) |>
  dplyr::tbl("gseid_srrid_variant") |>
  dplyr::collect()

conn_all_hetero_af |> DBI::dbListTables()

# for somatic plot
tbl_gseid_srrid_variant <- conn_all_hetero_af |>
  dplyr::tbl("gseid_srrid_variant")


tbl_gseid_srrid_variant |>
  dplyr::collect() |>
  dplyr::mutate(
    a = purrr::map(
      .x = variant_alltype,
      ~ {
        jsonlite::fromJSON(.x) |>
          purrr::pluck("heteroplasmic_variant") -> .v
        if (length(.v) == 0) {
          return(NULL)
        } else {
          return(tibble::tibble(variant = .v))
        }
      }
    )
  ) |>
  dplyr::select(-variant_alltype) |>
  tidyr::unnest(cols = c(a)) -> gseid_srrid_variant_hetero


gseid_srrid_variant_hetero |>
  # head(5) |>
  dplyr::mutate(
    forplot = parallel::mcmapply(
      FUN = \(
        thevariant,
        thesrrid
      ) {
        .d <- fn_forplot(thevariant, thesrrid)
        jsonlite::toJSON(.d, auto_unbox = TRUE, null = "null")
      },
      thevariant = variant,
      thesrrid = srrid,
      mc.cores = 20,
      SIMPLIFY = FALSE
    )
  ) -> gseid_srrid_variant_hetero_somatic_forplot


# gseid_srrid_variant_hetero_somatic_forplot |>
#   dplyr::select(forplot) |>
#   tidyr::unnest(
#     cols = c(forplot)
#   ) -> gseid_srrid_variant_hetero_somatic_forplot_

DBI::dbListTables(conn_all_hetero_af)
DBI::dbWriteTable(
  conn_all_hetero_af,
  "gseid_srrid_variant_hetero_somatic_forplot",
  gseid_srrid_variant_hetero_somatic_forplot,
  overwrite = TRUE,
  temporary = FALSE
)
DBI::dbListTables(conn_all_variant_cell)
conn_all_variant_cell |>
  dplyr::tbl("all_variant_cell")

DBI::dbListTables(conn_all_variant_cell)
# copy table all_variant_cell in  conn_all_variant_cell in to conn_all_hetero_af
DBI::dbWriteTable(
  conn_all_hetero_af,
  "all_variant_cell",
  dplyr::tbl(conn_all_variant_cell, "all_variant_cell") |> dplyr::collect(),
  overwrite = TRUE,
  temporary = FALSE
)

# DBI::dbListTables(conn_all_hetero_af)
# DBI::dbExecute(
#   conn_all_hetero_af,
#   "
#   ALTER TABLE allvariant_cell RENAME TO allvariants_cell
# "
# )

DBI::dbListTables(conn_all_hetero_af)
