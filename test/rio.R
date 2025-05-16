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
