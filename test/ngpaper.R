d2 <- readxl::read_excel(
  "/home/liuc9/github/scMOCHA-data/analysis/zzz/ngpaper/41588_2024_1838_MOESM4_ESM.xlsx",
  sheet = 2,
  skip = 3
)
d3 <- readxl::read_excel(
  "/home/liuc9/github/scMOCHA-data/analysis/zzz/ngpaper/41588_2024_1838_MOESM4_ESM.xlsx",
  sheet = 3,
  skip = 3
)

d4 <- readxl::read_excel(
  "/home/liuc9/github/scMOCHA-data/analysis/zzz/ngpaper/41588_2024_1838_MOESM4_ESM.xlsx",
  sheet = 4,
  skip = 3
)


d3 |> dplyr::arrange(-mean)
dplyr::filter(POS == 3243)
