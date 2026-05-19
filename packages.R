options(
  repos = c(CRAN = "https://cloud.r-project.org"),
  install.packages.check.source = "no"
)

Sys.setenv(
  R_REMOTES_NO_ERRORS_FROM_WARNINGS = "true",
  R_REMOTES_UPGRADE = "never"
)
remotes::install_github(
  "RcppCore/RcppParallel",
  upgrade = "never",
  dependencies = TRUE,
  force = TRUE
)


remotes::install_github(
  "JanMarvin/openxlsx2",
  upgrade = "never",
  dependencies = TRUE,
  force = TRUE
)


remotes::install_github(
  "qsbase/qs2",
  upgrade = "never",
  dependencies = TRUE,
  force = TRUE
)

remotes::install_github(
  "chunjie-sam-liu/jutils",
  upgrade = "never",
  dependencies = TRUE,
  force = TRUE
)

remotes::install_github(
  "satijalab/seurat-data",
  upgrade = "never",
  dependencies = TRUE,
  force = TRUE
)

remotes::install_github(
  "mojaveazure/seurat-disk",
  upgrade = "never",
  dependencies = TRUE,
  force = TRUE
)

remotes::install_github(
  "duckdb/duckdb-r",
  upgrade = "never",
  dependencies = TRUE,
  force = TRUE
)

# remotes::install_github(
#   "satijalab/azimuth",
#   upgrade = "never",
#   dependencies = TRUE,
#   force = TRUE
# )
