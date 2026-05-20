options(
  repos = c(CRAN = "https://cloud.r-project.org"),
  install.packages.check.source = "no"
)

Sys.setenv(
  R_REMOTES_NO_ERRORS_FROM_WARNINGS = "true",
  R_REMOTES_UPGRADE = "never"
)

# Install or verify required GitHub packages (only install when missing)
cli::cli_h1("GitHub Package Installation")
install_if_missing("RcppCore/RcppParallel")
install_if_missing("JanMarvin/openxlsx2")
install_if_missing("qsbase/qs2")
install_if_missing("chunjie-sam-liu/jutils")
install_if_missing("satijalab/seurat-data")
install_if_missing("mojaveazure/seurat-disk")
install_if_missing("duckdb/duckdb-r", pkg = "duckdb")
cli::cli_rule()

# Optional: azimuth (kept commented out)
# install_if_missing("satijalab/azimuth")
