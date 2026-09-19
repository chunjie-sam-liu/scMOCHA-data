# @AUTHOR: Chunjie Liu
# @CONTACT: chunjie.sam.liu.at.gmail.com
# @DATE: 2026-09-15
# @DESCRIPTION: Install the R packages that exist only on GitHub, after
#   `pixi install` has built the conda side. Run with `pixi run remote-install`.
#   Anything available on conda-forge or bioconda belongs in pixi.toml instead.
# @VERSION: v0.1.0

# Install GitHub package if not already present ---------------------------
install_if_missing <- function(repo, pkg = NULL) {
  pkg_name <- pkg %||% basename(repo)
  if (requireNamespace(pkg_name, quietly = TRUE)) {
    cli::cli_alert_success("Already installed: {.pkg {pkg_name}}")
    return(invisible(TRUE))
  }
  cli::cli_progress_step("Installing {.val {repo}}")
  remotes::install_github(repo, upgrade = "never", dependencies = TRUE)
  invisible(TRUE)
}
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
install_if_missing("davidsjoberg/ggsankey")
install_if_missing("duckdb/duckdb-r", pkg = "duckdb")

# Domain packages. Uncomment the ones this project actually loads.
# install_if_missing("satijalab/seurat-data")
# install_if_missing("mojaveazure/seurat-disk")
# install_if_missing("satijalab/azimuth")
# install_if_missing("WenjianBi/SPACox")
# install_if_missing("QuantGen/BEDMatrix")
# install_if_missing("boxiangliu/locuscomparer")
# install_if_missing("perishky/meffil")
# install_if_missing("hhhh5/ewastools")

cli::cli_rule()

# preprocessCore: rebuild single-threaded if broken -----------------------
# The conda build links pthreads and every normalize.quantiles() call dies with
# "return code from pthread_create() is 22" on some clusters. minfi's
# preprocessQuantile / preprocessFunnorm depend on it. Probe at runtime so a
# working install is never reinstalled, and so this is a no-op in a project
# that never installs preprocessCore.
fix_preprocesscore <- function() {
  if (!requireNamespace("preprocessCore", quietly = TRUE)) {
    cli::cli_alert_info("preprocessCore not installed, skipping threading fix")
    return(invisible(FALSE))
  }
  probe <- try(
    preprocessCore::normalize.quantiles(matrix(stats::rnorm(30), ncol = 3)),
    silent = TRUE
  )
  if (!inherits(probe, "try-error")) {
    cli::cli_alert_success("preprocessCore threading OK")
    return(invisible(TRUE))
  }
  cli::cli_progress_step("Rebuilding preprocessCore with --disable-threading")
  BiocManager::install(
    "preprocessCore",
    configure.args = c(preprocessCore = "--disable-threading"),
    type = "source",
    force = TRUE,
    update = FALSE,
    ask = FALSE
  )
  invisible(TRUE)
}

cli::cli_h1("preprocessCore threading check")
fix_preprocesscore()
cli::cli_rule()
