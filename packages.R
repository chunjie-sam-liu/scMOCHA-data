options(
  repos = c(CRAN = "https://cloud.r-project.org"),
  install.packages.check.source = "no"
)

Sys.setenv(
  R_REMOTES_NO_ERRORS_FROM_WARNINGS = "true",
  R_REMOTES_UPGRADE = "never"
)
is_installed_any <- function(pkgs) {
  for (p in pkgs) {
    if (suppressWarnings(requireNamespace(p, quietly = TRUE))) return(TRUE)
  }
  FALSE
}

pkg_candidates <- function(repo_name) {
  base <- sub('.*/', '', repo_name)
  parts <- strsplit(base, "-")[[1]]
  camel_lower <- if (length(parts) > 1) {
    paste0(
      parts[1],
      paste0(toupper(substring(parts[-1], 1, 1)), substring(parts[-1], 2))
    )
  } else {
    base
  }
  pascal <- paste0(
    toupper(substring(parts, 1, 1)),
    substring(parts, 2),
    collapse = ""
  )
  c(
    base,
    gsub("-", "", base),
    gsub("-", ".", base),
    camel_lower,
    pascal,
    paste0(toupper(substring(base, 1, 1)), substring(base, 2))
  ) |>
    unique()
}

install_if_missing <- function(repo, pkg = NULL) {
  if (is.null(pkg)) {
    candidates <- pkg_candidates(repo)
  } else {
    candidates <- c(pkg, pkg_candidates(repo)) |> unique()
  }

  if (is_installed_any(candidates)) {
    message(sprintf(
      "Already installed: %s (candidates: %s)",
      repo,
      paste(candidates, collapse = ", ")
    ))
    return(invisible(TRUE))
  }

  message(sprintf("Installing %s from GitHub...", repo))
  remotes::install_github(
    repo,
    upgrade = "never",
    dependencies = TRUE,
    force = TRUE
  )

  if (is_installed_any(candidates)) {
    message(sprintf("Installed: %s", repo))
    return(invisible(TRUE))
  }

  warning(sprintf(
    "Installation attempted but package not found for %s. Tried candidates: %s",
    repo,
    paste(candidates, collapse = ", ")
  ))
  invisible(FALSE)
}

# Install or verify required GitHub packages (only install when missing)
install_if_missing("RcppCore/RcppParallel")
install_if_missing("JanMarvin/openxlsx2")
install_if_missing("qsbase/qs2")
install_if_missing("chunjie-sam-liu/jutils")
install_if_missing("satijalab/seurat-data")
install_if_missing("mojaveazure/seurat-disk")
install_if_missing("duckdb/duckdb-r")

# Optional: azimuth (kept commented out)
# install_if_missing("satijalab/azimuth")
