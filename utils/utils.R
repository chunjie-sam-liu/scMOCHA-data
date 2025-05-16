#' Format Numbers for Human Readability
#'
#' This function formats numbers to be more easily read by humans,
#' adapting the format based on the magnitude of the number.
#'
#' @param .x A numeric value to be formatted
#'
#' @return A character string representing the formatted number
#'
#' @details
#' The function applies different formatting based on the value range:
#' - For values >= 0.1: Uses 2 significant digits
#' - For values between 0.001 and 0.1: Uses 2 significant digits
#' - For values < 0.001 but > 0: Uses scientific notation with 3 digits
#' - For 0: Returns "0"
#' - Negative numbers preserve their sign in the output
#'
#' @examples
#' human_read(123.456) # Returns "123"
#' human_read(0.0456) # Returns "0.046"
#' human_read(0.0000123) # Returns "1.23e-05"
#' human_read(-42.5) # Returns "-42.5"
#'
#' @export

human_read <- function(.x) {
  .sign = ifelse(.x < 0, TRUE, FALSE)
  .x <- abs(.x)

  if (.x >= 0.1) {
    .x %>%
      signif(digits = 2) %>%
      toString() -> .xx
  } else if (.x < 0.1 && .x >= 0.001) {
    .x %>%
      signif(digits = 2) %>%
      toString() -> .xx
  } else if (.x < 0.001 && .x > 0) {
    .x %>% format(digits = 3, scientific = TRUE) -> .xx
  } else {
    .xx <- "0"
  }

  ifelse(.sign, paste0("-", .xx), .xx)
}

#' Format P-values for LaTeX or Human Readable Output
#'
#' This function takes a P-value and formats it for human-readable presentation,
#' particularly in LaTeX format. It handles both standard decimal and scientific notation.
#'
#' @param .x Character. The P-value to format.
#' @param .s Character. Optional statistic value to include before the P-value. Default is NA.
#' @param .tex Logical. If TRUE (default), return a LaTeX expression object.
#'        If FALSE, return a character string.
#'
#' @return If .tex is TRUE, returns a LaTeX expression object that can be used
#'         in plotting functions. If .tex is FALSE, returns a character string.
#'
#' @details The function handles P-values in both decimal format (e.g., "0.05")
#'          and scientific notation (e.g., "1e-5"). For scientific notation,
#'          it formats the output as "P=1 × 10^{-5}" in LaTeX format.
#'
#'          If a statistic value .s is provided, it will be included before the P-value
#'          in the format "statistic, P=value".
#'
#' @examples
#' human_read_latex_pval("0.05") # Returns LaTeX expression for P=0.05
#' human_read_latex_pval("1e-5") # Returns LaTeX expression for P=1 × 10^{-5}
#' human_read_latex_pval("0.01", "t = 2.45") # Returns with statistic included
#' human_read_latex_pval("1e-5", .tex = FALSE) # Returns character string
#' human_read_latex_pval(
#'   .x = human_read(cor_test$p.value),
#'   .s = glue::glue("R={round(cor_test$estimate,3)}")
#' )

#' @importFrom glue glue
#' @importFrom latex2exp TeX
human_read_latex_pval <- function(.x, .s = NA, .tex = TRUE) {
  if (is.na(.s)) {
    if (grepl(pattern = "e", x = .x)) {
      sub("-0", "-", strsplit(split = "e", x = .x, fixed = TRUE)[[1]]) -> .xx
      thestr <- glue::glue("$\\textit{P}=<<.xx[1]>> \\times 10^{<<.xx[2]>>}$", .open = "<<", .close = ">>")
    } else {
      thestr <- glue::glue("$\\textit{P}=<<.x>>$", .open = "<<", .close = ">>")
    }
  } else {
    if (grepl(pattern = "e", x = .x)) {
      sub("-0", "-", strsplit(split = "e", x = .x, fixed = TRUE)[[1]]) -> .xx
      thestr <- glue::glue("<<.s>>, $\\textit{P}=<<.xx[1]>> \\times 10^{<<.xx[2]>>}$", .open = "<<", .close = ">>")
    } else {
      thestr <- glue::glue("<<.s>>, $\\textit{P}=<<.x>>$", .open = "<<", .close = ">>")
    }
  }
  if (isTRUE(.tex)) {
    latex2exp::TeX(thestr)
  } else {
    thestr
  }
}

#' Export Data to Various File Formats
#'
#' This function exports data to files in various formats, based on a specified format
#' or detected from the file extension.
#'
#' @param x Data object. The data to be exported.
#' @param file Character. The path to save the file.
#' @param format Character. The format to use for export. Options are "qs", "fst", "csv", "tsv", "xlsx",
#'        "json", "yaml", or "feather". If not specified, the format is detected from the file extension.
#'
#' @return Invisibly returns the base filename (without extension) used for export.
#'
#' @details
#' The function uses different packages based on the format:
#' - qs: Uses qs::qsave with 3 threads for performance
#' - csv: Uses data.table::fwrite for fast CSV export
#' - tsv: Uses data.table::fwrite with tab separator for TSV export
#' - fst: Uses fst::write_fst for binary CSV export
#' - xlsx: Uses writexl::write_xlsx for Excel export
#' - json: Uses jsonlite::toJSON for JSON export
#' - yaml: Uses yaml::write_yaml for YAML export
#' - feather: Uses arrow::write_feather for Feather export
#'
#' @examples
#' # Export based on file extension
#' export(mydata, "mydata.csv")
#' export(mydata, "mydata.tsv")
#' export(mydata, "mydata.qs")
#'
#' # Explicitly specify format
#' export(mydata, "mydata", format = "fst")
#'
export <- function(x, file, format = c("qs", "csv", "tsv", "fst", "xlsx", "json", "yaml", "feather")) {
  # Match and validate the format argument
  format <- match.arg(format)

  # Get the file extension
  ext <- tools::file_ext(file)

  # If extension is one of our supported formats and no format specified, use the extension
  if (ext %in% c("qs", "csv", "tsv", "fst", "xlsx", "json", "yaml", "feather") && missing(format)) {
    format <- ext
  }

  # Create new filename without extension
  new_file <- sub(paste0("\\.", ext, "$"), "", file)

  # Save based on the specified format
  if (format == "qs") {
    qs::qsave(x, file = paste0(new_file, ".qs"), nthreads = 3)
    message("Saved as ", paste0(new_file, ".qs"))
  }

  if (format == "csv") {
    data.table::fwrite(x, file = paste0(new_file, ".csv"))
    message("Saved as ", paste0(new_file, ".csv"))
  }

  if (format == "tsv") {
    data.table::fwrite(x, file = paste0(new_file, ".tsv"), sep = "\t")
    message("Saved as ", paste0(new_file, ".tsv"))
  }

  if (format == "fst") {
    fst::write_fst(x, path = paste0(new_file, ".fst"))
    message("Saved as ", paste0(new_file, ".fst"))
  }

  if (format == "xlsx") {
    writexl::write_xlsx(x, path = paste0(new_file, ".xlsx"))
    message("Saved as ", paste0(new_file, ".xlsx"))
  }

  if (format == "json") {
    jsonlite::write_json(x, path = paste0(new_file, ".json"))
    message("Saved as ", paste0(new_file, ".json"))
  }

  if (format == "yaml") {
    yaml::write_yaml(x, file = paste0(new_file, ".yaml"))
    message("Saved as ", paste0(new_file, ".yaml"))
  }

  if (format == "feather") {
    arrow::write_feather(x, sink = paste0(new_file, ".feather"))
    message("Saved as ", paste0(new_file, ".feather"))
  }

  # Return the base filename used (without extension)
  invisible(new_file)
}

#' Import Data from Various File Formats
#'
#' This function imports data from files in various formats, automatically detecting
#' the format based on file extension or using a specified format.
#'
#' @param file Character. The path to the file to import.
#' @param format Character. Optional. The format to use for import. If not specified,
#'        the format is detected from the file extension.
#'        Supported formats are: "qs", "rds", "csv", "tsv", "fst", "xlsx", "json", "yaml", "feather".
#'
#' @return The imported data
#'
#' @details
#' The function uses different packages based on the format:
#' - qs: Uses qs::qread
#' - rds: Uses readr::read_rds
#' - csv: Uses data.table::fread
#' - tsv: Uses data.table::fread with tab separator
#' - fst: Uses fst::read_fst
#' - xlsx: Uses readxl::read_xlsx
#' - json: Uses jsonlite::fromJSON
#' - yaml: Uses yaml::read_yaml
#' - feather: Uses arrow::read_feather
#'
#' @examples
#' # Import based on file extension
#' data <- import("mydata.csv")
#' data <- import("mydata.tsv")
#' data <- import("mydata.qs")
#' data <- import("mydata.rds")
#'
#' # Explicitly specify format
#' data <- import("mydata.txt", format = "csv")
#' data <- import("mydata.txt", format = "tsv")
#'
#' @export
import <- function(file, format = NULL) {
  # Check if file exists
  if (!file.exists(file)) {
    stop("File does not exist: ", file)
  }

  # Get the file extension if format is not specified
  if (is.null(format)) {
    format <- tools::file_ext(file)
  }

  # Convert to lowercase for case-insensitive matching
  format <- tolower(format)

  # Import based on the format
  switch(format,
    "qs" = {
      message("Importing qs file: ", file)
      qs::qread(file)
    },
    "rds" = {
      message("Importing rds file: ", file)
      readr::read_rds(file)
    },
    "csv" = {
      message("Importing csv file: ", file)
      data.table::fread(file)
    },
    "tsv" = {
      message("Importing tsv file: ", file)
      data.table::fread(file, sep = "\t")
    },
    "fst" = {
      message("Importing fst file: ", file)
      fst::read_fst(file)
    },
    "xlsx" = {
      message("Importing xlsx file: ", file)
      readxl::read_xlsx(file)
    },
    "json" = {
      message("Importing json file: ", file)
      jsonlite::fromJSON(file)
    },
    "yaml" = {
      message("Importing yaml file: ", file)
      yaml::read_yaml(file)
    },
    "feather" = {
      message("Importing feather file: ", file)
      arrow::read_feather(file)
    },
    {
      stop("Unsupported file format: ", format)
    }
  )
}
