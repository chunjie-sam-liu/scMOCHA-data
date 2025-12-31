# Function to find the peak of a density plot
find_density_peak <- function(x, bw = "nrd0", ...) {
  # Calculate density
  dx <- density(x, bw = bw, ...)

  # Find the index of maximum density
  peak_idx <- which.max(dx$y)

  # Return peak information
  peak_info <- list(
    x = dx$x[peak_idx], # x-coordinate of peak
    y = dx$y[peak_idx], # y-coordinate (density value) of peak
    density_obj = dx # full density object for plotting
  )

  return(peak_info)
}

# Function to find multiple peaks (local maxima)
find_density_peaks <- function(
  x,
  bw = "nrd0",
  min_height = 0.01,
  min_prominence = 0.1,
  min_distance = 0.05,
  ...
) {
  # Calculate density
  dx <- density(x, bw = bw, ...)

  # Find local maxima
  n <- length(dx$y)
  peaks <- c()

  # Check each point (excluding endpoints)
  for (i in 2:(n - 1)) {
    if (
      dx$y[i] > dx$y[i - 1] && dx$y[i] > dx$y[i + 1] && dx$y[i] >= min_height
    ) {
      peaks <- c(peaks, i)
    }
  }

  # Filter peaks by prominence and distance
  if (length(peaks) > 0) {
    peaks_info <- data.frame(
      idx = peaks,
      x = dx$x[peaks],
      y = dx$y[peaks]
    )

    # Calculate prominence for each peak
    peaks_info$prominence <- 0
    for (i in 1:nrow(peaks_info)) {
      peak_idx <- peaks_info$idx[i]
      peak_height <- peaks_info$y[i]

      # Find minimum heights on both sides
      left_min <- min(dx$y[1:peak_idx])
      right_min <- min(dx$y[peak_idx:n])
      baseline <- max(left_min, right_min)

      peaks_info$prominence[i] <- peak_height - baseline
    }

    # Filter by prominence
    peaks_info <- peaks_info[peaks_info$prominence >= min_prominence, ]

    # Filter by minimum distance between peaks
    if (nrow(peaks_info) > 1) {
      # Sort by height (keep strongest peaks when too close)
      peaks_info <- peaks_info[order(peaks_info$y, decreasing = TRUE), ]

      keep <- rep(TRUE, nrow(peaks_info))
      for (i in 1:(nrow(peaks_info) - 1)) {
        if (keep[i]) {
          for (j in (i + 1):nrow(peaks_info)) {
            if (abs(peaks_info$x[i] - peaks_info$x[j]) < min_distance) {
              keep[j] <- FALSE
            }
          }
        }
      }
      peaks_info <- peaks_info[keep, ]
    }

    # Sort by density value (highest first) and remove helper columns
    peaks_info <- peaks_info[
      order(peaks_info$y, decreasing = TRUE),
      c("x", "y", "prominence")
    ]
  } else {
    peaks_info <- data.frame(
      x = numeric(0),
      y = numeric(0),
      prominence = numeric(0)
    )
  }

  return(list(
    peaks = peaks_info,
    density_obj = dx
  ))
}

x <- tbl_thevariant_data |>
  dplyr::collect() |>
  dplyr::filter(celltype == "NK") |>
  dplyr::pull(af)

dx <- density(x)

plot(dx, main = "Density of x", xlab = "x", ylab = "Density")


# Find multiple peaks with adaptive filtering for bimodal distributions
max_density <- max(dx$y)
adaptive_min_height <- max_density * 0.1 # 10% of maximum density
adaptive_min_prominence <- max_density * 0.2 # 20% of maximum density

peaks_result <- find_density_peaks(
  x,
  min_height = adaptive_min_height,
  min_prominence = adaptive_min_prominence,
  min_distance = 0.2
)
cat("Number of significant peaks found:", nrow(peaks_result$peaks), "\n")
cat(
  "Adaptive thresholds - min_height:",
  round(adaptive_min_height, 3),
  ", min_prominence:",
  round(adaptive_min_prominence, 3),
  "\n"
)

if (nrow(peaks_result$peaks) > 0) {
  print(peaks_result$peaks)

  # Add all peaks to the plot
  points(
    peaks_result$peaks$x,
    peaks_result$peaks$y,
    col = "blue",
    pch = 17,
    cex = 1.2
  )
  for (i in 1:nrow(peaks_result$peaks)) {
    text(
      peaks_result$peaks$x[i],
      peaks_result$peaks$y[i] + 0.1,
      paste(
        "P",
        i,
        " (prom:",
        round(peaks_result$peaks$prominence[i], 2),
        ")",
        sep = ""
      ),
      col = "blue",
      pos = 3,
      cex = 0.8
    )
  }
}
