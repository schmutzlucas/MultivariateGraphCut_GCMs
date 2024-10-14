#' Compute N-dimensional Histogram
#'
#' This function computes an n-dimensional histogram from a given dataset, where each variable
#' is assigned to a specific bin based on the provided range and number of bins. The function
#' returns a flattened vector representing the histogram counts for each combination of bins
#' across all variables.
#'
#' @param data A matrix where each column represents a variable and each row is an observation.
#' @param range_var A matrix of dimension `[n_vars, 2]` that defines the minimum and maximum values
#' for each variable. The first column corresponds to the minimum values, and the second column
#' corresponds to the maximum values.
#' @param nbins A numeric value specifying the number of bins for each variable.
#'
#' @return A numeric vector representing the flattened n-dimensional histogram. The vector has length
#' equal to `nbins^n_vars`, where `n_vars` is the number of variables.
#'
#' @examples
#' # Example data with two variables
#' data <- matrix(rnorm(100), ncol = 2)  # 100 observations of 2 variables
#' range_var <- matrix(c(-3, 3, -3, 3), ncol = 2)  # Min/max for each variable
#' nbins <- 10  # Number of bins per variable
#'
#' # Compute the n-dimensional histogram
#' hist_vector <- compute_histND(data, range_var, nbins)
#'
#' @export
compute_histND <- function(data, range_var, nbins) {
  n_vars <- ncol(data)  # Number of variables
  n_obs <- nrow(data)   # Number of observations

  # Calculate the bin edges for each variable
  bin_edges <- lapply(1:n_vars, function(v) {
    seq(range_var[v, 1], range_var[v, 2], length.out = nbins + 1)
  })

  # Initialize a matrix to hold bin indices for each observation and variable
  bin_indices <- matrix(NA, nrow = n_obs, ncol = n_vars)

  # Assign each data point to a bin
  for (v in 1:n_vars) {
    bin_indices[, v] <- findInterval(data[, v], vec = bin_edges[[v]], rightmost.closed = TRUE)
    # Correct any indices that are zero or out of bounds
    bin_indices[, v][bin_indices[, v] < 1] <- 1
    bin_indices[, v][bin_indices[, v] > nbins] <- nbins
  }

  # Compute linear indices for the n-dimensional histogram
  dims <- rep(nbins, n_vars)
  multiplier <- cumprod(c(1, dims[-length(dims)]))
  idx_linear <- as.numeric((bin_indices - 1) %*% multiplier) + 1

  # Initialize the histogram vector
  hist_vector <- rep(0, prod(dims))

  # Count occurrences using table()
  counts <- table(idx_linear)

  # Update the histogram vector
  hist_vector[as.numeric(names(counts))] <- counts

  return(hist_vector)
}