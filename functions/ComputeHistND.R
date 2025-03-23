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
  n_vars <- ncol(data)   # number of variables
  n_obs  <- nrow(data)    # number of observations

  # Compute bin edges for each variable
  bin_edges <- lapply(1:n_vars, function(v) {
    seq(range_var[v, 1], range_var[v, 2], length.out = nbins + 1)
  })

  # For each variable, determine the bin index for each observation
  bin_indices <- sapply(1:n_vars, function(v) {
    inds <- findInterval(data[, v], vec = bin_edges[[v]], rightmost.closed = TRUE)
    inds[inds < 1] <- 1
    inds[inds > nbins] <- nbins
    return(inds)
  })

  # Compute linear indices for the n-dimensional histogram.
  # The overall number of bins is nbins^n_vars.
  dims <- rep(nbins, n_vars)
  multiplier <- cumprod(c(1, dims[-length(dims)]))
  idx_linear <- as.vector((as.matrix(bin_indices) - 1) %*% multiplier) + 1

  # Build histogram vector using tabulate; its length is exactly nbins^n_vars.
  hist_vector <- tabulate(idx_linear, nbins = prod(dims))

  return(hist_vector)
}

