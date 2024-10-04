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