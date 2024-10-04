compute_histND <- function(data, range_var, nbins) {
  n_vars <- dim(range_var)[1]  # Number of variables

  # Ensure `data` has observations as rows and variables as columns for `mapply`
  if (ncol(data) == n_vars) {
    data <- t(data)  # Transpose if necessary to match observations in rows and variables in columns
  }

  # Calculate bin edges for each variable using a single `nbins` value
  bin_edges <- lapply(seq_len(n_vars), function(v) {
    seq(range_var[v, 1], range_var[v, 2], length.out = nbins + 1)
  })

  # Use `mapply` to find the bin indices for each data point in each dimension
  indices <- mapply(findInterval, split(data, col(data)), bin_edges, SIMPLIFY = TRUE)

  # Ensure `indices` is a matrix with correct dimensions
  if (!is.matrix(indices)) {
    indices <- matrix(indices, ncol = n_vars)
  }

  # Check if data points are inside the defined bins
  if (nrow(indices) > 0) {
    inside <- apply(indices > 0 & indices <= nbins, 1, all)
  } else {
    inside <- logical(0)  # No valid rows, return empty logical vector
  }

  # Create an n-dimensional array for the histogram
  hist_array <- array(0, rep(nbins, n_vars))

  # Populate the histogram with counts if `inside` has valid entries
  if (length(inside) > 0) {
    for (k in which(inside)) {
      idx <- as.list(indices[k, ])
      hist_array[idx] <- hist_array[idx] + 1
    }
  }

  return(hist_array)
}
