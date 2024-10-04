compute_histND <- function(data, range_var, nbins) {
  n_vars <- dim(range_var)[1]  # Number of variables

  # Calculate the bin edges based on predefined ranges
  bin_edges <- lapply(seq_len(n_vars), function(v) {
    seq(range_var[v, 1], range_var[v, 2], length.out = nbins + 1)
  })

  # Initialize an empty n-dimensional histogram array
  hist_array <- array(0, dim = rep(nbins, n_vars))

  # For each observation, determine which bin it belongs to
  for (i in seq_len(nrow(data))) {
    # Determine bin indices for each variable
    indices <- mapply(function(value, edges) {
      # Find the interval for the value in the given bin edges
      which.min(value < edges) - 1
    }, data[i, ], bin_edges)

    # Ensure indices are within bounds (i.e., between 1 and `nbins`)
    indices[indices < 1] <- 1
    indices[indices > nbins] <- nbins

    # Increment the count in the corresponding bin
    hist_array[indices] <- hist_array[indices] + 1
  }

  # Flatten the n-dimensional histogram to a 1D vector
  hist_vector <- c(hist_array)

  return(hist_vector)
}
