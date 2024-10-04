# Generalized function for n-dimensional histogram
compute_histND <- function(data, range_var, nbins) {
  # Determine the number of variables from the first dimension of `range_var`
  n_vars <- dim(range_var)[1]  # Number of variables

  # Calculate bin edges for each variable using a single nbins value
  bin_edges <- lapply(seq_len(n_vars), function(v) {
    seq(range_var[v, 1], range_var[v, 2], length.out = nbins + 1)
  })

  # Use `findInterval` to place each data point in the correct bin for each dimension
  indices <- mapply(findInterval, data, bin_edges, SIMPLIFY = TRUE)

  # Check if data points are inside the defined bins
  inside <- apply(indices > 0 & indices <= nbins, 1, all)

  # Create an n-dimensional array for the histogram with equal bins for each variable
  hist_array <- array(0, rep(nbins, n_vars))

  # Populate the histogram with counts
  for (k in which(inside)) {
    idx <- as.list(indices[k, ])
    hist_array[idx] <- hist_array[idx] + 1
  }

  return(hist_array)
}
