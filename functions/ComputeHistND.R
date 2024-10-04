# Generalized function for n-dimensional histogram
compute_histND <- function(data, range_var, bins) {
  # Calculate bin edges for each variable using range_var
  bin_edges <- lapply(seq_along(data), function(v) seq(range_var[v, 1], range_var[v, 2], length.out = bins[v] + 1))

  # Use findInterval to place each data point in the correct bin for each dimension
  indices <- mapply(findInterval, data, bin_edges, SIMPLIFY = TRUE)

  # Check if data points are outside the defined bins
  inside <- apply(indices > 0 & indices <= bins, 1, all)

  # Create an n-dimensional array for the histogram
  hist_array <- array(0, bins)

  # Populate the histogram with counts
  for (k in which(inside)) {
    idx <- as.list(indices[k, ])
    hist_array[idx] <- hist_array[idx] + 1
  }

  return(hist_array)
}
