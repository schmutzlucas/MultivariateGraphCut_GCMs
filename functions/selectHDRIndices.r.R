select_hdr_indices <- function(pdf_vector, tau = 0.10) {
  # pdf_vector: a flat vector of probability mass values (assumed normalized).
  # tau: fraction of the mass considered as outliers (e.g., tau=0.10 means keep the central 90%).

  # Compute the target mass for the central region.
  target_mass <- 1 - tau

  # Order indices by decreasing probability
  sorted_indices_desc <- order(pdf_vector, decreasing = TRUE)

  # Compute the cumulative sum of the sorted PDF values.
  cumulative_sum_desc <- cumsum(pdf_vector[sorted_indices_desc])

  # Find the index in the sorted order at which the cumulative sum reaches target_mass.
  threshold_index <- which(cumulative_sum_desc >= target_mass)[1]

  # Determine the threshold density, kappa, from the sorted vector.
  kappa <- pdf_vector[sorted_indices_desc[threshold_index]]

  # Return the indices of bins whose probability is at least kappa.
  central_indices <- which(pdf_vector >= kappa)
  return(central_indices)
}