gradient_hdist <- function(hdist_map) {
  # Initialize the gradient error matrix
  gradient_error <- matrix(0, nrow(hdist_map), ncol(hdist_map))
  denom <- matrix(0, nrow(hdist_map), ncol(hdist_map))  # To count contributing neighbors

  # Define neighboring indices
  ileft <- 1:(nrow(hdist_map) - 1)
  iright <- 2:nrow(hdist_map)
  itop <- 1:(ncol(hdist_map) - 1)
  ibottom <- 2:ncol(hdist_map)

  # Update denominator for the number of neighbors
  denom[ileft, ] <- denom[ileft, ] + 1
  denom[iright, ] <- denom[iright, ] + 1
  denom[, itop] <- denom[, itop] + 1
  denom[, ibottom] <- denom[, ibottom] + 1

  # Compute gradient errors
  gradient_error[ileft, ] <- gradient_error[ileft, ] + abs(hdist_map[ileft, ] - hdist_map[iright, ])
  gradient_error[iright, ] <- gradient_error[iright, ] + abs(hdist_map[iright, ] - hdist_map[ileft, ])
  gradient_error[, itop] <- gradient_error[, itop] + abs(hdist_map[, itop] - hdist_map[, ibottom])
  gradient_error[, ibottom] <- gradient_error[, ibottom] + abs(hdist_map[, ibottom] - hdist_map[, itop])

  # Normalize by the number of neighbors
  return(gradient_error / denom)
}
