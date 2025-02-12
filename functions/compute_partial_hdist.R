compute_partial_hdist <- function(pdf_ref, pdf_model, selected_indices) {
  # Ensure pdf_model can be 1D (vector), 2D (bins x models), 3D (lon x lat x bins x models)
  dim_ref <- length(dim(pdf_ref))
  dim_model <- length(dim(pdf_model))

  # Case 1: Single grid point, single model (1D vector)
  if (is.null(dim(pdf_ref)) && is.null(dim(pdf_model))) {
    if (length(selected_indices) > 0) {
      return(sqrt(sum((sqrt(pdf_model[selected_indices]) - sqrt(pdf_ref[selected_indices]))^2)) / sqrt(2))
    } else {
      return(NA)
    }
  }

  # Case 2: Single grid point, multiple models (2D: bins x models)
  if (dim_ref == 1 && dim_model == 2) {
    n_models <- dim(pdf_model)[2]
    h_dist_partial <- rep(NA, n_models)

    for (m in 1:n_models) {
      if (length(selected_indices) > 0) {
        h_dist_partial[m] <- sqrt(sum((sqrt(pdf_model[selected_indices, m]) - sqrt(pdf_ref[selected_indices]))^2)) / sqrt(2)
      }
    }
    return(h_dist_partial)
  }

  # Case 3: Full grid (3D: lon x lat x bins) or (4D: lon x lat x bins x models)
  lon_size <- dim(pdf_ref)[1]
  lat_size <- dim(pdf_ref)[2]
  n_models <- ifelse(length(dim(pdf_model)) == 4, dim(pdf_model)[4], 1)

  # Initialize output
  h_dist_partial_future <- array(NA, dim = c(lon_size, lat_size, n_models))

  for (j in 1:lat_size) {
    for (i in 1:lon_size) {
      selected_indices_grid <- selected_indices[[i]][[j]]  # Retrieve bin indices for this grid point

      if (length(selected_indices_grid) > 0) {
        for (m in 1:n_models) {
          model_pdf <- if (n_models == 1) pdf_model[i, j, selected_indices_grid] else pdf_model[i, j, selected_indices_grid, m]

          h_dist_partial_future[i, j, m] <- sqrt(
            sum((sqrt(model_pdf) - sqrt(pdf_ref[i, j, selected_indices_grid]))^2)
          ) / sqrt(2)
        }
      }
    }
  }

  return(h_dist_partial_future)
}
