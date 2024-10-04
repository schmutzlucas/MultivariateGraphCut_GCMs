# Main function to compute n-dimensional PDFs
compute_nd_pdf <- function(variable_list, model_names, data_dir, year_interest, lon, lat, range_var, nbins) {
  n_var <- length(variable_list)  # Number of variables
  num_models <- length(model_names)
  pdf_matrix <- array(NA, dim = c(length(lon), length(lat), prod(nbins), num_models))  # Output array

  for (m in seq_along(model_names)) {
    model_name <- model_names[m]

    # Read data for all variables for the current model
    var_data_list <- list()
    for (v in seq_along(variable_list)) {
      file_path <- paste0(data_dir, model_name, '/', variable_list[v], '/', list.files(path = paste0(data_dir, model_name, '/', variable_list[v], '/'), pattern = glob2rx(paste0(variable_list[v], "_", model_name, "*.nc")))[1])
      nc_var <- nc_open(file_path)

      yyyy <- substr(as.character(nc.get.time.series(nc_var)), 1, 4)
      iyyyy <- which(yyyy %in% year_interest)
      lon_var <- ncvar_get(nc_var, "lon")
      lat_var <- ncvar_get(nc_var, "lat")

      lon_indices <- which(lon_var %in% lon)
      lat_indices <- which(lat_var %in% lat)
      start_lon <- min(lon_indices)
      start_lat <- min(lat_indices)

      # Extract data slice for the specified variable and close the file
      var_data <- ncvar_get(nc_var, variable_list[v], start = c(start_lon, start_lat, min(iyyyy)), count = c(length(lon_indices), length(lat_indices), length(iyyyy)))
      nc_close(nc_var)

      if (variable_list[v] == 'pr') var_data <- log(var_data + 1)  # Log transform if variable is precipitation
      var_data_list[[v]] <- var_data
    }

    # Compute the histogram for each pixel in the lon-lat grid
    for (i in seq_along(lon)) {
      for (j in seq_along(lat)) {
        # Collect data for each variable at the current pixel
        pixel_data <- sapply(var_data_list, function(var) var[i, j, ])

        # Transpose the matrix to match the expected input shape for `compute_histND()`
        pixel_data <- t(pixel_data)
        print(dim(pixel_data))

        # Compute n-dimensional histogram using compute_histND
        hist_tmp <- compute_histND(pixel_data, range_var[i, j, , ], nbins)

        # Normalize and store the histogram as a vector in pdf_matrix
        pdf_matrix[i, j, , m] <- c(hist_tmp / sum(hist_tmp))
      }
    }
  }

  return(pdf_matrix)
}
