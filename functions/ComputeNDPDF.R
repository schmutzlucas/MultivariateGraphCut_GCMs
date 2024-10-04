compute_nd_pdf <- function(variables, model_names, data_dir, year_interest, lon, lat, range_var, nbins) {
  n_var <- length(variables)  # Number of variables
  num_models <- length(model_names)
  pdf_matrix <- array(NA, dim = c(length(lon), length(lat), nbins^n_var, num_models))  # Output array

  for (m in seq_along(model_names)) {
    model_name <- model_names[m]
    cat(paste0("Processing model: ", model_name, "\n"))  # Debug print for current model

    # Read data for all variables for the current model
    var_data_list <- list()
    for (v in seq_along(variables)) {
      # Read the NetCDF file for the variable
      file_path <- paste0(data_dir, model_name, '/', variables[v], '/', list.files(path = paste0(data_dir, model_name, '/', variables[v], '/'), pattern = glob2rx(paste0(variables[v], "_", model_name, "*.nc")))[1])
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
      var_data <- ncvar_get(nc_var, variables[v], start = c(start_lon, start_lat, min(iyyyy)), count = c(length(lon_indices), length(lat_indices), length(iyyyy)))
      nc_close(nc_var)

      if (variables[v] == 'pr') var_data <- log(var_data + 1)  # Log transform if variable is precipitation
      var_data_list[[v]] <- var_data
    }

    # Compute the histogram for each pixel in the lon-lat grid
    for (i in seq_along(lon)) {
      for (j in seq_along(lat)) {
        # Collect data for each variable at the current pixel
        pixel_data <- sapply(var_data_list, function(var) var[i, j, ])

        # Compute n-dimensional histogram using compute_histND
        hist_tmp <- compute_histND(pixel_data, range_var[i, j, , ], nbins)

        # Check that the sum of counts equals the number of observations
        total_counts <- sum(hist_vector)
        cat("Total counts in histogram:", total_counts, "\n")
        cat("Number of observations:", nrow(pixel_data), "\n")

        # Normalize and store the histogram as a vector in pdf_matrix
        pdf_matrix[i, j, , m] <- c(hist_tmp / sum(hist_tmp))

        asfafs
      }
    }
  }

  return(pdf_matrix)
}
