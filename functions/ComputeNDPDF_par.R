# Main function to compute n-dimensional PDFs for two time periods
compute_nd_pdf_optimized <- function(variables, model_names, data_dir, year_present, year_future, lon, lat, range_var, nbins) {
  n_var <- length(variables)  # Number of variables
  num_models <- length(model_names)
  pdf_matrix_present <- array(NA, dim = c(length(lon), length(lat), nbins^n_var, num_models))  # PDF matrix for present period
  pdf_matrix_future <- array(NA, dim = c(length(lon), length(lat), nbins^n_var, num_models))   # PDF matrix for future period

  # Set up parallel processing for each model
  plan(multisession, workers = 4)  # Use a limited number of workers

  # Increase the maximum size of globals
  options(future.globals.maxSize = 8 * 1024^3)  # Allow up to 8 GiB for exporting globals

  # Process models in parallel for each time period
  pdf_matrix_list <- future_lapply(seq_along(model_names), function(m) {
    model_name <- model_names[m]
    cat(paste0("Processing model: ", model_name, " at ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n"))

    # Read data for all variables for the current model within each worker
    var_data_list_present <- list()
    var_data_list_future <- list()

    for (v in seq_along(variables)) {
      # Construct the file path for the current variable and model
      file_path <- paste0(data_dir, model_name, '/', variables[v], '/', list.files(path = paste0(data_dir, model_name, '/', variables[v], '/'), pattern = glob2rx(paste0(variables[v], "_", model_name, "*.nc")))[1])
      nc_var <- nc_open(file_path)

      # Extract the time range and data for the specified variables
      yyyy <- substr(as.character(nc.get.time.series(nc_var)), 1, 4)
      lon_var <- ncvar_get(nc_var, "lon")
      lat_var <- ncvar_get(nc_var, "lat")

      # Find indices matching the longitude and latitude ranges
      lon_indices <- which(lon_var %in% lon)
      lat_indices <- which(lat_var %in% lat)
      start_lon <- min(lon_indices)
      start_lat <- min(lat_indices)

      # Get the data slice for the present period
      iyyyy_present <- which(yyyy %in% year_present)
      var_data_present <- ncvar_get(nc_var, variables[v], start = c(start_lon, start_lat, min(iyyyy_present)), count = c(length(lon_indices), length(lat_indices), length(iyyyy_present)))

      # Get the data slice for the future period
      iyyyy_future <- which(yyyy %in% year_future)
      var_data_future <- ncvar_get(nc_var, variables[v], start = c(start_lon, start_lat, min(iyyyy_future)), count = c(length(lon_indices), length(lat_indices), length(iyyyy_future)))

      nc_close(nc_var)

      # Apply log transform for precipitation if required
      if (variables[v] == 'pr') {
        var_data_present <- log(var_data_present + 1)
        var_data_future <- log(var_data_future + 1)
      }

      var_data_list_present[[v]] <- var_data_present
      var_data_list_future[[v]] <- var_data_future
    }

    # Compute the PDF for each pixel sequentially for the present period
    model_pdf_matrix_present <- array(NA, dim = c(length(lon), length(lat), nbins^n_var))  # Temporary storage for present period's PDF
    for (i in seq_along(lon)) {
      for (j in seq_along(lat)) {
        # Collect data for each variable at the current pixel for present period
        pixel_data_present <- sapply(var_data_list_present, function(var) var[i, j, ])

        # Compute the n-dimensional histogram using `compute_histND`
        hist_tmp_present <- compute_histND(pixel_data_present, range_var[i, j, , ], nbins)

        # Normalize and store the histogram as a vector
        model_pdf_matrix_present[i, j, ] <- hist_tmp_present / sum(hist_tmp_present)
      }
    }

    # Compute the PDF for each pixel sequentially for the future period
    model_pdf_matrix_future <- array(NA, dim = c(length(lon), length(lat), nbins^n_var))  # Temporary storage for future period's PDF
    for (i in seq_along(lon)) {
      for (j in seq_along(lat)) {
        # Collect data for each variable at the current pixel for future period
        pixel_data_future <- sapply(var_data_list_future, function(var) var[i, j, ])

        # Compute the n-dimensional histogram using `compute_histND`
        hist_tmp_future <- compute_histND(pixel_data_future, range_var[i, j, , ], nbins)

        # Normalize and store the histogram as a vector
        model_pdf_matrix_future[i, j, ] <- hist_tmp_future / sum(hist_tmp_future)
      }
    }

    return(list(present = model_pdf_matrix_present, future = model_pdf_matrix_future))
  })

  # Combine results back into the main `pdf_matrix` for each period
  for (m in seq_along(model_names)) {
    pdf_matrix_present[, , , m] <- pdf_matrix_list[[m]]$present
    pdf_matrix_future[, , , m] <- pdf_matrix_list[[m]]$future
  }

  # Properly close parallel workers
  plan(sequential)  # Reset the plan to sequential to terminate workers

  return(list(present = pdf_matrix_present, future = pdf_matrix_future))
}
