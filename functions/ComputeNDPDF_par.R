# Main function to compute n-dimensional PDFs with efficient memory management
compute_nd_pdf_optimized <- function(variables, model_names, data_dir, year_interest, lon, lat, range_var, nbins) {
  n_var <- length(variables)  # Number of variables
  num_models <- length(model_names)
  pdf_matrix <- array(NA, dim = c(length(lon), length(lat), nbins^n_var, num_models))  # Output array

  # Set up parallel processing for each model
  plan(multisession, workers = 8)  # Use a limited number of workers
  
  # Increase the maximum size of globals
  options(future.globals.maxSize = 8 * 1024^3)  # Allow up to 8 GiB for exporting globals

  # Instead of passing large globals, we'll pass only the model name and read the data within the workers
  pdf_matrix_list <- future_lapply(seq_along(model_names), function(m) {
    model_name <- model_names[m]
    cat(paste0("Processing model: ", model_name, " at ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n"))

    # Read data for all variables for the current model within each worker
    var_data_list <- list()
    for (v in seq_along(variables)) {
      # Construct the file path for the current variable and model
      file_path <- paste0(data_dir, model_name, '/', variables[v], '/', list.files(path = paste0(data_dir, model_name, '/', variables[v], '/'), pattern = glob2rx(paste0(variables[v], "_", model_name, "*.nc")))[1])
      nc_var <- nc_open(file_path)

      # Extract the time range and data for the specified variables
      yyyy <- substr(as.character(nc.get.time.series(nc_var)), 1, 4)
      iyyyy <- which(yyyy %in% year_interest)
      lon_var <- ncvar_get(nc_var, "lon")
      lat_var <- ncvar_get(nc_var, "lat")

      # Find indices matching the longitude and latitude ranges
      lon_indices <- which(lon_var %in% lon)
      lat_indices <- which(lat_var %in% lat)
      start_lon <- min(lon_indices)
      start_lat <- min(lat_indices)

      # Extract the required data slice for the variable
      var_data <- ncvar_get(nc_var, variables[v], start = c(start_lon, start_lat, min(iyyyy)), count = c(length(lon_indices), length(lat_indices), length(iyyyy)))
      nc_close(nc_var)

      # Apply log transform for precipitation if required
      if (variables[v] == 'pr') var_data <- log(var_data + 1)
      var_data_list[[v]] <- var_data
    }

    # Compute the PDF for each pixel sequentially within each model
    model_pdf_matrix <- array(NA, dim = c(length(lon), length(lat), nbins^n_var))  # Temporary storage for each model's PDF
    for (i in seq_along(lon)) {
      for (j in seq_along(lat)) {
        # Collect data for each variable at the current pixel
        pixel_data <- sapply(var_data_list, function(var) var[i, j, ])

        # Compute the n-dimensional histogram using `compute_histND`
        hist_tmp <- compute_histND(pixel_data, range_var[i, j, , ], nbins)

        # Normalize and store the histogram as a vector
        model_pdf_matrix[i, j, ] <- hist_tmp / sum(hist_tmp)
      }
    }
    return(model_pdf_matrix)
  })

  # Combine results back into the main `pdf_matrix`
  for (m in seq_along(model_names)) {
    pdf_matrix[, , , m] <- pdf_matrix_list[[m]]
  }

  return(pdf_matrix)
}
