#' Compute N-dimensional PDFs for Two Time Periods
#'
#' This function computes n-dimensional Probability Density Functions (PDFs)
#' for a set of climate variables across multiple climate models. The PDFs are
#' computed for two time periods (present and future) and for each grid point
#' based on the provided longitude and latitude coordinates.
#'
#' @param variables A character vector of climate variables (e.g., 'pr', 'tas', 'psl').
#' @param model_names A character vector of climate model names to be processed.
#' @param data_dir A character string specifying the directory where the climate model data is stored.
#' @param year_present A numeric vector specifying the years representing the present time period.
#' @param year_future A numeric vector specifying the years representing the future time period.
#' @param lon A numeric vector of longitude coordinates.
#' @param lat A numeric vector of latitude coordinates.
#' @param range_var An array of dimension `[lon, lat, nvar, 2]` containing the range (min/max) for each variable.
#' @param nbins A numeric value specifying the number of bins to use when computing histograms for each variable.
#'
#' @return A list containing two arrays:
#' \item{present}{The n-dimensional PDF array for the present period.}
#' \item{future}{The n-dimensional PDF array for the future period.}
#' The arrays have dimensions `[lon, lat, nbins^nvar, num_models]`.
#'
#' @import ncdf4 future future.apply
#' @export

compute_nd_pdf_optimized_0centered <- function(variables, model_names, data_dir, year_present, year_future, lon, lat, range_var, nbins, workers) {
  n_var <- length(variables)
  num_models <- length(model_names)
  pdf_matrix_present <- array(NA, dim = c(length(lon), length(lat), nbins^n_var, num_models))
  pdf_matrix_future <- array(NA, dim = c(length(lon), length(lat), nbins^n_var, num_models))

  # Set up parallel processing
  plan(multisession, workers = workers)
  options(future.globals.maxSize = 8 * 1024^3)

  pdf_matrix_list <- future_lapply(seq_along(model_names), function(m) {
    model_name <- model_names[m]
    cat(paste0("Processing model: ", model_name, " at ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n"))

    var_data_list_present <- list()
    var_data_list_future <- list()

    for (v in seq_along(variables)) {
      # Construct the file path
      file_path <- paste0(data_dir, model_name, '/', variables[v], '/', list.files(
        path = paste0(data_dir, model_name, '/', variables[v], '/'),
        pattern = glob2rx(paste0(variables[v], "_", model_name, "*.nc"))
      )[1])

      nc_var <- nc_open(file_path)

      # Extract longitude & latitude from NetCDF file
      lon_var <- ncvar_get(nc_var, "lon")  # NetCDF longitude: 0 to 359
      lat_var <- ncvar_get(nc_var, "lat")

      # ✅ **Convert NetCDF longitude (0-359°) to expected format (-180 to 179°)**
      lon_var <- ifelse(lon_var >= 180, lon_var - 360, lon_var)  # Shift longitudes
      sorted_indices <- order(lon_var)  # Sorting ensures correct mapping
      lon_var <- lon_var[sorted_indices]  # Apply sorted order

      # ✅ **Find correct indices based on the updated longitude & latitude**
      lon_indices <- match(lon, lon_var)
      lat_indices <- match(lat, lat_var)

      # ✅ **Remove any missing indices (handle mismatches safely)**
      lon_indices <- lon_indices[!is.na(lon_indices)]
      lat_indices <- lat_indices[!is.na(lat_indices)]

      if (length(lon_indices) == 0 || length(lat_indices) == 0) {
        stop("Error: Longitude or latitude indices not found in NetCDF file.")
      }

      start_lon <- min(lon_indices)
      start_lat <- min(lat_indices)

      # Extract the time variable and match years
      yyyy <- extract_years_from_time(nc_var)

      # Get present period indices
      iyyyy_present <- which(yyyy %in% year_present)
      if (length(iyyyy_present) == 0) {
        stop(paste("Error: No matching present years found in NetCDF for", model_name))
      }

      var_data_present <- ncvar_get(nc_var, variables[v],
                                    start = c(start_lon, start_lat, min(iyyyy_present)),
                                    count = c(length(lon_indices), length(lat_indices), length(iyyyy_present)))

      # Get future period indices
      iyyyy_future <- which(yyyy %in% year_future)
      if (length(iyyyy_future) == 0) {
        stop(paste("Error: No matching future years found in NetCDF for", model_name))
      }

      var_data_future <- ncvar_get(nc_var, variables[v],
                                    start = c(start_lon, start_lat, min(iyyyy_future)),
                                    count = c(length(lon_indices), length(lat_indices), length(iyyyy_future)))

      nc_close(nc_var)

      # Print extracted dimensions for debugging
      cat("Extracted dimensions for", variables[v], "\n")
      print(dim(var_data_present))
      print(dim(var_data_future))

      # Apply log transform for precipitation
      if (variables[v] == 'pr') {
        var_data_present <- log(var_data_present + 1)
        var_data_future <- log(var_data_future + 1)
      }

      var_data_list_present[[v]] <- var_data_present
      var_data_list_future[[v]] <- var_data_future
    }

    # Compute PDFs for present period
    model_pdf_matrix_present <- array(NA, dim = c(length(lon), length(lat), nbins^n_var))
    for (i in seq_along(lon)) {
      for (j in seq_along(lat)) {
        pixel_data_present <- sapply(var_data_list_present, function(var) var[i, j, ])
        hist_tmp_present <- compute_histND(pixel_data_present, range_var[i, j, , ], nbins)
        model_pdf_matrix_present[i, j, ] <- hist_tmp_present / sum(hist_tmp_present)
      }
    }

    # Compute PDFs for future period
    model_pdf_matrix_future <- array(NA, dim = c(length(lon), length(lat), nbins^n_var))
    for (i in seq_along(lon)) {
      for (j in seq_along(lat)) {
        pixel_data_future <- sapply(var_data_list_future, function(var) var[i, j, ])
        hist_tmp_future <- compute_histND(pixel_data_future, range_var[i, j, , ], nbins)
        model_pdf_matrix_future[i, j, ] <- hist_tmp_future / sum(hist_tmp_future)
      }
    }

    return(list(present = model_pdf_matrix_present, future = model_pdf_matrix_future))
  })

  # Combine results into the main PDF matrices
  for (m in seq_along(model_names)) {
    pdf_matrix_present[, , , m] <- pdf_matrix_list[[m]]$present
    pdf_matrix_future[, , , m] <- pdf_matrix_list[[m]]$future
  }

  plan(sequential)

  return(list(present = pdf_matrix_present, future = pdf_matrix_future))
}
