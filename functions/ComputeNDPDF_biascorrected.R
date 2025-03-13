compute_nd_pdf_bias_corrected <- function(variables, reference_name, model_names, data_dir,
                                            year_present, year_future, lon, lat, nbins, workers,
                                            buffer = 0.05) {
  # Assumes that ncdf4, future, and future.apply packages (and the helper functions
  # extract_years_from_time and compute_histND) are already loaded.

  n_vars <- length(variables)
  nlon <- length(lon)
  nlat <- length(lat)

  ## ----------------------------
  ## SEGMENT 1: Process the reference data
  ## ----------------------------

  # Initialize arrays for joint PDFs (present and future)
  pdf_ref_present <- array(NA, dim = c(nlon, nlat, nbins^n_vars))
  pdf_ref_future  <- array(NA, dim = c(nlon, nlat, nbins^n_vars))

  # Lists to store per-variable reference statistics (each statistic is computed per grid point)
  reference_stats_present <- vector("list", n_vars)
  reference_stats_future  <- vector("list", n_vars)

  # Arrays to store the buffered range for PDF binning:
  # dimensions: [nlon, nlat, n_vars, 2] for present and future respectively.
  ref_range_present <- array(NA, dim = c(nlon, nlat, n_vars, 2))
  ref_range_future  <- array(NA, dim = c(nlon, nlat, n_vars, 2))

  # We also collect the raw reference data (for each variable and period) in a joint array for PDF computation.
  # Dimensions for present: [nlon, nlat, time_pres, n_vars]
  # and for future: [nlon, nlat, time_fut, n_vars]
  ref_data_present_all <- NULL
  ref_data_future_all  <- NULL

  for (v in seq_along(variables)) {
    var <- variables[v]

    # Construct file path for the reference
    ref_dir <- paste0(data_dir, reference_name, '/', var, '/')
    ref_file <- list.files(path = ref_dir, pattern = glob2rx(paste0(var, "_", reference_name, "*.nc")), full.names = TRUE)[1]
    if (is.na(ref_file)) stop("Reference file for ", var, " not found.")

    nc_ref <- nc_open(ref_file)
    # --- Reordering longitude ---
    lon_file <- ncvar_get(nc_ref, "lon")
    lat_file <- ncvar_get(nc_ref, "lat")
    if(any(lon_file >= 180)) {
      lon_file <- ifelse(lon_file >= 180, lon_file - 360, lon_file)
    }
    sorted_indices <- order(lon_file)
    lon_file_sorted <- lon_file[sorted_indices]
    # Match user-supplied lon and lat with file values (assumed lat is already in desired order)
    lon_indices <- match(lon, lon_file_sorted)
    lat_indices <- match(lat, lat_file)
    if(any(is.na(lon_indices)) || any(is.na(lat_indices)))
      stop("Grid indices for reference not found.")

    # Extract time information
    yyyy <- extract_years_from_time(nc_ref)
    iyear_pres <- which(yyyy %in% year_present)
    if (length(iyear_pres) == 0)
      stop("No present years found for reference ", reference_name)
    iyear_fut <- which(yyyy %in% year_future)
    if (length(iyear_fut) == 0)
      stop("No future years found for reference ", reference_name)

    # Read full spatial subset and reorder along longitude.
    full_data_pres <- ncvar_get(nc_ref, var,
                                start = c(1, min(lat_indices), min(iyear_pres)),
                                count = c(-1, length(lat_indices), length(iyear_pres)))
    full_data_pres <- full_data_pres[sorted_indices, , ]
    ref_data_pres <- full_data_pres[lon_indices, , ]

    full_data_fut <- ncvar_get(nc_ref, var,
                               start = c(1, min(lat_indices), min(iyear_fut)),
                               count = c(-1, length(lat_indices), length(iyear_fut)))
    full_data_fut <- full_data_fut[sorted_indices, , ]
    ref_data_fut <- full_data_fut[lon_indices, , ]

    nc_close(nc_ref)

    # Initialize matrices to hold per-grid point statistics for present period
    m_pres   <- matrix(NA, nlon, nlat)
    s_pres   <- matrix(NA, nlon, nlat)
    min_pres <- matrix(NA, nlon, nlat)
    max_pres <- matrix(NA, nlon, nlat)
    if (var == "pr") q90_pres <- matrix(NA, nlon, nlat)

    # And for future period:
    m_fut   <- matrix(NA, nlon, nlat)
    s_fut   <- matrix(NA, nlon, nlat)
    min_fut <- matrix(NA, nlon, nlat)
    max_fut <- matrix(NA, nlon, nlat)
    if (var == "pr") q90_fut <- matrix(NA, nlon, nlat)

    # Loop over each grid point and compute statistics (using the raw time series)
    for (i in seq_len(nlon)) {
      for (j in seq_len(nlat)) {
        ts_pres <- ref_data_pres[i, j, ]
        ts_fut  <- ref_data_fut[i, j, ]
        m_pres[i,j]   <- mean(ts_pres, na.rm = TRUE)
        s_pres[i,j]   <- sd(ts_pres, na.rm = TRUE)
        min_pres[i,j] <- min(ts_pres, na.rm = TRUE)
        max_pres[i,j] <- max(ts_pres, na.rm = TRUE)
        if (var == "pr")
          q90_pres[i,j] <- as.numeric(quantile(ts_pres, 0.90, na.rm = TRUE))

        m_fut[i,j]   <- mean(ts_fut, na.rm = TRUE)
        s_fut[i,j]   <- sd(ts_fut, na.rm = TRUE)
        min_fut[i,j] <- min(ts_fut, na.rm = TRUE)
        max_fut[i,j] <- max(ts_fut, na.rm = TRUE)
        if (var == "pr")
          q90_fut[i,j] <- as.numeric(quantile(ts_fut, 0.90, na.rm = TRUE))
      }
    }

    # Store computed statistics (for non-pr variables these are on raw data;
    # for pr we use raw values for Q90; later we transform the data for histogram computation)
    reference_stats_present[[v]] <- list(mean = m_pres, sd = s_pres, min = min_pres, max = max_pres)
    if (var == "pr")
      reference_stats_present[[v]]$q90 <- q90_pres
    reference_stats_future[[v]] <- list(mean = m_fut, sd = s_fut, min = min_fut, max = max_fut)
    if (var == "pr")
      reference_stats_future[[v]]$q90 <- q90_fut

    # For precipitation, apply log transform now (after statistics are computed on raw data)
    if (var == "pr") {
      ref_data_pres <- log(ref_data_pres + 1)
      ref_data_fut  <- log(ref_data_fut + 1)
    }

    # Compute buffered range for PDF binning at each grid point.
    # For each grid cell, range = [min - buffer*(max-min), max + buffer*(max-min)]
    # Note: For pr, we recompute the range on the log-transformed data.
    range_pres <- matrix(NA, nrow = nlon * nlat, ncol = 2)
    range_fut  <- matrix(NA, nrow = nlon * nlat, ncol = 2)
    for (idx in 1:(nlon * nlat)) {
      i <- ((idx - 1) %% nlon) + 1
      j <- ((idx - 1) %/% nlon) + 1
      if (var == "pr") {
        # Use the log-transformed data for range computation
        diff_pres <- max(ref_data_pres[i,j, ], na.rm = TRUE) - min(ref_data_pres[i,j, ], na.rm = TRUE)
        range_pres[idx, 1] <- min(ref_data_pres[i,j, ], na.rm = TRUE) - buffer * diff_pres
        range_pres[idx, 2] <- max(ref_data_pres[i,j, ], na.rm = TRUE) + buffer * diff_pres

        diff_fut <- max(ref_data_fut[i,j, ], na.rm = TRUE) - min(ref_data_fut[i,j, ], na.rm = TRUE)
        range_fut[idx, 1] <- min(ref_data_fut[i,j, ], na.rm = TRUE) - buffer * diff_fut
        range_fut[idx, 2] <- max(ref_data_fut[i,j, ], na.rm = TRUE) + buffer * diff_fut
      } else {
        diff_pres <- max_pres[i,j] - min_pres[i,j]
        range_pres[idx, 1] <- min_pres[i,j] - buffer * diff_pres
        range_pres[idx, 2] <- max_pres[i,j] + buffer * diff_pres
        diff_fut <- max_fut[i,j] - min_fut[i,j]
        range_fut[idx, 1] <- min_fut[i,j] - buffer * diff_fut
        range_fut[idx, 2] <- max_fut[i,j] + buffer * diff_fut
      }
    }
    # Reshape to [nlon, nlat, 2]
    range_pres_arr <- array(range_pres, dim = c(nlon, nlat, 2))
    range_fut_arr  <- array(range_fut, dim = c(nlon, nlat, 2))

    # Store the range for this variable.
    ref_range_present[,,v,] <- range_pres_arr
    ref_range_future[,,v,]  <- range_fut_arr

    # Collect the (transformed for pr) reference data across variables.
    if (v == 1) {
      ref_data_present_all <- array(NA, dim = c(nlon, nlat, length(iyear_pres), n_vars))
      ref_data_future_all  <- array(NA, dim = c(nlon, nlat, length(iyear_fut), n_vars))
    }
    ref_data_present_all[,,,v] <- ref_data_pres
    ref_data_future_all[,,,v]  <- ref_data_fut
  } # End loop over reference variables

  # Compute joint PDFs for the reference (gridpoint-wise) for both present and future.
  for (i in seq_len(nlon)) {
    for (j in seq_len(nlat)) {
      range_mat_pres <- matrix(NA, n_vars, 2)
      for (v in seq_len(n_vars))
        range_mat_pres[v, ] <- ref_range_present[i,j,v, ]
      pixel_data_pres <- sapply(1:n_vars, function(v) ref_data_present_all[i,j, , v])
      hist_pres <- compute_histND(pixel_data_pres, range_mat_pres, nbins)
      pdf_ref_present[i,j,] <- hist_pres / sum(hist_pres)

      range_mat_fut <- matrix(NA, n_vars, 2)
      for (v in seq_len(n_vars))
        range_mat_fut[v, ] <- ref_range_future[i,j,v, ]
      pixel_data_fut <- sapply(1:n_vars, function(v) ref_data_future_all[i,j, , v])
      hist_fut <- compute_histND(pixel_data_fut, range_mat_fut, nbins)
      pdf_ref_future[i,j,] <- hist_fut / sum(hist_fut)
    }
  }

  ## ----------------------------
  ## SEGMENT 2: Process each model (in parallel)
  ## ----------------------------

  num_models <- length(model_names)
  pdf_models_present <- array(NA, dim = c(nlon, nlat, nbins^n_vars, num_models))
  pdf_models_future  <- array(NA, dim = c(nlon, nlat, nbins^n_vars, num_models))

  plan(multisession, workers = workers)
  models_pdf_list <- future_lapply(seq_along(model_names), function(m_idx) {
    model_name <- model_names[m_idx]
    cat("Processing model:", model_name, "\n")

    # Lists to hold bias-corrected data (present and future) for each variable
    corrected_data_present_list <- vector("list", n_vars)
    corrected_data_future_list  <- vector("list", n_vars)

    for (v in seq_along(variables)) {
      var <- variables[v]
      mod_dir <- paste0(data_dir, model_name, '/', var, '/')
      mod_file <- list.files(path = mod_dir, pattern = glob2rx(paste0(var, "_", model_name, "*.nc")), full.names = TRUE)[1]
      if (is.na(mod_file)) stop("Model file for ", var, " not found for model ", model_name)

      nc_mod <- nc_open(mod_file)
      # --- Reordering longitude for model data ---
      lon_file <- ncvar_get(nc_mod, "lon")
      lat_file <- ncvar_get(nc_mod, "lat")
      if(any(lon_file >= 180)) {
        lon_file <- ifelse(lon_file >= 180, lon_file - 360, lon_file)
      }
      sorted_indices <- order(lon_file)
      lon_file_sorted <- lon_file[sorted_indices]
      lon_indices <- match(lon, lon_file_sorted)
      lat_indices <- match(lat, lat_file)
      if(any(is.na(lon_indices)) || any(is.na(lat_indices)))
        stop("Grid indices for model ", model_name, " not found.")

      yyyy <- extract_years_from_time(nc_mod)
      iyear_pres <- which(yyyy %in% year_present)
      if (length(iyear_pres)==0)
        stop("No present years for model ", model_name)
      iyear_fut <- which(yyyy %in% year_future)
      if (length(iyear_fut)==0)
        stop("No future years for model ", model_name)

      # Read full spatial subset and reorder along longitude
      full_data_pres <- ncvar_get(nc_mod, var,
                                  start = c(1, min(lat_indices), min(iyear_pres)),
                                  count = c(-1, length(lat_indices), length(iyear_pres)))
      full_data_pres <- full_data_pres[sorted_indices, , ]
      mod_data_pres <- full_data_pres[lon_indices, , ]

      full_data_fut <- ncvar_get(nc_mod, var,
                                 start = c(1, min(lat_indices), min(iyear_fut)),
                                 count = c(-1, length(lat_indices), length(iyear_fut)))
      full_data_fut <- full_data_fut[sorted_indices, , ]
      mod_data_fut <- full_data_fut[lon_indices, , ]

      nc_close(nc_mod)

      # Create arrays for corrected data (same dimensions as mod_data)
      corr_pres <- array(NA, dim = dim(mod_data_pres))
      corr_fut  <- array(NA, dim = dim(mod_data_fut))

      # Loop over each grid point and apply the bias correction.
      for (i in seq_len(nlon)) {
        for (j in seq_len(nlat)) {
          ts_mod_pres <- mod_data_pres[i,j, ]
          ts_mod_fut  <- mod_data_fut[i,j, ]
          if (var != "pr") {
            # Z-score bias correction: (x - mean_mod)/sd_mod scaled by reference sd and shifted by reference mean
            ref_mean_pres <- reference_stats_present[[v]]$mean[i,j]
            ref_sd_pres   <- reference_stats_present[[v]]$sd[i,j]
            mod_mean_pres <- mean(ts_mod_pres, na.rm = TRUE)
            mod_sd_pres   <- sd(ts_mod_pres, na.rm = TRUE)
            corr_pres[i,j,] <- ((ts_mod_pres - mod_mean_pres) / mod_sd_pres) * ref_sd_pres + ref_mean_pres

            ref_mean_fut <- reference_stats_future[[v]]$mean[i,j]
            ref_sd_fut   <- reference_stats_future[[v]]$sd[i,j]
            mod_mean_fut <- mean(ts_mod_fut, na.rm = TRUE)
            mod_sd_fut   <- sd(ts_mod_fut, na.rm = TRUE)
            corr_fut[i,j,] <- ((ts_mod_fut - mod_mean_fut) / mod_sd_fut) * ref_sd_fut + ref_mean_fut
          } else {
            # For precipitation: rescale so that the 90th percentile of raw data matches the reference raw Q90.
            ref_q90_pres <- reference_stats_present[[v]]$q90[i,j]
            mod_q90_pres <- as.numeric(quantile(ts_mod_pres, 0.90, na.rm = TRUE))
            corr_pres[i,j,] <- ts_mod_pres * (ref_q90_pres / mod_q90_pres)

            ref_q90_fut <- reference_stats_future[[v]]$q90[i,j]
            mod_q90_fut <- as.numeric(quantile(ts_mod_fut, 0.90, na.rm = TRUE))
            corr_fut[i,j,] <- ts_mod_fut * (ref_q90_fut / mod_q90_fut)

            # After bias correction, apply the log transform for precipitation.
            corr_pres[i,j,] <- log(corr_pres[i,j,] + 1)
            corr_fut[i,j,]  <- log(corr_fut[i,j,] + 1)
          }
        }
      }
      corrected_data_present_list[[v]] <- corr_pres
      corrected_data_future_list[[v]]  <- corr_fut
    } # End loop over variables for current model

    # Compute joint PDFs from the bias–corrected data for this model.
    model_pdf_pres <- array(NA, dim = c(nlon, nlat, nbins^n_vars))
    model_pdf_fut  <- array(NA, dim = c(nlon, nlat, nbins^n_vars))
    for (i in seq_len(nlon)) {
      for (j in seq_len(nlat)) {
        pixel_data_mod_pres <- sapply(1:n_vars, function(v) corrected_data_present_list[[v]][i,j, ])
        # Use the reference range (from ref_range_present) for binning consistency.
        range_mat <- matrix(NA, n_vars, 2)
        for (v in seq_len(n_vars))
          range_mat[v,] <- ref_range_present[i,j,v,]
        hist_mod_pres <- compute_histND(pixel_data_mod_pres, range_mat, nbins)
        model_pdf_pres[i,j,] <- hist_mod_pres / sum(hist_mod_pres)

        pixel_data_mod_fut <- sapply(1:n_vars, function(v) corrected_data_future_list[[v]][i,j, ])
        range_mat_fut <- matrix(NA, n_vars, 2)
        for (v in seq_len(n_vars))
          range_mat_fut[v,] <- ref_range_future[i,j,v,]
        hist_mod_fut <- compute_histND(pixel_data_mod_fut, range_mat_fut, nbins)
        model_pdf_fut[i,j,] <- hist_mod_fut / sum(hist_mod_fut)
      }
    }
    return(list(present = model_pdf_pres, future = model_pdf_fut))
  })  # End future_lapply over models
  plan(sequential)

  # Combine model results into overall arrays.
  for (m_idx in seq_along(model_names)) {
    pdf_models_present[,,,m_idx] <- models_pdf_list[[m_idx]]$present
    pdf_models_future[,,,m_idx]  <- models_pdf_list[[m_idx]]$future
  }

  return(list(
    pdf_reference = list(present = pdf_ref_present, future = pdf_ref_future),
    pdf_models    = list(present = pdf_models_present, future = pdf_models_future),
    reference_stats = list(present = reference_stats_present, future = reference_stats_future)
  ))
}
