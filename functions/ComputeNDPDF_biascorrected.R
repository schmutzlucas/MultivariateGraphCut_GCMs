compute_nd_pdf_bias_corrected <- function(variables, reference_name, model_names, data_dir,
                                            year_present, year_future, lon, lat, nbins, workers,
                                            buffer = 0.05) {
  # Load required libraries (assumes ncdf4, future, future.apply are already loaded)

  n_vars <- length(variables)
  nlon <- length(lon)
  nlat <- length(lat)

  ## ----------------------------
  ## SEGMENT 1: Process the reference data
  ## ----------------------------

  # Initialize arrays to store joint PDF for reference (present & future)
  pdf_ref_present <- array(NA, dim = c(nlon, nlat, nbins^n_vars))
  pdf_ref_future  <- array(NA, dim = c(nlon, nlat, nbins^n_vars))

  # We also store, for each variable, the reference statistics (mean, sd, min, max, and for pr, Q90)
  reference_stats_present <- vector("list", n_vars)
  reference_stats_future  <- vector("list", n_vars)

  # Also store the range to be used for PDF computation.
  # For each variable at each grid point, the range is [min - buffer*(max-min), max + buffer*(max-min)]
  ref_range_present <- array(NA, dim = c(nlon, nlat, n_vars, 2))
  ref_range_future  <- array(NA, dim = c(nlon, nlat, n_vars, 2))

  # We will also collect the raw reference data (for each variable, period) in a joint array for PDF computation.
  # For present and future separately, dimensions: [lon, lat, time, variable]
  ref_data_present_all <- NULL
  ref_data_future_all  <- NULL

  for (v in seq_along(variables)) {
    var <- variables[v]

    # Construct file path for the reference
    ref_dir <- paste0(data_dir, reference_name, '/', var, '/')
    ref_file <- list.files(path = ref_dir, pattern = glob2rx(paste0(var, "_", reference_name, "*.nc")), full.names = TRUE)[1]
    if (is.na(ref_file)) stop("Reference file for ", var, " not found.")

    nc_ref <- nc_open(ref_file)
    lon_var <- ncvar_get(nc_ref, "lon")
    lat_var <- ncvar_get(nc_ref, "lat")
    # (Adjust longitudes if needed – here we assume a matching grid)
    lon_indices <- match(lon, lon_var)
    lat_indices <- match(lat, lat_var)
    lon_indices <- lon_indices[!is.na(lon_indices)]
    lat_indices <- lat_indices[!is.na(lat_indices)]
    if(length(lon_indices) == 0 || length(lat_indices) == 0)
      stop("Grid indices for reference not found.")
    start_lon <- min(lon_indices)
    start_lat <- min(lat_indices)

    yyyy <- extract_years_from_time(nc_ref)

    # Extract present and future data slices
    iyear_pres <- which(yyyy %in% year_present)
    if (length(iyear_pres) == 0)
      stop("No present years found for reference ", reference_name)
    ref_data_pres <- ncvar_get(nc_ref, var,
                               start = c(start_lon, start_lat, min(iyear_pres)),
                               count = c(length(lon_indices), length(lat_indices), length(iyear_pres)))

    iyear_fut <- which(yyyy %in% year_future)
    if (length(iyear_fut) == 0)
      stop("No future years found for reference ", reference_name)
    ref_data_fut <- ncvar_get(nc_ref, var,
                              start = c(start_lon, start_lat, min(iyear_fut)),
                              count = c(length(lon_indices), length(lat_indices), length(iyear_fut)))

    nc_close(nc_ref)

    # For consistency with the original code, one might reorder data by sorted longitude;
    # here we assume lon and lat are already aligned.

    # Initialize matrices for stats
    m_pres <- matrix(NA, nlon, nlat)
    s_pres <- matrix(NA, nlon, nlat)
    min_pres <- matrix(NA, nlon, nlat)
    max_pres <- matrix(NA, nlon, nlat)
    if (var == "pr") q90_pres <- matrix(NA, nlon, nlat)

    m_fut <- matrix(NA, nlon, nlat)
    s_fut <- matrix(NA, nlon, nlat)
    min_fut <- matrix(NA, nlon, nlat)
    max_fut <- matrix(NA, nlon, nlat)
    if (var == "pr") q90_fut <- matrix(NA, nlon, nlat)

    # Loop over grid points (could be vectorized if needed)
    for (i in seq_len(nlon)) {
      for (j in seq_len(nlat)) {
        ts_pres <- ref_data_pres[i, j, ]
        ts_fut  <- ref_data_fut[i, j, ]
        m_pres[i,j] <- mean(ts_pres, na.rm = TRUE)
        s_pres[i,j] <- sd(ts_pres, na.rm = TRUE)
        min_pres[i,j] <- min(ts_pres, na.rm = TRUE)
        max_pres[i,j] <- max(ts_pres, na.rm = TRUE)
        if (var == "pr")
          q90_pres[i,j] <- as.numeric(quantile(ts_pres, 0.90, na.rm = TRUE))

        m_fut[i,j] <- mean(ts_fut, na.rm = TRUE)
        s_fut[i,j] <- sd(ts_fut, na.rm = TRUE)
        min_fut[i,j] <- min(ts_fut, na.rm = TRUE)
        max_fut[i,j] <- max(ts_fut, na.rm = TRUE)
        if (var == "pr")
          q90_fut[i,j] <- as.numeric(quantile(ts_fut, 0.90, na.rm = TRUE))
      }
    }

    reference_stats_present[[v]] <- list(mean = m_pres, sd = s_pres, min = min_pres, max = max_pres)
    if (var == "pr")
      reference_stats_present[[v]]$q90 <- q90_pres
    reference_stats_future[[v]] <- list(mean = m_fut, sd = s_fut, min = min_fut, max = max_fut)
    if (var == "pr")
      reference_stats_future[[v]]$q90 <- q90_fut

    # Compute buffered range for PDF binning
    range_pres <- matrix(NA, nrow = nlon * nlat, ncol = 2)
    range_fut  <- matrix(NA, nrow = nlon * nlat, ncol = 2)
    for (idx in 1:(nlon * nlat)) {
      i <- ((idx - 1) %% nlon) + 1
      j <- ((idx - 1) %/% nlon) + 1
      diff_pres <- max_pres[i,j] - min_pres[i,j]
      range_pres[idx,1] <- min_pres[i,j] - buffer * diff_pres
      range_pres[idx,2] <- max_pres[i,j] + buffer * diff_pres
      diff_fut <- max_fut[i,j] - min_fut[i,j]
      range_fut[idx,1] <- min_fut[i,j] - buffer * diff_fut
      range_fut[idx,2] <- max_fut[i,j] + buffer * diff_fut
    }
    # Reshape to [nlon, nlat, 2]
    range_pres_arr <- array(range_pres, dim = c(nlon, nlat, 2))
    range_fut_arr  <- array(range_fut, dim = c(nlon, nlat, 2))

    # Store each variable's range in the joint reference range arrays.
    ref_range_present[,,v,] <- range_pres_arr
    ref_range_future[,,v,]  <- range_fut_arr

    # Collect raw data across variables (assuming same time dimension length across variables)
    if (v == 1) {
      ref_data_present_all <- array(NA, dim = c(nlon, nlat, length(iyear_pres), n_vars))
      ref_data_future_all  <- array(NA, dim = c(nlon, nlat, length(iyear_fut), n_vars))
    }
    ref_data_present_all[,,,v] <- ref_data_pres
    ref_data_future_all[,,,v]  <- ref_data_fut
  } # end loop over reference variables

  # Now compute joint PDFs for the reference for both periods (gridpoint–wise).
  for (i in seq_len(nlon)) {
    for (j in seq_len(nlat)) {
      # For present, build a range matrix [n_vars x 2] from ref_range_present
      range_mat_pres <- matrix(NA, n_vars, 2)
      for (v in seq_len(n_vars))
        range_mat_pres[v,] <- ref_range_present[i,j,v,]
      pixel_data_pres <- sapply(1:n_vars, function(v) ref_data_present_all[i,j, , v])
      hist_pres <- compute_histND(pixel_data_pres, range_mat_pres, nbins)
      pdf_ref_present[i,j,] <- hist_pres / sum(hist_pres)

      # For future:
      range_mat_fut <- matrix(NA, n_vars, 2)
      for (v in seq_len(n_vars))
        range_mat_fut[v,] <- ref_range_future[i,j,v,]
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

    # For each variable, load the model data and apply bias correction
    corrected_data_present_list <- vector("list", n_vars)
    corrected_data_future_list  <- vector("list", n_vars)

    for (v in seq_along(variables)) {
      var <- variables[v]
      mod_dir <- paste0(data_dir, model_name, '/', var, '/')
      mod_file <- list.files(path = mod_dir, pattern = glob2rx(paste0(var, "_", model_name, "*.nc")), full.names = TRUE)[1]
      if (is.na(mod_file)) stop("Model file for ", var, " not found for model ", model_name)

      nc_mod <- nc_open(mod_file)
      lon_var <- ncvar_get(nc_mod, "lon")
      lat_var <- ncvar_get(nc_mod, "lat")
      lon_indices <- match(lon, lon_var)
      lat_indices <- match(lat, lat_var)
      lon_indices <- lon_indices[!is.na(lon_indices)]
      lat_indices <- lat_indices[!is.na(lat_indices)]
      if(length(lon_indices)==0 || length(lat_indices)==0)
        stop("Grid indices for model ", model_name, " not found.")
      start_lon <- min(lon_indices)
      start_lat <- min(lat_indices)

      yyyy <- extract_years_from_time(nc_mod)

      # Present data
      iyear_pres <- which(yyyy %in% year_present)
      if (length(iyear_pres)==0)
        stop("No present years for model ", model_name)
      mod_data_pres <- ncvar_get(nc_mod, var,
                                 start = c(start_lon, start_lat, min(iyear_pres)),
                                 count = c(length(lon_indices), length(lat_indices), length(iyear_pres)))
      # Future data
      iyear_fut <- which(yyyy %in% year_future)
      if (length(iyear_fut)==0)
        stop("No future years for model ", model_name)
      mod_data_fut <- ncvar_get(nc_mod, var,
                                start = c(start_lon, start_lat, min(iyear_fut)),
                                count = c(length(lon_indices), length(lat_indices), length(iyear_fut)))
      nc_close(nc_mod)

      # Create arrays to hold the bias-corrected values
      corr_pres <- array(NA, dim = dim(mod_data_pres))
      corr_fut  <- array(NA, dim = dim(mod_data_fut))

      # Loop over each grid point and perform the correction
      for (i in seq_len(nlon)) {
        for (j in seq_len(nlat)) {
          ts_mod_pres <- mod_data_pres[i,j, ]
          ts_mod_fut  <- mod_data_fut[i,j, ]

          # For non-precipitation: use z-score correction:
          if (var != "pr") {
            ref_mean_pres <- reference_stats_present[[v]]$mean[i,j]
            ref_sd_pres   <- reference_stats_present[[v]]$sd[i,j]
            mod_mean_pres <- mean(ts_mod_pres, na.rm = TRUE)
            mod_sd_pres   <- sd(ts_mod_pres, na.rm = TRUE)
            # Bias correction: (x - μ_mod)/σ_mod then scaled and shifted by reference stats.
            corr_pres[i,j,] <- ((ts_mod_pres - mod_mean_pres) / mod_sd_pres) * ref_sd_pres + ref_mean_pres

            ref_mean_fut <- reference_stats_future[[v]]$mean[i,j]
            ref_sd_fut   <- reference_stats_future[[v]]$sd[i,j]
            mod_mean_fut <- mean(ts_mod_fut, na.rm = TRUE)
            mod_sd_fut   <- sd(ts_mod_fut, na.rm = TRUE)
            corr_fut[i,j,] <- ((ts_mod_fut - mod_mean_fut) / mod_sd_fut) * ref_sd_fut + ref_mean_fut
          } else {
            # For precipitation: use Q90 scaling
            ref_q90_pres <- reference_stats_present[[v]]$q90[i,j]
            mod_q90_pres <- as.numeric(quantile(ts_mod_pres, 0.90, na.rm = TRUE))
            corr_pres[i,j,] <- ts_mod_pres * (ref_q90_pres / mod_q90_pres)

            ref_q90_fut <- reference_stats_future[[v]]$q90[i,j]
            mod_q90_fut <- as.numeric(quantile(ts_mod_fut, 0.90, na.rm = TRUE))
            corr_fut[i,j,] <- ts_mod_fut * (ref_q90_fut / mod_q90_fut)
          }
        }
      }
      corrected_data_present_list[[v]] <- corr_pres
      corrected_data_future_list[[v]]  <- corr_fut
    } # end loop over variables for this model

    # Now compute joint PDFs from the bias–corrected model data for this model.
    model_pdf_pres <- array(NA, dim = c(nlon, nlat, nbins^n_vars))
    model_pdf_fut  <- array(NA, dim = c(nlon, nlat, nbins^n_vars))
    for (i in seq_len(nlon)) {
      for (j in seq_len(nlat)) {
        # Build the matrix of corrected data for each variable (present)
        pixel_data_mod_pres <- sapply(1:n_vars, function(v) corrected_data_present_list[[v]][i,j, ])
        # Use the reference range (from ref_range_present) for consistency
        range_mat <- matrix(NA, n_vars, 2)
        for (v in seq_len(n_vars))
          range_mat[v,] <- ref_range_present[i,j,v,]
        hist_mod_pres <- compute_histND(pixel_data_mod_pres, range_mat, nbins)
        model_pdf_pres[i,j,] <- hist_mod_pres / sum(hist_mod_pres)

        # Similarly for future:
        pixel_data_mod_fut <- sapply(1:n_vars, function(v) corrected_data_future_list[[v]][i,j, ])
        range_mat_fut <- matrix(NA, n_vars, 2)
        for (v in seq_len(n_vars))
          range_mat_fut[v,] <- ref_range_future[i,j,v,]
        hist_mod_fut <- compute_histND(pixel_data_mod_fut, range_mat_fut, nbins)
        model_pdf_fut[i,j,] <- hist_mod_fut / sum(hist_mod_fut)
      }
    }
    return(list(present = model_pdf_pres, future = model_pdf_fut))
  })  # end future_lapply over models
  plan(sequential)

  # Combine the results into a single array for models.
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
