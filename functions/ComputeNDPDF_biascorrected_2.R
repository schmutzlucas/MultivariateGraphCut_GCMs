# Set up parallel processing option (64 GiB maximum globals)
options(future.globals.maxSize = 64 * 1024^3L)  # 64 GiB

#-------------------------------
# Helper function for processing a single model.
# This function is defined at the top level to reduce the size of globals captured.
#-------------------------------
process_model_pdf <- function(m_idx, model_names, nlon, nlat, variables, data_dir,
                              year_present, year_future, reference_stats_present,
                              ref_range_present, ref_range_future, nbins) {
  model_name <- as.character(model_names[[m_idx]])
  cat("Processing model:", model_name, "\n")

  n_vars <- length(variables)
  corrected_data_present_list <- vector("list", n_vars)
  corrected_data_future_list  <- vector("list", n_vars)
  raw_means_pres_list <- vector("list", n_vars)
  raw_means_fut_list  <- vector("list", n_vars)
  out_range_pres <- array(0, dim = c(nlon, nlat, n_vars))
  out_range_fut  <- array(0, dim = c(nlon, nlat, n_vars))

  # Read grid info from the first variable
  var0 <- variables[1]
  mod_dir0 <- paste0(data_dir, model_name, '/', var0, '/')
  mod_file0 <- list.files(path = mod_dir0, pattern = glob2rx(paste0(var0, "_", model_name, "*.nc")), full.names = TRUE)[1]
  if (is.na(mod_file0)) stop("Model file for ", var0, " not found for model ", model_name)
  nc_mod0 <- nc_open(mod_file0)

  lon_file <- ncvar_get(nc_mod0, "lon")
  lat_file <- ncvar_get(nc_mod0, "lat")

  # Normalize longitude to [-180, 180] if needed
  lon_file_adjusted <- ifelse(lon_file >= 180, lon_file - 360, lon_file)
  lon_user_adjusted <- ifelse(lon >= 180, lon - 360, lon)


  # Sort file longitudes and keep original indices
  lon_order <- order(lon_file_adjusted)
  lon_file_sorted <- lon_file_adjusted[lon_order]
  lon_file_original_sorted <- lon_file[lon_order]  # For NetCDF indexing

  # Match user-specified longitudes and latitudes to file grid
  lon_idx_unsorted <- match(lon_user_adjusted, lon_file_sorted)
  if (any(is.na(lon_idx_unsorted))) stop("Some user-specified longitudes not found in NetCDF.")
  lon_idx_in_file <- lon_order[lon_idx_unsorted]  # indices in original NetCDF file

  lat_idx_in_file <- match(lat, lat_file)
  if (any(is.na(lat_idx_in_file))) stop("Some user-specified latitudes not found in NetCDF.")

  nc_close(nc_mod0)
  gc()

  # Loop over variables
  for (v in seq_along(variables)) {
    var <- variables[v]
    mod_dir <- paste0(data_dir, model_name, '/', var, '/')
    mod_file <- list.files(path = mod_dir, pattern = glob2rx(paste0(var, "_", model_name, "*.nc")), full.names = TRUE)[1]
    if (is.na(mod_file)) stop("Model file for ", var, " not found for model ", model_name)

    nc_mod <- nc_open(mod_file)
    yyyy <- extract_years_from_time(nc_mod)
    iyear_pres <- which(yyyy %in% year_present)
    if (length(iyear_pres) == 0) stop("No present years for model ", model_name)
    iyear_fut <- which(yyyy %in% year_future)
    if (length(iyear_fut) == 0) stop("No future years for model ", model_name)

    # Compute NetCDF read ranges
    start_lon <- min(lon_idx_in_file)
    count_lon <- max(lon_idx_in_file) - start_lon + 1
    start_lat <- min(lat_idx_in_file)
    count_lat <- max(lat_idx_in_file) - start_lat + 1

    # ---- Present data ----
    start_time_pres <- min(iyear_pres)
    count_time_pres <- max(iyear_pres) - start_time_pres + 1

    full_data_pres <- ncvar_get(nc_mod, var,
                                start = c(start_lon, start_lat, start_time_pres),
                                count = c(count_lon, count_lat, count_time_pres))

    lon_file_subset <- lon_file[start_lon:(start_lon + count_lon - 1)]
    lat_file_subset <- lat_file[start_lat:(start_lat + count_lat - 1)]
    lon_idx_local <- match(lon_file_original_sorted[lon_idx_unsorted], lon_file_subset)
    lat_idx_local <- match(lat, lat_file_subset)

    mod_data_pres <- full_data_pres[lon_idx_local, lat_idx_local, ]

    # ---- Future data ----
    start_time_fut <- min(iyear_fut)
    count_time_fut <- max(iyear_fut) - start_time_fut + 1

    full_data_fut <- ncvar_get(nc_mod, var,
                               start = c(start_lon, start_lat, start_time_fut),
                               count = c(count_lon, count_lat, count_time_fut))

    mod_data_fut <- full_data_fut[lon_idx_local, lat_idx_local, ]

    nc_close(nc_mod)
    gc()


    # Compute raw means (without correction) for each pixel for present and future.
    # These are matrices of size [nlon, nlat] for the current variable.
    raw_means_pres_list[[v]] <- apply(mod_data_pres, c(1,2), mean, na.rm = TRUE)
    raw_means_fut_list[[v]]  <- apply(mod_data_fut, c(1,2), mean, na.rm = TRUE)

    # Pre-allocate corrected arrays
    corr_pres <- array(NA, dim = dim(mod_data_pres))
    corr_fut  <- array(NA, dim = dim(mod_data_fut))

    for (i in seq_len(nlon)) {
      for (j in seq_len(nlat)) {
        ts_mod_pres <- mod_data_pres[i, j, ]
        ts_mod_fut  <- mod_data_fut[i, j, ]
        if (var != "pr") {
          ref_mean_pres <- reference_stats_present[[v]]$mean[i, j]
          ref_sd_pres   <- reference_stats_present[[v]]$sd[i, j]

          mod_mean_pres <- mean(ts_mod_pres, na.rm = TRUE)
          mod_sd_pres   <- sd(ts_mod_pres, na.rm = TRUE)
          corr_pres[i, j, ] <- ((ts_mod_pres - mod_mean_pres) / mod_sd_pres) * ref_sd_pres + ref_mean_pres

          mod_mean_fut <- mean(ts_mod_fut, na.rm = TRUE)
          mod_sd_fut   <- sd(ts_mod_fut, na.rm = TRUE)
          corr_fut[i, j, ] <- ((ts_mod_fut - mod_mean_pres) / mod_sd_pres) * ref_sd_pres + ref_mean_pres

        } else {
          ref_q90_pres <- reference_stats_present[[v]]$q90[i, j]
          mod_q90_pres <- as.numeric(quantile(ts_mod_pres, 0.90, na.rm = TRUE))
          corr_pres[i, j, ] <- ts_mod_pres * (ref_q90_pres / mod_q90_pres)

          mod_q90_fut <- as.numeric(quantile(ts_mod_fut, 0.90, na.rm = TRUE))
          corr_fut[i, j, ] <- ts_mod_fut * (ref_q90_pres / mod_q90_pres)

          # Force negative values to 0 before applying the log transform
          corr_pres[i, j, ][corr_pres[i, j, ] < 0] <- 0
          corr_fut[i, j, ][corr_fut[i, j, ] < 0] <- 0

          corr_pres[i, j, ] <- log(corr_pres[i, j, ] + 1)
          corr_fut[i, j, ]  <- log(corr_fut[i, j, ] + 1)
        }
      }
    }
    corrected_data_present_list[[v]] <- corr_pres
    corrected_data_future_list[[v]]  <- corr_fut
    rm(corr_pres, corr_fut)
    gc()
  } # End loop over variables

  # Compute PDFs for bias-corrected data (as before)
  model_pdf_pres <- array(NA, dim = c(nlon, nlat, nbins^n_vars))
  model_pdf_fut  <- array(NA, dim = c(nlon, nlat, nbins^n_vars))
  for (i in seq_len(nlon)) {
    for (j in seq_len(nlat)) {
      # --- Compute PDF for the present period ---
      pixel_data_mod_pres <- sapply(1:n_vars, function(v) corrected_data_present_list[[v]][i, j, ])
      range_mat_pres <- matrix(NA, n_vars, 2)
      for (v in seq_len(n_vars)) {
        range_mat_pres[v, ] <- ref_range_present[i, j, v, ]
      }
      count_vec_pres <- sapply(1:n_vars, function(v)
        sum(pixel_data_mod_pres[, v] < range_mat_pres[v, 1] | pixel_data_mod_pres[, v] > range_mat_pres[v, 2]))
      out_range_pres[i, j, ] <- count_vec_pres
      hist_mod_pres <- compute_histND(pixel_data_mod_pres, range_mat_pres, nbins)
      model_pdf_pres[i, j, ] <- hist_mod_pres / sum(hist_mod_pres)

      # --- Compute PDF for the future period ---
      # Use present reference range for bias-corrected future data.
      pixel_data_mod_fut <- sapply(1:n_vars, function(v) corrected_data_future_list[[v]][i, j, ])
      range_mat_fut <- matrix(NA, n_vars, 2)
      for (v in seq_len(n_vars)) {
        range_mat_fut[v, ] <- ref_range_future[i, j, v, ]
      }
      count_vec_fut <- sapply(1:n_vars, function(v)
        sum(pixel_data_mod_fut[, v] < range_mat_fut[v, 1] | pixel_data_mod_fut[, v] > range_mat_fut[v, 2]))
      out_range_fut[i, j, ] <- count_vec_fut
      hist_mod_fut <- compute_histND(pixel_data_mod_fut, range_mat_fut, nbins)
      model_pdf_fut[i, j, ] <- hist_mod_fut / sum(hist_mod_fut)
    }
  }

  rm(mod_data_pres, mod_data_fut)
  gc()

  # Assemble raw means arrays (uncorrected) for the model.
  model_means_pres <- array(NA, dim = c(nlon, nlat, n_vars))
  model_means_fut  <- array(NA, dim = c(nlon, nlat, n_vars))
  for (v in seq_len(n_vars)) {
    model_means_pres[,,v] <- raw_means_pres_list[[v]]
    model_means_fut[,,v]  <- raw_means_fut_list[[v]]
  }

  return(list(present = model_pdf_pres, future = model_pdf_fut,
              out_range = list(present = out_range_pres, future = out_range_fut),
              raw_means = list(present = model_means_pres, future = model_means_fut)))
}



#-------------------------------
# Main function
#-------------------------------
compute_nd_pdf_bias_corrected_2 <- function(variables, reference_name, model_names, data_dir,
                                            year_present, year_future, lon, lat, nbins, workers,
                                            buffer, verbose) {

  n_vars <- length(variables)
  nlon <- length(lon)
  nlat <- length(lat)

  pdf_ref_present <- array(NA, dim = c(nlon, nlat, nbins^n_vars))
  pdf_ref_future  <- array(NA, dim = c(nlon, nlat, nbins^n_vars))

  reference_stats_present <- vector("list", n_vars)
  reference_stats_future  <- vector("list", n_vars)

  ref_range_present <- array(NA, dim = c(nlon, nlat, n_vars, 2))
  ref_range_future  <- array(NA, dim = c(nlon, nlat, n_vars, 2))

  ref_data_present_all <- NULL
  ref_data_future_all  <- NULL

  for (v in seq_along(variables)) {
    var <- variables[v]
    ref_dir <- paste0(data_dir, reference_name, '/', var, '/')
    ref_file <- list.files(path = ref_dir, pattern = glob2rx(paste0(var, "_", reference_name, "*.nc")), full.names = TRUE)[1]
    if (is.na(ref_file)) stop("Reference file for ", var, " not found.")

    nc_ref <- nc_open(ref_file)

    lon_file <- ncvar_get(nc_ref, "lon")
    lat_file <- ncvar_get(nc_ref, "lat")

    # Normalize longitudes to [-180, 180] range for both file and user inputs
    lon_file_adjusted <- ifelse(lon_file >= 180, lon_file - 360, lon_file)
    lon_user_adjusted <- ifelse(lon >= 180, lon - 360, lon)

    lon_order <- order(lon_file_adjusted)
    lon_file_sorted <- lon_file_adjusted[lon_order]
    lon_file_original_sorted <- lon_file[lon_order]

    lon_idx_unsorted <- match(lon_user_adjusted, lon_file_sorted)
    if (any(is.na(lon_idx_unsorted))) stop("Some user-specified longitudes not found in reference NetCDF.")
    lon_idx_in_file <- lon_order[lon_idx_unsorted]

    lat_idx_in_file <- match(lat, lat_file)
    if (any(is.na(lat_idx_in_file))) stop("Some user-specified latitudes not found in reference NetCDF.")

    yyyy <- extract_years_from_time(nc_ref)
    iyear_pres <- which(yyyy %in% year_present)
    if (length(iyear_pres)==0) stop("No present years found for reference ", reference_name)
    iyear_fut <- which(yyyy %in% year_future)
    if (length(iyear_fut)==0) stop("No future years found for reference ", reference_name)

    start_lon <- min(lon_idx_in_file)
    count_lon <- max(lon_idx_in_file) - start_lon + 1
    start_lat <- min(lat_idx_in_file)
    count_lat <- max(lat_idx_in_file) - start_lat + 1

    start_time_pres <- min(iyear_pres)
    count_time_pres <- max(iyear_pres) - start_time_pres + 1

    full_data_pres <- ncvar_get(nc_ref, var,
                                start = c(start_lon, start_lat, start_time_pres),
                                count = c(count_lon, count_lat, count_time_pres))

    lon_file_subset <- lon_file[start_lon:(start_lon + count_lon - 1)]
    lat_file_subset <- lat_file[start_lat:(start_lat + count_lat - 1)]
    lon_idx_local <- match(lon_file_original_sorted[lon_idx_unsorted], lon_file_subset)
    lat_idx_local <- match(lat, lat_file_subset)

    ref_data_pres <- full_data_pres[lon_idx_local, lat_idx_local, , drop = FALSE]
    rm(full_data_pres)

    start_time_fut <- min(iyear_fut)
    count_time_fut <- max(iyear_fut) - start_time_fut + 1
    full_data_fut <- ncvar_get(nc_ref, var,
                               start = c(start_lon, start_lat, start_time_fut),
                               count = c(count_lon, count_lat, count_time_fut))
    ref_data_fut <- full_data_fut[lon_idx_local, lat_idx_local, , drop = FALSE]
    rm(full_data_fut)

    nc_close(nc_ref)
    gc()

    m_pres   <- matrix(NA, nlon, nlat)
    s_pres   <- matrix(NA, nlon, nlat)
    min_pres <- matrix(NA, nlon, nlat)
    max_pres <- matrix(NA, nlon, nlat)
    if (var == "pr") q90_pres <- matrix(NA, nlon, nlat)

    m_fut   <- matrix(NA, nlon, nlat)
    s_fut   <- matrix(NA, nlon, nlat)
    min_fut <- matrix(NA, nlon, nlat)
    max_fut <- matrix(NA, nlon, nlat)
    if (var == "pr") q90_fut <- matrix(NA, nlon, nlat)

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

    reference_stats_present[[v]] <- list(mean = m_pres, sd = s_pres, min = min_pres, max = max_pres)
    if (var == "pr")
      reference_stats_present[[v]]$q90 <- q90_pres
    reference_stats_future[[v]] <- list(mean = m_fut, sd = s_fut, min = min_fut, max = max_fut)
    if (var == "pr")
      reference_stats_future[[v]]$q90 <- q90_fut

    if (var == "pr") {
      ref_data_pres <- log(ref_data_pres + 1)
      ref_data_fut  <- log(ref_data_fut + 1)
    }

    range_pres <- matrix(NA, nrow = nlon * nlat, ncol = 2)
    range_fut  <- matrix(NA, nrow = nlon * nlat, ncol = 2)
    for (idx in 1:(nlon * nlat)) {
      i <- ((idx - 1) %% nlon) + 1
      j <- ((idx - 1) %/% nlon) + 1
      if (var == "pr") {
        diff_pres <- max(ref_data_pres[i,j, ], na.rm = TRUE)
        range_pres[idx, 1] <- 0
        range_pres[idx, 2] <- max(ref_data_pres[i,j, ], na.rm = TRUE) + buffer * diff_pres

        diff_fut <- max(ref_data_fut[i,j, ], na.rm = TRUE)
        range_fut[idx, 1] <- 0
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
    range_pres_arr <- array(range_pres, dim = c(nlon, nlat, 2))
    range_fut_arr  <- array(range_fut, dim = c(nlon, nlat, 2))

    ref_range_present[,,v,] <- range_pres_arr
    ref_range_future[,,v,]  <- range_fut_arr

    if (v == 1) {
      ref_data_present_all <- array(NA, dim = c(nlon, nlat, length(iyear_pres), n_vars))
      ref_data_future_all  <- array(NA, dim = c(nlon, nlat, length(iyear_fut), n_vars))
    }
    ref_data_present_all[,,,v] <- ref_data_pres
    ref_data_future_all[,,,v]  <- ref_data_fut
  } # End loop over reference variables

  for (i in seq_len(nlon)) {
    if (verbose && (i %% 5 == 0)) {
      cat(sprintf("[%s] Processing reference grid row %d of %d\n",
                  format(Sys.time(), "%Y-%m-%d %H:%M:%S"), i, nlon))
    }
    for (j in seq_len(nlat)) {
      range_mat_pres <- matrix(NA, n_vars, 2)
      for (v in seq_len(n_vars)){
        range_mat_pres[v, ] <- ref_range_present[i, j, v, ]
      }
      pixel_data_pres <- sapply(1:n_vars, function(v) ref_data_present_all[i, j, , v])
      hist_pres <- compute_histND(pixel_data_pres, range_mat_pres, nbins)

      pdf_ref_present[i, j, ] <- hist_pres / sum(hist_pres)

      range_mat_fut <- matrix(NA, n_vars, 2)
      for (v in seq_len(n_vars)){
        range_mat_fut[v, ] <- ref_range_future[i, j, v, ]
      }
      pixel_data_fut <- sapply(1:n_vars, function(v) ref_data_future_all[i, j, , v])
      hist_fut <- compute_histND(pixel_data_fut, range_mat_fut, nbins)

      pdf_ref_future[i, j, ] <- hist_fut / sum(hist_fut)
    }
  }

  rm(ref_data_present_all, ref_data_future_all)
  gc()

  ## ----------------------------
  ## SEGMENT 2: Process each model (in parallel)
  ## ----------------------------

  num_models <- length(model_names)
  pdf_models_present <- array(NA, dim = c(nlon, nlat, nbins^n_vars, num_models))
  pdf_models_future  <- array(NA, dim = c(nlon, nlat, nbins^n_vars, num_models))
  # New: Allocate arrays for raw (uncorrected) model means.
  means_models_present <- array(NA, dim = c(nlon, nlat, n_vars, num_models))
  means_models_future  <- array(NA, dim = c(nlon, nlat, n_vars, num_models))

  # Arrays to record out-of-range counts [nlon, nlat, n_vars, num_models]
  out_range_pres_all <- array(0, dim = c(nlon, nlat, n_vars, num_models))
  out_range_fut_all  <- array(0, dim = c(nlon, nlat, n_vars, num_models))

  plan(multisession, workers = workers)
  library(future)
  globals <- future::globalsOf(
    function() future_lapply(1, function(x) NULL),
    substitute = FALSE
  )
  print(names(globals))
  sapply(globals, object.size)
  models_pdf_list <- future_lapply(seq_along(model_names),
                                   FUN = function(m_idx) process_model_pdf(m_idx, model_names, nlon, nlat, variables, data_dir,
                                                                           year_present, year_future, reference_stats_present,
                                                                           ref_range_present, ref_range_future, nbins))
  plan(sequential)

  for (m_idx in seq_along(model_names)) {
    pdf_models_present[,,,m_idx] <- models_pdf_list[[m_idx]]$present
    pdf_models_future[,,,m_idx]  <- models_pdf_list[[m_idx]]$future
    means_models_present[,,,m_idx] <- models_pdf_list[[m_idx]]$raw_means$present
    means_models_future[,,,m_idx]  <- models_pdf_list[[m_idx]]$raw_means$future
  }

  for (m_idx in seq_along(model_names)) {
    out_range_pres_all[,,,m_idx] <- models_pdf_list[[m_idx]]$out_range$present
    out_range_fut_all[,,,m_idx]  <- models_pdf_list[[m_idx]]$out_range$future
  }

  return(list(
    pdf_reference = list(present = pdf_ref_present, future = pdf_ref_future),
    pdf_models    = list(present = pdf_models_present, future = pdf_models_future),
    raw_means_models = list(present = means_models_present, future = means_models_future),
    reference_stats = list(present = reference_stats_present, future = reference_stats_future),
    out_of_range_counts = list(present = out_range_pres_all, future = out_range_fut_all),
    ref_range_present = ref_range_present,
    ref_range_future = ref_range_future
  ))
}

