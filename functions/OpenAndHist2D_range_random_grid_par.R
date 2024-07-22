library(ncdf4)
library(parallel)

OpenAndHist2D_range_lonlat_2_par <- function(model_names, variables, year_interest, range_var, lon_vec, lat_vec, nbins1d) {
  # Initialize data structures
  pdf_matrix <- array(0, c(length(lon_vec), length(lat_vec), nbins1d^2, length(model_names)))
  x_breaks <- array(0, c(length(lon_vec), length(lat_vec), nbins1d+1, length(model_names)))
  y_breaks <- array(0, c(length(lon_vec), length(lat_vec), nbins1d+1, length(model_names)))

  process_model <- function(model_name) {
    result <- list(pdf_matrix = array(0, c(length(lon_vec), length(lat_vec), nbins1d^2)),
                   x_breaks = array(0, c(length(lon_vec), length(lat_vec), nbins1d+1)),
                   y_breaks = array(0, c(length(lon_vec), length(lat_vec), nbins1d+1)))

    for (var_idx in 1:2) {
      dir_path <- paste0(data_dir, model_name, '/', variables[var_idx], '/')
      pattern <- glob2rx(paste0(variables[var_idx], "_", model_name, "*.nc"))
      file_name <- list.files(path = dir_path, pattern = pattern, full.names = TRUE)
      stopifnot(length(file_name) == 1)
      nc <- nc_open(file_name)
      if (var_idx == 1) nc1 <- nc else nc2 <- nc
    }

    yyyy2 <- substr(as.character(nc.get.time.series(nc2)), 1, 4)
    iyyyy2 <- which(yyyy2 %in% year_interest)
    lon_dim1 <- length(ncvar_get(nc1, "lon"))
    lat_dim1 <- length(ncvar_get(nc1, "lat"))
    lon_dim2 <- length(ncvar_get(nc2, "lon"))
    lat_dim2 <- length(ncvar_get(nc2, "lat"))

    for (i in seq_along(lon_vec)) {
      for (j in seq_along(lat_vec)) {
        if (lon_vec[i] > lon_dim1 || lat_vec[j] > lat_dim1 || lon_vec[i] < 1 || lat_vec[j] < 1) {
          stop("Error: Index exceeds dimension bound for lon/lat")
        }

        tmp1 <- ncvar_get(nc1, variables[1], start = c(lon_vec[i], lat_vec[j], min(iyyyy2)), count = c(1, 1, length(iyyyy2)))
        tmp2 <- ncvar_get(nc2, variables[2], start = c(lon_vec[i], lat_vec[j], min(iyyyy2)), count = c(1, 1, length(iyyyy2)))

        h2d <- hist2(x = tmp1, y = tmp2, nbins = nbins1d, xlim = range_var[[1]], ylim = range_var[[2]])
        result$pdf_matrix[i, j, ] <- as.vector(h2d$counts)
        result$x_breaks[i, j, ] <- h2d$breaks[[1]]
        result$y_breaks[i, j, ] <- h2d$breaks[[2]]
      }
    }
    nc_close(nc1)
    nc_close(nc2)

    return(result)
  }

  cl <- makeCluster(detectCores() - 1)
  clusterExport(cl, list("data_dir", "variables", "year_interest", "range_var", "lon_vec", "lat_vec", "nbins1d", "nc.get.time.series", "hist2"))
  results <- parLapply(cl, model_names, function(model_name) process_model(model_name))
  stopCluster(cl)

  for (m in seq_along(model_names)) {
    pdf_matrix[,,,m] <- results[[m]]$pdf_matrix
    x_breaks[,,,m] <- results[[m]]$x_breaks
    y_breaks[,,,m] <- results[[m]]$y_breaks
  }

  return(list(pdf_matrix, x_breaks, y_breaks))
}
