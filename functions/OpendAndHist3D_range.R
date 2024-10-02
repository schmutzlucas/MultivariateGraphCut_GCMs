OpenAndHist3D_range <- function(model_names, variables, year_interest, range_var) {
  # Initialize data structures for 3D histograms
  pdf_matrix <- array(0, c(length(lon), length(lat), nbins1d^3, length(model_names)))
  x_breaks <- array(0, c(length(lon), length(lat), nbins1d+1, length(model_names)))
  y_breaks <- array(0, c(length(lon), length(lat), nbins1d+1, length(model_names)))
  z_breaks <- array(0, c(length(lon), length(lat), nbins1d+1, length(model_names)))

  m <- 1
  for(model_name in model_names) {
    # Open data for three variables
    # Load variable 1
    file_path1 <- paste0(data_dir, model_name, '/', variables[1], '/', list.files(path = paste0(data_dir, model_name, '/', variables[1], '/'), pattern = glob2rx(paste0(variables[1], "_", model_name, "*.nc"))))
    nc1 <<- nc_open(file_path1)
    tmp_grid_1 <- ncvar_get(nc1, variables[1], start = c(1, 1, min(iyyyy1)), count = c(-1, -1, length(iyyyy1)))

    # Load variable 2
    file_path2 <- paste0(data_dir, model_name, '/', variables[2], '/', list.files(path = paste0(data_dir, model_name, '/', variables[2], '/'), pattern = glob2rx(paste0(variables[2], "_", model_name, "*.nc"))))
    nc2 <<- nc_open(file_path2)
    tmp_grid_2 <- ncvar_get(nc2, variables[2], start = c(1, 1, min(iyyyy2)), count = c(-1, -1, length(iyyyy2)))

    # Load variable 3
    file_path3 <- paste0(data_dir, model_name, '/', variables[3], '/', list.files(path = paste0(data_dir, model_name, '/', variables[3], '/'), pattern = glob2rx(paste0(variables[3], "_", model_name, "*.nc"))))
    nc3 <<- nc_open(file_path3)
    tmp_grid_3 <- ncvar_get(nc3, variables[3], start = c(1, 1, min(iyyyy3)), count = c(-1, -1, length(iyyyy3)))

    # Compute histograms at each grid point
    for (i in seq_along(lon)) {
      for (j in seq_along(lat)) {
        tmp1 <- tmp_grid_1[i, j, ]
        tmp2 <- tmp_grid_2[i, j, ]
        tmp3 <- tmp_grid_3[i, j, ]

        # Define fixed breaks for 3D histogram
        breaks1 <- seq(from = range_var[[variables[1]]][i, j, 1], to = range_var[[variables[1]]][i, j, 2], length.out = nbins1d + 1)
        breaks2 <- seq(from = range_var[[variables[2]]][i, j, 1], to = range_var[[variables[2]]][i, j, 2], length.out = nbins1d + 1)
        breaks3 <- seq(from = range_var[[variables[3]]][i, j, 1], to = range_var[[variables[3]]][i, j, 2], length.out = nbins1d + 1)

        # Compute 3D histogram
        hist_tmp <- get_histNd(cbind(tmp1, tmp2, tmp3), bw = c(bin_width1, bin_width2, bin_width3), bo = c(bin_origin1, bin_origin2, bin_origin3), limvar = limvar)
        pdf_matrix[i, j, , m] <- c(hist_tmp / sum(hist_tmp))  # Normalize
      }
    }
    nc_close(nc1)
    nc_close(nc2)
    nc_close(nc3)
    m <- m + 1
  }
  return(list(pdf_matrix, range_var, x_breaks, y_breaks, z_breaks))
}
