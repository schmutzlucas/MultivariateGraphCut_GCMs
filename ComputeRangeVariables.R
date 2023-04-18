range_var <- list()
range_var[[variables[1]]] <- array(data = NA, dim = c(length(lon), length(lat), 2, length(model_names)))
range_var[[variables[2]]] <- array(data = NA, dim = c(length(lon), length(lat), 2, length(model_names)))

range_var_final <- list()
range_var_final[[variables[1]]] <- array(data = NA, dim = c(length(lon), length(lat), 2))
range_var_final[[variables[2]]] <- array(data = NA, dim = c(length(lon), length(lat), 2))

year_interest <- 1975:2014

# Loop through variables and models
m <- 1
for(model_name in model_names){
  # Variable 1
  dir_path <- paste0(data_dir, model_name, '/', variables[1], '/')
  # Create the pattern
  pattern <- glob2rx(paste0(variables[1], "_", model_name, "_", period, "*.nc"))

  # Get the filepath
  file_name <- list.files(path = dir_path,
                          pattern = pattern)
  file_path <- paste0(dir_path, file_name)
  print(file_path)

  # Check that there is only one matching file
  nc1 <<- nc_open(file_path)

  # Extract and average data for present
  yyyy1 <- substr(as.character(nc.get.time.series(nc1)), 1, 4)
  iyyyy1 <- which(yyyy1 %in% year_interest)

  # Get the entire 2D-time model as array
  tmp_grid_1 <- ncvar_get(nc1, variables[1], start = c(1, 1, min(iyyyy1)), count = c(-1, -1, length(iyyyy1)))


  # Variable 2
  dir_path <- paste0(data_dir, model_name, '/', variables[2], '/')
  # Create the pattern
  pattern <- glob2rx(paste0(variables[2], "_", model_name, "_", period, "*.nc"))

  # Get the filepath
  file_name <- list.files(path = dir_path,
                          pattern = pattern)
  file_path <- paste0(dir_path, file_name)
  print(file_path)

  # Check that there is only one matching file
  nc2 <<- nc_open(file_path)

  # Extract and average data for present
  yyyy2 <- substr(as.character(nc.get.time.series(nc2)), 1, 4)
  iyyyy2 <- which(yyyy2 %in% year_interest)

  # Get the entire 2D-time model as array
  tmp_grid_2 <- ncvar_get(nc2, variables[2], start = c(1, 1, min(iyyyy2)), count = c(-1, -1, length(iyyyy2)))

  # For each grid point...
  for (i in seq_along(lon)) {
    for (j in seq_along(lat)) {
      if (i%%10 == 0 && j %% 100 == 0) {
        print(c(m, model_name, i, j))
      }
      tmp1 <- tmp_grid_1[i, j, ]
      tmp2 <- tmp_grid_2[i, j, ]

      if (variables[1] == 'pr') {
        tmp1 <- log2((tmp1 * 86400) + 1)
      }
      if (variables[2] == 'pr') {
        tmp2 <- log2((tmp2 * 86400) + 1)
      }

      if (variables[1] == 'pr') {
        range_var[[variables[1]]][i, j, 1, m] <- 0
        range_var[[variables[1]]][i, j, 2, m] <- range(tmp1)[2]

        range_var[[variables[2]]][i, j, 1, m] <- range(tmp2)[1]
        range_var[[variables[2]]][i, j, 2, m] <- range(tmp2)[2]
      }
      else if (variables[2] == 'pr') {
        range_var[[variables[2]]][i, j, 1] <- 0
        range_var[[variables[2]]][i, j, 2] <- range(tmp2)[2]

        range_var[[variables[1]]][i, j, 1, m] <- range(tmp1)[1]
        range_var[[variables[1]]][i, j, 2, m] <- range(tmp1)[2]
      }
      else{
        range_var[[variables[1]]][i, j, 1, m] <- range(tmp1)[1]
        range_var[[variables[1]]][i, j, 2, m] <- range(tmp1)[2]

        range_var[[variables[2]]][i, j, 1, m] <- range(tmp2)[1]
        range_var[[variables[2]]][i, j, 2, m] <- range(tmp2)[2]
      }
    }
  }

  # Close the file
  nc_close(nc1)
  nc_close(nc2)
  # Update counter for models
  m <- m + 1
}

# For each grid point...
for (i in seq_along(lon)) {
  for (j in seq_along(lat)) {
    range_var_final[[variables[1]]][i, j, 1] <- min(range_var[[variables[1]]][i,j,1,])
    range_var_final[[variables[1]]][i, j, 2] <- max(range_var[[variables[1]]][i,j,1,])

    range_var_final[[variables[2]]][i, j, 1] <- min(range_var[[variables[2]]][i,j,1,])
    range_var_final[[variables[2]]][i, j, 2] <- max(range_var[[variables[2]]][i,j,1,])
  }
}