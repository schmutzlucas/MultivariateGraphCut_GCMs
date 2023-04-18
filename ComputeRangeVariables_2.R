# Function to calculate ranges for a single variable
calculate_ranges <- function(variable, model_names, data_dir, period, year_interest, lon, lat) {
  range_var <- array(data = NA, dim = c(length(lon), length(lat), 2, length(model_names)))

  # Loop through models
  m <- 1
  for (model_name in model_names) {
    # Variable
    dir_path <- paste0(data_dir, model_name, '/', variable, '/')
    # Create the pattern
    pattern <- glob2rx(paste0(variable, "_", model_name, "_", period, "*.nc"))

    # Get the filepath
    file_name <- list.files(path = dir_path,
                            pattern = pattern)
    file_path <- paste0(dir_path, file_name)
    print(file_path)

    # Check that there is only one matching file
    nc_var <<- nc_open(file_path)

    # Extract and average data for present
    yyyy <- substr(as.character(nc.get.time.series(nc_var)), 1, 4)
    iyyyy <- which(yyyy %in% year_interest)

    # Get the entire 2D-time model as array
    tmp_grid_var <- ncvar_get(nc_var, variable, start = c(1, 1, min(iyyyy)), count = c(-1, -1, length(iyyyy)))

    # For each grid point...
    for (i in seq_along(lon)) {
      for (j in seq_along(lat)) {
        if (i %% 10 == 0 && j %% 100 == 0) {
          print(c(m, model_name, i, j))
        }
        tmp_var <- tmp_grid_var[i, j, ]

        if (variable == 'pr') {
          tmp_var <- log2((tmp_var * 86400) + 1)
        }

        range_var[i, j, 1, m] <- min(tmp_var)
        range_var[i, j, 2, m] <- max(tmp_var)
      }
    }

    # Close the file
    nc_close(nc_var)
    # Update counter for models
    m <- m + 1
  }

  return(range_var)
}

year_interest <- 1975:2014

# Calculate ranges for the first variable
range_var_1 <- calculate_ranges(variables[1], model_names, data_dir, period, year_interest, lon, lat)
# Calculate ranges for the second variable
range_var_2 <- calculate_ranges(variables[2], model_names, data_dir, period, year_interest, lon, lat)

saveRDS(range_var_1, 'ranges/pr_log_range_7Models_1975-2014.rds', compress = FALSE)
saveRDS(range_var_2, 'ranges/tas_range_7Models_1975-2014.rds', compress = FALSE)
range_var_final <- list()
range_var_final[[variables[1]]] <- array(data = NA, dim = c(length(lon), length(lat), 2))
range_var_final[[variables[2]]] <- array(data = NA, dim = c(length(lon), length(lat), 2))

# For each grid point...
for (i in seq_along(lon)) {
  for (j in seq_along(lat)) {
    range_var_final[[variables[1]]][i, j, 1] <- min(range_var_1[i, j, 1, ])
    range_var_final[[variables[1]]][i, j, 2] <- max(range_var_1[i, j, 2, ])

    range_var_final[[variables[2]]][i, j, 1] <- min(range_var_2[i, j, 1, ])
    range_var_final[[variables[2]]][i, j, 2] <- max(range_var_2[i, j, 2, ])
  }
}

saveRDS(range_var_final, 'ranges/range_var_final_7Models_1975-2014.rds', compress = FALSE)
