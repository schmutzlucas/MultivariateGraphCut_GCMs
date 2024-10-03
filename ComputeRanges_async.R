# Install and load necessary libraries
list_of_packages <- read.table("package_list.txt", sep = "\n")$V1
new.packages <- list_of_packages[!(list_of_packages %in% installed.packages()[, "Package"])]
if (length(new.packages)) install.packages(new.packages, repos = "https://cloud.r-project.org")

library(future)
library(future.apply)
library(ncdf4)
library(devtools)
lapply(list_of_packages, library, character.only = TRUE)
install_github("thaos/gcoWrapR")

# Loading local functions
source_code_dir <- 'functions/'
file_paths <- list.files(source_code_dir, full.names = TRUE)
for (path in file_paths) source(path)

# Global variables
lon <- 0:359
lat <- -70:70
year_interest <- 1950:2022
data_dir <- 'data/CMIP6_merged_all/'
variables <- c('pr', 'tas', 'psl')
model_names <- read.table('model_names_pr_tas_psl.txt')$V1

# Set up parallel backend using `future`
plan(multisession, workers = 3)

# Function to calculate ranges for each variable across all models
calculate_ranges <- function(variable, model_names, data_dir, year_interest, lon, lat) {
  range_var <- array(NA, dim = c(length(lon), length(lat), 2, length(model_names)))

  for (m in seq_along(model_names)) {
    model_name <- model_names[m]
    file_path <- paste0(data_dir, model_name, '/', variable, '/', list.files(path = paste0(data_dir, model_name, '/', variable, '/'), pattern = glob2rx(paste0(variable, "_", model_name, "*.nc")))[1])
    nc_var <- nc_open(file_path)
    yyyy <- substr(as.character(nc.get.time.series(nc_var)), 1, 4)
    iyyyy <- which(yyyy %in% year_interest)
    lon_var <- ncvar_get(nc_var, "lon")
    lat_var <- ncvar_get(nc_var, "lat")

    lon_indices <- which(lon_var %in% lon)
    lat_indices <- which(lat_var %in% lat)
    start_lon <- min(lon_indices)
    start_lat <- min(lat_indices)

    tmp_grid_var <- ncvar_get(nc_var, variable, start = c(start_lon, start_lat, min(iyyyy)), count = c(length(lon_indices), length(lat_indices), length(iyyyy)))
    if (variable == 'pr') tmp_grid_var <- log(tmp_grid_var + 1)

    range_var[ , , 1, m] <- apply(tmp_grid_var, c(1, 2), min, na.rm = TRUE)
    range_var[ , , 2, m] <- apply(tmp_grid_var, c(1, 2), max, na.rm = TRUE)

    nc_close(nc_var)
  }

  return(range_var)
}

# Calculate ranges for each variable asynchronously
range_results <- future_lapply(variables, function(var) calculate_ranges(var, model_names, data_dir, year_interest, lon, lat))

# Assign the results to respective variables
range_var_1 <- range_results[[1]]
range_var_2 <- range_results[[2]]
range_var_3 <- range_results[[3]]

# Save each variable's range separately
saveRDS(range_var_1, 'ranges/pr_log_range_AllModelsPar_1950-2022_70deg.rds')
saveRDS(range_var_2, 'ranges/tas_range_AllModelsPar_1950-2022_70deg.rds')
saveRDS(range_var_3, 'ranges/psl_range_AllModelsPar_1950-2022_70deg.rds')

# Final global range computation by merging results
range_var_final <- list()
for (var in variables) range_var_final[[var]] <- array(NA, dim = c(length(lon), length(lat), 2))

for (i in seq_along(lon)) {
  for (j in seq_along(lat)) {
    range_var_final[[variables[1]]][i, j, 1] <- min(range_var_1[i, j, 1, ], na.rm = TRUE)
    range_var_final[[variables[1]]][i, j, 2] <- max(range_var_1[i, j, 2, ], na.rm = TRUE)
    range_var_final[[variables[2]]][i, j, 1] <- min(range_var_2[i, j, 1, ], na.rm = TRUE)
    range_var_final[[variables[2]]][i, j, 2] <- max(range_var_2[i, j, 2, ], na.rm = TRUE)
    range_var_final[[variables[3]]][i, j, 1] <- min(range_var_3[i, j, 1, ], na.rm = TRUE)
    range_var_final[[variables[3]]][i, j, 2] <- max(range_var_3[i, j, 2, ], na.rm = TRUE)
  }
}

# Save the final merged ranges
saveRDS(range_var_final, 'ranges/range_var_final_allModelsPar_1950-2022_70deg.rds')
