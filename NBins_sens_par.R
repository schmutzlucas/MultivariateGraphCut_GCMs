# Install and load necessary libraries
list_of_packages <- read.table("package_list.txt", sep="\n")$V1
new.packages <- list_of_packages[!(list_of_packages %in% installed.packages()[,"Package"])]
if(length(new.packages))
  install.packages(new.packages, repos = "https://cloud.r-project.org")

library(devtools)
lapply(list_of_packages, library, character.only = TRUE)
install_github("thaos/gcoWrapR")

# Loading local functions
source_code_dir <- 'functions/' # The directory where all functions are saved.
file_paths <- list.files(source_code_dir, full.names = TRUE)
for(path in file_paths){source(path)}

# Method 3
model_names <- read.table('model_names_long.txt')
model_names <- as.list(model_names[['V1']])

# Loading the ranges
range_var_final <- readRDS('ranges/range_var_final_allModelsPar_1950-2022_new.rds')

# Temporal ranges
year_present <- 1990:2021

# Data directory
data_dir <- 'data/CMIP6_merged_all/'

# List of the variable used
variables <- c('pr', 'tas')

nbins1d_values <- 2^(5:12)

# Define the longitude indices with step size (1 to 360, should be within 1 to 360)
lon <- seq(1, 360, length.out = 10)

# Define symmetric latitude indices around the equator (30 to 150, should be within 1 to 181)
lat <- seq(30, 150, length.out = 10)

# Ensure they are integers
lon <- as.integer(lon)
lat <- as.integer(lat)

# Initialize list to store h_dist for each nbins1d
h_dist_list_long <- list()

# Load necessary libraries for parallel processing
library(foreach)
library(doParallel)
library(ncdf4)
library(ncdf4.helpers)

# Set up parallel backend
num_cores <- 8
cl <- makeCluster(num_cores)
registerDoParallel(cl)

# Export necessary objects and functions to cluster
clusterExport(cl, c("model_names", "variables", "year_present", "range_var_final", "lon", "lat", "OpenAndHist2D_range_lonlat_2"))
clusterEvalQ(cl, {
  library(ncdf4)
  library(ncdf4.helpers)
  library(gcoWrapR)
  library(squash)
})

# Parallel loop using foreach
h_dist_list_long <- foreach(nbins1d = nbins1d_values, .combine = 'c', .packages = c('gcoWrapR', 'foreach', 'doParallel', 'ncdf4', 'ncdf4.helpers'), .export = c('model_names', 'variables', 'year_present', 'range_var_final', 'lon', 'lat')) %dopar% {
  # Explicitly run garbage collection
  gc()

  # Model names
  model_names <- read.table('model_names_long.txt')
  model_names <- as.list(model_names[['V1']])
  # Index of the reference
  ref_index <- 1

  # Call the function with the current nbins1d
  tmp <- OpenAndHist2D_range_lonlat_2(model_names, variables, year_present, range_var_final, lon, lat, nbins1d)

  pdf_matrix <- tmp[[1]]
  kde_models <- pdf_matrix[ , , , -ref_index]
  kde_models <- array(kde_models, dim = c(length(lon), length(lat), nbins1d^2, length(model_names)-1))
  kde_ref <- pdf_matrix[ , , , ref_index]
  kde_ref <- array(kde_ref, dim = c(length(lon), length(lat), nbins1d^2))
  rm(pdf_matrix)
  range_matrix <- tmp[[2]]
  x_breaks <- tmp[[3]]
  y_breaks <- tmp[[4]]

  # Choose the reference in the models
  reference_name <- model_names[ref_index]
  model_names <- model_names[-ref_index]

  h_dist <- array(NaN, c(length(lon), length(lat),  length(model_names)))

  # Loop through variables and models
  m <- 1
  for (model_name in model_names) {
    for (i in 1:length(lon)) {
      for (j in 1:length(lat)){
        # Compute Hellinger distance
        h_dist_unchecked <- sqrt(sum((sqrt(kde_models[i, j, , m]) - sqrt(kde_ref[i, j, ]))^2)) / sqrt(2)

        # Replace NaN with 0
        h_dist[i, j, m] <- replace(h_dist_unchecked, is.nan(h_dist_unchecked), 0)
      }
    }
    m <- m + 1
  }

  # Return the h_dist for the current nbins1d as a named list element
  list(setNames(list(h_dist), as.character(nbins1d)))
}

# Stop parallel backend
stopCluster(cl)

# Combine the list of lists into a single list
h_dist_list_combined <- list()
for (item in h_dist_list_long) {
  h_dist_list_combined <- c(h_dist_list_combined, item)
}

# Get the current date and time
current_time <- Sys.time()

# Format the date and time as a string in the format 'yyyymmddhhmm'
formatted_time <- format(current_time, "%Y%m%d%H%M")

# Concatenate the formatted time string with your desired filename
filename <- paste0(formatted_time, "h_dist_list.rds")

# Save the results to plot later
saveRDS(h_dist_list_combined, file = filename)

h_dist_values <- h_dist_list_combined

# Initialize plot_data as an empty data frame
plot_data <- data.frame(
  nbins1d = integer(),
  h_dist = numeric()
)

# Loop through the nbins1d_values and prepare the data for plotting
for (nbins in nbins1d_values) {
  h_dist <- h_dist_list_combined[[as.character(nbins)]]
  h_dist_values <- as.vector(h_dist[!is.nan(h_dist)])

  # Check if h_dist_values is not empty before appending to plot_data
  if (length(h_dist_values) > 0) {
    plot_data <- rbind(plot_data, data.frame(nbins1d = nbins, h_dist = h_dist_values))
  }
}

# Create the box plots
ggplot(plot_data, aes(x = factor(nbins1d), y = h_dist)) +
  geom_boxplot() +
  labs(title = "Hellinger Distance Distribution for Different nbins1d",
       x = "Number of Bins (nbins1d)",
       y = "Hellinger Distance") +
  theme_minimal()

# Get the current date and time
current_time <- Sys.time()

# Format the date and time as a string in the format 'yyyymmddhhmm'
formatted_time <- format(current_time, "%Y%m%d%H%M")

# Concatenate the formatted time string with your desired filename
filename <- paste0(formatted_time, "_my_workspace_NBins_sens.RData")

# Save the workspace using the generated filename
save.image(file = filename, compress = FALSE)
