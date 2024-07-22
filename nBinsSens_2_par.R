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
lapply(file_paths, source)

# Method 3
model_names <- as.list(read.table('model_names_long.txt')$V1)

# Loading the ranges
range_var_final <- readRDS('ranges/range_var_final_allModelsPar_1950-2022_new.rds')

# Temporal ranges
year_present <- 1990:2021

# Data directory
data_dir <- 'data/CMIP6_merged_all/'

# List of the variables used
variables <- c('pr', 'tas')

nbins1d_values <- 2^(5:12)

# Define the longitude indices with step size (1 to 360, should be within 1 to 360)
lon <- as.integer(seq(1, 360, length.out = 10))

# Define symmetric latitude indices around the equator (30 to 150, should be within 1 to 181)
lat <- as.integer(seq(30, 150, length.out = 10))

# Initialize list to store h_dist for each nbins1d
h_dist_list_long <- list()

for (nbins1d in nbins1d_values) {
  # Call the refactored function with the current nbins1d
  tmp <- OpenAndHist2D_range_lonlat_2_par(model_names, variables, year_present, range_var_final, lon, lat, nbins1d)

  pdf_matrix <- tmp[[1]]
  kde_models <- array(pdf_matrix[ , , , -1], dim = c(length(lon), length(lat), nbins1d^2, length(model_names)-1))
  kde_ref <- array(pdf_matrix[ , , , 1], dim = c(length(lon), length(lat), nbins1d^2))

  # Placeholder: Further processing to compute h_dist
  # Example calculation for h_dist: Mean squared error between model PDFs and reference PDF
  h_dist <- array(0, c(length(lon), length(lat), length(model_names)-1))
  for (m in 1:(length(model_names)-1)) {
    h_dist[ , , m] <- rowMeans((kde_models[ , , , m] - kde_ref)^2)
  }

  h_dist_list_long[[as.character(nbins1d)]] <- h_dist
}

# Save the h_dist_list_long to a file for later use
saveRDS(h_dist_list_long, file = 'h_dist_list_long.rds')

# Optionally, output summary statistics or visualizations
# Example: Mean h_dist across all models and grid points for each nbins1d
mean_h_dist <- sapply(h_dist_list_long, function(h_dist) mean(h_dist, na.rm = TRUE))
print(mean_h_dist)
