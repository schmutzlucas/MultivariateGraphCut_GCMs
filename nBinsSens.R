
# Install and load necessary libraries
list_of_packages <- read.table("package_list.txt", sep="\n")$V1
new.packages <- list_of_packages[!(list_of_packages %in% installed.packages()[,"Package"])]
if(length(new.packages))
  install.packages(new.packages, repos = "https://cloud.r-project.org")

library(devtools)
lapply(list_of_packages, library, character.only = TRUE)
install_github("thaos/gcoWrapR")


# Loading local functions
source_code_dir <- 'functions/' #The directory where all functions are saved.
file_paths <- list.files(source_code_dir, full.names = T)
for(path in file_paths){source(path)}


# Method 3
model_names <- read.table('model_names_long.txt')
model_names <- as.list(model_names[['V1']])


# Loading the ranges
range_var_final <- readRDS('ranges/range_var_final_allModelsPar_1950-2022_new.rds')


# Setting global variables
lon <- c(150)
lat <- c(55)
# Temporal ranges
year_present <<- 1990:2021

# data directory
data_dir <<- 'data/CMIP6_merged_all/'


# List of the variable used
variables <- c('pr', 'tas')

nbins1d <<- c(3)


# Model names
model_names <- read.table('model_names_long.txt')
model_names <- as.list(model_names[['V1']])
# Index of the reference
ref_index <<- 1



tmp <- OpenAndHist2D_range_lonlat(model_names, variables, year_present, range_var_final, lon, lat)

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
reference_name <<- model_names[ref_index]
model_names <<- model_names[-ref_index]


h_dist <- array(NaN, c(length(lon), length(lat),  length(model_names)))

# Loop through variables and models
m <- 1
for (model_name in model_names) {
  for (i in 1:length(lon)) {

      # Compute Hellinger distance
      h_dist_unchecked <- sqrt(sum((sqrt(kde_models[i, i, , m]) - sqrt(kde_ref[i, i, ]))^2)) / sqrt(2)

      # Replace NaN with 0
      h_dist[i, i, m] <- replace(h_dist_unchecked, is.nan(h_dist_unchecked), 0)
    }
  m <- m + 1
}
