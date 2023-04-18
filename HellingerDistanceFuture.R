# Install and load necessary libraries
list_of_packages <- read.table("package_list.txt", sep="\n")$V1
new.packages <- list_of_packages[!(list_of_packages %in% installed.packages()[,"Package"])]
if(length(new.packages))
  install.packages(new.packages, repos = "https://cran.us.r-project.org")

library(devtools)
lapply(list_of_packages, library, character.only = TRUE)
install_github("thaos/gcoWrapR")


# Loading local functions

source_code_dir <- 'functions/' #The directory where all functions are saved.
file_paths <- list.files(source_code_dir, full.names = T)
for(path in file_paths){source(path)}


# Setting global variables
lon <<- 0:359
lat <<- -90:90
# Temporal ranges
year_present <<- 1985:1999
year_future <<- 2000:2014
# data directory
data_dir <<- 'data/CMIP6/'

# Bins for the kde
nbins1d <<- 32

# Period
period <<- 'historical'

# List of the variable used
# variables <- c('tas', 'tasmax',
#                'tasmin', 'pr',
#                'psl', 'hur',
#                'huss')
# variables <- c('pr', 'tas', 'tasmin', 'tasmax')
variables <- c('pr', 'tas')

# Index of the reference
ref_index <<- 1
# Obtains the list of models from the model names or from a file
# Method 1
# # TODO This needs ajustements to remove prefixes and suffixes
# dir_path <- paste0('data/CMIP6/')
# model_names <- list.dirs(dir_path, recursive = FALSE)
# model_names <- basename(model_names)
model_names <- c('GFDL-ESM4',
                 'FGOALS-f3-L',
                 'MPI-ESM1-2-HR',
                 'ACCESS-CM2',
                 'ACCESS-ESM1-5',
                 'INM-CM5-0',
                 'MIROC6')


tmp <- OpenAndHist2D_range(model_names, variables, year_future , period, range_var_final)

pdf_matrix <- tmp[[1]]
kde_models <- pdf_matrix[ , , , -ref_index]
kde_ref <- pdf_matrix[ , , , ref_index]
range_matrix <- tmp[[2]]
x_breaks <- tmp[[3]]
y_breaks <- tmp[[4]]

# Choose the reference in the models
reference_name <<- model_names[ref_index]
model_names <<- model_names[-ref_index]


# Computing the sum of hellinger distances between models and reference --> used as datacost
h_dist_future <- array(data = 0, dim = c(length(lon), length(lat),
                                         length(model_names)))

# Loop through variables and models
m <- 1
for (model_name in model_names) {
  for (i in seq_along(lon)) {
    for (j in seq_along(lat)) {
      # Compute Hellinger distance
      h_dist_unchecked <- sqrt(sum((sqrt(kde_models[i, j, , m]) - sqrt(kde_ref[i, j, ]))^2)) / sqrt(2)

      # Replace NaN with 0
      h_dist_future[i, j, m] <- replace(h_dist_unchecked, is.nan(h_dist_unchecked), 0)
    }
  }
  m <- m + 1
}

# Get the current date and time
current_time <- Sys.time()

# Format the date and time as a string in the format 'yyyymmddhhmm'
formatted_time <- format(current_time, "%Y%m%d%H%M")

# Concatenate the formatted time string with your desired filename
filename <- paste0(formatted_time, "_my_workspace_with_hdist_future.RData")

# Save the workspace using the generated filename
save.image(file = filename, compress = FALSE)

