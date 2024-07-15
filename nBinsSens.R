
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
lon <- c(-150, 0, 37)
lat <- c(-55, 3, 22)
# Temporal ranges
year_present <<- 1990:2021
year_future <<- 2000:2022
# data directory
data_dir <<- 'data/CMIP6_merged_all/'



# List of the variable used
variables <- c('pr', 'tas')

# Model names
model_names <- read.table('model_names_long.txt')
model_names <- as.list(model_names[['V1']])
# Index of the reference
ref_index <<- 1



tmp <- OpenAndHist2D_range(model_names, variables, year_present, range_var_final)


pdf_matrix <- tmp[[1]]
kde_models <- pdf_matrix[ , , , -ref_index]
kde_ref <- pdf_matrix[ , , , ref_index]
rm(pdf_matrix)
range_matrix <- tmp[[2]]
x_breaks <- tmp[[3]]
y_breaks <- tmp[[4]]





# Loop through variables and models
m <- 1
for (model_name in model_names) {
  for (i in seq_along(lon)) {
    for (j in seq_along(lat)) {
      # Compute Hellinger distance
      h_dist_unchecked <- sqrt(sum((sqrt(kde_models_future[i, j, , m]) - sqrt(kde_ref_future[i, j, ]))^2)) / sqrt(2)

      # Replace NaN with 0
      h_dist_future[i, j, m] <- replace(h_dist_unchecked, is.nan(h_dist_unchecked), 0)
    }
  }
  m <- m + 1
}

