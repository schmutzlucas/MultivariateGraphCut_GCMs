# Install and load necessary libraries
list_of_packages <- read.table("package_list.txt", sep="\n")$V1
new.packages <- list_of_packages[!(list_of_packages %in% installed.packages()[,"Package"])]
if(length(new.packages))
  install.packages(new.packages, repos = "https://cran.us.r-project.org")

library(devtools)
lapply(list_of_packages, library, character.only = TRUE)
install_github("thaos/gcoWrapR")

# Import libraries
#library(c(list_of_packages))


# Loading local functions

source_code_dir <- 'functions/' #The directory where all functions are saved.
file_paths <- list.files(source_code_dir, full.names = T)
for(path in file_paths){source(path)}

# Temporal ranges
year_present <- 1979:1998
year_future <- 1999:2019

# List of the variable used
# variables <- c('tas', 'tasmax',
#                'tasmin', 'pr',
#                'psl', 'hur',
#                'huss')
variables <- c('tas', 'pr', 'tas')


# TODO choose the method for the list of model names
# Selected models
selected_models <- c(
  2,34,16,9,8,28
)
# Obtains the list of models from the model names or from a file
# Method 1
# TODO This needs ajustements to remove prefixes and suffixes
tmp <- 'pr'
model_names <- list.files(paste0('data/CMIP5/', tmp))[selected_models]
model_names <- sub(paste0(tmp, '_'), '', model_names)
model_names <- sub('.nc', '', model_names)
rm(tmp)
reference_names <- c('era5')

