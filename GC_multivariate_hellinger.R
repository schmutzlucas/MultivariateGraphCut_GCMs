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

# Setting global variables
lon <<- 0:359
lat <<- -90:90
# Temporal ranges
year_present <<- 1979:1998
year_future <<- 1999:2019
# data directory
data_dir <<- 'data/CMIP6/'

# Bins for the kde
nbins1d <<- 512

# Period
period <<- 'historical'

range_var <<- list()

# List of the variable used
# variables <- c('tas', 'tasmax',
#                'tasmin', 'pr',
#                'psl', 'hur',
#                'huss')
variables <- c('pr', 'tas', 'tasmax')


# Obtains the list of models from the model names or from a file
# Method 1
# TODO This needs ajustements to remove prefixes and suffixes
dir_path <- paste0('data/CMIP6/')
model_names <- list.dirs(dir_path, recursive = FALSE)
model_names <- basename(model_names)

tmp <- OpenAndKDE1D_new(
  model_names, variables, year_present, year_future, period
)

saveRDS(tmp, file = 'kde1d_cmip6.rds')

kde_models <- tmp$kde_matrix[ , , , -1, ]
kde_ref <- tmp$kde_matrix[ , , , 1, ]

h_dist <- array(data = NA, dim = c(length(lon), length(lat),
                                   length(model_names), length(variables)))

# Loop through variables and models
v <- 1
for(var in variables){
  m <- 1
  for(model_name in model_names){
    for (i in seq_along(lon)){
      for (j in seq_along(lat)){
        h_dist[i, j, m, v] <- sqrt(sum((sqrt(kde_models[i, j, , m, v]) -
                                    sqrt(kde_ref[i, j, , 1, v]))^2)) / sqrt(2)

      }
    }
  }
}





