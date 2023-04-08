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



{
  readRDS(kde1d_cmip6.rds)


  kde_models <- tmp[[1]][ , , , -1, ]
  kde_ref <- tmp[[1]][ , , , 1, ]
}



# Obtains the list of models from the model names or from a file
# Method 1
# TODO This needs ajustements to remove prefixes and suffixes
dir_path <- paste0('data/CMIP6/')
model_names <- list.dirs(dir_path, recursive = FALSE)
model_names <- basename(model_names)

# Choose the reference in the models
reference_name <- model_names[1]
model_names <- model_names[-1]


# Open and average the models for the selected time periods
tmp <- OpenAndAverageCMIP6(
  model_names, variables, year_present, year_future, period
)
models_list <- tmp[[1]]
models_matrix <- tmp[[2]]
rm(tmp)

tmp <- OpenAndAverageCMIP6(
  reference_name, variables, year_present, year_future, period
)
reference_list <- tmp[[1]]
reference_matrix <- tmp[[2]]
rm(tmp)

models_matrix_nrm <- list()
models_matrix_nrm <- NormalizeVariables(models_matrix, variables, 'StdSc')

reference_matrix_nrm <- list()
reference_matrix_nrm <- NormalizeVariables(reference_matrix, variables, 'StdSc')

# MinBias labelling
MinBias_labels <- MinBiasOptimization(reference_matrix_nrm$present,
                                      models_matrix_nrm$present)



# Graphcut labelling
GC_result <- list()
GC_result <- GraphCutHellinger(kde_ref = kde_ref,
                                  kde_models = kde_models,
                                  models_smoothcost = models_matrix_nrm$future,
                                  weight_data = 1,
                                  weight_smooth = 1,
                                  verbose = TRUE)

saveRDS(GC_result, file = 'GC_result.rds')




