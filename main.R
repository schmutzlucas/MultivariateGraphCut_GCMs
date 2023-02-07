# Install and load necessary libraries
{
  list.of.packages <- c('RcppXPtrUtils','devtools', 'Rcpp', 'ncdf4',
                        'ncdf4.helpers', 'abind')
  new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
  if(length(new.packages))
    install.packages(new.packages, repos = "https://cran.us.r-project.org")
  library(devtools)
  lapply(list.of.packages, library, character.only = TRUE)
  install_github("thaos/gcoWrapR")
}
{
  library(gcoWrapR)
  library(ncdf4)
  library(ncdf4.helpers)
  library(abind)
  library(roxygen2)
}

# Loading local functions
{
  source_code_dir <- 'functions/' #The directory where all functions are saved.
  file_paths <- list.files(source_code_dir, full.names = T)
  for(path in file_paths){source(path)}
}


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
                      28,34,16,9,8,2
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


# #Method2
# # TODO : test if it works
# model_names <- read.table("model_names.txt", sep="\n")$V1

# Open and average the models for the selected time periods
# TODO : Add path as argument
tmp <- OpenAndAverageModels(
  model_names, variables, year_present, year_future
  )
models_list <- tmp[[1]]
models_matrix <- tmp[[2]]
rm(tmp)

# Open and average the refeences for the selected time periods
tmp <- OpenAndAverageModels(
  reference_names, variables, year_present, year_future
  )
reference_list <- tmp[[1]]
reference_matrix <- tmp[[2]]

models_matrix_nrm <- list()
models_matrix_nrm <- NormalizeVariables(models_matrix, variables, 'StdSc')

reference_matrix_nrm <- list()
reference_matrix_nrm <- NormalizeVariables(reference_matrix, variables, 'StdSc')

# MinBias labelling
# TODO Not working!!!
# MinBias_label <- MinBiasOptimization(reference_list, models_list)

# TODO implement bias_var from old code


# Graphcut labelling
GC_result <- GraphCutOptimization(reference = reference_matrix_nrm$present,
                                  models_datacost = models_matrix_nrm$present,
                                  models_smoothcost = models_matrix_nrm$future,
                                  weight_data = 1,
                                  weight_smooth = 1)
