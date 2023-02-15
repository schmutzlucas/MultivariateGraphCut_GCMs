# Install and load necessary libraries
packages <- c('RcppXPtrUtils', 'devtools', 'Rcpp', 'ncdf4', 'ncdf4.helpers', 'abind', 'gcoWrapR')
new_packages <- packages[!(packages %in% installed.packages()[,"Package"])]
if(length(new_packages)) {
  install.packages(new_packages, repos = "https://cran.us.r-project.org")
}
lapply(packages, library, character.only = TRUE)

# Load local functions
source_code_dir <- 'functions/' # The directory where all functions are saved.
file_paths <- list.files(source_code_dir, full.names = TRUE)
lapply(file_paths, source)

# Temporal ranges
year_present <- 1979:1998
year_future <- 1999:2019

# List of variables used
variables <- c('tas', 'pr', 'tas')

# Selected models
selected_models <- c(28, 34, 16, 9, 8, 2)

# Obtain the list of model names from the files
model_names <- list.files(paste0('data/CMIP5/', variables[1]))[selected_models]
model_names <- sub(paste0(variables[1], '_'), '', model_names)
model_names <- sub('.nc', '', model_names)

reference_names <- c('era5')

# Open and average the models for the selected time periods
OpenAndAverageModels <- function(model_names, variables, year_present, year_future) {
  # TODO: Add path as argument
  models_list <- lapply(model_names, function(model) {
    lapply(variables, function(var) {
      opennc(paste0('data/CMIP5/', var, '/', var, '_', model, '.nc'))[, year_present:year_future, , ]
    })
  })
  models_matrix <- lapply(models_list, function(model) {
    abind(model, along = 4)
  })
  return(list(models_list = models_list, models_matrix = models_matrix))
}

tmp <- OpenAndAverageModels(model_names, variables, year_present, year_future)
models_list <- tmp$models_list
models_matrix <- tmp$models_matrix

tmp <- OpenAndAverageModels(reference_names, variables, year_present, year_future)
reference_list <- tmp$models_list
reference_matrix <- tmp$models_matrix

# Normalize variables
NormalizeVariables <- function(matrix_list, variables, method) {
  normalized_list <- lapply(matrix_list, function(matrix) {
    if (method == 'StdSc') {
      lapply(1:ncol(matrix), function(i) {
        matrix[, i, , ] <- (matrix[, i, , ] - mean(matrix[, i, , ])) / sd(matrix[, i, , ])
        matrix[, i, , ]
      })
    } else {
      matrix
    }
  })
  names(normalized_list) <- variables
  return(normalized_list)
}

models_matrix_nrm <- NormalizeVariables(models_matrix, variables, 'StdSc')
reference_matrix_nrm <- NormalizeVariables(reference_matrix, variables, 'StdSc')

# Graphcut labelling
GC_result <- GraphCutOptimization(
  reference = reference_matrix_nrm$present,
  models_datacost = models_matrix_nrm$present,
  models_smoothcost = models_matrix_nrm$future,
  weight_data = 1,
  weight_smooth = 1,
  verbose = TRUE
)
saveRDS(GC_result, file = 'GC_result.rds')

print('Test finished line 63')