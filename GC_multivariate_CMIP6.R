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
year_present <- 1971:1990
year_future <- 1991:2010

# List of the variable used
# variables <- c('tas', 'tasmax',
#                'tasmin', 'pr',
#                'psl', 'hur',
#                'huss')
variables <- c('tas', 'pr', 'tasmax')

period <- c('historical')

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

# TODO implement bias_var from old code

# Graphcut labelling
GC_result <- list()
GC_result <- GraphCutOptimization(reference = reference_matrix_nrm$present,
                                  models_datacost = models_matrix_nrm$present,
                                  models_smoothcost = models_matrix_nrm$future,
                                  weight_data = 1,
                                  weight_smooth = 1,
                                  verbose = TRUE)
saveRDS(GC_result, file = 'GC_result.rds')

GC_labels <- GC_result$label_attribution

