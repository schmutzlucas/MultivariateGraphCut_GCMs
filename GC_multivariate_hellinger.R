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
year_present <<- 1979:1994
year_future <<- 1999:2014
# data directory
data_dir <<- 'data/CMIP6/'

# Bins for the kde
nbins1d <<- 512

# Period
period <<- 'historical'



# List of the variable used
# variables <- c('tas', 'tasmax',
#                'tasmin', 'pr',
#                'psl', 'hur',
#                'huss')
variables <- c('pr', 'tas', 'tasmin', 'tasmax')
#variables <- c('pr', 'tas')

# Obtains the list of models from the model names or from a file
# Method 1
# TODO This needs ajustements to remove prefixes and suffixes
dir_path <- paste0('data/CMIP6/')
model_names <- list.dirs(dir_path, recursive = FALSE)
model_names <- basename(model_names)

# Initialize data structures
pdf_matrix <<- array(0, c(length(lon), length(lat), nbins1d,
                          length(model_names), length(variables)))
range_var <<- list()

tmp <- OpenAndKDE1D_new(
  model_names, variables, year_present, year_future, period
)

# # Get the current date and time
# current_time <- Sys.time()
#
# # Format the date and time as yyyymmddhhmm
# formatted_time <- format(current_time, "%Y%m%d%H%M")
#
# # Create the filename with the formatted timestamp and "GC_results" at the end
# filename <- paste0(formatted_time, "_pdf1d_cmip6.rds")
# saveRDS(tmp, filename)



pdf_matrix <- tmp[[1]]
kde_models <- pdf_matrix[ , , , -2, ]
kde_ref <- pdf_matrix[ , , , 2, ]

# Choose the reference in the models
reference_name <<- model_names[2]
model_names <<- model_names[-2]


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

# # Get the current date and time
# current_time <- Sys.time()
#
# # Format the date and time as a string in the format 'yyyymmddhhmm'
# formatted_time <- format(current_time, "%Y%m%d%H%M")
#
# # Concatenate the formatted time string with your desired filename
# filename <- paste0(formatted_time, "_my_workspace.RData")
#
# # Save the workspace using the generated filename
# save.image(file = filename)



# Computing the sum of hellinger distances between models and reference --> used as datacost
h_dist <- array(data = NA, dim = c(length(lon), length(lat),
                                   length(model_names), length(variables)))

# # Loop through variables and models
# v <- 1
# for(var in variables){
#   m <- 1
#   for(model_name in model_names){
#     for (i in seq_along(lon)){
#       for (j in seq_along(lat)){
#         h_dist[i, j, m, v] <- sqrt(sum((sqrt(kde_models[i, j, , m, v]) -
#           sqrt(kde_ref[i, j, , v]))^2)) / sqrt(2)
#
#       }
#     }
#     m <- m + 1
#   }
#   v <- v + 1
# }
# remove(v,m)

# Graphcut hellinger labelling
GC_result_hellinger <- list()
GC_result_hellinger <- GraphCutHellinger(kde_ref = kde_ref,
                               kde_models = kde_models,
                               models_smoothcost = models_matrix_nrm$future,
                               weight_data = 1,
                               weight_smooth = 1,
                               verbose = TRUE)


# Graphcut labelling
GC_result <- list()
GC_result <- GraphCutOptimization(reference = reference_matrix_nrm$present,
                                  models_datacost = models_matrix_nrm$present,
                                  models_smoothcost = models_matrix_nrm$future,
                                  weight_data = 1,
                                  weight_smooth = 1,
                                  verbose = TRUE)


# Get the current date and time
current_time <- Sys.time()

# Format the date and time as yyyymmddhhmm
formatted_time <- format(current_time, "%Y%m%d%H%M")

# Create the filename with the formatted timestamp and "GC_results" at the end
filename <- paste0(formatted_time, "_GC_results_hellinger.rds")

# Save the RDS file with the timestamped filename
#saveRDS(GC_result, file = filename)


# Get the current date and time
current_time <- Sys.time()

# Format the date and time as a string in the format 'yyyymmddhhmm'
formatted_time <- format(current_time, "%Y%m%d%H%M")

# Concatenate the formatted time string with your desired filename
filename <- paste0(formatted_time, "_my_workspace.RData")

# Save the workspace using the generated filename
save.image(file = filename, compression = FALSE)


#
#
# idx <- which.max(h_dist)
# subscripts <- arrayInd(idx, dim(h_dist))
#
# i <- subscripts[1, 1]
# j <- subscripts[1, 2]
# m <- subscripts[1, 3]
# v <- subscripts[1, 4]
#
# pdf_max <- kde_models[i, j, , m, v]
#
#


