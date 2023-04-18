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
# dir_path <- paste0('data/CMIP6/')
# model_names <- list.dirs(dir_path, recursive = FALSE)
# model_names <- basename(model_names)

# Method 2
model_names <- c('GFDL-ESM4',
                 'FGOALS-f3-L',
                 'MPI-ESM1-2-HR',
                 'ACCESS-CM2',
                 'ACCESS-ESM1-5',
                 'INM-CM5-0',
                 'MIROC6')


tmp <- OpenAndHist2D_range(model_names, variables, year_present , period, range_var_final)

pdf_matrix <- tmp[[1]]
kde_models <- pdf_matrix[ , , , -ref_index]
kde_ref <- pdf_matrix[ , , , ref_index]
range_matrix <- tmp[[2]]
x_breaks <- tmp[[3]]
y_breaks <- tmp[[4]]

# Choose the reference in the models
reference_name <<- model_names[ref_index]
model_names <<- model_names[-ref_index]


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
GC_result_hellinger <- GraphCutHellinger2D(kde_ref = kde_ref,
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

# Graphcut labelling
MinBiasHellinger <- list()
MinBiasHellinger <- GraphCutHellinger2D(kde_ref = kde_ref,
                                        kde_models = kde_models,
                                        models_smoothcost = models_matrix_nrm$future,
                                        weight_data = 1,
                                        weight_smooth = 0,
                                        verbose = TRUE)



MMM <- list()
MMM$tas <- apply(models_matrix$future[,,,2], c(1, 2), mean)

MMM$pr <- apply(models_matrix$future[,,,1]*86400, c(1, 2), mean)

GC_hellinger_projections <- list()
j <- 1
for(var in variables){
  for(l in 0:(length(model_names))){
    islabel <- which(GC_result_hellinger$label_attribution == l)
    GC_hellinger_projections[[var]][islabel] <- models_matrix$future[,,(l),j][islabel]
  }
  j <- j + 1
}
GC_hellinger_projections$tas <- matrix(GC_hellinger_projections$tas, nrow = 360)
GC_hellinger_projections$pr <- matrix(GC_hellinger_projections$pr * 86400, nrow = 360)




# # Get the current date and time
# current_time <- Sys.time()
#
# # Format the date and time as yyyymmddhhmm
# formatted_time <- format(current_time, "%Y%m%d%H%M")
#
# # Create the filename with the formatted timestamp and "GC_results" at the end
# filename <- paste0(formatted_time, "_GC_results_hellinger.rds", compress = FALSE)
#
# # Save the RDS file with the timestamped filename
# #saveRDS(GC_result, file = filename)


# Get the current date and time
current_time <- Sys.time()

# Format the date and time as a string in the format 'yyyymmddhhmm'
formatted_time <- format(current_time, "%Y%m%d%H%M")

# Concatenate the formatted time string with your desired filename
filename <- paste0(formatted_time, "_my_workspace.RData")

# Save the workspace using the generated filename
save.image(file = filename, compress = FALSE)


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