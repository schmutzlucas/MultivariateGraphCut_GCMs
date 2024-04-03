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

range_var_final <- readRDS('ranges/range_var_final_allModelsPar_1950-2022.rds')

# Setting global variables
lon <<- 0:359
lat <<- -90:90
# Temporal ranges
year_present <<- 1977:1999
year_future <<- 2000:2022
# data directory
data_dir <<- 'data/CMIP6_merged_all/'

# Bins for the kde
nbins1d <<- 32



# List of the variable used
# variables <- c('tas', 'tasmax',
#                'tasmin', 'pr',
#                'psl', 'hur',
#                'huss')
# variables <- c('pr', 'tas', 'tasmin', 'tasmax')
variables <- c('pr', 'tas')

# Obtains the list of models from the model names or from a file
# # Method 1
# dir_path <- paste0('data/CMIP6_merged_all/')
# model_names <- list.dirs(dir_path, recursive = FALSE)
# model_names <- basename(model_names)

# Method 2
# model_names <- c(  'ERA5',
#                                   'MIROC6',
#                                   'IPSL-CM6A-LR',
#                                   'NorESM2-MM',
#                                   'UKESM1-0-LL')

# Method 3
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


tmp <- OpenAndHist2D_range(model_names, variables, year_future , range_var_final)

pdf_matrix <- tmp[[1]]
kde_models_future <- pdf_matrix[ , , , -ref_index]
kde_ref_future <- pdf_matrix[ , , , ref_index]
range_matrix_future <- tmp[[2]]
x_breaks_future <- tmp[[3]]
y_breaks_future <- tmp[[4]]
rm(pdf_matrix)
rm(tmp)

# Choose the reference in the models
reference_name <<- model_names[ref_index]
model_names <<- model_names[-ref_index]


# Open and average the models for the selected time periods
tmp <- OpenAndAverageCMIP6(
  model_names, variables, year_present, year_future
)
models_list <- tmp[[1]]
models_matrix <- tmp[[2]]
rm(tmp)


tmp <- OpenAndAverageCMIP6(
  reference_name, variables, year_present, year_future
)
reference_list <- tmp[[1]]
reference_matrix <- tmp[[2]]
rm(tmp)

models_matrix_nrm <- list()
models_matrix_nrm <- NormalizeVariables(models_matrix, variables, 'StdSc')

reference_matrix_nrm <- list()
reference_matrix_nrm <- NormalizeVariables(reference_matrix, variables, 'StdSc')


# Get the current date and time

current_time <- Sys.time()

# Format the date and time as a string in the format 'yyyymmddhhmm'
formatted_time <- format(current_time, "%Y%m%d%H%M")

# Concatenate the formatted time string with your desired filename
filename <- paste0(formatted_time, "_my_workspace_ERA5_allModels_beforeOptim.RData")

# Save the workspace using the generated filename
save.image(file = filename, compress = FALSE)



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
                                           weight_smooth = 0.5,
                                           verbose = TRUE)

# Graphcut hellinger labelling
GC_result_hellinger_new <- list()
GC_result_hellinger_new <- GraphCutHellinger2D_new2(kde_ref = kde_ref,
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

# Format the date and time as a string in the format 'yyyymmddhhmm'
formatted_time <- format(current_time, "%Y%m%d%H%M")

# Concatenate the formatted time string with your desired filename
filename <- paste0(formatted_time, "_my_workspace_ERA5_allmodels.RData")

# Save the workspace using the generated filename
save.image(file = filename, compress = FALSE)



# Graphcut labelling
MinBiasHellinger <- list()
MinBiasHellinger <- GraphCutHellinger2D(kde_ref = kde_ref,
                                        kde_models = kde_models,
                                        models_smoothcost = models_matrix_nrm$future,
                                        weight_data = 1,
                                        weight_smooth = 0,
                                        verbose = TRUE)



# Get the current date and time

current_time <- Sys.time()

# Format the date and time as a string in the format 'yyyymmddhhmm'
formatted_time <- format(current_time, "%Y%m%d%H%M")

# Concatenate the formatted time string with your desired filename
filename <- paste0(formatted_time, "_my_workspace_ERA5_allmodelsKDE_present-future.RData")

# Save the workspace using the generated filename
save.image(file = filename, compress = FALSE)



# Computing the sum of hellinger distances between models and reference --> used as datacost
h_dist_future <- array(data = 0, dim = c(length(lon), length(lat),
                                         length(model_names)))

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

n1 <- dim(kde_models)[1]
n2 <- dim(kde_models)[2]
n3 <- dim(kde_models)[3]
n4 <- dim(kde_models)[4]

MMM_KDE_future <- array(0, dim = c(n1, n2, n3))

for (i in 1:n1) {
  for (j in 1:n2) {
    for (k in 1:n3) {
      MMM_KDE_future[i,j,k] <- mean(kde_models_future[i,j,k,])
    }
  }
}



# Computing the sum of hellinger distances between models and reference --> used as datacost
MMM_h_dist_future <- array(data = 0, dim = dim(MinBias_labels))

# Loop through grip-points
  for (i in 1:360) {
    for (j in 1:181) {
      # Compute Hellinger distance
      h_dist_unchecked <- sqrt(sum((sqrt(MMM_KDE_future[i,j, ]) - sqrt(kde_ref[i, j, ]))^2)) / sqrt(2)

      # Replace NaN with 0
      MMM_h_dist_future[i, j] <- replace(h_dist_unchecked, is.nan(h_dist_unchecked), 0)
    }
  }



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



# Get the current date and time
current_time <- Sys.time()

# Format the date and time as a string in the format 'yyyymmddhhmm'
formatted_time <- format(current_time, "%Y%m%d%H%M")

# Concatenate the formatted time string with your desired filename
filename <- paste0(formatted_time, "_my_workspace_ERA5_final_allModels.RData")

# Save the workspace using the generated filename
save.image(file = filename, compress = FALSE)

