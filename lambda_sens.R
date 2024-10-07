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

load('202407231828_my_workspace_ERA5_allModels_beforeOptim.RData')

lambda_values <- c(1, 5, 10, 25, 50, 100)

# Pre-computing H
h_dist <- array(0, c(length(lon), length(lat), length(model_names)))
m <- 1
for(model_name in model_names){
  for (i in seq_along(lon)){
    for (j in seq_along(lat)){
      h_dist[i, j, m] <- sqrt(sum((sqrt(kde_models[i, j, , m]) -
        sqrt(kde_ref[i, j, ]))^2)) / sqrt(2)

    }
  }
  m <- m + 1
}
remove(m)

# Graphcut hellinger labelling
GC_result_hellinger_new <- list()
n <- 1
for (lambda in lambda_values) {
  GC_result_hellinger_new[[n]] <-
    GraphCutHellinger2D_new2(pdf_ref = kde_ref,
                             pdf_models_future = kde_models_future,
                             h_dist = h_dist,
                             weight_data = lambda,
                             weight_smooth = 10,
                             nBins = nbins1d^2,
                             seed = 123,
                             verbose = TRUE,
                             rebuild = TRUE)
  n <- n + 1
}

GC_result_hellinger_new <- setNames(GC_result_hellinger_new, lambda_values)


# Get the current date and time

current_time <- Sys.time()

# Format the date and time as a string in the format 'yyyymmddhhmm'
formatted_time <- format(current_time, "%Y%m%d%H%M")

# Concatenate the formatted time string with your desired filename
filename <- paste0(formatted_time, "_my_workspace_ERA5_allModels_lambda_sens.RData")

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

# Loop through variables and models
GC_results_lambda_h_dist_future <- list()
GC_results_lambda_h_dist_mean_future <- list()
i <- 1
for (lambda in lambda_values) {
  tmp <- array(NA, dim = dim(GC_result_hellinger_new[[i]]$label_attribution))

  for(l in 0:(length(model_names))){
    islabel <- which(GC_result_hellinger_new[[i]]$label_attribution== l)
    tmp[islabel] <- h_dist_future[,,(l)][islabel]
  }
  GC_results_lambda_h_dist_future[[i]] <- tmp
  GC_results_lambda_h_dist_mean_future[[i]] <- mean(tmp)
  i <- i + 1
}



