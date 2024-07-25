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

  GC_result_hellinger_new[[n]] <- GraphCutHellinger2D_new2(kde_ref = kde_ref,
                                                                kde_models_future = kde_models_future,
                                                                h_dist = h_dist,
                                                                weight_data = lambda,
                                                                weight_smooth = 10,
                                                                nBins = nbins1d^2,
                                                                verbose = TRUE,
                                                                rebuild = TRUE)
  n <- n + 1
}


