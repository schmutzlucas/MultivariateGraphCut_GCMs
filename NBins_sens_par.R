# Load necessary libraries
library(parallel)
library(devtools)
library(ggplot2)  # For plotting

# Install and load necessary packages
list_of_packages <- read.table("package_list.txt", sep="\n")$V1
new.packages <- list_of_packages[!(list_of_packages %in% installed.packages()[,"Package"])]
if(length(new.packages) > 0) {
  install.packages(new.packages, repos = "https://cloud.r-project.org")
}
lapply(list_of_packages, library, character.only = TRUE)

# Install GitHub package
devtools::install_github("thaos/gcoWrapR")

# Load local functions
source_code_dir <- 'functions/'
file_paths <- list.files(source_code_dir, full.names = TRUE)
lapply(file_paths, source)

# Initialize parallel cluster with 2 cores
cl <- makeCluster(2)
clusterExport(cl, "file_paths")
clusterEvalQ(cl, {
    lapply(list_of_packages, require, character.only = TRUE)
    file_paths <- list.files('functions/', full.names = TRUE)
    lapply(file_paths, source)
})

# Prepare data
range_var_final <- readRDS('ranges/range_var_final_allModelsPar_1950-2022_new.rds')
model_names <- read.table('model_names_long.txt')$V1
variables <- c('pr', 'tas')
year_present <- 1990:2021
lon <- as.integer(seq(1, 360, length.out = 10))
lat <- as.integer(seq(30, 150, length.out = 10))
nbins1d_values <- 2^(5:6)

# Export necessary objects to cluster nodes
clusterExport(cl, list("model_names", "variables", "year_present", "range_var_final", "lon", "lat"))

# Reference model index
ref_index <- 1

# Parallel computation for each nbins1d
h_dist_results <- parLapply(cl, nbins1d_values, function(nbins1d) {
    tmp <- OpenAndHist2D_range_lonlat_2(model_names, variables, year_present, range_var_final, lon, lat, nbins1d)
    pdf_matrix <- tmp[[1]]

    # Extract KDEs: Reference and models
    kde_ref <- pdf_matrix[, , , ref_index]
    kde_models <- pdf_matrix[, , , -ref_index]

    h_dist <- array(NaN, c(length(lon), length(lat), length(model_names)-1))
    for (m in 1:(length(model_names)-1)) {
        for (i in 1:length(lon)) {
            for (j in 1:length(lat)) {
                h_dist_unchecked <- sqrt(sum((sqrt(kde_models[i, j, , m]) - sqrt(kde_ref[i, j, ]))^2)) / sqrt(2)
                h_dist[i, j, m] <- replace(h_dist_unchecked, is.nan(h_dist_unchecked), 0)
            }
        }
    }
    mean_h_dist <- mean(h_dist, na.rm = TRUE)
    return(list(nbins1d = nbins1d, mean_h_dist = mean_h_dist))
})

# Stop the cluster
stopCluster(cl)

# Convert results to a data frame for plotting
results_df <- do.call(rbind, lapply(h_dist_results, function(x) data.frame(nbins1d = x$nbins1d, mean_h_dist = x$mean_h_dist)))
results_df$nbins1d <- factor(results_df$nbins1d, levels = as.character(nbins1d_values))

# Plotting
p <- ggplot(results_df, aes(x = nbins1d, y = mean_h_dist)) +
    geom_bar(stat = "identity") +
    theme_minimal() +
    labs(title = "Mean Hellinger Distance by nbins1d", x = "Number of Bins (nbins1d)", y = "Mean Hellinger Distance")

# Display plot
print(p)

# Optionally, save the plot
ggsave("Hellinger_Distances.png", plot = p, width = 10, height = 8)
