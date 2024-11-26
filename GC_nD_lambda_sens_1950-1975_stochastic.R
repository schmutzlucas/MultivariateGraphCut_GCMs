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

range_var_final <- readRDS('ranges/range_var_final_allModelsPar_1950-2023_90deg_3v.rds')

# Setting global variables
lon <- 0:359
lat <- -90:90
# Temporal ranges
year_present <<- 1950:1975
year_future <<- 1998:2023
# data directory
data_dir <<- 'data/CMIP6_merged_all/'

# Bins for the pdfs
nbins1d <<- 8


# List of the variable used
variables <- c('pr', 'tas', 'psl')

# Obtains the list of models from the model names or from a file
model_names <- read.table('model_names_pr_tas_psl.txt')
model_names <- as.list(model_names[['V1']])
# Index of the reference
ref_index <<- 1
# Custom function to format time into human-readable format
format_time <- function(time_seconds) {
  hours <- floor(time_seconds / 3600)
  minutes <- floor((time_seconds %% 3600) / 60)
  seconds <- round(time_seconds %% 60, 2)
  if (hours > 0) {
    return(paste(hours, "hours", minutes, "minutes", seconds, "seconds"))
  } else if (minutes > 0) {
    return(paste(minutes, "minutes", seconds, "seconds"))
  } else {
    return(paste(seconds, "seconds"))
  }
}

# Time the execution of the optimized function
time_optimized <- system.time({
  # todo add number of workers as argument
  tmp <- compute_nd_pdf_optimized(variables, model_names, data_dir, year_present, year_future,
                                  lon, lat, aperm(abind(range_var_final, along = 4), c(1, 2, 4, 3)), nbins1d)
})
cat("Time taken for compute_nd_pdf_optimized: ", format_time(time_optimized["elapsed"]), "\n")

# Choose the reference in the models
reference_name <<- model_names[ref_index]
model_names <<- model_names[-ref_index]


# Get the current date and time
current_time <- Sys.time()

# Format the date and time as a string in the format 'yyyymmddhhmm'
formatted_time <- format(current_time, "%Y%m%d%H%M")

# Concatenate the formatted time string with your desired filename
filename <- paste0(formatted_time, "_my_workspace_ERA5_allModels_beforeOptim_3v.RData")

# Save the workspace using the generated filename
save.image(file = filename, compress = FALSE)

pdf_present <- tmp$present
pdf_future <- tmp$future

pdf_ref_present <- pdf_present[ , , , 1]
pdf_models_present <- pdf_present[ , , , -1]

pdf_ref_future <- pdf_future[ , , , 1]
pdf_models_future <- pdf_future[ , , , -1]

rm(pdf_present, pdf_future)


# Computing the sum of hellinger distances between models and reference --> used as datacost
h_dist <- array(data = 0, dim = c(length(lon), length(lat),
                                  length(model_names)))
h_dist_unchecked <- array(data = 0, dim = c(length(lon), length(lat),
                                            length(model_names)))

h_dist_future <- array(data = 0, dim = c(length(lon), length(lat),
                                         length(model_names)))
h_dist_unchecked_future <- array(data = 0, dim = c(length(lon), length(lat),
                                                   length(model_names)))

# Loop through variables and models
m <- 1
for (model_name in model_names) {
  for (i in seq_along(lon)) {
    for (j in seq_along(lat)) {
      # Compute Hellinger distance
      h_dist_unchecked[i, j, m] <- sqrt(sum((sqrt(pdf_models_present[i, j, , m]) - sqrt(pdf_ref_present[i, j, ]))^2)) / sqrt(2)
      h_dist_unchecked_future[i, j, m] <- sqrt(sum((sqrt(pdf_models_future[i, j, , m]) - sqrt(pdf_ref_future[i, j, ]))^2)) / sqrt(2)
    }
  }
  m <- m + 1
}

hist(h_dist_unchecked)
# Replace NaN with 0
h_dist[,,] <- replace(h_dist_unchecked[,,], is.nan(h_dist_unchecked), 0)
h_dist_future[,,] <- replace(h_dist_unchecked_future[,,], is.nan(h_dist_unchecked), 0)
hist(h_dist)
rm(h_dist_unchecked, h_dist_unchecked_future)


# Get the current date and time
current_time <- Sys.time()

# Format the date and time as a string in the format 'yyyymmddhhmm'
formatted_time <- format(current_time, "%Y%m%d%H%M")

# Concatenate the formatted time string with your desired filename
filename <- paste0(formatted_time, "_my_workspace_ERA5_allModels_beforeOptim_3v.RData")

# Save the workspace using the generated filename
save.image(file = filename, compress = FALSE)


# Initialize lists to store results and seeds
GC_results_stoch <- list()

# Loop through lambda values from 0 to 0.8 in increments of 0.1
for (lambda in seq(0, 0.8, by = 0.1)) {
  # Initialize a sub-list to store results for each lambda
  GC_results_stoch[[paste0("lambda_", lambda)]] <- list()

  # Run 10 iterations for each lambda with different seeds
  for (i in 1:10) {
    # Wrap each iteration in tryCatch to handle errors gracefully
    tryCatch({
      # Run Graph Cut with the current lambda and seed
      GC_result_hellinger <- GraphCutHellinger_nD(
        pdf_models_future = pdf_models_present,
        h_dist = h_dist,
        weight_data = 1,               # Fixed data weight
        weight_smooth = lambda,        # Varying lambda
        nBins = nbins1d^3,
        seed = i,               # Use the pre-generated seed
        verbose = TRUE,
        rebuild = FALSE
      )

      # Store the result in the sub-list for this lambda
      GC_results_stoch[[paste0("lambda_", lambda)]][[paste0("iteration_", i)]] <- GC_result_hellinger

      # Save results after each iteration to ensure progress is not lost
      save(GC_results_stoch, file = "GC_result_hellinger_lambda.RData", compress = FALSE)

    }, error = function(e) {
      cat("Error encountered with lambda =", lambda, "and iteration =", i, ": ", e$message, "\n")
    })
  }
}

# Get the current date and time
current_time <- Sys.time()

# Format the date and time as a string in the format 'yyyymmddhhmm'
formatted_time <- format(current_time, "%Y%m%d%H%M")

# Concatenate the formatted time string with your desired filename
filename <- paste0(formatted_time, "_my_workspace_ERA5_allModels_beforeOptim_3v.RData")

# Save the workspace using the generated filename
save.image(file = filename, compress = FALSE)


library(ggplot2)

# Initialize lists to store costs across iterations for each lambda
data_costs_all <- list()
smooth_costs_all <- list()

# Loop through each lambda and extract the costs across iterations
for (i in seq_along(GC_results_stoch)) {
  # Extract all iterations for the current lambda
  iterations <- GC_results_stoch[[i]]

  # Initialize vectors to store costs for each iteration
  data_costs_iter <- numeric(length(iterations))
  smooth_costs_iter <- numeric(length(iterations))

  # Extract costs for each iteration
  for (j in seq_along(iterations)) {
    data_costs_iter[j] <- iterations[[j]]$`Data and smooth cost`$`Data cost`
    smooth_costs_iter[j] <- iterations[[j]]$`Data and smooth cost`$`Smooth cost`
  }

  # Store the costs for the current lambda
  data_costs_all[[i]] <- data_costs_iter
  smooth_costs_all[[i]] <- smooth_costs_iter
}

# Compute mean and standard deviation for data costs and smooth costs
lambda_values <- as.numeric(sub("lambda_", "", names(GC_results_stoch)))
data_costs_mean <- sapply(data_costs_all, mean)
data_costs_sd <- sapply(data_costs_all, sd)
smooth_costs_mean <- sapply(smooth_costs_all, mean)
smooth_costs_sd <- sapply(smooth_costs_all, sd)

# Adjust data costs by subtracting the mean at lambda = 0
data_cost_at_zero <- data_costs_mean[which(lambda_values == 0)]
data_costs_mean <- data_costs_mean - data_cost_at_zero

# Normalize smooth costs by dividing by lambda values
normalized_smooth_costs_mean <- smooth_costs_mean / lambda_values
normalized_smooth_costs_sd <- smooth_costs_sd / lambda_values

# Handle cases where lambda_values is 0 to avoid division by zero
normalized_smooth_costs_mean[is.nan(normalized_smooth_costs_mean) | is.infinite(normalized_smooth_costs_mean)] <- 0
normalized_smooth_costs_sd[is.nan(normalized_smooth_costs_sd) | is.infinite(normalized_smooth_costs_sd)] <- 0

# Create a data frame for ggplot
plot_data <- data.frame(
  lambda = lambda_values,
  data_cost_mean = data_costs_mean,
  data_cost_lower = data_costs_mean - data_costs_sd,
  data_cost_upper = data_costs_mean + data_costs_sd,
  smooth_cost_mean = normalized_smooth_costs_mean,
  smooth_cost_lower = normalized_smooth_costs_mean - normalized_smooth_costs_sd,
  smooth_cost_upper = normalized_smooth_costs_mean + normalized_smooth_costs_sd
)

# Plot the data with confidence intervals using ggplot2
ggplot(plot_data, aes(x = lambda)) +
  geom_line(aes(y = data_cost_mean, color = "Data Cost")) +
  geom_ribbon(aes(ymin = data_cost_lower, ymax = data_cost_upper, fill = "Data Cost"), alpha = 0.2) +
  geom_line(aes(y = smooth_cost_mean, color = "Normalized Smooth Cost")) +
  geom_ribbon(aes(ymin = smooth_cost_lower, ymax = smooth_cost_upper, fill = "Normalized Smooth Cost"), alpha = 0.2) +
  labs(
    title = "Data Cost and Normalized Smooth Cost with Confidence Intervals",
    x = "Lambda",
    y = "Cost",
    color = "Metric",
    fill = "Metric"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 16, hjust = 0.5),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 12),
    legend.position = "bottom"
  )
