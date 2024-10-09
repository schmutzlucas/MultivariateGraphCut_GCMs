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

range_var_final <- readRDS('ranges/range_var_final_allModelsPar_1950-2022_70deg.rds')

# Setting global variables
lon <- 0:359
lat <- -70:70
# Temporal ranges
year_present <<- 1977:1999
year_future <<- 2000:2022
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


# Step 1: Extract the PDF for the desired pixel (i = 54, j = 56) and model 1
pdf_present_pixel <- tmp$present[54, 56, , 1]  # Present period for model 1
pdf_future_pixel <- tmp$future[54, 56, , 1]    # Future period for model 1

# Step 2: Reshape the 1D PDF vector into a 3D array
# Assuming we have 8 bins per variable (as set in nbins1d)
nbins <- nbins1d  # Number of bins
pdf_present_3d <- array(pdf_present_pixel, dim = c(nbins, nbins, nbins))
pdf_future_3d <- array(pdf_future_pixel, dim = c(nbins, nbins, nbins))

# Step 3: Prepare data for 3D plotting
# Convert the 3D array into a format suitable for plotly
plot_data_present <- melt(pdf_present_3d)
plot_data_present <- subset(plot_data_present, value > 0)  # Filter non-zero values


# Example pixel and model to plot
lon_idx <- 54
lat_idx <- 56
model_idx <- 1

# Extract the 3D histogram values for the specific pixel and model
hist_values <- tmp$present[lon_idx, lat_idx, , model_idx]

# Extract the range values for the specific pixel from the list structure
pr_range_pixel <- range_var_final$pr[lon_idx, lat_idx, ]
tas_range_pixel <- range_var_final$tas[lon_idx, lat_idx, ]
psl_range_pixel <- range_var_final$psl[lon_idx, lat_idx, ]

# Compute the bin edges using the extracted ranges
pr_bins <- seq(pr_range_pixel[1], pr_range_pixel[2], length.out = 9)  # nbins + 1
tas_bins <- seq(tas_range_pixel[1], tas_range_pixel[2], length.out = 9)
psl_bins <- seq(psl_range_pixel[1], psl_range_pixel[2], length.out = 9)

# Create a data frame for the 3D scatter plot based on the histogram values
plot_data <- expand.grid(Pr = pr_bins[-length(pr_bins)],
                         Tas = tas_bins[-length(tas_bins)],
                         Psl = psl_bins[-length(psl_bins)])
plot_data$Value <- as.vector(hist_values)

# Remove zero values to focus on non-empty bins
plot_data <- subset(plot_data, Value > 0)

# Create a 3D scatter plot with appropriate axes and color scale
plot_ly(data = plot_data, x = ~Pr, y = ~Tas, z = ~Psl, size = ~Value, color = ~Value, colors = c("blue", "red")) %>%
  add_markers(sizemode = "diameter", marker = list(sizeref = 0.1)) %>%  # Adjust sizeref for better scaling
  layout(scene = list(xaxis = list(title = 'Precipitation (pr)', range = c(min(pr_bins), max(pr_bins))),
                      yaxis = list(title = 'Temperature (tas)', range = c(min(tas_bins), max(tas_bins))),
                      zaxis = list(title = 'Pressure (psl)', range = c(min(psl_bins), max(psl_bins)))),
         title = paste("3D Histogram at Pixel (", lon_idx, ",", lat_idx, ") for Model", model_idx))



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

# Loop through variables and models
m <- 1
for (model_name in model_names) {
  for (i in seq_along(lon)) {
    for (j in seq_along(lat)) {
      # Compute Hellinger distance
      h_dist_unchecked[i, j, m] <- sqrt(sum((sqrt(pdf_models_present[i, j, , m]) - sqrt(pdf_ref_present[i, j, ]))^2)) / sqrt(2)
    }
  }
  m <- m + 1
}

hist(h_dist_unchecked)
# Replace NaN with 0
h_dist[,,] <- replace(h_dist_unchecked[,,], is.nan(h_dist_unchecked), 0)
hist(h_dist)
rm(h_dist_unchecked)

# Graphcut hellinger labelling
GC_result_hellinger_new <- list()
GC_result_hellinger_new <- GraphCutHellinger2D_new3(pdf_models_future = pdf_models_future[,,, ] ,
                                                    h_dist = h_dist,
                                                    weight_data = 1,
                                                    weight_smooth = 1,
                                                    nBins = nbins1d^3,
                                                    seed = 2,
                                                    verbose = TRUE,
                                                    rebuild = TRUE)