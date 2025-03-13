# Install and load necessary libraries
list_of_packages <- read.table("package_list.txt", sep="\n")$V1
new.packages <- list_of_packages[!(list_of_packages %in% installed.packages()[,"Package"])]
if(length(new.packages))
  install.packages(new.packages, repos = "https://cloud.r-project.org")

library(devtools)
lapply(list_of_packages, library, character.only = TRUE)
install_github("schmutzlucas/gcoWrapR")

# Loading local functions
source_code_dir <- 'functions/'  # The directory where all functions are saved.
file_paths <- list.files(source_code_dir, full.names = TRUE)
for(path in file_paths){
  source(path)
}


# Setting global variables
lon <- -180:179
lat <- -90:90
lon_size <- length(lon)
lat_size <- length(lat)

# Temporal ranges
year_present <<- 1950:1975
year_future <<- 1998:2023

# Data directory
data_dir <<- 'data/CMIP6_merged_all/'

# List of variables used
variables <- c('pr', 'tas', 'psl')

# Bins for the PDFs (nbins1d is the number of bins per variable)
nbins1d <<- 8
# (The joint PDF will have nbins1d^n_vars bins)

# Obtain the list of models from a file
model_names <- read.table('model_names_pr_tas_psl_short.txt')
model_names <- as.list(model_names[['V1']])
ref_index <<- 1

# Separate the reference model from the rest.
reference_name <<- model_names[[ref_index]]
model_names <<- model_names[-ref_index]

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

# Time the execution of the new bias-corrected function.
time_optimized <- system.time({
  results <- compute_nd_pdf_bias_corrected_2(variables, reference_name, model_names, data_dir,
                                             year_present, year_future, lon, lat, nbins1d,
                                             workers = 4, buffer = 0.15, verbose = TRUE)
})
cat("Time taken for compute_nd_pdf_bias_corrected: ",
    format_time(time_optimized["elapsed"]), "\n")

# Store the returned components.
pdf_ref    <- results$pdf_reference      # List with $present and $future PDFs for the reference.
pdf_models <- results$pdf_models           # List with $present and $future PDFs for the models.
ref_stats  <- results$reference_stats      # List with computed reference statistics (present and future).

# Get the current date and time
current_time <- Sys.time()

# Format the date and time as a string in the format 'yyyymmddhhmm'
formatted_time <- format(current_time, "%Y%m%d%H%M")

# Concatenate the formatted time string with your desired filename
filename <- paste0(formatted_time, "_my_workspace_ERA5_bias_corrected_3models.RData")

# Save the workspace using the generated filename
save.image(file = filename, compress = FALSE)

# Plotting the 3d pdf of one gridpoint
{
  # Example indices (adjust as needed)
lon_index <- 10   # selected longitude index
lat_index <- 10   # selected latitude index
model_idx <- 3     # select one model (from the pdf_models array)
nbins <- 8       # number of bins per variable

# Extract the PDF vector for the chosen grid point and model.
# Here pdf_models is from results$pdf_models$present and has dimensions:
# [lon, lat, nbins^n_vars, num_models].
pdf_vector <- results$pdf_models$present[lon_index, lat_index, , model_idx]

# Reshape the 1D PDF vector into a 3D array.
pdf_3d <- array(pdf_vector, dim = c(nbins, nbins, nbins))

# Extract the reference range for that grid point.
# We assume results$ref_range_present has dimensions: [lon, lat, n_vars, 2]
range_mat <- results$ref_range_present[lon_index, lat_index, , ]  # dimensions: [3, 2]

# Compute bin edges and centers for each variable.
centers <- list()
for(v in 1:3) {
  bin_edges <- seq(range_mat[v, 1], range_mat[v, 2], length.out = nbins + 1)
  centers[[v]] <- (bin_edges[-1] + bin_edges[-length(bin_edges)])/2
}

# Create a grid of bin centers.
grid <- expand.grid(x = centers[[1]], y = centers[[2]], z = centers[[3]])

# Flatten the 3D PDF into a vector.
pdf_flat <- as.vector(pdf_3d)
# Normalize PDF values for marker size.
normalized_pdf <- pdf_flat / max(pdf_flat, na.rm = TRUE)

library(plotly)
fig <- plot_ly(
  data = grid,
  x = ~x,
  y = ~y,
  z = ~z,
  type = "scatter3d",
  mode = "markers",
  marker = list(
    size = ~normalized_pdf * 75,  # Adjust scaling factor as needed
    color = ~pdf_flat,
    colorscale = "Viridis",
    showscale = TRUE
  ),
  text = ~paste("PDF Value:", round(pdf_flat, 4))
) %>% layout(
  scene = list(
    xaxis = list(title = "pr"),
    yaxis = list(title = "tas"),
    zaxis = list(title = "psl")
  ),
  title = paste("3D PDF for Model", model_idx, "Grid Point (Lon:", lon_index, ", Lat:", lat_index, ")")
)

fig

}

# --- Begin Complete Hellinger Distance Computation using compute_partial_hdist ---

# Extract PDFs from the results
pdf_ref_present <- pdf_ref$present
pdf_models_present <- pdf_models$present
pdf_ref_future <- pdf_ref$future
pdf_models_future <- pdf_models$future

# Total number of bins (joint PDF) per grid point
n_bins_total <- nbins1d^(length(variables))

# Create a nested list of selected indices for every grid point (using all bins)
all_bins <- 1:n_bins_total
selected_indices_all <- vector("list", length(lon))
for(i in seq_along(lon)) {
  selected_indices_all[[i]] <- vector("list", length(lat))
  for(j in seq_along(lat)) {
    selected_indices_all[[i]][[j]] <- all_bins
  }
}

# Compute complete Hellinger distances for present and future using compute_partial_hdist:
h_dist_present <- compute_partial_hdist(pdf_ref_present, pdf_models_present, selected_indices_all)
h_dist_future  <- compute_partial_hdist(pdf_ref_future, pdf_models_future, selected_indices_all)

# Plot histograms for debugging
hist(h_dist_present, main = "Complete Hellinger Distance (Present)")
hist(h_dist_future, main = "Complete Hellinger Distance (Future)")

# --- End Complete Hellinger Distance Computation ---

smooth_cost <- 1
tryCatch({
  GC_result11 <- GraphCutHellinger_nD_lat(
    pdf_models_present = pdf_models_present,
    h_dist = h_dist_future,
    weight_data = 1,               # Fixed data weight
    weight_smooth = smooth_cost,   # Varying smooth cost
    nBins = nbins1d^3,
    lat = lat,
    seed = 1,
    verbose = TRUE,
    rebuild = TRUE
  )
}, error = function(e) {
  cat("Error encountered with smooth cost =", smooth_cost, ": ", e$message, "\n")
})
