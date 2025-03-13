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

# Read the reference range (if needed elsewhere)
range_var_final_0 <- readRDS('ranges/range_var_final_GreenwichCentered_1950-2023_90deg_3v.rds')

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
model_names <- read.table('model_names_pr_tas_psl_3mod.txt')
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
  results <- compute_nd_pdf_bias_corrected(variables, reference_name, model_names, data_dir,
                                             year_present, year_future, lon, lat, nbins1d,
                                             workers = 3, buffer = 0.05)
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


