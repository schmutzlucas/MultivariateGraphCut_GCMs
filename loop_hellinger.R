  for(var in variables){
    i <- 1
    for(model_name in model_names){
      print(var)
      print(model_name)
      # Set the directory path and file name pattern
      dir_path <- paste0('data/CMIP6/', model_name, '/', var, '/')
      print(dir_path)
      file_pattern <- paste0(var, '_', model_name, '_', exp, '_','*.nc')
      print(file_pattern)

      # Get the full file path of the matching file
      file_path <- list.files(path = dir_path, pattern = file_pattern, full.names = TRUE)
      print(file_path)

      # Check that there is only one matching file
      if (length(file_path) == 1) {
        # Open the file using nc_open
        nc <- nc_open(file_path)
        # do something with the netCDF file
        # ...
      } else {
        # Handle the case where there are multiple or no matching files
        print("Error: Found multiple or no matching files")
      }

      # # Create dimensions and initialize matrices for present and future data if this is the first model and variable
      # if(j == 1 && i == 1) {
      #   lat <- ncvar_get(nc, "lat")
      #   lon <- ncvar_get(nc, "lon")
      #   data_matrix$present <- data_matrix$future <-
      #     array(0, c(length(lon), length(lat),
      #                     length(model_names), length(variables)))
      # }

      # # Extract and average data for present
      # yyyy <- substr(as.character(nc.get.time.series(nc)), 1, 4)
      # iyyyy <- which(yyyy %in% year_present)
      # tmp <- ncvar_get(nc, var, start = c(1, 1, min(iyyyy)), count = c(-1, -1, length(iyyyy)))

    }
  }

library(stringr)

dir_path <- "data/"

subdirs <- dir(dir_path, full.names = TRUE, recursive = FALSE)

model_names <- unique(str_extract(subdirs, "(?<=/)[A-Za-z0-9-]+"))

model_names <- str_replace_all(cmip6_models, "^data/", "")


library(ncdf4)

# Define the directory containing the netCDF files
dir_path <- "data"

# Get the list of model names from the folder names in dir_path
model_list <- list.dirs(dir_path, full.names = TRUE, recursive = FALSE)
model_list <- sub(paste0("^data/", dir_path), "", model_list)

# Define a list of variable names to open
var_list <- c("tas", "tasmax")

# Loop through each combination of variable and model name
for (var_name in var_list) {
  for (model_name in model_list) {
    # Define the file path for the current variable and model name
    file_path <- paste0(dir_path, "/", model_name, "/", var_name, "/", var_name, "_", model_name, "_historical_r1i1p1f1_19500101-20141230_merged_regridded_v20190614.nc")
    # Open the file using nc_open() from the ncdf4 package
    nc <- nc_open(file_path)
    # Print some information about the file
    print(paste0("Opened file: ", file_path))
    print(paste0("Variable: ", var_name))
    print(paste0("Model: ", model_name))
    print(paste0("Dimensions: ", dimnames(nc)))
    print(paste0("Variables: ", names(nc$var)))
    print(paste0("Attributes: ", attributes(nc)))
    print("------------------------")
    # Close the file connection
    nc_close(nc)
  }
}

