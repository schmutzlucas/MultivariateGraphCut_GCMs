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


tmp <- compute_nd_pdf(variables, model_names, data_dir, year_present,
                      lon, lat, aperm(abind(range_var_final, along = 4), c(1, 2, 4, 3)), nbins1d)

