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


# Case 1
# Ref : Horizontal stripes
# Model 1 : Uniform
# Model 2 : Uniform

# Bins for the pdfs
nbins1d <<- 8

# Create a list of models
model_names <- c("model1", "model2")


# Initialize dimensions
lon <- 45
lat <- 22
pdf_bins <- 512

# Create an empty array to hold the data
reference_matrix <- array(0, dim = c(lon, lat, pdf_bins))

# Define the bins for the pdf (10 bins at 0.1 for pair/impair stripes)
even_bins <- 100:109
odd_bins <- 410:419

# Fill the matrix with horizontal stripe pattern
for (j in 1:lat) {
  if (j %% 2 == 0) {
    # Pair stripes
    reference_matrix[, j, even_bins] <- 0.1
  } else {
    # Impair stripes
    reference_matrix[, j, odd_bins] <- 0.1
  }
}


# Create empty arrays for the two models
model1_matrix <- array(0, dim = c(lon, lat, pdf_bins))
model2_matrix <- array(0, dim = c(lon, lat, pdf_bins))

# Define the bins for the pdf (10 bins at 0.1 for each model's case)
even_bins <- 100:109
odd_bins <- 410:419

# Fill Model 1 uniformly with the pair stripe pdf (bins 100-109 with 0.1)
model1_matrix[, , even_bins] <- 0.1

# Fill Model 2 uniformly with the impair stripe pdf (bins 410-419 with 0.1)
model2_matrix[, , odd_bins] <- 0.1

# Create an empty 4D array to hold both models in the fourth dimension
pdf_models <- array(0, dim = c(lon, lat, pdf_bins, length(model_names)))

# Assign model1 to the first "slice" (m = 1) and model2 to the second "slice" (m = 2)
pdf_models[,,,1] <- model1_matrix
pdf_models[,,,2] <- model2_matrix


# Computing the sum of hellinger distances between models and reference --> used as datacost
h_dist <- array(data = 0, dim = c(lon, lat,
                                  length(model_names)))
h_dist_unchecked <- array(data = 0, dim = c(lon, lat,
                                            length(model_names)))

# Loop through variables and models
m <- 1
for (model_name in model_names) {
  for (i in seq_along(1:lon)) {
    for (j in seq_along(1:lat)) {
      # Compute Hellinger distance
      h_dist_unchecked[i, j, m] <- sqrt(sum((sqrt(pdf_models[i, j, , m]) - sqrt(reference_matrix[i, j, ]))^2)) / sqrt(2)
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
GC_result_hellinger_new <- GraphCutHellinger(pdf_models_future = pdf_models,
                                                    h_dist = h_dist,
                                                    weight_data = 1,
                                                    weight_smooth = 0.1,
                                                    nBins = nbins1d^3,
                                                    seed = 1,
                                                    verbose = TRUE,
                                                    rebuild = TRUE)



image(GC_result_hellinger_new$label_attribution)



# Case 2
# Ref : Vertical stripes
# Model 1 : Uniform
# Model 2 : Uniform


# Initialize dimensions
lon <- 45
lat <- 22
pdf_bins <- 512

# Define the bins for the pdf (10 bins at 0.1 for even/odd stripes)
even_bins <- 100:109
odd_bins <- 410:419

# Create an empty array to hold the reference matrix with vertical stripes
reference_matrix <- array(0, dim = c(lon, lat, pdf_bins))

# Fill the matrix with a vertical stripe pattern
for (i in 1:lon) {
  if (i %% 2 == 0) {
    # Even stripes
    reference_matrix[i, , even_bins] <- 0.1
  } else {
    # Odd stripes
    reference_matrix[i, , odd_bins] <- 0.1
  }
}

# Create empty arrays for the two models
model1_matrix <- array(0, dim = c(lon, lat, pdf_bins))
model2_matrix <- array(0, dim = c(lon, lat, pdf_bins))

# Fill Model 1 uniformly with the even stripe pdf (bins 100-109 with 0.1)
model1_matrix[, , even_bins] <- 0.1

# Fill Model 2 uniformly with the odd stripe pdf (bins 410-419 with 0.1)
model2_matrix[, , odd_bins] <- 0.1

# Create an empty 4D array to hold both models in the fourth dimension
pdf_models <- array(0, dim = c(lon, lat, pdf_bins, length(model_names)))

# Assign model1 to the first "slice" (m = 1) and model2 to the second "slice" (m = 2)
pdf_models[,,,1] <- model1_matrix
pdf_models[,,,2] <- model2_matrix

# Computing the sum of Hellinger distances between models and reference --> used as datacost
h_dist <- array(data = 0, dim = c(lon, lat, length(model_names)))
h_dist_unchecked <- array(data = 0, dim = c(lon, lat, length(model_names)))

# Loop through variables and models to calculate the Hellinger distance
m <- 1
for (model_name in model_names) {
  for (i in 1:lon) {
    for (j in 1:lat) {
      # Compute Hellinger distance for each lon-lat cell and model
      h_dist_unchecked[i, j, m] <- sqrt(sum((sqrt(pdf_models[i, j, , m]) - sqrt(reference_matrix[i, j, ]))^2)) / sqrt(2)
    }
  }
  m <- m + 1
}

# Plotting Hellinger distances
hist(h_dist_unchecked)

# Replace NaN with 0 in Hellinger distance calculations
h_dist[,,] <- replace(h_dist_unchecked[,,], is.nan(h_dist_unchecked), 0)
hist(h_dist)
rm(h_dist_unchecked)

# GraphCut Hellinger labelling
GC_result_hellinger_new <- list()
GC_result_hellinger_new <- GraphCutHellinger(pdf_models_future = pdf_models,
                                                    h_dist = h_dist,
                                                    weight_data = 1,
                                                    weight_smooth = 0.1,
                                                    nBins = nbins1d^3,
                                                    seed = 11,
                                                    verbose = TRUE,
                                                    rebuild = TRUE
)

# Visualize label attribution result
image(GC_result_hellinger_new$label_attribution)





# Case 3: Chessboard Pattern with Larger Squares
# Ref: Chessboard pattern where each "square" is larger than 1x1 (e.g., 3x3),
#      alternating between two different PDF patterns for even and odd squares.
# Model 1: Uniform PDF pattern matching the even "squares" of the reference (bins 100-109)
# Model 2: Uniform PDF pattern matching the odd "squares" of the reference (bins 410-419)

# Bins for the pdfs
nbins1d <<- 8

# Create a list of models
model_names <- c("model1", "model2")

# Initialize dimensions
lon <- 2
lat <- 2
pdf_bins <- 512
square_size <- 1  # Size of each "chessboard" square (3x3)

# Define the bins for the pdf (10 bins at 0.1 for even/odd squares)
even_bins <- 100:109
odd_bins <- 410:419

# Create an empty array to hold the reference matrix with a chessboard pattern
reference_matrix <- array(0, dim = c(lon, lat, pdf_bins))

# Fill the matrix with a chessboard pattern of larger squares
for (i in seq(1, lon, by = square_size)) {
  for (j in seq(1, lat, by = square_size)) {
    # Determine if the current square should use even or odd bins
    if ((i %/% square_size + j %/% square_size) %% 2 == 0) {
      # Even "chessboard" square
      reference_matrix[i:(min(i + square_size - 1, lon)), j:(min(j + square_size - 1, lat)), even_bins] <- 0.1
    } else {
      # Odd "chessboard" square
      reference_matrix[i:(min(i + square_size - 1, lon)), j:(min(j + square_size - 1, lat)), odd_bins] <- 0.1
    }
  }
}


# Create empty arrays for the two models
model1_matrix <- array(0, dim = c(lon, lat, pdf_bins))
model2_matrix <- array(0, dim = c(lon, lat, pdf_bins))

# Fill Model 1 uniformly with the even stripe pdf (bins 100-109 with 0.1)
model1_matrix[, , even_bins] <- 0.1

# Fill Model 2 uniformly with the odd stripe pdf (bins 410-419 with 0.1)
model2_matrix[, , odd_bins] <- 0.1

# Create an empty 4D array to hold both models in the fourth dimension
pdf_models <- array(0, dim = c(lon, lat, pdf_bins, length(model_names)))

# Assign model1 to the first "slice" (m = 1) and model2 to the second "slice" (m = 2)
pdf_models[,,,1] <- model1_matrix
pdf_models[,,,2] <- model2_matrix

# Computing the sum of Hellinger distances between models and reference --> used as datacost
h_dist <- array(data = 0, dim = c(lon, lat, length(model_names)))
h_dist_unchecked <- array(data = 0, dim = c(lon, lat, length(model_names)))

# Loop through variables and models to calculate the Hellinger distance
m <- 1
for (model_name in model_names) {
  for (i in 1:lon) {
    for (j in 1:lat) {
      # Compute Hellinger distance for each lon-lat cell and model
      h_dist_unchecked[i, j, m] <- sqrt(sum((sqrt(pdf_models[i, j, , m]) - sqrt(reference_matrix[i, j, ]))^2)) / sqrt(2)
    }
  }
  m <- m + 1
}

# Plotting Hellinger distances
hist(h_dist_unchecked)

# Replace NaN with 0 in Hellinger distance calculations
h_dist[,,] <- replace(h_dist_unchecked[,,], is.nan(h_dist_unchecked), 0)
hist(h_dist)
rm(h_dist_unchecked)


# GraphCut Hellinger labelling
GC_result_hellinger_new <- list()
GC_result_hellinger_new <- GraphCutHellinger_xD(
  pdf_models_future = pdf_models,
  h_dist = h_dist,
  weight_data = 1,
  weight_smooth = 0.1,
  nBins = nbins1d^3,
  seed = 10,
  verbose = TRUE,
  rebuild = TRUE
)

# Visualize label attribution result
image(GC_result_hellinger_new$label_attribution)



#  Now systematically testing
# Create a directory for saving images if it doesn't exist
if (!dir.exists("figure")) dir.create("figure")

# Initialize a list to store results if needed (consider reducing stored data if memory is limited)
GC_results <- list()

# Loop through smooth cost values from 0 to 1 in increments of 0.05
for (smooth_cost in seq(0, 1, by = 0.05)) {
  # Wrap each iteration in tryCatch to handle errors gracefully
  tryCatch({
    # Run Graph Cut with the varying smooth cost
    GC_result_hellinger_new <- GraphCutHellinger(
      pdf_models_future = pdf_models,
      h_dist = h_dist,
      weight_data = 1,              # Fixed data weight
      weight_smooth = smooth_cost,   # Varying smoothness cost
      nBins = nbins1d^3,
      seed = 1,
      verbose = TRUE,
      rebuild = TRUE
    )

    # Store only essential results if memory is limited (optional)
    GC_results[[paste0("smooth_", smooth_cost)]] <- GC_result_hellinger_new$label_attribution

    # Plot using image with a title and save
    image(GC_result_hellinger_new$label_attribution,
          main = paste("Label Attribution (Smooth Cost =", smooth_cost, ")"))

    # Save the plot as a PNG
    file_name <- paste0("figure/LabelAttribution_smooth_", smooth_cost, ".png")
    dev.copy(png, filename = file_name, width = 35 * 96, height = 25 * 96, res = 300)
    dev.off()

    # Run garbage collection to free memory
    gc()

  }, error = function(e) {
    cat("Error encountered with smooth cost =", smooth_cost, ": ", e$message, "\n")
  })
}

# If not required, the GC_results list can be reduced or stored externally if memory usage is an issue
