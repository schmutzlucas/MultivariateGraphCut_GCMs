# Install and load necessary libraries
list_of_packages <- read.table("package_list.txt", sep="\n")$V1
new.packages <- list_of_packages[!(list_of_packages %in% installed.packages()[,"Package"])]
if(length(new.packages))
  install.packages(new.packages, repos = "https://cloud.r-project.org")

library(devtools)
lapply(list_of_packages, library, character.only = TRUE)
# install_github("schmutzlucas/gcoWrapR")

# Loading local functions
source_code_dir <- 'functions/'  # The directory where all functions are saved.
file_paths <- list.files(source_code_dir, full.names = TRUE)
for(path in file_paths){
  source(path)
}


# Setting global variables
lon <- -10:120
lat <- 0:75
lon_size <- length(lon)
lat_size <- length(lat)

# Temporal ranges
year_present <<- 1950:1975
year_future <<- 2075:2100



# # Setting global variables
# lon <- -180:179
# lat <- -90:90
# lon_size <- length(lon)
# lat_size <- length(lat)

# # Temporal ranges
# year_present <<- 1950:1975
# year_future <<- 1998:2023

workers <- 8

# Data directory
data_dir <<- 'data/CMIP6_merged_all/'

# List of variables used
variables <- c('pr', 'tas', 'psl')

# Bins for the PDFs (nbins1d is the number of bins per variable)
nbins1d <<- 8
nbins <<- nbins1d^(length(variables))
# (The joint PDF will have nbins1d^n_vars bins)

# Obtain the list of models from a file
model_names <- read.table('model_names_pr_tas_psl_perfect_model.txt')
model_names <- as.list(model_names[['V1']])


GC06_result_list <- list()
GC06_present_list <- list()
GC06_future_list <- list()
GC06_hdist_present_list <- list()
GC06_hdist_future_list <- list()
GC06_partial_hdist_future_list <- list()
MMM_hdist_present_list <- list()
MMM_hdist_future_list <- list()
MMM_partial_hdist_future_list <- list()
MMM_present_list <- list()
MMM_future_list <- list()



for (m in seq_along(model_names)) {

  # Obtain the list of models from a file
  model_names <- read.table('model_names_pr_tas_psl_perfect_model.txt')
  model_names <- as.list(model_names[['V1']])

  # Separate the reference model from the rest.
  reference_name <- model_names[[m]]
  model_names <- model_names[-m]

  cat("Processing model", reference_name, "as reference\n")


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

  # # Time the execution of the new bias-corrected function.
  # time_optimized <- system.time({
  #   results <- compute_nd_pdf_bias_corrected(
  #     variables, reference_name, model_names, data_dir,
  #     year_present, year_future, lon, lat, nbins1d,
  #     range_var = aperm(abind(range_var_final$ranges, along = 4), c(1, 2, 4, 3)), verbose = TRUE
  #   )
  #
  # })
  # cat("Time taken for compute_nd_pdf_bias_corrected: ",
  #     format_time(time_optimized["elapsed"]), "\n")


  # Time the execution of the new bias-corrected function.
  time_optimized <- system.time({
    results <- compute_nd_pdf_bias_corrected_2(
      variables,
      reference_name,
      model_names,
      data_dir,
      year_present,
      year_future,
      lon,
      lat,
      nbins1d,
      workers = workers,    # Adjust the number of workers as needed
      buffer = 0.10,
      verbose = TRUE
    )
  })
  cat("Time taken for compute_nd_pdf_bias_corrected: ",
      format_time(time_optimized["elapsed"]), "\n")



  # Extract PDFs from the unified function results
  pdf_ref_present <- results$pdf_ref$present
  pdf_models_present <- results$pdf_models$present
  pdf_ref_future <- results$pdf_ref$future
  pdf_models_future <- results$pdf_models$future


  # Create a nested list of selected indices for every grid point (using all bins)
  all_bins <- 1:nbins
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

  # GC
  smooth_cost <- 0.6
  tryCatch({
    GC06_result <- GraphCutHellinger_nD_lat(
      pdf_models_future = pdf_models_future,
      h_dist = h_dist_present,
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
  gc()

  # H Dist GC
  # Initialize a lon x lat matrix for each smooth cost
  GC06_hdist_present <- matrix(NA, nrow = length(lon), ncol = length(lat))
  GC06_hdist_future <- matrix(NA, nrow = length(lon), ncol = length(lat))

  GC06_present <- array(NA, dim = c(length(lon), length(lat), nbins))
  GC06_future <- array(NA, dim = c(length(lon), length(lat), nbins))

  for(l in 1:(length(model_names))){  # Ensure that indexing aligns with model names
    islabel <- which(GC06_result$label_attribution == l)
    GC06_hdist_present[islabel] <- h_dist_present[,,l][islabel]
    GC06_hdist_future[islabel] <- h_dist_future[,,l][islabel]

    GC06_present[islabel] <- pdf_models_present[,,,l][islabel]
    GC06_future[islabel] <- pdf_models_future[,,,l][islabel]
  }

  # MMM
  # Compute Multi-Model Mean for Present
  MMM_present <- apply(pdf_models_present, c(1, 2, 3), mean)

  # Compute Multi-Model Mean for Future
  MMM_future <- apply(pdf_models_future, c(1, 2, 3), mean)


  # Initialize arrays to store the Hellinger distance for MMM
  MMM_hdist_present <- array(NA, dim = c(length(lon), length(lat)))
  MMM_hdist_future <- array(NA, dim = c(length(lon), length(lat)))

  # Compute Hellinger distance for the Multi-Model Mean
  for (i in seq_along(lon)) {
    for (j in seq_along(lat)) {
      # Compute Hellinger distance for present
      MMM_hdist_present[i, j] <- sqrt(sum((sqrt(MMM_present[i, j, ]) - sqrt(pdf_ref_present[i, j, ]))^2)) / sqrt(2)

      # Compute Hellinger distance for future
      MMM_hdist_future[i, j] <- sqrt(sum((sqrt(MMM_future[i, j, ]) - sqrt(pdf_ref_future[i, j, ]))^2)) / sqrt(2)
    }
  }
  gc()

  # Partial H_dist computation

  # --- Step 1: Compute ldr_indices for each grid point ---
  # The idea is to compute, for each grid point, the indices of the bins that are NOT in the "central" region.
  # The central region is defined using select_hdr_indices(), which returns the bins containing the central (1-tau) mass.
  # Here, tau=0.1 means we keep the central 90% and treat the remaining 10% as "low density" (i.e., extreme) bins.
  #
  # Assuming:
  #   pdf_ref_future has dimensions [lon_size, lat_size, nbins]
  #   nbins is defined, and lon_size = length(lon), lat_size = length(lat)
  ldr_indices <- vector("list", lon_size)
  for (i in seq_len(lon_size)) {
    ldr_indices[[i]] <- vector("list", lat_size)
    for (j in seq_len(lat_size)) {
      # Compute the central region indices using tau=0.1
      central_indices <- select_hdr_indices(pdf_ref_future[i, j, ], tau = 0.10)
      # Then, define the low density indices as those not in the central region
      ldr_indices[[i]][[j]] <- setdiff(seq_len(nbins), central_indices)
    }
  }

  # --- Step 2: Compute partial Hellinger distances on the ldr regions ---
  # We use compute_partial_hdist() to compute, for each grid point and for each model,
  # the Hellinger distance between the reference and model PDFs, but only over the bins defined by the ldr_indices.
  # This returns an array with dimensions [lon_size, lat_size, n_models]
  partial_hdist_future <- compute_partial_hdist(pdf_ref_future, pdf_models_future, ldr_indices)

  # --- Step 3: (Optional) Compute Mean Partial Hellinger Distances per model ---
  # This gives a summary (scalar) for each model.
  n_models <- length(model_names)
  mean_partial_hdist <- sapply(seq_len(n_models), function(m) {
    mean(partial_hdist_future[,, m], na.rm = TRUE)
  })

  # --- Step 4: (Optional) Compute Full Hellinger Distances for Comparison ---
  # For example, if you have previously computed the full Hellinger distances in an array h_dist_future
  mean_h_dist_future <- sapply(seq_len(n_models), function(m) {
    mean(h_dist_future[,, m], na.rm = TRUE)
  })


  # Compute partial H dist on the GC results :
  # H Dist GC
  # Initialize a lon x lat matrix for each smooth cost
  GC06_partial_hdist_future <- matrix(NA, nrow = length(lon), ncol = length(lat))

  for(l in 1:(length(model_names))){  # Ensure that indexing aligns with model names
    islabel <- which(GC06_result$label_attribution == l)
    GC06_partial_hdist_future[islabel] <- partial_hdist_future[,,l][islabel]
  }


  # --- Compute Partial Hellinger Distance for MMM (Future) on Low Density Regions ---

  # Initialize the output matrix for partial Hellinger distances (MMM)

  MMM_partial_hdist_future <- array(NA, dim = c(lon_size, lat_size))

  for (i in seq_len(lon_size)) {
    for (j in seq_len(lat_size)) {
      # Extract the reference PDF vector at grid point (i, j)
      pdf_ref_vec <- pdf_ref_future[i, j, ]
      # Get the selected low density indices for this grid point
      selected_bins <- ldr_indices[[i]][[j]]

      if (length(selected_bins) > 0) {
        MMM_partial_hdist_future[i, j] <- sqrt(
          sum(( sqrt(MMM_future[i, j, selected_bins]) - sqrt(pdf_ref_vec[selected_bins]) )^2)
        ) / sqrt(2)
      } else {
        MMM_partial_hdist_future[i, j] <- NA
      }
    }
  }
  # Print summary statistics
  cat("Mean Partial Hellinger Distance for MMM (future):", mean(MMM_partial_hdist_future, na.rm = TRUE), "\n")


  # Store results under the current reference model name
  GC06_result_list[[reference_name]] <-                 GC06_result
  GC06_present_list[[reference_name]] <-                GC06_present
  GC06_future_list[[reference_name]] <-                 GC06_future
  GC06_hdist_present_list[[reference_name]] <-          GC06_hdist_present
  GC06_hdist_future_list[[reference_name]] <-           GC06_hdist_future
  GC06_partial_hdist_future_list[[reference_name]]  <-  GC06_partial_hdist_future
  MMM_hdist_present_list[[reference_name]] <-           MMM_hdist_present
  MMM_hdist_future_list[[reference_name]] <-            MMM_hdist_future
  MMM_partial_hdist_future_list[[reference_name]] <-    MMM_partial_hdist_future
  MMM_present_list[[reference_name]] <-                 MMM_present
  MMM_future_list[[reference_name]] <-                  MMM_future

}


# Get the current date and time
current_time <- Sys.time()
# Format the date and time as a string in the format 'yyyymmddhhmm'
formatted_time <- format(current_time, "%Y%m%d%H%M")
# Concatenate the formatted time string with your desired filename
filename <- paste0(formatted_time, "_my_workspace_PME_bias_corrected_22models_10-10.RData")
# Save the workspace using the generated filename
save.image(file = filename, compress = FALSE)



# Boxplot of the H dist projection by ref
# Get the name of the first reference model
ref_name <- names(GC06_hdist_future_list)[1]

# Extract the corresponding Hellinger matrices
gc_h <- GC06_hdist_future_list[[ref_name]]
mmm_h <- MMM_hdist_future_list[[ref_name]]

# Flatten the matrices into vectors
gc_values <- as.vector(gc_h)
mmm_values <- as.vector(mmm_h)

# Combine into a data frame for ggplot
df <- data.frame(
  Hellinger = c(gc_values, mmm_values),
  Method = factor(rep(c("GraphCut", "MMM"), each = length(gc_values)))
)

# Remove NAs if any
df <- na.omit(df)

# Plot
p <- ggplot(df, aes(x = Method, y = Hellinger, fill = Method)) +
  geom_boxplot(outlier.size = 0.8, outlier.alpha = 0.3) +
  scale_fill_manual(values = c("GraphCut" = "#0072B2", "MMM" = "#D55E00")) +
  labs(
    title = paste("Hellinger Distance (Future) - Reference:", ref_name),  # <-- fixed here
    y = "Hellinger Distance",
    x = ""
  ) +
  theme_minimal(base_size = 14) +
  theme(legend.position = "none")

print(p)


# Violin plot of the H dist by ref
{
  # Boxplot of the H dist projection by ref
  # Get the name of the first reference model
  ref_name <- names(GC06_hdist_future_list)[1]

  # Extract the corresponding Hellinger matrices
  gc_h <- GC06_hdist_future_list[[ref_name]]
  mmm_h <- MMM_hdist_future_list[[ref_name]]

  # Flatten the matrices into vectors
  gc_values <- as.vector(gc_h)
  mmm_values <- as.vector(mmm_h)

  # Combine into a data frame for ggplot
  df <- data.frame(
    Hellinger = c(gc_values, mmm_values),
    Method = factor(rep(c("GraphCut", "MMM"), each = length(gc_values)))
  )

  # Remove NAs if any
  df <- na.omit(df)
  library(ggplot2)
  library(dplyr)

  df_all$Reference <- factor(df_all$Reference, levels = unique(df_all$Reference))

  ggplot(df_all, aes(x = Reference, y = Hellinger, fill = Method)) +
    geom_violin(position = position_dodge(width = 0.6),
                width = 2,      # Plus large horizontalement
                adjust = 1,     # Lissage plus doux (tu peux tester 1, 1.5, 2)
                alpha = 0.79) +
    scale_fill_manual(values = c("GraphCut" = "#0072B2", "MMM" = "#D55E00")) +
    labs(
      title = "Distribution of Hellinger Distance (Future) by Reference Model",
      x = "Reference Model",
      y = "Hellinger Distance"
    ) +
    theme_minimal(base_size = 14) +
    theme(
      axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 8),
      legend.position = "top",
      plot.title = element_text(size = 16, face = "bold")
    )

}

# Violin plot of the partial H dist by ref
{

  # Boxplot of the H dist projection by ref
  # Get the name of the first reference model
  ref_name <- names(GC06_hdist_future_list)[1]

  # Extract the corresponding Hellinger matrices
  gc_h <- GC06_partial_hdist_future_list[[ref_name]]
  mmm_h <- MMM_partial_hdist_future_list[[ref_name]]

  # Flatten the matrices into vectors
  gc_values <- as.vector(gc_h)
  mmm_values <- as.vector(mmm_h)

  # Combine into a data frame for ggplot
  df <- data.frame(
    Hellinger = c(gc_values, mmm_values),
    Method = factor(rep(c("GraphCut", "MMM"), each = length(gc_values)))
  )

  # Remove NAs if any
  df <- na.omit(df)
  library(ggplot2)
  library(dplyr)

  df_all$Reference <- factor(df_all$Reference, levels = unique(df_all$Reference))

  ggplot(df_all, aes(x = Reference, y = Hellinger, fill = Method)) +
    geom_violin(position = position_dodge(width = 0.6),
                width = 2,      # Plus large horizontalement
                adjust = 1,     # Lissage plus doux (tu peux tester 1, 1.5, 2)
                alpha = 0.79) +
    scale_fill_manual(values = c("GraphCut" = "#0072B2", "MMM" = "#D55E00")) +
    labs(
      title = "Distribution of Hellinger Distance (Future) by Reference Model",
      x = "Reference Model",
      y = "Hellinger Distance"
    ) +
    theme_minimal(base_size = 14) +
    theme(
      axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 8),
      legend.position = "top",
      plot.title = element_text(size = 16, face = "bold")
    )

}

# Aggregated violin by method
{
  library(ggplot2)
  library(dplyr)

  # Flatten both lists into a long data.frame with method and value
  gc_data <- do.call(c, lapply(GC06_hdist_future_list, as.vector))
  mmm_data <- do.call(c, lapply(MMM_hdist_future_list, as.vector))

  # Remove NA values
  gc_data <- gc_data[!is.na(gc_data)]
  mmm_data <- mmm_data[!is.na(mmm_data)]

  # Combine into a data.frame
  df_all_methods <- data.frame(
    Method = factor(rep(c("GraphCut", "MMM"), times = c(length(gc_data), length(mmm_data))),
                    levels = c("GraphCut", "MMM")),
    Hellinger = c(gc_data, mmm_data)
  )

  # Plot
  ggplot(df_all_methods, aes(x = Method, y = Hellinger, fill = Method)) +
    geom_violin(scale = "area", trim = TRUE, adjust = 1.5) +
    scale_fill_manual(values = c("GraphCut" = "#1f78b4", "MMM" = "#e66101")) +
    theme_minimal(base_size = 14) +
    labs(
      title = "Distribution of Hellinger Distance (Future)",
      subtitle = "Aggregated across all reference models",
      y = "Hellinger Distance",
      x = NULL
    )

}