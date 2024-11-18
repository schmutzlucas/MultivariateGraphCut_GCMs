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

GC_results <- list()
GC_result_hellinger <- list()
# Loop through smooth cost values from 0 to 1 in increments of 0.05
for (smooth_cost in seq(0, 1.2, by = 0.1)) {
  # Wrap each iteration in tryCatch to handle errors gracefully
  tryCatch({
    # Run Graph Cut with the varying smooth cost
    GC_result_hellinger <- GraphCutHellinger_nD(
      pdf_models_future = pdf_models_present,
      h_dist = h_dist,
      weight_data = 1,               # Fixed data weight
      weight_smooth = smooth_cost,   # Varying smooth cost
      nBins = nbins1d^3,
      seed = 1,
      verbose = TRUE,
      rebuild = FALSE
    )

    # Store only essential results if memory is limited (optional)
    GC_results[[paste0("smooth_", smooth_cost)]] <- GC_result_hellinger

    save(GC_results, file = "GC_result_hellinger_lambda.RData", compress = FALSE)

  }, error = function(e) {
    cat("Error encountered with smooth cost =", smooth_cost, ": ", e$message, "\n")
  })
}

# Get the current date and time
current_time <- Sys.time()

# Format the date and time as a string in the format 'yyyymmddhhmm'
formatted_time <- format(current_time, "%Y%m%d%H%M")

# Concatenate the formatted time string with your desired filename
filename <- paste0(formatted_time, "_my_workspace_ERA5_allModels_beforeOptim_3v.RData")

# Save the workspace using the generated filename
save.image(file = filename, compress = FALSE)



# Initialize empty vectors to store the costs
data_costs <- numeric(length(GC_results))
smooth_costs <- numeric(length(GC_results))

# Loop through each result in GC_results and extract the costs
for (i in seq_along(GC_results)) {
  data_costs[i] <- GC_results[[i]]$`Data and smooth cost`$`Data cost`
  smooth_costs[i] <- GC_results[[i]]$`Data and smooth cost`$`Smooth cost`
}

# Print the vectors to check the results
print(data_costs)
print(smooth_costs)


# Extract the smooth cost values from the names
smooth_cost_names <- names(GC_results)
smooth_cost_values <- as.numeric(sub("smooth_", "", smooth_cost_names))

# Get the order of the smooth cost values
order_indices <- order(smooth_cost_values)

# Reorder GC_results, data_costs, and smooth_costs based on the order_indices
GC_results <- GC_results[order_indices]
data_costs <- sapply(GC_results, function(x) x$`Data and smooth cost`$`Data cost`)
smooth_costs <- sapply(GC_results, function(x) x$`Data and smooth cost`$`Smooth cost`)

# Adjust data costs by subtracting the data cost at smooth cost = 0
data_cost_at_zero <- data_costs[which(smooth_cost_values == 0)]
data_costs <- data_costs - data_cost_at_zero

# Normalize smooth costs by dividing by smooth cost values
smooth_costs <- smooth_costs / smooth_cost_values

# Handle cases where smooth_cost_values is 0 to avoid division by zero
smooth_costs[is.nan(smooth_costs) | is.infinite(smooth_costs)] <- 0

# Reorder smooth cost values
smooth_cost_values <- smooth_cost_values[order_indices]

# Plot data costs
plot(smooth_cost_values, data_costs, type = "o", col = "blue",
     xlab = "Smooth Cost", ylab = "Cost",
     main = "Data Cost and Normalized Smooth Cost vs Smoothness Weight",
     pch = 16, ylim = range(c(data_costs, smooth_costs)))

# Add normalized smooth costs to the same plot
lines(smooth_cost_values, smooth_costs, type = "o", col = "red", pch = 16)

# Add a legend
legend("topright", legend = c("Data Cost", "Normalized Smooth Cost"),
       col = c("blue", "red"), pch = 16, lty = 1)




# Load necessary packages
library(ggplot2)
library(pals)
library(reshape2)

# Generate the polychrome color palette and create a named color mapping
color_palette <- pals::glasbey(length(model_names))
names(color_palette) <- model_names  # Associate each color with a model name

# Loop through each smooth cost result in GC_results
for (smooth_cost in names(GC_results)) {
  # Extract the label attribution for the current smooth cost
  GC_labels <- GC_results[[smooth_cost]]$label_attribution

  # Convert the label matrix to a data frame for plotting
  label_df <- melt(GC_labels, c("lon", "lat"), value.name = "label_attribution")
  label_df$lat <- label_df$lat - 90  # Adjust latitudes if necessary

  # Convert label_attribution to a factor with ALL model names as levels
  label_df$label_attribution <- factor(label_df$label_attribution, levels = seq_along(model_names), labels = model_names)

  # Create the plot
  p <- ggplot() +
    geom_tile(data = label_df, aes(x = lon, y = lat, fill = label_attribution)) +
    scale_fill_manual(values = color_palette, na.value = "white", guide = guide_legend(title = "Model Names", ncol = 2)) +  # Keep all model names in the legend
    ggtitle(paste("Label GC Hellinger - Smooth Weight:", smooth_cost)) +
    borders("world2", colour = 'black', lwd = 0.12) +
    scale_x_continuous(expand = c(0, 0)) +
    scale_y_continuous(limits = c(-90, 90), expand = c(0, 0)) +  # Set y-axis limits
    theme(legend.position = 'bottom') +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
    theme(panel.background = element_blank()) +
    xlab('Longitude') +
    ylab('Latitude') +
    theme_bw() +
    theme(
      legend.key.size = unit(0.5, 'cm'),        # Reduce legend key size
      legend.key.height = unit(0.5, 'cm'),      # Reduce legend key height
      legend.key.width = unit(0.5, 'cm'),       # Reduce legend key width
      legend.title = element_text(size = 10),   # Reduce legend title font size
      legend.text = element_text(size = 8),     # Reduce legend text font size
      plot.title = element_text(size = 16),
      plot.subtitle = element_text(size = 12, hjust = 0.5),
      axis.text = element_text(size = 10),
      axis.title = element_text(size = 12)
    ) +
    easy_center_title()

  # Generate file name based on the smooth cost
  name <- paste0("figure/Labels_GC_Hellinger_smooth_1950-1975_", smooth_cost)

  # Save the plot as both PDF and PNG
  ggsave(paste0(name, ".pdf"), plot = p, width = 20, height = 15, units = "cm", dpi = 300)
  ggsave(paste0(name, ".png"), plot = p, width = 20, height = 15, units = "cm", dpi = 300)
}







GC_hdist_future <- list()
GC_hdist <- list()

for (smooth_cost in names(GC_results)) {
  # Initialize a lon x lat matrix for each smooth cost
  GC_hdist[[smooth_cost]] <- matrix(NA, nrow = length(lon), ncol = length(lat))
  GC_hdist_future[[smooth_cost]] <- matrix(NA, nrow = length(lon), ncol = length(lat))

  for(l in 1:(length(model_names))){  # Ensure that indexing aligns with model names
    islabel <- which(GC_results[[smooth_cost]]$label_attribution == l)
    GC_hdist[[smooth_cost]][islabel] <- h_dist[,,l][islabel]
    GC_hdist_future[[smooth_cost]][islabel] <- h_dist_future[,,l][islabel]
  }
}

library(ggplot2)
library(reshape2)
library(dplyr)

# Prepare a data frame for plotting
plot_data <- data.frame()

for (smooth_cost in names(GC_results)) {
  # Extract present and future distances as vectors
  h_dist_present_vect <- as.vector(GC_hdist[[smooth_cost]])
  h_dist_future_vect <- as.vector(GC_hdist_future[[smooth_cost]])

  # Create a data frame with both present and future distances, along with smooth cost and type labels
  smooth_cost_data <- data.frame(
    h_dist = c(h_dist_present_vect, h_dist_future_vect),  # Present first, then Future
    Type = factor(rep(c("Present", "Future"), each = length(h_dist_present_vect)), levels = c("Present", "Future")),  # Order "Present" first
    Smooth_Cost = as.numeric(sub("smooth_", "", smooth_cost))  # Extract numerical value of smooth cost
  )

  # Bind this to the main data frame
  plot_data <- rbind(plot_data, smooth_cost_data)
}

# Convert smooth cost to a factor for ordering in the plot
plot_data$Smooth_Cost <- factor(plot_data$Smooth_Cost, levels = sort(unique(plot_data$Smooth_Cost)))

# Plotting the data with ggplot2
ggplot(plot_data, aes(x = Smooth_Cost, y = h_dist, fill = Type)) +
  geom_boxplot(position = position_dodge(width = 0.75), outlier.shape = NA) +  # Remove outliers from boxplot
  labs(
    title = "Hellinger Distance Comparison for Present and Future",
    x = "Smooth Cost",
    y = "Hellinger Distance",
    fill = "Type"
  ) +
  scale_y_continuous(limits = c(0, 0.4)) +  # Limit y-axis to 0.4
  scale_fill_manual(values = c("Present" = "#00BFC4", "Future" = "#F8766D")) +  # Ensure "Present" color first
  theme_minimal() +
  theme(
    plot.title = element_text(size = 16, hjust = 0.5),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 10),
    legend.position = "bottom"
  )




library(ggplot2)
library(reshape2)


for (smooth_cost in names(GC_results)) {
  # Plot for present Hellinger distance
  h_dist_map <- GC_hdist[[smooth_cost]]
  present_df <- melt(h_dist_map, c("lon", "lat"), value.name = "Bias")

  p_present <- ggplot() +
    geom_tile(data = present_df, aes(x = lon, y = lat - 90, fill = Bias)) +
    labs(subtitle = '',
         title = paste0('GraphCut Hellinger - Present - Smooth Weight: ', smooth_cost,
                        '\nMean Hellinger distance = ', round(mean(h_dist_map, na.rm = TRUE), 3))) +
    scale_fill_gradient(low = "white", high = "#015a8c", limits = c(0.1, 0.70), oob = scales::squish) +
    borders("world2", colour = 'black', lwd = 0.12) +
    scale_x_continuous(expand = c(0, 0)) +
    scale_y_continuous(expand = c(0, 0)) +
    theme(legend.position = 'bottom',
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank()) +
    xlab('Longitude') +
    ylab('Latitude') +
    labs(fill = 'Hellinger \nDistance') +
    theme_bw() +
    theme(legend.key.size = unit(1, 'cm'),
          legend.key.height = unit(1.4, 'cm'),
          legend.key.width = unit(0.4, 'cm'),
          legend.title = element_text(size = 16),
          legend.text = element_text(size = 12),
          plot.title = element_text(size = 24),
          plot.subtitle = element_text(size = 20, hjust = 0.5),
          axis.text = element_text(size = 14),
          axis.title = element_text(size = 16)) +
    easy_center_title()

  # Save present plot
  name_present <- paste0('figure/H_dist_present_GC_hellinger_1950-1975_smooth_', smooth_cost)
  ggsave(paste0(name_present, '.pdf'), plot = p_present, width = 35, height = 25, units = "cm", dpi = 300)
  ggsave(paste0(name_present, '.png'), plot = p_present, width = 35, height = 25, units = "cm", dpi = 300)

  # Plot for future Hellinger distance
  h_dist_map_future <- GC_hdist_future[[smooth_cost]]
  future_df <- melt(h_dist_map_future, c("lon", "lat"), value.name = "Bias")

  p_future <- ggplot() +
    geom_tile(data = future_df, aes(x = lon, y = lat - 90, fill = Bias)) +
    labs(subtitle = 'Projection period: 1999 - 2014',
         title = paste0('GraphCut Hellinger - Future - Smooth Weight: ', smooth_cost,
                        '\nMean Hellinger distance = ', round(mean(h_dist_map_future, na.rm = TRUE), 3))) +
    scale_fill_gradient(low = "white", high = "#015a8c", limits = c(0.1, 0.70), oob = scales::squish) +
    borders("world2", colour = 'black', lwd = 0.12) +
    scale_x_continuous(expand = c(0, 0)) +
    scale_y_continuous(expand = c(0, 0)) +
    theme(legend.position = 'bottom',
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank()) +
    xlab('Longitude') +
    ylab('Latitude') +
    labs(fill = 'Hellinger \nDistance') +
    theme_bw() +
    theme(legend.key.size = unit(1, 'cm'),
          legend.key.height = unit(1.4, 'cm'),
          legend.key.width = unit(0.4, 'cm'),
          legend.title = element_text(size = 16),
          legend.text = element_text(size = 12),
          plot.title = element_text(size = 24),
          plot.subtitle = element_text(size = 20, hjust = 0.5),
          axis.text = element_text(size = 14),
          axis.title = element_text(size = 16)) +
    easy_center_title()

  # Save future plot
  name_future <- paste0('figure/H_dist_future_GC_hellinger_1950-1975_smooth_', smooth_cost)
  ggsave(paste0(name_future, '.pdf'), plot = p_future, width = 35, height = 25, units = "cm", dpi = 300)
  ggsave(paste0(name_future, '.png'), plot = p_future, width = 35, height = 25, units = "cm", dpi = 300)
}
