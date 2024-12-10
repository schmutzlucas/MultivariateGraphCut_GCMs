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

range_var_final <- readRDS('ranges/range_var_final_allModelsPar_1950-2100_90deg_3v.rds')

# Setting global variables
lon <- 0:359
lat <- -90:90
# Temporal ranges
year_present <<- 1950:1975
year_future <<- 2076:2100
# data directory
data_dir <<- 'data/CMIP6_merged_all/'

# Bins for the pdfs
nbins1d <<- 8


# List of the variable used
variables <- c('pr', 'tas', 'psl')

# Obtains the list of models from the model names or from a file
model_names <- read.table('model_names_pr_tas_psl_perfect_model.txt')
model_names <- as.list(model_names[['V1']])
# Index of the reference
ref_index <<- 8
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
filename <- paste0(formatted_time, "_my_workspace_PerfectModel_pdf.RData")

# Save the workspace using the generated filename
save.image(file = filename, compress = FALSE)

pdf_present <- tmp$present
pdf_future <- tmp$future

pdf_ref_present <- pdf_present[ , , , ref_index]
pdf_models_present <- pdf_present[ , , , -ref_index]

pdf_ref_future <- pdf_future[ , , , ref_index]
pdf_models_future <- pdf_future[ , , , -ref_index]

rm(pdf_present, pdf_future)


# Initialize arrays for Hellinger distances
h_dist <- array(NA, dim = c(length(lon), length(lat), length(model_names)))
h_dist_future <- array(NA, dim = c(length(lon), length(lat), length(model_names)))

# Compute Hellinger distances for each model
m <- 1
for (model_name in model_names) {
  for (i in seq_along(lon)) {
    for (j in seq_along(lat)) {
      # Compute Hellinger distance for the present
      h_dist[i, j, m] <- sqrt(sum((sqrt(pdf_models_present[i, j, , m]) - sqrt(pdf_ref_present[i, j, ]))^2)) / sqrt(2)

      # Compute Hellinger distance for the future
      h_dist_future[i, j, m] <- sqrt(sum((sqrt(pdf_models_future[i, j, , m]) - sqrt(pdf_ref_future[i, j, ]))^2)) / sqrt(2)
    }
  }
  m <- m + 1
}

# Replace NaN values with 0 in Hellinger distance arrays
h_dist <- replace(h_dist, is.nan(h_dist), 0)
h_dist_future <- replace(h_dist_future, is.nan(h_dist_future), 0)

# Check distributions
hist(h_dist)
hist(h_dist_future)



# Get the current date and time
current_time <- Sys.time()

# Format the date and time as a string in the format 'yyyymmddhhmm'
formatted_time <- format(current_time, "%Y%m%d%H%M")

# Concatenate the formatted time string with your desired filename
filename <- paste0(formatted_time, "_my_workspace_PerfectModel_pdf_hdist.RData")

# Save the workspace using the generated filename
save.image(file = filename, compress = FALSE)


# Initialize lists to store results and seeds
GC_results_stoch <- list()
lambdas_loop <- c(0, 0.1, 0.2, 0.4, 0.6, 0.8, 1.0, 2)

# Loop through the specified lambda values
for (lambda in lambdas_loop) {
  # Initialize a sub-list to store results for each lambda
  GC_results_stoch[[paste0("lambda_", lambda)]] <- list()

  # Run 10 iterations for each lambda with different seeds
  for (i in 1:1) {
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
filename <- paste0(formatted_time, "_my_workspace_PerfectModel_stochResults.RData")

# Save the workspace using the generated filename
save.image(file = filename, compress = FALSE)



# Creating figure for the total cost (Data + Smooth) with min-max range
{
  library(ggplot2)

  # Initialize lists to store total costs across iterations for each lambda
  total_costs_all <- list()

  # Loop through each lambda and extract the costs across iterations
  for (i in seq_along(GC_results_stoch)) {
    # Extract all iterations for the current lambda
    iterations <- GC_results_stoch[[i]]

    # Initialize vector to store total costs for each iteration
    total_costs_iter <- numeric(length(iterations))

    # Extract costs for each iteration
    for (j in seq_along(iterations)) {
      data_cost <- iterations[[j]]$`Data and smooth cost`$`Data cost`
      smooth_cost <- iterations[[j]]$`Data and smooth cost`$`Smooth cost`
      total_costs_iter[j] <- data_cost + smooth_cost  # Sum raw data and smooth costs
    }

    # Store the total costs for the current lambda
    total_costs_all[[i]] <- total_costs_iter
  }

  # Compute mean, min, and max for total costs
  lambda_values <- as.numeric(sub("lambda_", "", names(GC_results_stoch)))
  total_costs_mean <- sapply(total_costs_all, mean)
  total_costs_min <- sapply(total_costs_all, min)
  total_costs_max <- sapply(total_costs_all, max)

  # Create a data frame for ggplot
  plot_data <- data.frame(
    lambda = lambda_values,
    total_cost_mean = total_costs_mean,
    total_cost_min = total_costs_min,
    total_cost_max = total_costs_max
  )

  # Plot the total cost with min-max range using ggplot2
  p <- ggplot(plot_data, aes(x = lambda)) +
    geom_line(aes(y = total_cost_mean, color = "Total Cost"), size = 1) +
    geom_point(aes(y = total_cost_mean, color = "Total Cost"), size = 2) +
    geom_ribbon(aes(ymin = total_cost_min, ymax = total_cost_max, fill = "Total Cost"), alpha = 0.2) +
    labs(
      title = "Total Cost (Data + Smooth) with Min-Max Range",
      x = "Lambda",
      y = "Total Cost",
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

  # Save the plot in the figure folder
  output_dir <- "figure/Total_Cost"
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  file_name_pdf <- file.path(output_dir, "Total_Cost_Lambda.pdf")
  ggsave(file_name_pdf, plot = p, width = 25, height = 20, units = "cm", dpi = 300)
  file_name_png <- file.path(output_dir, "Total_Cost_Lambda.png")
  ggsave(file_name_png, plot = p, width = 25, height = 20, units = "cm", dpi = 300)

  # Print the plot
  print(p)
}

# Creating figure of data (standardized) and smooth (normalized) cost with min max interval
{
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
  p <- ggplot(plot_data, aes(x = lambda)) +
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
  p
  # Save the plot in the figure folder
  output_dir <- "figure/Total_Cost"
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  file_name_pdf <- file.path(output_dir, "data_smooth_cost_Lambda.pdf")
  ggsave(file_name_pdf, plot = p, width = 25, height = 20, units = "cm", dpi = 300)
  file_name_png <- file.path(output_dir, "data_smooth_cost_Lambda.png")
  ggsave(file_name_png, plot = p, width = 25, height = 20, units = "cm", dpi = 300)
}

# Figure of the labelling maps
{
  # Load necessary packages
  library(ggplot2)
  library(pals)
  library(reshape2)
  library(fs)  # To handle folder creation
  library(ggeasy)

  # Generate the polychrome color palette and create a named color mapping
  color_palette <- pals::glasbey(length(model_names))
  names(color_palette) <- model_names  # Associate each color with a model name

  # Loop through each lambda and its iterations in GC_results_stoch
  for (lambda in names(GC_results_stoch)) {
    # Create a folder for the current lambda
    lambda_folder <- paste0("figure/labelling/", lambda)
    dir_create(lambda_folder)  # Create the lambda folder if it doesn't exist

    # Loop through each iteration for the current lambda
    for (iteration in names(GC_results_stoch[[lambda]])) {
      # Extract the label attribution for the current iteration
      GC_labels <- GC_results_stoch[[lambda]][[iteration]]$label_attribution

      # Convert the label matrix to a data frame for plotting
      label_df <- melt(GC_labels, c("lon", "lat"), value.name = "label_attribution")
      label_df$lat <- label_df$lat - 90  # Adjust latitudes if necessary

      # Convert label_attribution to a factor with ALL model names as levels
      label_df$label_attribution <- factor(label_df$label_attribution, levels = seq_along(model_names), labels = model_names)

      # Create the plot
      p <- ggplot() +
        geom_tile(data = label_df, aes(x = lon, y = lat, fill = label_attribution)) +
        scale_fill_manual(values = color_palette, na.value = "white", guide = guide_legend(title = "Model Names", ncol = 2)) +  # Keep all model names in the legend
        ggtitle(paste("Label GC Hellinger - Lambda:", lambda)) +
        labs(subtitle = paste("Seed:", iteration)) +
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

      # Generate file name based on the lambda and iteration
      name <- paste0(lambda_folder, "/Labels_GC_Hellinger_lambda_", lambda, "_seed_", iteration)

      # Save the plot as both PDF and PNG
      # ggsave(paste0(name, ".pdf"), plot = p, width = 20, height = 15, units = "cm", dpi = 300)
      ggsave(paste0(name, ".png"), plot = p, width = 20, height = 15, units = "cm", dpi = 300)
    }
  }
}

# Computing the hellinger distance grid for each iteration
{
  GC_hdist_future <- list()
  GC_hdist <- list()

  for (lambda in names(GC_results_stoch)) {
    # Initialize sub-lists to store gradient errors for each iteration under the current lambda
    GC_hdist[[lambda]] <- list()
    GC_hdist_future[[lambda]] <- list()

    # Loop through each iteration for the current lambda
    for (iteration in names(GC_results_stoch[[lambda]])) {
      # Initialize a lon x lat matrix for the current iteration
      GC_hdist[[lambda]][[iteration]] <- matrix(NA, nrow = length(lon), ncol = length(lat))
      GC_hdist_future[[lambda]][[iteration]] <- matrix(NA, nrow = length(lon), ncol = length(lat))

      # Compute the Hellinger distances for present and future for each model
      for (l in 1:(length(model_names))) {  # Ensure that indexing aligns with model names
        islabel <- which(GC_results_stoch[[lambda]][[iteration]]$label_attribution == l)
        GC_hdist[[lambda]][[iteration]][islabel] <- h_dist[,,l][islabel]
        GC_hdist_future[[lambda]][[iteration]][islabel] <- h_dist_future[,,l][islabel]
      }
    }
  }
}

# Maps of Hellinger distance for each iteration and lambda
{# Initialize separate data frames for present and future average Hellinger distances
  average_hdist_present <- data.frame(Lambda = numeric(), Seed = numeric(), Average_Hellinger = numeric())
  average_hdist_future <- data.frame(Lambda = numeric(), Seed = numeric(), Average_Hellinger = numeric())

  # Loop through lambdas
  for (lambda in names(GC_hdist)) {
    # Extract the numeric value of lambda
    lambda_value <- as.numeric(sub("lambda_", "", lambda))

    # Create subfolders for present and future
    dir.create(file.path("figure", "H_dist_present", lambda), recursive = TRUE, showWarnings = FALSE)
    dir.create(file.path("figure", "H_dist_future", lambda), recursive = TRUE, showWarnings = FALSE)

    # Loop through iterations
    for (seed in names(GC_hdist[[lambda]])) {
      # Process present Hellinger distances
      h_dist_map <- GC_hdist[[lambda]][[seed]]
      avg_hdist_present <- mean(h_dist_map, na.rm = TRUE)

      # Store the average Hellinger distance for present
      average_hdist_present <- rbind(average_hdist_present, data.frame(
        Lambda = lambda_value,
        Seed = as.numeric(sub("iteration_", "", seed)),
        Average_Hellinger = avg_hdist_present
      ))

      # Prepare data for present visualization
      present_df <- melt(h_dist_map, varnames = c("lon", "lat"), value.name = "Hellinger_Distance")
      present_df$lat <- present_df$lat - 90  # Adjust latitude for plotting

      # Create the plot for present Hellinger distance
      p_present <- ggplot() +
        geom_tile(data = present_df, aes(x = lon, y = lat, fill = Hellinger_Distance)) +
        ggtitle("GraphCut Hellinger Distance - Present") +
        labs(
          subtitle = paste0("Lambda: ", lambda_value, " | Seed: ", seed,
                            "\nMean Hellinger Distance: ", round(avg_hdist_present, 4))
        ) +
        scale_fill_gradient(low = "white", high = "#015a8c", oob = scales::squish) +
        borders("world2", colour = "black", lwd = 0.12) +
        scale_x_continuous(expand = c(0, 0)) +
        scale_y_continuous(expand = c(0, 0)) +
        theme_minimal() +
        xlab("Longitude") +
        ylab("Latitude") +
        labs(fill = "Hellinger \nDistance") +
        theme(
          plot.title = element_text(size = 18),
          plot.subtitle = element_text(size = 12),
          axis.title = element_text(size = 14),
          axis.text = element_text(size = 10),
          legend.title = element_text(size = 12),
          legend.text = element_text(size = 10)
        )

      # Save present plot
      name_present <- file.path("figure", "H_dist_present", lambda, paste0("seed_", seed))
      # ggsave(paste0(name_present, ".pdf"), plot = p_present, width = 20, height = 15, units = "cm", dpi = 300)
      ggsave(paste0(name_present, ".png"), plot = p_present, width = 20, height = 15, units = "cm", dpi = 300)

      # Process future Hellinger distances
      h_dist_map_future <- GC_hdist_future[[lambda]][[seed]]
      avg_hdist_future <- mean(h_dist_map_future, na.rm = TRUE)

      # Store the average Hellinger distance for future
      average_hdist_future <- rbind(average_hdist_future, data.frame(
        Lambda = lambda_value,
        Seed = as.numeric(sub("iteration_", "", seed)),
        Average_Hellinger = avg_hdist_future
      ))

      # Prepare data for future visualization
      future_df <- melt(h_dist_map_future, varnames = c("lon", "lat"), value.name = "Hellinger_Distance")
      future_df$lat <- future_df$lat - 90  # Adjust latitude for plotting

      # Create the plot for future Hellinger distance
      p_future <- ggplot() +
        geom_tile(data = future_df, aes(x = lon, y = lat, fill = Hellinger_Distance)) +
        ggtitle("GraphCut Hellinger Distance - Future") +
        labs(
          subtitle = paste0("Lambda: ", lambda_value, " | Seed: ", seed,
                            "\nMean Hellinger Distance: ", round(avg_hdist_future, 4))
        ) +
        scale_fill_gradient(low = "white", high = "#015a8c", oob = scales::squish) +
        borders("world2", colour = "black", lwd = 0.12) +
        scale_x_continuous(expand = c(0, 0)) +
        scale_y_continuous(expand = c(0, 0)) +
        theme_minimal() +
        xlab("Longitude") +
        ylab("Latitude") +
        labs(fill = "Hellinger \nDistance") +
        theme(
          plot.title = element_text(size = 18),
          plot.subtitle = element_text(size = 12),
          axis.title = element_text(size = 14),
          axis.text = element_text(size = 10),
          legend.title = element_text(size = 12),
          legend.text = element_text(size = 10)
        )

      # Save future plot
      name_future <- file.path("figure", "H_dist_future", lambda, paste0("seed_", seed))
      # ggsave(paste0(name_future, ".pdf"), plot = p_future, width = 20, height = 15, units = "cm", dpi = 300)
      ggsave(paste0(name_future, ".png"), plot = p_future, width = 20, height = 15, units = "cm", dpi = 300)
    }
  }

  # Save average Hellinger distances
  save(average_hdist_present, file = "average_hdist_present.RData")
  save(average_hdist_future, file = "average_hdist_future.RData")
}

# Average Hellinger distance with min max interval
{
  library(ggplot2)

  # Compute min, max, and mean for present and future
  summary_hdist_present <- aggregate(Average_Hellinger ~ Lambda, data = average_hdist_present,
                                     FUN = function(x) c(mean = mean(x), min = min(x), max = max(x)))
  summary_hdist_future <- aggregate(Average_Hellinger ~ Lambda, data = average_hdist_future,
                                    FUN = function(x) c(mean = mean(x), min = min(x), max = max(x)))

  # Expand the results into separate columns
  summary_hdist_present <- do.call(data.frame, summary_hdist_present)
  names(summary_hdist_present) <- c("Lambda", "Mean", "Min", "Max")

  summary_hdist_future <- do.call(data.frame, summary_hdist_future)
  names(summary_hdist_future) <- c("Lambda", "Mean", "Min", "Max")

  # Create the plot with ggplot2
  p <- ggplot() +
    # Add the observed range shaded area for present
    geom_ribbon(
      data = summary_hdist_present,
      aes(x = Lambda, ymin = Min, ymax = Max),
      fill = "blue", alpha = 0.2
    ) +
    # Add the line for present
    geom_line(data = summary_hdist_present, aes(x = Lambda, y = Mean), color = "blue", linewidth = 1) +
    geom_point(data = summary_hdist_present, aes(x = Lambda, y = Mean), color = "blue", size = 2) +

    # Add the observed range shaded area for future
    geom_ribbon(
      data = summary_hdist_future,
      aes(x = Lambda, ymin = Min, ymax = Max),
      fill = "red", alpha = 0.2
    ) +
    # Add the line for future
    geom_line(data = summary_hdist_future, aes(x = Lambda, y = Mean), color = "red", linewidth = 1) +
    geom_point(data = summary_hdist_future, aes(x = Lambda, y = Mean), color = "red", size = 2) +

    labs(
      title = "Observed Range of Hellinger Distance Across Iterations",
      x = "Lambda",
      y = "Average Hellinger Distance",
      fill = "Legend"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(size = 16, hjust = 0.5),
      axis.title = element_text(size = 14),
      axis.text = element_text(size = 12),
      legend.position = "top"
    )

  # Save the plot as both PDF and PNG
  name <- "figure/average_hdist_vs_lambda"
  # ggsave(paste0(name, ".pdf"), plot = p, width = 35, height = 25, units = "cm", dpi = 300)
  ggsave(paste0(name, ".png"), plot = p, width = 35, height = 25, units = "cm", dpi = 300)
  p
}

# Maps of the gradients of Hellinger distance
{
  # Initialize a data frame to store average gradients for present and future
  average_gradients <- data.frame(Lambda = numeric(), Average_Gradient_Present = numeric(), Average_Gradient_Future = numeric())

  # Define fixed color scale limits
  limits <- c(0, 0.4)
  v_limits <- signif(seq(limits[1], limits[2], length.out = 5), 2)  # Legend ticks with 2 significant figures

  # Loop through each lambda
  for (lambda in names(GC_hdist_future)) {
    # Create folders for the current lambda
    lambda_value <- as.numeric(sub("lambda_", "", lambda))
    dir.create(file.path("figure", "Gradient_Hellinger", "Present", paste0("lambda_", lambda_value)), recursive = TRUE, showWarnings = FALSE)
    dir.create(file.path("figure", "Gradient_Hellinger", "Future", paste0("lambda_", lambda_value)), recursive = TRUE, showWarnings = FALSE)

    # Initialize data frames to store iteration-wise average gradients for present and future
    iteration_gradients_present <- data.frame(Iteration = integer(), Average_Gradient = numeric())
    iteration_gradients_future <- data.frame(Iteration = integer(), Average_Gradient = numeric())

    # Loop through iterations for the current lambda
    for (iteration in names(GC_hdist_future[[lambda]])) {
      # Compute the gradient of the Hellinger distance for present
      gradient_map_present <- gradient_hdist(GC_hdist[[lambda]][[iteration]])
      avg_gradient_present <- mean(abs(gradient_map_present), na.rm = TRUE)

      # Compute the gradient of the Hellinger distance for future
      gradient_map_future <- gradient_hdist(GC_hdist_future[[lambda]][[iteration]])
      avg_gradient_future <- mean(abs(gradient_map_future), na.rm = TRUE)

      # Store the average gradients for the iteration
      iteration_gradients_present <- rbind(iteration_gradients_present, data.frame(
        Iteration = as.integer(sub("iteration_", "", iteration)),
        Average_Gradient = avg_gradient_present
      ))

      iteration_gradients_future <- rbind(iteration_gradients_future, data.frame(
        Iteration = as.integer(sub("iteration_", "", iteration)),
        Average_Gradient = avg_gradient_future
      ))

      # Prepare data for present visualization
      plot_df_present <- melt(gradient_map_present, varnames = c("lon", "lat"), value.name = "Gradient")
      plot_df_present$lat <- plot_df_present$lat - 90  # Adjust latitude for plotting
      plot_df_present$Gradient[plot_df_present$Gradient > limits[2]] <- limits[2]  # Cap values at the upper limit

      # Create the plot for present
      p_present <- ggplot() +
        geom_tile(data = plot_df_present, aes(x = lon, y = lat, fill = Gradient)) +
        ggtitle("GraphCut Gradient of Hellinger Distance - Present") +
        labs(
          subtitle = paste0(
            "Lambda: ", lambda_value,
            " | Iteration: ", sub("iteration_", "", iteration),
            " | Average Gradient: ", round(avg_gradient_present, 4)
          )
        ) +
        scale_fill_gradientn(
          colors = c("white", "red"),  # Gradient from white to red
          breaks = v_limits,          # Legend breaks with 2 significant figures
          limits = limits             # Fixed color scale limits
        ) +
        borders("world2", colour = "black", lwd = 0.12) +
        scale_x_continuous(expand = c(0, 0)) +
        scale_y_continuous(expand = c(0, 0)) +
        theme(legend.position = "bottom") +
        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
        theme(panel.background = element_blank()) +
        xlab("Longitude") +
        ylab("Latitude") +
        labs(fill = "Gradient \n") +
        theme_bw() +
        theme(
          legend.key.size = unit(1, "cm"),
          legend.key.height = unit(1.4, "cm"),
          legend.key.width = unit(0.4, "cm"),
          legend.title = element_text(size = 16),
          legend.text = element_text(size = 12),
          plot.title = element_text(size = 24),
          plot.subtitle = element_text(size = 20, hjust = 0.5, margin = margin(b = 10)),
          axis.text = element_text(size = 14),
          axis.title = element_text(size = 16)
        ) +
        easy_center_title()

      # Save the plot for present iteration
      file_name_present <- file.path("figure", "Gradient_Hellinger", "Present", paste0("lambda_", lambda_value), paste0("Gradient_Hellinger_Present_Lambda_", lambda_value, "_Iteration_", sub("iteration_", "", iteration)))
      ggsave(paste0(file_name_present, ".png"), plot = p_present, width = 35, height = 25, units = "cm", dpi = 300)

      # Prepare data for future visualization
      plot_df_future <- melt(gradient_map_future, varnames = c("lon", "lat"), value.name = "Gradient")
      plot_df_future$lat <- plot_df_future$lat - 90  # Adjust latitude for plotting
      plot_df_future$Gradient[plot_df_future$Gradient > limits[2]] <- limits[2]  # Cap values at the upper limit

      # Create the plot for future
      p_future <- ggplot() +
        geom_tile(data = plot_df_future, aes(x = lon, y = lat, fill = Gradient)) +
        ggtitle("GraphCut Gradient of Hellinger Distance - Future") +
        labs(
          subtitle = paste0(
            "Lambda: ", lambda_value,
            " | Iteration: ", sub("iteration_", "", iteration),
            " | Average Gradient: ", round(avg_gradient_future, 4)
          )
        ) +
        scale_fill_gradientn(
          colors = c("white", "red"),  # Gradient from white to red
          breaks = v_limits,          # Legend breaks with 2 significant figures
          limits = limits             # Fixed color scale limits
        ) +
        borders("world2", colour = "black", lwd = 0.12) +
        scale_x_continuous(expand = c(0, 0)) +
        scale_y_continuous(expand = c(0, 0)) +
        theme(legend.position = "bottom") +
        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
        theme(panel.background = element_blank()) +
        xlab("Longitude") +
        ylab("Latitude") +
        labs(fill = "Gradient \n") +
        theme_bw() +
        theme(
          legend.key.size = unit(1, "cm"),
          legend.key.height = unit(1.4, "cm"),
          legend.key.width = unit(0.4, "cm"),
          legend.title = element_text(size = 16),
          legend.text = element_text(size = 12),
          plot.title = element_text(size = 24),
          plot.subtitle = element_text(size = 20, hjust = 0.5, margin = margin(b = 10)),
          axis.text = element_text(size = 14),
          axis.title = element_text(size = 16)
        ) +
        easy_center_title()

      # Save the plot for future iteration
      file_name_future <- file.path("figure", "Gradient_Hellinger", "Future", paste0("lambda_", lambda_value), paste0("Gradient_Hellinger_Future_Lambda_", lambda_value, "_Iteration_", sub("iteration_", "", iteration)))
      ggsave(paste0(file_name_future, ".png"), plot = p_future, width = 35, height = 25, units = "cm", dpi = 300)
    }

    # Compute the overall average gradients across iterations for the current lambda
    lambda_avg_gradient_present <- mean(iteration_gradients_present$Average_Gradient)
    lambda_avg_gradient_future <- mean(iteration_gradients_future$Average_Gradient)

    # Store the overall average gradients in the results table
    average_gradients <- rbind(average_gradients, data.frame(
      Lambda = lambda_value,
      Average_Gradient_Present = lambda_avg_gradient_present,
      Average_Gradient_Future = lambda_avg_gradient_future
    ))
  }

  # Save the overall average gradients for later analysis
  write.csv(average_gradients, file = "figure/Gradient_Hellinger/Average_Gradients.csv", row.names = FALSE)

}

# Plot of the average gradients
{
  # Initialize a data frame to store min, max, and mean gradients for present and future
  gradient_stats <- data.frame(
    Lambda = numeric(),
    Mean_Gradient_Present = numeric(),
    Min_Gradient_Present = numeric(),
    Max_Gradient_Present = numeric(),
    Mean_Gradient_Future = numeric(),
    Min_Gradient_Future = numeric(),
    Max_Gradient_Future = numeric()
  )

  # Loop through each lambda to compute statistics for present and future
  for (lambda in names(GC_hdist_future)) {
    # Extract lambda value
    lambda_value <- as.numeric(sub("lambda_", "", lambda))

    # Compute gradients for all iterations for present
    iteration_gradients_present <- sapply(GC_hdist[[lambda]], function(x) mean(abs(gradient_hdist(x)), na.rm = TRUE))

    # Compute gradients for all iterations for future
    iteration_gradients_future <- sapply(GC_hdist_future[[lambda]], function(x) mean(abs(gradient_hdist(x)), na.rm = TRUE))

    # Add statistics to the data frame
    gradient_stats <- rbind(gradient_stats, data.frame(
      Lambda = lambda_value,
      Mean_Gradient_Present = mean(iteration_gradients_present),
      Min_Gradient_Present = min(iteration_gradients_present),
      Max_Gradient_Present = max(iteration_gradients_present),
      Mean_Gradient_Future = mean(iteration_gradients_future),
      Min_Gradient_Future = min(iteration_gradients_future),
      Max_Gradient_Future = max(iteration_gradients_future)
    ))
  }

  # Sort by lambda for consistent plotting
  gradient_stats <- gradient_stats[order(gradient_stats$Lambda), ]

  # Plot average Hellinger gradients for present and future with ranges
  p <- ggplot() +
    # Present
    geom_line(data = gradient_stats, aes(x = Lambda, y = Mean_Gradient_Present), color = "blue", size = 1) +
    geom_point(data = gradient_stats, aes(x = Lambda, y = Mean_Gradient_Present), color = "blue", size = 2) +
    geom_ribbon(data = gradient_stats, aes(x = Lambda, ymin = Min_Gradient_Present, ymax = Max_Gradient_Present), fill = "blue", alpha = 0.2) +
    # Future
    geom_line(data = gradient_stats, aes(x = Lambda, y = Mean_Gradient_Future), color = "red", size = 1) +
    geom_point(data = gradient_stats, aes(x = Lambda, y = Mean_Gradient_Future), color = "red", size = 2) +
    geom_ribbon(data = gradient_stats, aes(x = Lambda, ymin = Min_Gradient_Future, ymax = Max_Gradient_Future), fill = "red", alpha = 0.2) +
    # Labels and theme
    labs(
      title = "Average Hellinger Gradient with Range",
      subtitle = "Blue: Present | Red: Future",
      x = "Lambda",
      y = "Average Hellinger Gradient"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(size = 16, hjust = 0.5),
      plot.subtitle = element_text(size = 14, hjust = 0.5),
      axis.title = element_text(size = 14),
      axis.text = element_text(size = 12)
    )

  # Save the plot
  output_dir <- "figure/Gradient_Hellinger/Average"
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  file_name_pdf <- file.path(output_dir, "Average_Gradient_Hellinger_Present_Future.pdf")
  ggsave(file_name_pdf, plot = p, width = 25, height = 20, units = "cm", dpi = 300)
  file_name_png <- file.path(output_dir, "Average_Gradient_Hellinger_Present_Future.png")
  ggsave(file_name_png, plot = p, width = 25, height = 20, units = "cm", dpi = 300)

  # Print the plot
  print(p)
}

# Plot of the average gradients
{
  # Initialize a data frame to store min, max, and mean gradients for present and future
  gradient_stats <- data.frame(
    Lambda = numeric(),
    Mean_Gradient_Present = numeric(),
    Min_Gradient_Present = numeric(),
    Max_Gradient_Present = numeric(),
    Mean_Gradient_Future = numeric(),
    Min_Gradient_Future = numeric(),
    Max_Gradient_Future = numeric()
  )

  # Initialize a data frame to store model-level gradient stats
  model_gradient_stats <- data.frame(
    Mean_Model_Gradient = numeric(),
    Min_Model_Gradient = numeric(),
    Max_Model_Gradient = numeric()
  )

  # Compute model-level gradients for each model
  for (model_idx in seq_len(dim(h_dist)[3])) {
    present_gradients <- apply(h_dist[, , model_idx], 1:2, gradient_hdist)
    future_gradients <- apply(h_dist_future[, , model_idx], 1:2, gradient_hdist)

    avg_present_gradient <- mean(abs(present_gradients), na.rm = TRUE)
    avg_future_gradient <- mean(abs(future_gradients), na.rm = TRUE)

    # Store min, max, and mean across the model
    model_gradient_stats <- rbind(
      model_gradient_stats,
      data.frame(
        Mean_Model_Gradient = (avg_present_gradient + avg_future_gradient) / 2,
        Min_Model_Gradient = min(c(avg_present_gradient, avg_future_gradient), na.rm = TRUE),
        Max_Model_Gradient = max(c(avg_present_gradient, avg_future_gradient), na.rm = TRUE)
      )
    )
  }

  # Overall stats for the gradient floor
  floor_mean <- mean(model_gradient_stats$Mean_Model_Gradient)
  floor_min <- min(model_gradient_stats$Min_Model_Gradient)
  floor_max <- max(model_gradient_stats$Max_Model_Gradient)

  # Loop through each lambda to compute statistics for present and future
  for (lambda in names(GC_hdist_future)) {
    # Extract lambda value
    lambda_value <- as.numeric(sub("lambda_", "", lambda))

    # Compute gradients for all iterations for present
    iteration_gradients_present <- sapply(GC_hdist[[lambda]], function(x) mean(abs(gradient_hdist(x)), na.rm = TRUE))

    # Compute gradients for all iterations for future
    iteration_gradients_future <- sapply(GC_hdist_future[[lambda]], function(x) mean(abs(gradient_hdist(x)), na.rm = TRUE))

    # Add statistics to the data frame
    gradient_stats <- rbind(gradient_stats, data.frame(
      Lambda = lambda_value,
      Mean_Gradient_Present = mean(iteration_gradients_present),
      Min_Gradient_Present = min(iteration_gradients_present),
      Max_Gradient_Present = max(iteration_gradients_present),
      Mean_Gradient_Future = mean(iteration_gradients_future),
      Min_Gradient_Future = min(iteration_gradients_future),
      Max_Gradient_Future = max(iteration_gradients_future)
    ))
  }

  # Sort by lambda for consistent plotting
  gradient_stats <- gradient_stats[order(gradient_stats$Lambda), ]

  # Plot average Hellinger gradients for present and future with ranges and model-level floor
  p <- ggplot() +
    # Present
    geom_line(data = gradient_stats, aes(x = Lambda, y = Mean_Gradient_Present), color = "blue", size = 1) +
    geom_point(data = gradient_stats, aes(x = Lambda, y = Mean_Gradient_Present), color = "blue", size = 2) +
    geom_ribbon(data = gradient_stats, aes(x = Lambda, ymin = Min_Gradient_Present, ymax = Max_Gradient_Present), fill = "blue", alpha = 0.2) +
    # Future
    geom_line(data = gradient_stats, aes(x = Lambda, y = Mean_Gradient_Future), color = "red", size = 1) +
    geom_point(data = gradient_stats, aes(x = Lambda, y = Mean_Gradient_Future), color = "red", size = 2) +
    geom_ribbon(data = gradient_stats, aes(x = Lambda, ymin = Min_Gradient_Future, ymax = Max_Gradient_Future), fill = "red", alpha = 0.2) +
    # Model-Level Floor Gradient
    geom_line(aes(x = gradient_stats$Lambda, y = floor_mean), color = "black", linetype = "dashed", size = 1) +
    geom_ribbon(aes(x = gradient_stats$Lambda, ymin = floor_min, ymax = floor_max), fill = "grey", alpha = 0.3) +
    # Labels and theme
    labs(
      title = "Average Hellinger Gradient with Range and Floor",
      subtitle = "Blue: Present | Red: Future | Grey: Model-Level Gradient Floor (Min-Max)",
      x = "Lambda",
      y = "Average Hellinger Gradient"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(size = 16, hjust = 0.5),
      plot.subtitle = element_text(size = 14, hjust = 0.5),
      axis.title = element_text(size = 14),
      axis.text = element_text(size = 12)
    )

  # Save the plot
  output_dir <- "figure/Gradient_Hellinger/Average"
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  file_name_pdf <- file.path(output_dir, "Average_Gradient_Hellinger_Present_Future_Model_Floor.pdf")
  ggsave(file_name_pdf, plot = p, width = 25, height = 20, units = "cm", dpi = 300)
  file_name_png <- file.path(output_dir, "Average_Gradient_Hellinger_Present_Future_Model_Floor.png")
  ggsave(file_name_png, plot = p, width = 25, height = 20, units = "cm", dpi = 300)

  # Print the plot
  print(p)
}
