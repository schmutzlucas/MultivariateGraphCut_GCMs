# Install and load necessary libraries
list_of_packages <- read.table("package_list.txt", sep="\n")$V1
new.packages <- list_of_packages[!(list_of_packages %in% installed.packages()[,"Package"])]
if(length(new.packages))
  install.packages(new.packages, repos = "https://cloud.r-project.org")

library(devtools)
lapply(list_of_packages, library, character.only = TRUE)
install_github("schmutzlucas/gcoWrapR")


library(future)
plan(sequential)  # no R process parallelism, just native threads



# Loading local functions
source_code_dir <- 'functions/' #The directory where all functions are saved.
file_paths <- list.files(source_code_dir, full.names = T)
for(path in file_paths){source(path)}

load("workspaces/202511221241_workspace_multiRes_3v_gc_results.RData")

# Setting global variables
lon <- -180:179
lat <- -90:90
lon_size <- length(lon)
lat_size <- length(lat)
# Temporal ranges
year_present <<- 1950:1975
year_future <<- 1999:2024
# data directory
data_dir <<- 'data/CMIP6_summer_Apr15-Oct14'


# List of the variable used
variables <- c('pr', 'tas', 'psl')

## 1.  bin-resolution choices
nbins3d <- 8     # 3-D joint
nbins2d <- 16    # 2-D pairs
nbins1d <- 32    # 1-D marginals

# Vector of lambda values (smooth costs)
smooth_vals <- c(
  0.025, 0.05, 0.075, 0.10, 0.125,
  0.15,  0.20, 0.25, 0.30, 0.40, 0.50
)

# List to store GC_result for each lambda
GC_results_lambda <- vector("list", length(smooth_vals))
names(GC_results_lambda) <- paste0("lambda_", smooth_vals)

# (optional) simple status table you can fill later
GC_lambda_status <- data.frame(
  lambda       = smooth_vals,
  success      = NA,        # TRUE/FALSE
  message      = NA_character_,
  stringsAsFactors = FALSE
)

# Fixed seed
seed_gc <- 1L

for (i in seq_along(smooth_vals)) {
  smooth_cost <- smooth_vals[i]
  cat("\n----------------------------------------------------\n")
  cat("Running GraphCut for lambda (smooth_cost) =", smooth_cost,
      "with seed =", seed_gc, "\n")

  GC_result_i <- tryCatch({
    GraphCutHellinger_nD_lat(
      pdf_models_future = pdf3_models_fut,  # 3-D PDFs for labeling
      h_dist            = h_dist_pres,      # datacost (present Hellinger)
      weight_data       = 1,
      weight_smooth     = smooth_cost,
      nBins             = nbins_total3d,
      lat               = lat,
      seed              = seed_gc,
      verbose           = TRUE,
      rebuild           = FALSE  # IMPORTANT: assume C++ already compiled once
    )
  }, error = function(e) {
    cat("  ⚠ GraphCut failed at lambda =", smooth_cost, ":\n",
        e$message, "\n")
    GC_lambda_status$success[i] <- FALSE
    GC_lambda_status$message[i] <- e$message
    return(NULL)
  })

  if (!is.null(GC_result_i)) {
    GC_results_lambda[[i]]    <- GC_result_i
    GC_lambda_status$success[i] <- TRUE
    GC_lambda_status$message[i] <- ""
  }

  gc()
}


# ======================================================================
# SAVE all GC_results after lambda sweep
# ======================================================================

saveRDS(
  GC_results_lambda,
  file = "GC_results_lambda_seed1.rds"
)

cat("\nSaved GC_results_lambda to: GC_results_lambda_seed1.rds\n")

# =====================================================================
#  Build all h_dist maps (future + present) for all lambda values
#  from stored GraphCut results (GC_results_lambda)
# =====================================================================

GC_maps <- vector("list", length(smooth_vals))
names(GC_maps) <- paste0("lambda_", smooth_vals)

for (i in seq_along(smooth_vals)) {

  lam <- smooth_vals[i]
  GC_result <- GC_results_lambda[[i]]

  cat("\nComputing H-dist maps for lambda =", lam, "\n")

  if (is.null(GC_result)) {
    GC_maps[[i]] <- NULL
    next
  }

  # ---------------------------------------------------------------
  # Initialize lon × lat matrices
  # ---------------------------------------------------------------
  GC_hdist_pres <- matrix(NA_real_, nrow = length(lon), ncol = length(lat))
  GC_hdist_fut  <- matrix(NA_real_, nrow = length(lon), ncol = length(lat))

  # ---------------------------------------------------------------
  # Fill values based on which model was selected at each pixel
  # ---------------------------------------------------------------
  for (l in seq_along(model_names)) {

    islabel <- which(GC_result$label_attribution == l)

    if (length(islabel) == 0L) next

    GC_hdist_pres[islabel] <- h_dist_pres[ , , l][islabel]
    GC_hdist_fut [islabel] <- h_dist_fut [ , , l][islabel]
  }

  # ---------------------------------------------------------------
  # Store maps + summary metrics
  # ---------------------------------------------------------------
  GC_maps[[i]] <- list(
    lambda         = lam,
    GC_hdist_pres  = GC_hdist_pres,
    GC_hdist_fut   = GC_hdist_fut,
    mean_pres      = mean(GC_hdist_pres, na.rm = TRUE),
    mean_fut       = mean(GC_hdist_fut, na.rm = TRUE)
  )

  gc()
}

cat("\nAll H-dist maps computed and stored in GC_maps\n")



# ======================================================================
#  Loop over all lambda values and save labelling maps to a folder
# ======================================================================

# Folder where to store all maps
out_dir <- "figure/Labelling/lambda_sweep_seed1"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

for (i in seq_along(smooth_vals)) {

  lam <- smooth_vals[i]
  GC_result <- GC_results_lambda[[i]]

  if (is.null(GC_result)) {
    message("Skipping lambda = ", lam, " (no GC_result)")
    next
  }

  message("Producing map for lambda = ", lam)

  # --------------------------------------------------------------------
  # Labelling Map
  # --------------------------------------------------------------------
  GC_labels <- GC_result$label_attribution  # [lon_idx, lat_idx] matrix

  # Convert the label matrix to a data frame for plotting
  label_df <- reshape2::melt(
    GC_labels,
    varnames  = c("lon_idx", "lat_idx"),
    value.name = "label_attribution"
  )

  # Map indices to actual lon/lat values
  label_df$lon <- lon[label_df$lon_idx]
  label_df$lat <- lat[label_df$lat_idx]

  # Ensure labels are factors aligned with model_names
  label_df$label_attribution <- factor(
    label_df$label_attribution,
    levels = seq_along(model_names),
    labels = model_names
  )

  # Color palette (named by model_names)
  color_palette <- c(
    "#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd",
    "#8c564b", "#e377c2", "#7f7f7f", "#bcbd22", "#17becf",
    "#aec7e8", "#ffbb78", "#98df8a", "#ff9896", "#c5b0d5",
    "#c49c94", "#f7b6d2", "#c7c7c7", "#dbdb8d", "#9edae5",
    "#393b79", "#5254a3", "#6b6ecf"
  )
  names(color_palette) <- model_names

  p6 <- ggplot() +
    geom_tile(
      data = label_df,
      aes(x = lon, y = lat, fill = label_attribution)
    ) +
    scale_fill_manual(
      values = color_palette,
      na.value = "white",
      guide = guide_legend(title = "Model Names", ncol = 1)
    ) +
    ggtitle(sprintf("Label GC Hellinger - Lambda: %.3f (seed = 1)", lam)) +
    borders("world", colour = "black", size = 0.12) +
    theme_bw() +
    coord_fixed(
      ratio = 1.3,
      xlim  = c(-180, 180),
      ylim  = c(-84, 90),
      expand = FALSE
    ) +
    theme(
      legend.position   = "right",
      panel.grid.major  = element_blank(),
      panel.grid.minor  = element_blank(),
      panel.background  = element_blank(),
      legend.key.size   = unit(0.5, "cm"),
      legend.key.height = unit(0.5, "cm"),
      legend.key.width  = unit(0.5, "cm"),
      legend.title      = element_text(size = 10),
      legend.text       = element_text(size = 8),
      plot.title        = element_text(size = 16),
      plot.subtitle     = element_text(size = 12, hjust = 0.5),
      axis.text         = element_text(size = 10),
      axis.title        = element_text(size = 12)
    ) +
    xlab("Longitude") +
    ylab("Latitude") +
    easy_center_title()

  # Save to PDF + PNG in the target folder
  base_name <- sprintf("%s/Labelling_GC_lambda%.3f_seed1", out_dir, lam)

  ggsave(paste0(base_name, ".pdf"),
         plot   = p6,
         width  = 20,
         height = 15,
         units  = "cm",
         dpi    = 300)

  ggsave(paste0(base_name, ".png"),
         plot   = p6,
         width  = 20,
         height = 15,
         units  = "cm",
         dpi    = 300)
}


# ======================================================================
# Produce Hellinger (hdist3) maps for all lambda values
# GC_maps[[i]] must contain:
#   - GC_hdist_fut   (lon × lat matrix)
#   - lambda         (scalar)
# ======================================================================

out_dir <- "figure/Hellinger3_lambda_sweep_seed1"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

for (i in seq_along(GC_maps)) {

  entry <- GC_maps[[i]]

  if (is.null(entry)) {
    message("Skipping lambda index ", i, " (NULL GC_maps entry)")
    next
  }

  lam <- entry$lambda
  GC_hdist_fut <- entry$GC_hdist_fut

  message("Producing H-dist3 map for λ = ", lam)

  # --------------------------------------------------------------------
  # 1) Melt to dataframe
  # --------------------------------------------------------------------
  test_df <- reshape2::melt(
    GC_hdist_fut,
    varnames = c("lon_idx", "lat_idx"),
    value.name = "H_dist"
  )

  test_df$lon <- lon[test_df$lon_idx]
  test_df$lat <- lat[test_df$lat_idx]

  # Wrap longitude to [-180, 180]
  test_df$lon_wrapped <- ifelse(test_df$lon > 180, test_df$lon - 360, test_df$lon)

  # --------------------------------------------------------------------
  # 2) Compute global weighted mean H-dist
  # --------------------------------------------------------------------
  cos_weights <- cos(lat * pi/180)
  weight_matrix <- matrix(rep(cos_weights, each = length(lon)),
                          nrow = length(lon), byrow = FALSE)

  global_h <- sum(GC_hdist_fut * weight_matrix, na.rm = TRUE) /
    sum(weight_matrix, na.rm = TRUE)

  # --------------------------------------------------------------------
  # 3) Plot settings
  # --------------------------------------------------------------------
  limit <- 0.5
  limits <- c(0.1, limit)
  v_limits <- seq(limits[1], limits[2], length.out = 3)

  p_hdist <- ggplot() +
    geom_tile(
      data = test_df,
      aes(x = lon_wrapped, y = lat, fill = H_dist)
    ) +
    labs(subtitle = "Projection period: 1998–2023") +
    ggtitle(sprintf("GraphCut (λ = %.3f) — Mean H = %.3f", lam, global_h)) +
    scale_fill_gradient(
      low = "white",
      high = "#015a8c",
      limits = limits,
      breaks = v_limits,
      oob = scales::squish
    ) +
    borders("world", colour = "black", linewidth = 0.12) +
    scale_x_continuous(expand = c(0, 0)) +
    scale_y_continuous(expand = c(0, 0)) +
    coord_fixed(ratio = 1.3,
                xlim = c(-180, 180), ylim = c(-84, 90), expand = FALSE) +
    theme_bw() +
    theme(
      legend.position   = "right",
      panel.grid.major  = element_blank(),
      panel.grid.minor  = element_blank(),
      panel.background  = element_blank(),
      legend.key.size   = unit(1, "cm"),
      legend.key.height = unit(1.4, "cm"),
      legend.key.width  = unit(0.4, "cm"),
      legend.title      = element_text(size = 16),
      legend.text       = element_text(size = 12),
      plot.title        = element_text(size = 24),
      plot.subtitle     = element_text(size = 20, hjust = 0.5),
      axis.text         = element_text(size = 14),
      axis.title        = element_text(size = 16)
    ) +
    xlab("Longitude") +
    ylab("Latitude") +
    labs(fill = "H") +
    easy_center_title()

  print(p_hdist)

  # --------------------------------------------------------------------
  # 4) Save for this λ
  # --------------------------------------------------------------------
  base <- sprintf("%s/hdist3_GC_lambda%.3f_seed1", out_dir, lam)

  ggsave(paste0(base, ".pdf"),
         plot = p_hdist, width = 20, height = 15, units = "cm", dpi = 300)

  ggsave(paste0(base, ".png"),
         plot = p_hdist, width = 20, height = 15, units = "cm", dpi = 300)
}



# ======================================================================
# Mean H_dist (future) vs lambda — two curves:
#   1) mean H_dist
#   2) latitude-weighted mean H_dist
# ======================================================================

# Select only valid entries
valid_idx <- which(!sapply(GC_maps, is.null))

# Extract lambda, mean_fut
lambda_vec <- vapply(GC_maps[valid_idx], function(x) x$lambda, numeric(1))
mean_fut   <- vapply(GC_maps[valid_idx], function(x) x$mean_fut, numeric(1))

# ----------------------------------------------------------------------
# Compute latitude-weighted H for each lambda (as in your hdist3 plots)
# ----------------------------------------------------------------------

weighted_fut <- numeric(length(valid_idx))

cos_weights <- cos(lat * pi/180)
weight_matrix <- matrix(
  rep(cos_weights, each = length(lon)),
  nrow = length(lon),
  byrow = FALSE
)

for (i in seq_along(valid_idx)) {

  GC_hdist_fut <- GC_maps[[ valid_idx[i] ]]$GC_hdist_fut

  weighted_fut[i] <-
    sum(GC_hdist_fut * weight_matrix, na.rm = TRUE) /
      sum(weight_matrix, na.rm = TRUE)
}

# ----------------------------------------------------------------------
# Combine into a single dataframe for ggplot
# ----------------------------------------------------------------------

df_lambda <- data.frame(
  lambda          = lambda_vec,
  mean_hdist_fut  = mean_fut,
  weighed_hdist_f = weighted_fut
)

# Sort by lambda
df_lambda <- df_lambda[order(df_lambda$lambda), ]

# ----------------------------------------------------------------------
# Plot
# ----------------------------------------------------------------------

out_dir <- "figure/Labelling/lambda_sweep_seed1"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

p_lambda <- ggplot(df_lambda, aes(x = lambda)) +
  geom_line(aes(y = mean_hdist_fut, colour = "Mean H_dist")) +
  geom_point(aes(y = mean_hdist_fut, colour = "Mean H_dist"), size = 2) +
  geom_line(aes(y = weighed_hdist_f, colour = "Weighted mean H_dist")) +
  geom_point(aes(y = weighed_hdist_f, colour = "Weighted mean H_dist"), size = 2) +
  scale_colour_manual(
    values = c("Mean H_dist" = "black", "Weighted mean H_dist" = "red"),
    name   = NULL
  ) +
  scale_x_continuous(
    breaks = df_lambda$lambda,
    labels = format(df_lambda$lambda, digits = 3)
  ) +
  xlab("Lambda (smoothness weight)") +
  ylab("Hellinger distance (future)") +
  ggtitle("Mean and Latitude-weighted Hellinger distance vs lambda\n(seed = 1)") +
  theme_bw() +
  theme(
    panel.grid.minor = element_blank(),
    axis.text.x      = element_text(angle = 45, hjust = 1),
    plot.title       = element_text(size = 14),
    legend.position  = "right"
  )

print(p_lambda)

# Save
base_name <- file.path(out_dir, "Mean_and_Weighted_Hdist_vs_lambda_seed1")
ggsave(paste0(base_name, ".pdf"), plot = p_lambda,
       width = 18, height = 10, units = "cm", dpi = 300)
ggsave(paste0(base_name, ".png"), plot = p_lambda,
       width = 18, height = 10, units = "cm", dpi = 300)

