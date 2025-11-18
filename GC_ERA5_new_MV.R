# Install and load necessary libraries
list_of_packages <- read.table("package_list.txt", sep="\n")$V1
new.packages <- list_of_packages[!(list_of_packages %in% installed.packages()[,"Package"])]
if(length(new.packages))
  install.packages(new.packages, repos = "https://cloud.r-project.org")

library(devtools)
lapply(list_of_packages, library, character.only = TRUE)
install_github("schmutzlucas/gcoWrapR")


# Loading local functions
source_code_dir <- 'functions/' #The directory where all functions are saved.
file_paths <- list.files(source_code_dir, full.names = T)
for(path in file_paths){source(path)}


range_var_final <- readRDS('ranges/range_var_final_ERA5_1950-2023_3v.rds')

# ------------------------------------------------------------------
# A. build permutation that converts 0…359 → -180…+179 order
lon_file        <- 0:359
lon_adj         <- ifelse(lon_file >= 180, lon_file - 360, lon_file)
lon_order_adj   <- order(lon_adj)                # length 360

# ------------------------------------------------------------------
# B. re-order every variable’s range matrix
rng_list <- lapply(range_var_final$ranges, \(mat) {
  mat[ lon_order_adj, , , drop = FALSE ]         # keep 3-dim structure
})

# ------------------------------------------------------------------
# C. bind into one 4-D [lon, lat, var, min/max] array
vars <- names(rng_list)               # "pr" "tas" "psl"
range_var_gc <- array(
  NA_real_,
  dim = c(360, 181, length(vars), 2),
  dimnames = list(lon = -180:179, lat = -90:90, var = vars, bound = c("min","max"))
)


for (k in seq_along(vars)) {
  range_var_gc[ , , k, ] <- rng_list[[k]]
}

# Setting global variables
lon <- -180:179
lat <- -90:90
lon_size <- length(lon)
lat_size <- length(lat)
# Temporal ranges
year_present <<- 1950:1975
year_future <<- 1998:2024
# data directory
data_dir <<- 'data/CMIP6_summer_Apr15-Oct14'


# List of the variable used
variables <- c('pr', 'tas', 'psl')

## 1.  bin-resolution choices
nbins3d <- 8     # 3-D joint
nbins2d <- 16    # 2-D pairs
nbins1d <- 32    # 1-D marginals

## 2.  model list
model_names <- scan("model_names_pr_tas_psl.txt", what = "", quiet = TRUE)

## 3.  number of parallel workers
workers <- 4   # adapt to your machine

## 4.  call the multi-resolution histogram builder
cat("→ building PDFs and means …\n")
t0 <- Sys.time()
tmp <- compute_nd_pdf_multi(
  variables, model_names,
  data_dir,
  year_present, year_future,
  lon, lat,
  range_var_gc,
  nbins3d = nbins3d,
  nbins2d = nbins2d,
  nbins1d = nbins1d,
  workers  = workers)
cat("   done in", round(difftime(Sys.time(), t0, units = "mins"), 2), "minutes\n")
gc()

# --------------------------------------------------------------------
#  (optional) choose a reference model now, just like before
# --------------------------------------------------------------------
ref_index      <- 1                       # or any other index
reference_name <- model_names[ref_index]
model_names    <- model_names[-ref_index]

pdf3_present <- tmp$pdf3$present          # 3-D joint PDFs
pdf3_future  <- tmp$pdf3$future

pdf2_present <- tmp$pdf2$present          # list of 3 two-D vectors
pdf2_future  <- tmp$pdf2$future

pdf1_present <- tmp$pdf1$present          # list of 3 one-D vectors
pdf1_future  <- tmp$pdf1$future

mean_present <- tmp$mean$present          # cell-wise means
mean_future  <- tmp$mean$future

# drop the tmp list if memory is tight
rm(tmp)
gc()

# --------------------------------------------------------------------
#  Save the workspace
# --------------------------------------------------------------------
stamp    <- format(Sys.time(), "%Y%m%d%H%M")
filename <- file.path("workspaces", paste0(stamp, "_workspace_multiRes_3v.RData"))
save.image(file = filename, compress = FALSE)
cat("✓ workspace saved to", filename, "\n")


# --------------------------------------------------------------------
#  1) Split out reference vs. models
# --------------------------------------------------------------------
# assume: pdf3_present/future and ref_index are in scope
pdf3_ref_pres    <- pdf3_present[,,, ref_index]
pdf3_models_pres <- pdf3_present[,,,-ref_index]

pdf3_ref_fut    <- pdf3_future[,,, ref_index]
pdf3_models_fut <- pdf3_future[,,,-ref_index]

rm(pdf3_future, pdf3_present)

# --------------------------------------------------------------------
#  2) build “all‐bins” index list for full Hellinger
# --------------------------------------------------------------------
nbins_total3d <- nbins3d^3
selected_indices_all <- lapply(seq_len(lon_size), function(i) {
  lapply(seq_len(lat_size), function(j) {
    seq_len(nbins_total3d)
  })
})

# --------------------------------------------------------------------
#  3) compute full Hellinger distance maps
# --------------------------------------------------------------------
h_dist_pres <- compute_partial_hdist(
  pdf_ref      = pdf3_ref_pres,
  pdf_model    = pdf3_models_pres,
  selected_indices = selected_indices_all
)

h_dist_fut  <- compute_partial_hdist(
  pdf_ref      = pdf3_ref_fut,
  pdf_model    = pdf3_models_fut,
  selected_indices = selected_indices_all
)

# --------------------------------------------------------------------
#  4) replace NaN with zero
# --------------------------------------------------------------------
h_dist_pres[is.nan(h_dist_pres)] <- 0
h_dist_fut [is.nan(h_dist_fut )] <- 0

# quick check
hist(h_dist_pres, main="H-dist Present (all bins)")
hist(h_dist_fut,  main="H-dist Future  (all bins)")

# --------------------------------------------------------------------
#  5) run GraphCut on the 3-D cost map
# --------------------------------------------------------------------
smooth_cost <- 0.1

GC_result <- tryCatch({
  GraphCutHellinger_nD_lat(
    pdf_models_future = pdf3_models_fut,   # using “present” PDFs for labeling
    h_dist            = h_dist_pres,        # datacost = Hellinger(pres)
    weight_data       = 1,
    weight_smooth     = smooth_cost,
    nBins             = nbins_total3d,
    lat               = lat,
    seed              = 1,
    verbose           = TRUE,
    rebuild           = TRUE
  )
}, error = function(e) {
  cat("⚠️  GraphCut failed at smooth_cost =", smooth_cost, ":\n", e$message, "\n")
  NULL
})
gc()

# Compute the pdf3 for MMM and GC | Compute the hdist_3
{
  # Initialize arrays for GraphCut PDFs
  pdf3_GC_pres <- array(NA, dim = c(length(lon), length(lat), dim(pdf3_models_pres)[3]))
  pdf3_GC_fut <- array(NA, dim = c(length(lon), length(lat), dim(pdf3_models_fut)[3]))

  # Assign PDFs based on GraphCut label attribution
  for (i in seq_along(lon)) {
    for (j in seq_along(lat)) {
      # Get the GraphCut-selected label (indexing in R is 1-based)
      label <- GC_result$label_attribution[i, j]  # [lat, lon] ordering!

      # Assign the PDF of the selected model
      pdf3_GC_pres[i, j, ] <- pdf3_models_pres[i, j, , label]
      pdf3_GC_fut[i, j, ] <- pdf3_models_fut[i, j, , label]
    }
  }


  # Initialize a lon x lat matrix for each smooth cost
  GC_hdist_pres <- matrix(NA, nrow = length(lon), ncol = length(lat))
  GC_hdist_fut <- matrix(NA, nrow = length(lon), ncol = length(lat))

  for(l in 1:(length(model_names))){  # Ensure that indexing aligns with model names
    islabel <- which(GC_result$label_attribution == l)
    GC_hdist_pres[islabel] <- h_dist_pres[,,l][islabel]
    GC_hdist_fut[islabel] <- h_dist_fut[,, l][islabel]
  }


  hist(GC_hdist_fut)


  mean(GC_hdist_fut)

  # Compute Multi-Model Mean for Present
  pdf3_MMM_pres <- apply(pdf3_models_pres, c(1, 2, 3), mean)

  # Compute Multi-Model Mean for Future
  pdf3_MMM_fut <- apply(pdf3_models_fut, c(1, 2, 3), mean)


  # Initialize arrays to store the Hellinger distance for MMM
  MMM_hdist_pres <- array(NA, dim = c(length(lon), length(lat)))
  MMM_hdist_fut <- array(NA, dim = c(length(lon), length(lat)))

  # Compute Hellinger distance for the Multi-Model Mean
  for (i in seq_along(lon)) {
    for (j in seq_along(lat)) {
      # Compute Hellinger distance for present
      MMM_hdist_pres[i, j] <- sqrt(sum((sqrt(pdf3_MMM_pres[i, j, ]) - sqrt(pdf3_ref_pres[i, j, ]))^2)) / sqrt(2)

      # Compute Hellinger distance for future
      MMM_hdist_fut[i, j] <- sqrt(sum((sqrt(pdf3_MMM_fut[i, j, ]) - sqrt(pdf3_ref_fut[i, j, ]))^2)) / sqrt(2)
    }
  }

  # Replace NaN values with 0
  MMM_hdist_pres <- replace(MMM_hdist_pres, is.nan(MMM_hdist_pres), 0)
  MMM_hdist_fut <- replace(MMM_hdist_fut, is.nan(MMM_hdist_fut), 0)

  # Visualize or analyze the results
  hist(MMM_hdist_pres)
  hist(GC_hdist_pres)
  hist(MMM_hdist_fut)
  hist(GC_hdist_fut)

  mean(MMM_hdist_fut)
  mean(GC_hdist_fut)
  gc()

}
# --------------------------------------------------------------------
#  Save the workspace
# --------------------------------------------------------------------
stamp    <- format(Sys.time(), "%Y%m%d%H%M")
filename <- paste0(stamp, "_workspace_multiRes_3v_gc_results.RData")
save.image(file = filename, compress = FALSE)
cat("✓ workspace saved to", filename, "\n")
#
# # Labelling Map
# {
#   # Extract the label attribution for the current smooth cost
#   GC_labels <- GC_result$label_attribution
#
#   # Convert the label matrix to a data frame for plotting
#   label_df <- reshape2::melt(GC_labels, varnames = c("lon_idx", "lat_idx"), value.name = "label_attribution")
#
#   # Explicitly assign correct longitude and latitude values
#   label_df$lon <- lon[label_df$lon_idx]  # Map longitude indices to values
#   label_df$lat <- lat[label_df$lat_idx]    # Map latitude indices to values
#
#   # Force label_attribution to be a factor with levels in the exact order of model_names.
#   label_df$label_attribution <- factor(label_df$label_attribution,
#                                        levels = seq_along(model_names),
#                                        labels = model_names)
#   # Create a named color palette: each model name is explicitly mapped to its color.
#   color_palette <- c(
#     "#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd",
#     "#8c564b", "#e377c2", "#7f7f7f", "#bcbd22", "#17becf",
#     "#aec7e8", "#ffbb78", "#98df8a", "#ff9896", "#c5b0d5",
#     "#c49c94", "#f7b6d2", "#c7c7c7", "#dbdb8d", "#9edae5",
#     "#393b79", "#5254a3", "#6b6ecf"
#   )
#   # Name the palette vector with model_names (in the same order)
#   names(color_palette) <- model_names
#
#   p6 <- ggplot() +
#     geom_tile(data = label_df, aes(x = lon, y = lat, fill = label_attribution)) +
#     scale_fill_manual(
#       values = color_palette,
#       na.value = "white",
#       guide = guide_legend(title = "Model Names", ncol = 1)
#     ) +
#     ggtitle(paste("Label GC Hellinger - Lambda:", 0.1)) +
#     borders("world", colour = 'black', size = 0.12) +
#     theme_bw() +
#     coord_fixed(ratio = 1.3, xlim = c(-180, 180), ylim = c(-84, 90), expand = FALSE)+
#     theme(
#       legend.position = 'right',
#       panel.grid.major = element_blank(),
#       panel.grid.minor = element_blank(),
#       panel.background = element_blank(),
#       legend.key.size = unit(0.5, 'cm'),
#       legend.key.height = unit(0.5, 'cm'),
#       legend.key.width = unit(0.5, 'cm'),
#       legend.title = element_text(size = 10),
#       legend.text = element_text(size = 8),
#       plot.title = element_text(size = 16),
#       plot.subtitle = element_text(size = 12, hjust = 0.5),
#       axis.text = element_text(size = 10),
#       axis.title = element_text(size = 12)
#     ) +
#     xlab('Longitude') +
#     ylab('Latitude') +
#     easy_center_title()
#
#   p6
#
#   # 9) Save
#   name <- paste0("figure/Labelling/Labelling_GC_0.1_seed1_2")
#   ggsave(paste0(name, ".pdf"), plot = p6, width = 20, height = 15, units = "cm", dpi = 300)
#   ggsave(paste0(name, ".png"), plot = p6, width = 20, height = 15, units = "cm", dpi = 300)
# }
#
#
# # Single point PDF 3d
# {
#   library(plotly)
#
#   # Example grid point indices
#   lon_index <- 180 + 64 # Example longitude index
#   lat_index <- 90 # Example latitude index
#
#   # Extract the 512-bin PDF vector for the specific grid point
#   pdf_vector <- pdf3_models_fut[lon_index, lat_index, ,2]
#
#   # Reshape the PDF vector into a 3D array of dimensions [8, 8, 8]
#   nbins <- 8
#   pdf_3d <- array(pdf_vector, dim = c(nbins, nbins, nbins))
#
#   # Extract the ranges for the variables from range_var_final
#   range_var <- range_var_final$ranges
#   var1_min <- range_var$pr[lon_index, lat_index, 1]
#   var1_max <- range_var$pr[lon_index, lat_index, 2]
#   var2_min <- range_var$tas[lon_index, lat_index, 1]
#   var2_max <- range_var$tas[lon_index, lat_index, 2]
#   var3_min <- range_var$psl[lon_index, lat_index, 1]
#   var3_max <- range_var$psl[lon_index, lat_index, 2]
#
#   # Create bin edges for each variable
#   x_bins <- seq(var1_min, var1_max, length.out = nbins + 1)
#   y_bins <- seq(var2_min, var2_max, length.out = nbins + 1)
#   z_bins <- seq(var3_min, var3_max, length.out = nbins + 1)
#
#   # Create the coordinates for the centers of the bins
#   x_centers <- (x_bins[-1] + x_bins[-length(x_bins)]) / 2
#   y_centers <- (y_bins[-1] + y_bins[-length(y_bins)]) / 2
#   z_centers <- (z_bins[-1] + z_bins[-length(z_bins)]) / 2
#
#   # Expand the grid of coordinates
#   grid <- expand.grid(x = x_centers, y = y_centers, z = z_centers)
#
#   # Flatten the PDF array into a vector
#   pdf_flat <- as.vector(pdf_3d)
#
#   # Combine the coordinates with the PDF values
#   plot_data <- data.frame(
#     x = grid$x,
#     y = grid$y,
#     z = grid$z,
#     value = pdf_flat
#   )
#
#   # Normalize PDF values for marker size
#   normalized_pdf <- pdf_flat / max(pdf_flat, na.rm = TRUE)
#
#   # Plot the 3D histogram
#   fig <- plot_ly(
#     data = plot_data,
#     x = ~x,
#     y = ~y,
#     z = ~z,
#     type = "scatter3d",
#     mode = "markers",
#     marker = list(
#       size = ~normalized_pdf * 75,  # Adjust size scaling factor
#       color = ~value,
#       colorscale = "Viridis",
#       showscale = TRUE
#     ),
#     text = ~paste("PDF Value:", round(value, 4))
#   ) %>%
#     layout(
#       scene = list(
#         xaxis = list(title = "Variable 1 (pr)"),
#         yaxis = list(title = "Variable 2 (tas)"),
#         zaxis = list(title = "Variable 3 (psl)")
#       ),
#       title = paste("3D PDF for Grid Point (Lon:", lon_index, ", Lat:", lat_index, ")")
#     )
#
#   # Show the plot
#   fig
# }
#
#
# # Map of tas bias for all models :
# for (i in seq_along(model_names)) {
#
#   # 1) Compute bias
#   bias_tas <- mean_future$tas[,,i+1] - mean_future$tas[,,1]  # because index 1 is ref, models start at 2
#
#   # 2) Melt to dataframe
#   test_df <- melt(bias_tas, varnames = c("lon_idx", "lat_idx"), value.name = "Bias")
#   test_df$lon <- lon[test_df$lon_idx]
#   test_df$lat <- lat[test_df$lat_idx]
#
#   # 3) Apply limits
#   limit <- 8
#   limits <- c(-8, 8)
#   v_limits <- as.numeric(format(seq(-limit, limit, len=6), digits = 3))
#   test_df$Bias[test_df$Bias < -limit] <- -limit
#   test_df$Bias[test_df$Bias > limit] <- limit
#
#   # 4) Wrap longitude
#   test_df$lon_wrapped <- ifelse(test_df$lon > 180, test_df$lon - 360, test_df$lon)
#
#   # 5) Cosine latitude weights for *absolute* global mean bias
#   cos_weights <- cos(lat * pi/180)
#   weight_matrix <- matrix(rep(cos_weights, each = length(lon)), nrow = length(lon), byrow = FALSE)
#   global_bias <- sum(abs(bias_tas) * weight_matrix, na.rm = TRUE) / sum(weight_matrix, na.rm = TRUE)
#
#   # 6) Model name for the title
#   model_name <- model_names[i]
#
#   # 7) Plot
#   p_bias <- ggplot() +
#     geom_tile(data = test_df, aes(x = lon_wrapped, y = lat, fill = Bias)) +
#     labs(subtitle = 'Projection period: 1998 - 2023') +
#     ggtitle(paste0(model_name, ': Mean |Bias| = ', round(global_bias, 3), ' K')) +
#     scale_fill_gradientn(
#       colours = rev(brewer.pal(11, "RdBu")),
#       breaks = v_limits,
#       limits = limits,
#       oob = scales::squish
#     ) +
#     borders("world", colour = 'black', lwd = 0.12) +
#     scale_x_continuous(expand = c(0, 0)) +
#     scale_y_continuous(expand = c(0, 0)) +
#     coord_fixed(ratio = 1.3, xlim = c(-180, 180), ylim = c(-84, 90), expand = FALSE) +
#     theme_bw() +
#     theme(
#       legend.position = 'right',  # Legend to the right
#       panel.grid.major = element_blank(),
#       panel.grid.minor = element_blank(),
#       panel.background = element_blank(),
#       legend.key.size = unit(1, 'cm'),
#       legend.key.height = unit(1.4, 'cm'),
#       legend.key.width = unit(0.4, 'cm'),
#       legend.title = element_text(size = 16),
#       legend.text = element_text(size = 12),
#       plot.title = element_text(size = 24),
#       plot.subtitle = element_text(size = 20, hjust = 0.5),
#       axis.text = element_text(size = 14),
#       axis.title = element_text(size = 16)
#     ) +
#     xlab('Longitude') +
#     ylab('Latitude') +
#     labs(fill = '[K]') +
#     easy_center_title()
#
#
#   # 8) Save
#   name <- paste0("figure/Bias2/tas/Bias_tas_", model_name)
#   ggsave(paste0(name, ".pdf"), plot = p_bias, width = 20, height = 15, units = "cm", dpi = 300)
#   ggsave(paste0(name, ".png"), plot = p_bias, width = 20, height = 15, units = "cm", dpi = 300)
#
# }
#
#
# # Map of psl bias for all models
# for (i in seq_along(model_names)) {
#
#   # 1) Compute bias
#   bias_psl <- mean_future$psl[,,i+1] - mean_future$psl[,,1]  # index 1 is ref, models start at 2
#
#   # 2) Melt to dataframe
#   test_df <- melt(bias_psl, varnames = c("lon_idx", "lat_idx"), value.name = "Bias")
#   test_df$lon <- lon[test_df$lon_idx]
#   test_df$lat <- lat[test_df$lat_idx]
#
#   # 3) Apply limits
#   limit <- 1500  # Based on your max (2230) → rounded up to 2500 Pa
#   limits <- c(-limit, limit)
#   v_limits <- as.numeric(format(seq(-limit, limit, len=6), digits = 3))
#   test_df$Bias[test_df$Bias < -limit] <- -limit
#   test_df$Bias[test_df$Bias > limit] <- limit
#
#   # 4) Wrap longitude
#   test_df$lon_wrapped <- ifelse(test_df$lon > 180, test_df$lon - 360, test_df$lon)
#
#   # 5) Cosine latitude weights for *absolute* global mean bias
#   cos_weights <- cos(lat * pi/180)
#   weight_matrix <- matrix(rep(cos_weights, each = length(lon)), nrow = length(lon), byrow = FALSE)
#   global_bias <- sum(abs(bias_psl) * weight_matrix, na.rm = TRUE) / sum(weight_matrix, na.rm = TRUE)
#
#   # 6) Model name for the title
#   model_name <- model_names[i]
#
#   # 7) Plot
#   p_bias <- ggplot() +
#     geom_tile(data = test_df, aes(x = lon_wrapped, y = lat, fill = Bias)) +
#     labs(subtitle = 'Projection period: 1998 - 2023') +
#     ggtitle(paste0(model_name, ': Mean |Bias| = ', round(global_bias, 1), ' Pa')) +
#     scale_fill_gradientn(
#       colours = rev(brewer.pal(11, "RdBu")),
#       breaks = v_limits,
#       limits = limits,
#       oob = scales::squish
#     ) +
#     borders("world", colour = 'black', lwd = 0.12) +
#     scale_x_continuous(expand = c(0, 0)) +
#     scale_y_continuous(expand = c(0, 0)) +
#     coord_fixed(ratio = 1.3, xlim = c(-180, 180), ylim = c(-84, 90), expand = FALSE) +
#     theme_bw() +
#     theme(
#       legend.position = 'right',  # Legend to the right
#       panel.grid.major = element_blank(),
#       panel.grid.minor = element_blank(),
#       panel.background = element_blank(),
#       legend.key.size = unit(1, 'cm'),
#       legend.key.height = unit(1.4, 'cm'),
#       legend.key.width = unit(0.4, 'cm'),
#       legend.title = element_text(size = 16),
#       legend.text = element_text(size = 12),
#       plot.title = element_text(size = 24),
#       plot.subtitle = element_text(size = 20, hjust = 0.5),
#       axis.text = element_text(size = 14),
#       axis.title = element_text(size = 16)
#     ) +
#     xlab('Longitude') +
#     ylab('Latitude') +
#     labs(fill = '[Pa]') +
#     easy_center_title()
#
#
#   # 8) Save
#   name <- paste0("figure/Bias2/psl/Bias_psl_", model_name)
#   ggsave(paste0(name, ".pdf"), plot = p_bias, width = 20, height = 15, units = "cm", dpi = 300)
#   ggsave(paste0(name, ".png"), plot = p_bias, width = 20, height = 15, units = "cm", dpi = 300)
#
# }
#
# # Map of pr bias for all models
# for (i in seq_along(model_names)) {
#
#   # 1) Back-transform precipitation data (from log(pr + 1))
#   pr_model_back <- exp(mean_future$pr[,,i+1]) - 1
#   pr_ref_back <- exp(mean_future$pr[,,1]) - 1
#
#   # 2) Compute bias in mm/day
#   bias_pr <- pr_model_back - pr_ref_back
#
#   # 3) Melt to dataframe
#   test_df <- melt(bias_pr, varnames = c("lon_idx", "lat_idx"), value.name = "Bias")
#   test_df$lon <- lon[test_df$lon_idx]
#   test_df$lat <- lat[test_df$lat_idx]
#
#   # 4) Apply limits
#   limit <- 4  # Precipitation in mm/day
#   limits <- c(-limit, limit)
#   v_limits <- as.numeric(format(seq(-limit, limit, len=6), digits = 3))
#   test_df$Bias[test_df$Bias < -limit] <- -limit
#   test_df$Bias[test_df$Bias > limit] <- limit
#
#   # 5) Wrap longitude
#   test_df$lon_wrapped <- ifelse(test_df$lon > 180, test_df$lon - 360, test_df$lon)
#
#   # 6) Cosine latitude weights for *absolute* global mean bias
#   cos_weights <- cos(lat * pi/180)
#   weight_matrix <- matrix(rep(cos_weights, each = length(lon)), nrow = length(lon), byrow = FALSE)
#   global_bias <- sum(abs(bias_pr) * weight_matrix, na.rm = TRUE) / sum(weight_matrix, na.rm = TRUE)
#
#   # 7) Model name for the title
#   model_name <- model_names[i]
#
#   # 8) Plot
#   p_bias <- ggplot() +
#     geom_tile(data = test_df, aes(x = lon_wrapped, y = lat, fill = Bias)) +
#     labs(subtitle = 'Projection period: 1998 - 2023') +
#     ggtitle(paste0(model_name, ': Mean |Bias| = ', round(global_bias, 2), 'mm/day')) +
#     scale_fill_gradientn(
#       colours = rev(brewer.pal(11, "RdBu")),
#       breaks = v_limits,
#       limits = limits,
#       oob = scales::squish
#     ) +
#     borders("world", colour = 'black', lwd = 0.12) +
#     scale_x_continuous(expand = c(0, 0)) +
#     scale_y_continuous(expand = c(0, 0)) +
#     coord_fixed(ratio = 1.3, xlim = c(-180, 180), ylim = c(-84, 90), expand = FALSE) +
#     theme_bw() +
#     theme(
#       legend.position = 'right',  # Legend to the right
#       panel.grid.major = element_blank(),
#       panel.grid.minor = element_blank(),
#       panel.background = element_blank(),
#       legend.key.size = unit(1, 'cm'),
#       legend.key.height = unit(1.4, 'cm'),
#       legend.key.width = unit(0.4, 'cm'),
#       legend.title = element_text(size = 16),
#       legend.text = element_text(size = 12),
#       plot.title = element_text(size = 24),
#       plot.subtitle = element_text(size = 20, hjust = 0.5),
#       axis.text = element_text(size = 14),
#       axis.title = element_text(size = 16)
#     ) +
#     xlab('Longitude') +
#     ylab('Latitude') +
#     labs(fill = '[]') +
#     easy_center_title()
#
#
#   # 9) Save
#   name <- paste0("figure/Bias2/pr/Bias_pr_", model_name)
#   ggsave(paste0(name, ".pdf"), plot = p_bias, width = 20, height = 15, units = "cm", dpi = 300)
#   ggsave(paste0(name, ".png"), plot = p_bias, width = 20, height = 15, units = "cm", dpi = 300)
# }
#
# # MMM bias maps :
# {
#   # Clear up memory
#   gc()
#
#   # Load RColorBrewer if needed
#   library(RColorBrewer)
#
#   # Variables to loop over
#   variables <- c("tas", "pr", "psl")
#   limits_list <- list(
#     tas = 8,
#     pr  = 4,
#     psl = 500
#   )
#   unit_list <- list(
#     tas = "[K]",
#     pr  = "[mm]",
#     psl = "[Pa]"
#   )
#
#   for (var in variables) {
#     # 1) Compute MMM
#     MMM <- apply(mean_future[[var]][,,2:23], c(1,2), mean, na.rm=TRUE)
#
#     # 2) If var is pr, back-transform (log(pr+1) to mm/day)
#     if (var == "pr") {
#       MMM <- exp(MMM) - 1
#       ref_back <- exp(mean_future[[var]][,,1]) - 1
#     } else {
#       ref_back <- mean_future[[var]][,,1]
#     }
#
#     # 3) Compute bias with reference
#     bias <- MMM - ref_back
#
#     # 4) Melt to dataframe
#     test_df <- melt(bias, varnames = c("lon_idx", "lat_idx"), value.name = "Bias")
#     test_df$lon <- lon[test_df$lon_idx]
#     test_df$lat <- lat[test_df$lat_idx]
#
#     # 5) Apply limits
#     limit <- limits_list[[var]]
#     limits <- c(-limit, limit)
#     v_limits <- as.numeric(format(seq(-limit, limit, len=6), digits = 3))
#     test_df$Bias[test_df$Bias < -limit] <- -limit
#     test_df$Bias[test_df$Bias > limit] <- limit
#
#     # 6) Wrap longitude
#     test_df$lon_wrapped <- ifelse(test_df$lon > 180, test_df$lon - 360, test_df$lon)
#
#     # 7) Compute global mean absolute bias
#     cos_weights <- cos(lat * pi/180)
#     weight_matrix <- matrix(rep(cos_weights, each = length(lon)), nrow = length(lon), byrow = FALSE)
#     global_bias <- sum(abs(bias) * weight_matrix, na.rm = TRUE) / sum(weight_matrix, na.rm = TRUE)
#
#     # 8) Plot
#     p_bias <- ggplot() +
#       geom_tile(data = test_df, aes(x = lon_wrapped, y = lat, fill = Bias)) +
#       labs(subtitle = 'Projection period: 1998 - 2023') +
#       ggtitle(paste0('MMM: Mean |Bias| = ', round(global_bias, 2), ' ', unit_list[[var]])) +
#       scale_fill_gradientn(
#         colours = rev(brewer.pal(11, "RdBu")),
#         breaks = v_limits,
#         limits = limits,
#         oob = scales::squish
#       ) +
#       borders("world", colour = 'black', lwd = 0.12) +
#       scale_x_continuous(expand = c(0, 0)) +
#       scale_y_continuous(expand = c(0, 0)) +
#       coord_fixed(ratio = 1.3, xlim = c(-180, 180), ylim = c(-84, 90), expand = FALSE) +
#       theme_bw() +
#       theme(
#         legend.position = 'right',
#         panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank(),
#         panel.background = element_blank(),
#         legend.key.size = unit(1, 'cm'),
#         legend.key.height = unit(1.4, 'cm'),
#         legend.key.width = unit(0.4, 'cm'),
#         legend.title = element_text(size = 16),
#         legend.text = element_text(size = 12),
#         plot.title = element_text(size = 24),
#         plot.subtitle = element_text(size = 20, hjust = 0.5),
#         axis.text = element_text(size = 14),
#         axis.title = element_text(size = 16)
#       ) +
#       xlab('Longitude') +
#       ylab('Latitude') +
#       labs(fill = unit_list[[var]]) +
#       easy_center_title()
#
#     print(p_bias)
#
#     # 9) Save
#     name <- paste0("figure/Bias2/", var, "/Bias_", var, "_MMM")
#     ggsave(paste0(name, ".pdf"), plot = p_bias, width = 20, height = 15, units = "cm", dpi = 300)
#     ggsave(paste0(name, ".png"), plot = p_bias, width = 20, height = 15, units = "cm", dpi = 300)
#   }
#
# }
#
# # Map bias GC
# {
#   # Load RColorBrewer if needed
#   library(RColorBrewer)
#
#   # Variables to loop over
#   variables <- c("tas", "pr", "psl")
#   limits_list <- list(
#     tas = 8,
#     pr  = 4,    # Adjusted limit for precipitation
#     psl = 500
#   )
#   unit_list <- list(
#     tas = "[K]",
#     pr  = "[mm]",
#     psl = "[Pa]"
#   )
#
#   for (var in variables) {
#
#     # 1) Build GraphCut map of the variable
#     GC_map <- array(NA, dim = dim(GC_result$label_attribution))
#
#     for (l in seq_along(model_names)) {
#       islabel <- which(GC_result$label_attribution == l)
#       GC_map[islabel] <- mean_future[[var]][,, l+1][islabel]  # because index 1 is ref
#     }
#
#     # 2) If var is pr, back-transform (log(pr+1) to mm/day)
#     if (var == "pr") {
#       GC_map <- exp(GC_map) - 1
#       ref_back <- exp(mean_future[[var]][,,1]) - 1
#     } else {
#       ref_back <- mean_future[[var]][,,1]
#     }
#
#     # 3) Compute bias with reference
#     bias <- GC_map - ref_back
#
#     # 4) Melt to dataframe
#     test_df <- melt(bias, varnames = c("lon_idx", "lat_idx"), value.name = "Bias")
#     test_df$lon <- lon[test_df$lon_idx]
#     test_df$lat <- lat[test_df$lat_idx]
#
#     # 5) Apply limits
#     limit <- limits_list[[var]]
#     limits <- c(-limit, limit)
#     v_limits <- as.numeric(format(seq(-limit, limit, len=6), digits = 3))
#     test_df$Bias[test_df$Bias < -limit] <- -limit
#     test_df$Bias[test_df$Bias > limit] <- limit
#
#     # 6) Wrap longitude
#     test_df$lon_wrapped <- ifelse(test_df$lon > 180, test_df$lon - 360, test_df$lon)
#
#     # 7) Compute global mean absolute bias
#     cos_weights <- cos(lat * pi/180)
#     weight_matrix <- matrix(rep(cos_weights, each = length(lon)), nrow = length(lon), byrow = FALSE)
#     global_bias <- sum(abs(bias) * weight_matrix, na.rm = TRUE) / sum(weight_matrix, na.rm = TRUE)
#
#     # 8) Plot
#     p_bias <- ggplot() +
#       geom_tile(data = test_df, aes(x = lon_wrapped, y = lat, fill = Bias)) +
#       labs(subtitle = 'Projection period: 1998 - 2023') +
#       ggtitle(paste0('GC MMM: Mean |Bias| = ', round(global_bias, 2), ' ', unit_list[[var]])) +
#       scale_fill_gradientn(
#         colours = rev(brewer.pal(11, "RdBu")),
#         breaks = v_limits,
#         limits = limits,
#         oob = scales::squish
#       ) +
#       borders("world", colour = 'black', lwd = 0.12) +
#       scale_x_continuous(expand = c(0, 0)) +
#       scale_y_continuous(expand = c(0, 0)) +
#       coord_fixed(ratio = 1.3, xlim = c(-180, 180), ylim = c(-84, 90), expand = FALSE) +
#       theme_bw() +
#       theme(
#         legend.position = 'right',
#         panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank(),
#         panel.background = element_blank(),
#         legend.key.size = unit(1, 'cm'),
#         legend.key.height = unit(1.4, 'cm'),
#         legend.key.width = unit(0.4, 'cm'),
#         legend.title = element_text(size = 16),
#         legend.text = element_text(size = 12),
#         plot.title = element_text(size = 24),
#         plot.subtitle = element_text(size = 20, hjust = 0.5),
#         axis.text = element_text(size = 14),
#         axis.title = element_text(size = 16)
#       ) +
#       xlab('Longitude') +
#       ylab('Latitude') +
#       labs(fill = unit_list[[var]]) +
#       easy_center_title()
#
#     print(p_bias)
#
#     # 9) Save
#     name <- paste0("figure/Bias2/", var, "/Bias_", var, "_GC_MMM")
#     ggsave(paste0(name, ".pdf"), plot = p_bias, width = 20, height = 15, units = "cm", dpi = 300)
#     ggsave(paste0(name, ".png"), plot = p_bias, width = 20, height = 15, units = "cm", dpi = 300)
#
#   }
#
# }
#
# # Map of H_dist3 for all models
# {
#   # Define limits for Hellinger distances
#   limit <- 0.7
#   limits <- c(0.1, limit)
#   v_limits <- seq(limits[1], limits[2], length.out = 4)
#
#
#   # Loop over each model
#   for (i in seq_along(model_names)) {
#
#     # 1) Extract Hellinger distance map for the current model
#     hdist_model <- h_dist_fut[,,i]
#
#     # 2) Melt to dataframe
#     test_df <- melt(hdist_model, varnames = c("lon_idx", "lat_idx"), value.name = "H_dist")
#     test_df$lon <- lon[test_df$lon_idx]
#     test_df$lat <- lat[test_df$lat_idx]
#
#     # 3) Wrap longitude to [-180, 180]
#     test_df$lon_wrapped <- ifelse(test_df$lon > 180, test_df$lon - 360, test_df$lon)
#
#     # 4) Compute global weighted mean Hellinger distance
#     cos_weights <- cos(lat * pi/180)
#     weight_matrix <- matrix(rep(cos_weights, each = length(lon)), nrow = length(lon), byrow = FALSE)
#     global_h <- sum(hdist_model * weight_matrix, na.rm = TRUE) / sum(weight_matrix, na.rm = TRUE)
#
#     # 5) Model name for title
#     model_name <- model_names[i]
#
#     # 6) Plot
#     p_hdist <- ggplot() +
#       geom_tile(data = test_df, aes(x = lon_wrapped, y = lat, fill = H_dist)) +
#       labs(subtitle = 'Projection period: 1998 - 2023') +
#       ggtitle(paste0(model_name, ': Mean H = ', round(global_h, 3))) +
#       scale_fill_gradient(
#         low = "white",
#         high = "#015a8c",
#         limits = limits,
#         breaks = v_limits,
#         oob = scales::squish
#       ) +
#       borders("world", colour = 'black', lwd = 0.12) +
#       scale_x_continuous(expand = c(0, 0)) +
#       scale_y_continuous(expand = c(0, 0)) +
#       coord_fixed(ratio = 1.3, xlim = c(-180, 180), ylim = c(-84, 90), expand = FALSE) +
#       theme_bw() +
#       theme(
#         legend.position = 'right',
#         panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank(),
#         panel.background = element_blank(),
#         legend.key.size = unit(1, 'cm'),
#         legend.key.height = unit(1.4, 'cm'),
#         legend.key.width = unit(0.4, 'cm'),
#         legend.title = element_text(size = 16),
#         legend.text = element_text(size = 12),
#         plot.title = element_text(size = 24),
#         plot.subtitle = element_text(size = 20, hjust = 0.5),
#         axis.text = element_text(size = 14),
#         axis.title = element_text(size = 16)
#       ) +
#       xlab('Longitude') +
#       ylab('Latitude') +
#       labs(fill = 'H') +
#       easy_center_title()
#
#
#     # 7) Save
#     name <- paste0("figure/Hellinger3_2/hdist3_", model_name)
#     ggsave(paste0(name, ".pdf"), plot = p_hdist, width = 20, height = 15, units = "cm", dpi = 300)
#     ggsave(paste0(name, ".png"), plot = p_hdist, width = 20, height = 15, units = "cm", dpi = 300)
#   }
# }
#
# # MMM map of hdist3
# {
#
#   # Define limits for Hellinger distances
#   limit <- 0.5
#   limits <- c(0.1, limit)
#   v_limits <- seq(limits[1], limits[2], length.out = 3)
#
#   # 1) Melt to dataframe
#   test_df <- melt(MMM_hdist_fut, varnames = c("lon_idx", "lat_idx"), value.name = "H_dist")
#   test_df$lon <- lon[test_df$lon_idx]
#   test_df$lat <- lat[test_df$lat_idx]
#
#   # 2) Wrap longitude to [-180, 180]
#   test_df$lon_wrapped <- ifelse(test_df$lon > 180, test_df$lon - 360, test_df$lon)
#
#   # 3) Compute global weighted mean Hellinger distance
#   cos_weights <- cos(lat * pi/180)
#   weight_matrix <- matrix(rep(cos_weights, each = length(lon)), nrow = length(lon), byrow = FALSE)
#   global_h <- sum(MMM_hdist_fut * weight_matrix, na.rm = TRUE) / sum(weight_matrix, na.rm = TRUE)
#
#   # 4) Plot
#   p_hdist <- ggplot() +
#     geom_tile(data = test_df, aes(x = lon_wrapped, y = lat, fill = H_dist)) +
#     labs(subtitle = 'Projection period: 1998 - 2023') +
#     ggtitle(paste0('MMM: Mean H = ', round(global_h, 3))) +
#     scale_fill_gradient(
#       low = "white",
#       high = "#015a8c",
#       limits = limits,
#       breaks = v_limits,
#       oob = scales::squish
#     ) +
#     borders("world", colour = 'black', lwd = 0.12) +
#     scale_x_continuous(expand = c(0, 0)) +
#     scale_y_continuous(expand = c(0, 0)) +
#     coord_fixed(ratio = 1.3, xlim = c(-180, 180), ylim = c(-84, 90), expand = FALSE) +
#     theme_bw() +
#     theme(
#       legend.position = 'right',
#       panel.grid.major = element_blank(),
#       panel.grid.minor = element_blank(),
#       panel.background = element_blank(),
#       legend.key.size = unit(1, 'cm'),
#       legend.key.height = unit(1.4, 'cm'),
#       legend.key.width = unit(0.4, 'cm'),
#       legend.title = element_text(size = 16),
#       legend.text = element_text(size = 12),
#       plot.title = element_text(size = 24),
#       plot.subtitle = element_text(size = 20, hjust = 0.5),
#       axis.text = element_text(size = 14),
#       axis.title = element_text(size = 16)
#     ) +
#     xlab('Longitude') +
#     ylab('Latitude') +
#     labs(fill = 'H') +
#     easy_center_title()
#
#   print(p_hdist)
#
#   # 5) Save
#   name <- "figure/Hellinger3_2/hdist3_MMM"
#   ggsave(paste0(name, ".pdf"), plot = p_hdist, width = 20, height = 15, units = "cm", dpi = 300)
#   ggsave(paste0(name, ".png"), plot = p_hdist, width = 20, height = 15, units = "cm", dpi = 300)
#
# }
#
# # GC Map hdist3
# {
#   # Define limits for Hellinger distances
#   limit <- 0.5
#   limits <- c(0.1, limit)
#   v_limits <- seq(limits[1], limits[2], length.out = 3)
#
#   # 1) Melt to dataframe
#   test_df <- melt(GC_hdist_fut, varnames = c("lon_idx", "lat_idx"), value.name = "H_dist")
#   test_df$lon <- lon[test_df$lon_idx]
#   test_df$lat <- lat[test_df$lat_idx]
#
#   # 2) Wrap longitude to [-180, 180]
#   test_df$lon_wrapped <- ifelse(test_df$lon > 180, test_df$lon - 360, test_df$lon)
#
#   # 3) Compute global weighted mean Hellinger distance
#   cos_weights <- cos(lat * pi/180)
#   weight_matrix <- matrix(rep(cos_weights, each = length(lon)), nrow = length(lon), byrow = FALSE)
#   global_h <- sum(GC_hdist_fut * weight_matrix, na.rm = TRUE) / sum(weight_matrix, na.rm = TRUE)
#
#   # 4) Plot
#   p_hdist <- ggplot() +
#     geom_tile(data = test_df, aes(x = lon_wrapped, y = lat, fill = H_dist)) +
#     labs(subtitle = 'Projection period: 1998 - 2023') +
#     ggtitle(paste0('GraphCut: Mean H = ', round(global_h, 3))) +
#     scale_fill_gradient(
#       low = "white",
#       high = "#015a8c",
#       limits = limits,
#       breaks = v_limits,
#       oob = scales::squish
#     ) +
#     borders("world", colour = 'black', lwd = 0.12) +
#     scale_x_continuous(expand = c(0, 0)) +
#     scale_y_continuous(expand = c(0, 0)) +
#     coord_fixed(ratio = 1.3, xlim = c(-180, 180), ylim = c(-84, 90), expand = FALSE) +
#     theme_bw() +
#     theme(
#       legend.position = 'right',
#       panel.grid.major = element_blank(),
#       panel.grid.minor = element_blank(),
#       panel.background = element_blank(),
#       legend.key.size = unit(1, 'cm'),
#       legend.key.height = unit(1.4, 'cm'),
#       legend.key.width = unit(0.4, 'cm'),
#       legend.title = element_text(size = 16),
#       legend.text = element_text(size = 12),
#       plot.title = element_text(size = 24),
#       plot.subtitle = element_text(size = 20, hjust = 0.5),
#       axis.text = element_text(size = 14),
#       axis.title = element_text(size = 16)
#     ) +
#     xlab('Longitude') +
#     ylab('Latitude') +
#     labs(fill = 'H') +
#     easy_center_title()
#
#   print(p_hdist)
#
#   # 5) Save
#   name <- "figure/Hellinger3_2/hdist3_GC"
#   ggsave(paste0(name, ".pdf"), plot = p_hdist, width = 20, height = 15, units = "cm", dpi = 300)
#   ggsave(paste0(name, ".png"), plot = p_hdist, width = 20, height = 15, units = "cm", dpi = 300)
# }
#
#
# # MMM Map hdist2
# {
#   # Initialize matrix for Hellinger distance
#   mmm_hdist_prtas_fut <- matrix(NA, nrow = length(lon), ncol = length(lat))
#
#   # Compute for each grid point
#   for (i in seq_along(lon)) {
#     for (j in seq_along(lat)) {
#       # Reference PDF
#       pdf_ref <- pdf2_future$pr_tas[i, j, , 1]
#
#       # MMM PDF (mean over models 2:23)
#       pdf_models <- pdf2_future$pr_tas[i, j, , 2:23]
#       mmm_pdf <- rowMeans(pdf_models, na.rm = TRUE)
#
#       # Hellinger distance
#       hd <- sqrt(sum((sqrt(mmm_pdf) - sqrt(pdf_ref))^2)) / sqrt(2)
#       mmm_hdist_prtas_fut[i, j] <- hd
#     }
#   }
#
#   # Replace any NaN with 0
#   mmm_hdist_prtas_fut[is.nan(mmm_hdist_prtas_fut)] <- 0
#
#   # Basic visualization to check
#   hist(mmm_hdist_prtas_fut, breaks = 50, main = "Hellinger Distance for MMM (pr_tas)")
#
#   # Compute global mean
#   cos_weights <- cos(lat * pi/180)
#   weight_matrix <- matrix(rep(cos_weights, each = length(lon)), nrow = length(lon), byrow = FALSE)
#   global_mean_h <- sum(mmm_hdist_prtas_fut * weight_matrix, na.rm = TRUE) / sum(weight_matrix, na.rm = TRUE)
#   cat("Global weighted mean H for MMM (pr_tas):", round(global_mean_h, 3), "\n")
#
#   # Remove temporary variables to save memory
#   rm(pdf_ref, pdf_models, mmm_pdf)
#   gc()
#
#   # Plot for MMM Hellinger Distance (pr_tas)
# {
#   # Define limits for Hellinger distances
#   limit <- 0.5
#   limits <- c(0.1, limit)
#   v_limits <- seq(limits[1], limits[2], length.out = 3)
#
#   # 1) Melt to dataframe
#   test_df <- melt(mmm_hdist_prtas_fut, varnames = c("lon_idx", "lat_idx"), value.name = "H_dist")
#   test_df$lon <- lon[test_df$lon_idx]
#   test_df$lat <- lat[test_df$lat_idx]
#
#   # 2) Wrap longitude to [-180, 180]
#   test_df$lon_wrapped <- ifelse(test_df$lon > 180, test_df$lon - 360, test_df$lon)
#
#   # 3) Compute global weighted mean Hellinger distance
#   cos_weights <- cos(lat * pi/180)
#   weight_matrix <- matrix(rep(cos_weights, each = length(lon)), nrow = length(lon), byrow = FALSE)
#   global_h <- sum(mmm_hdist_prtas_fut * weight_matrix, na.rm = TRUE) / sum(weight_matrix, na.rm = TRUE)
#
#   # 4) Plot
#   p_hdist <- ggplot() +
#     geom_tile(data = test_df, aes(x = lon_wrapped, y = lat, fill = H_dist)) +
#     labs(subtitle = 'Projection period: 1998 - 2023') +
#     ggtitle(paste0('MMM pr_tas: Mean H = ', round(global_h, 3))) +
#     scale_fill_gradient(
#       low = "white",
#       high = "#015a8c",
#       limits = limits,
#       breaks = v_limits,
#       oob = scales::squish
#     ) +
#     borders("world", colour = 'black', lwd = 0.12) +
#     scale_x_continuous(expand = c(0, 0)) +
#     scale_y_continuous(expand = c(0, 0)) +
#     coord_fixed(ratio = 1.3, xlim = c(-180, 180), ylim = c(-84, 90), expand = FALSE) +
#     theme_bw() +
#     theme(
#       legend.position = 'right',
#       panel.grid.major = element_blank(),
#       panel.grid.minor = element_blank(),
#       panel.background = element_blank(),
#       legend.key.size = unit(1, 'cm'),
#       legend.key.height = unit(1.4, 'cm'),
#       legend.key.width = unit(0.4, 'cm'),
#       legend.title = element_text(size = 16),
#       legend.text = element_text(size = 12),
#       plot.title = element_text(size = 24),
#       plot.subtitle = element_text(size = 20, hjust = 0.5),
#       axis.text = element_text(size = 14),
#       axis.title = element_text(size = 16)
#     ) +
#     xlab('Longitude') +
#     ylab('Latitude') +
#     labs(fill = 'H') +
#     easy_center_title()
#
#   print(p_hdist)
#
#   # 5) Save
#   name <- "figure/Hellinger_2/hdist2_MMM_prtas"
#   ggsave(paste0(name, ".pdf"), plot = p_hdist, width = 20, height = 15, units = "cm", dpi = 300)
#   ggsave(paste0(name, ".png"), plot = p_hdist, width = 20, height = 15, units = "cm", dpi = 300)
# }
#
# }
#
#
# # GC Map hdist pr_tas
# {
#   # Initialize matrix for Hellinger distance
#   gc_hdist_prtas_fut <- matrix(NA, nrow = length(lon), ncol = length(lat))
#
#   # Compute for each grid point
#   for (i in seq_along(lon)) {
#     for (j in seq_along(lat)) {
#       # Reference PDF
#       pdf_ref <- pdf2_future$pr_tas[i, j, , 1]
#
#       # GraphCut-selected model index (labels start at 1)
#       label <- GC_result$label_attribution[i, j]
#
#       # GraphCut model PDF at this point
#       pdf_gc <- pdf2_future$pr_tas[i, j, , label+1]
#
#       # Hellinger distance
#       hd <- sqrt(sum((sqrt(pdf_gc) - sqrt(pdf_ref))^2)) / sqrt(2)
#       gc_hdist_prtas_fut[i, j] <- hd
#     }
#   }
#
#   # Replace any NaN with 0
#   gc_hdist_prtas_fut[is.nan(gc_hdist_prtas_fut)] <- 0
#
#   # Basic visualization to check
#   hist(gc_hdist_prtas_fut, breaks = 50, main = "Hellinger Distance for GC (pr_tas)")
#
#   # Compute global mean
#   cos_weights <- cos(lat * pi/180)
#   weight_matrix <- matrix(rep(cos_weights, each = length(lon)), nrow = length(lon), byrow = FALSE)
#   global_mean_h <- sum(gc_hdist_prtas_fut * weight_matrix, na.rm = TRUE) / sum(weight_matrix, na.rm = TRUE)
#   cat("Global weighted mean H for GC (pr_tas):", round(global_mean_h, 3), "\n")
#
#   # Remove temporary variables to save memory
#   rm(pdf_ref, pdf_gc)
#   gc()
#
#   # Plot for GC Hellinger Distance (pr_tas)
# {
#   # Define limits for Hellinger distances
#   limit <- 0.5
#   limits <- c(0.1, limit)
#   v_limits <- seq(limits[1], limits[2], length.out = 3)
#
#   # 1) Melt to dataframe
#   test_df <- melt(gc_hdist_prtas_fut, varnames = c("lon_idx", "lat_idx"), value.name = "H_dist")
#   test_df$lon <- lon[test_df$lon_idx]
#   test_df$lat <- lat[test_df$lat_idx]
#
#   # 2) Wrap longitude to [-180, 180]
#   test_df$lon_wrapped <- ifelse(test_df$lon > 180, test_df$lon - 360, test_df$lon)
#
#   # 3) Compute global weighted mean Hellinger distance
#   cos_weights <- cos(lat * pi/180)
#   weight_matrix <- matrix(rep(cos_weights, each = length(lon)), nrow = length(lon), byrow = FALSE)
#   global_h <- sum(gc_hdist_prtas_fut * weight_matrix, na.rm = TRUE) / sum(weight_matrix, na.rm = TRUE)
#
#   # 4) Plot
#   p_hdist <- ggplot() +
#     geom_tile(data = test_df, aes(x = lon_wrapped, y = lat, fill = H_dist)) +
#     labs(subtitle = 'Projection period: 1998 - 2023') +
#     ggtitle(paste0('GraphCut pr_tas: Mean H = ', round(global_h, 3))) +
#     scale_fill_gradient(
#       low = "white",
#       high = "#015a8c",
#       limits = limits,
#       breaks = v_limits,
#       oob = scales::squish
#     ) +
#     borders("world", colour = 'black', lwd = 0.12) +
#     scale_x_continuous(expand = c(0, 0)) +
#     scale_y_continuous(expand = c(0, 0)) +
#     coord_fixed(ratio = 1.3, xlim = c(-180, 180), ylim = c(-84, 90), expand = FALSE) +
#     theme_bw() +
#     theme(
#       legend.position = 'right',
#       panel.grid.major = element_blank(),
#       panel.grid.minor = element_blank(),
#       panel.background = element_blank(),
#       legend.key.size = unit(1, 'cm'),
#       legend.key.height = unit(1.4, 'cm'),
#       legend.key.width = unit(0.4, 'cm'),
#       legend.title = element_text(size = 16),
#       legend.text = element_text(size = 12),
#       plot.title = element_text(size = 24),
#       plot.subtitle = element_text(size = 20, hjust = 0.5),
#       axis.text = element_text(size = 14),
#       axis.title = element_text(size = 16)
#     ) +
#     xlab('Longitude') +
#     ylab('Latitude') +
#     labs(fill = 'H') +
#     easy_center_title()
#
#   print(p_hdist)
#
#   # 5) Save
#   name <- "figure/Hellinger_2/hdist2_GC_prtas"
#   ggsave(paste0(name, ".pdf"), plot = p_hdist, width = 20, height = 15, units = "cm", dpi = 300)
#   ggsave(paste0(name, ".png"), plot = p_hdist, width = 20, height = 15, units = "cm", dpi = 300)
# }
# }
#
#
# # Single point PDF 2d
# {
#   library(plotly)
#
#   # Example grid point indices
#   lon_index <- 180+6 # Example longitude index
#   lat_index <- 90+46  # Example latitude index
#
#   # Extract the 256-bin PDF vector for the specific grid point (2D PDF: pr_tas)
#   pdf_vector <- pdf2_future$pr_tas[lon_index, lat_index, , 3]
#
#   # Reshape the PDF vector into a 2D array of dimensions [16, 16]
#   nbins <- 16
#   pdf_2d <- matrix(pdf_vector, nrow = nbins, ncol = nbins)
#
#   # Extract the ranges for the variables from range_var_final
#   range_var <- range_var_final$ranges
#   var1_min <- range_var$pr[lon_index, lat_index, 1]
#   var1_max <- range_var$pr[lon_index, lat_index, 2]
#   var2_min <- range_var$tas[lon_index, lat_index, 1]
#   var2_max <- range_var$tas[lon_index, lat_index, 2]
#
#   # Create bin edges for each variable
#   x_bins <- seq(var1_min, var1_max, length.out = nbins + 1)
#   y_bins <- seq(var2_min, var2_max, length.out = nbins + 1)
#
#   # Create the coordinates for the centers of the bins
#   x_centers <- (x_bins[-1] + x_bins[-length(x_bins)]) / 2
#   y_centers <- (y_bins[-1] + y_bins[-length(y_bins)]) / 2
#
#   # Create a meshgrid for the bin centers
#   grid <- expand.grid(x = x_centers, y = y_centers)
#
#   # Flatten the PDF matrix into a vector
#   pdf_flat <- as.vector(pdf_2d)
#
#   # Combine the coordinates with the PDF values
#   plot_data <- data.frame(
#     x = grid$x,
#     y = grid$y,
#     value = pdf_flat
#   )
#
#   # Plot the 2D histogram as a heatmap
#   fig <- plot_ly(
#     data = plot_data,
#     x = ~x,
#     y = ~y,
#     z = ~value,
#     type = "heatmap",
#     colorscale = "Viridis",
#     colorbar = list(title = "PDF Value"),
#     text = ~paste("PDF Value:", round(value, 4))
#   ) %>%
#     layout(
#       xaxis = list(title = "Variable 1 (pr)"),
#       yaxis = list(title = "Variable 2 (tas)"),
#       title = paste("2D PDF for Grid Point (Lon:", lon_index, ", Lat:", lat_index, ")")
#     )
#
#   # Show the plot
#   fig
# }
#
#
# # HDR and lDR on 2d
# {
#
#   # Example grid point
#   lon <-  101
#   lat <- 90 + 4
#
#   # Fixed axis limits
#   tas_min <- 292
#   tas_max <- 305
#   tas_limit <- c(tas_min, tas_max)
#   pr_max <- max(log2(180))
#   pr_limit <- c(0, pr_max)
#
#   # Load time series data
#   nc <- nc_open('data/CMIP6_merged_all/EC-Earth3-CC/tas/tas_EC-Earth3-CC_19500101-21001230.nc')
#   var <- 'tas'
#   year_present <- 1970:2014
#   yyyy <- substr(as.character(nc.get.time.series(nc)), 1, 4)
#   iyyyy <- which(yyyy %in% year_present)
#   tas <- ncvar_get(nc, var, start = c(lon, lat, min(iyyyy)), count = c(1, 1, length(iyyyy)))
#
#   var <- 'pr'
#   nc <- nc_open('data/CMIP6_merged_all/EC-Earth3-CC/pr/pr_EC-Earth3-CC_19500101-21001230.nc')
#   pr <- ncvar_get(nc, var, start = c(lon, lat, min(iyyyy)), count = c(1, 1, length(iyyyy)))
#   pr <- log2(pr + 1)
#
#   # Create data frame
#   df <- data.frame(x = tas, y = pr)
#
#   # Bin definitions
#   tas_bins <- seq(tas_min, tas_max, length.out = 65)
#   pr_bins <- seq(0, pr_max, length.out = 65)
#
#   df$x_bin <- cut(df$x, breaks = tas_bins, include.lowest = TRUE)
#   df$y_bin <- cut(df$y, breaks = pr_bins, include.lowest = TRUE)
#
#   counts <- as.data.frame(table(df$x_bin, df$y_bin))
#   names(counts) <- c("x_bin", "y_bin", "count")
#
#   x_centers <- (tas_bins[-1] + tas_bins[-length(tas_bins)]) / 2
#   y_centers <- (pr_bins[-1] + pr_bins[-length(pr_bins)]) / 2
#   counts$x <- x_centers[as.numeric(counts$x_bin)]
#   counts$y <- y_centers[as.numeric(counts$y_bin)]
#
#   # Replace zeros with NA for full PDF plot to show as white
#   counts$value <- ifelse(counts$count == 0, NA, counts$count)
#
#   # Plot full PDF
#   p_full <- ggplot(counts, aes(x = x, y = y, fill = value)) +
#     geom_tile() +
#     scale_fill_gradientn(colors = viridis(64), name = "Count",
#                          limits = c(0, 80),
#                          oob = scales::squish,
#                          na.value = "white") +
#     labs(title = "Full PDF", x = "Temperature [K]", y = "Precipitation [log2(mm/day+1)]") +
#     scale_x_continuous(limits = tas_limit, expand = c(0, 0)) +
#     scale_y_continuous(limits = pr_limit, expand = c(0, 0)) +
#     theme_bw() +
#     easy_center_title()
#   print(p_full)
#
#   # Create a vector of probabilities
#   pdf_vector <- counts$count / sum(counts$count)
#
#   # Compute HDR indices
#   hdr_indices <- select_hdr_indices(pdf_vector, tau = 0.10)
#
#   # Compute LDR indices (complement)
#   all_indices <- seq_along(pdf_vector)
#   ldr_indices <- setdiff(all_indices, hdr_indices)
#
#   # Create HDR and LDR "value" columns separately for plotting
#   counts$hdr_value <- NA
#   counts$ldr_value <- NA
#   counts$hdr_value[hdr_indices] <- counts$count[hdr_indices]
#   counts$ldr_value[ldr_indices] <- counts$count[ldr_indices]
#
#   # Plot HDR
#   p_hdr <- ggplot(counts, aes(x = x, y = y, fill = hdr_value)) +
#     geom_tile() +
#     scale_fill_gradientn(colors = viridis(64), name = "Count",
#                          limits = c(0, 80),
#                          oob = scales::squish,
#                          na.value = "white") +
#     labs(title = "HDR (Top 90%)", x = "Temperature [K]", y = "Precipitation [log2(mm/day+1)]") +
#     scale_x_continuous(limits = tas_limit, expand = c(0, 0)) +
#     scale_y_continuous(limits = pr_limit, expand = c(0, 0)) +
#     theme_bw() +
#     easy_center_title()
#   print(p_hdr)
#
#   # Filter to bins that are part of LDR and have a count > 0
#   ldr_df <- counts[ldr_indices, ]
#   ldr_df <- ldr_df[ldr_df$count > 0, ]
#
#   p_ldr <- ggplot(ldr_df, aes(x = x, y = y, fill = count)) +
#     geom_tile() +
#     scale_fill_gradientn(colors = viridis(64), name = 'Count', limits = c(0, 80), oob = scales::squish) +
#     labs(title = 'LDR (Lowest 10%)', x = 'Temperature [K]', y = 'Precipitation [log2(mm/day+1)]') +
#     scale_x_continuous(limits = tas_limit, expand = c(0, 0)) +
#     scale_y_continuous(limits = pr_limit, expand = c(0, 0)) +
#     theme_bw() +
#     easy_center_title()
#   print(p_ldr)
#
#
#   # Create a new column for LDR with full-grid style
#   counts$ldr_full_value <- NA
#   counts$ldr_full_value[ldr_indices] <- counts$count[ldr_indices]
#
#   # Plot LDR with full grid: empty LDR bins are still shown in color
#   p_ldr_full <- ggplot(counts, aes(x = x, y = y, fill = ldr_full_value)) +
#     geom_tile() +
#     scale_fill_gradientn(colors = viridis(64), name = "Count",
#                          limits = c(0, 80),
#                          oob = scales::squish,
#                          na.value = "white") +
#     labs(title = "LDR (Lowest 10%)", x = "Temperature [K]", y = "Precipitation [log2(mm/day+1)]") +
#     scale_x_continuous(limits = tas_limit, expand = c(0, 0)) +
#     scale_y_continuous(limits = pr_limit, expand = c(0, 0)) +
#     theme_bw() +
#     easy_center_title()
#   print(p_ldr_full)
#
#   # Save this fourth plot
#   ggsave('figure/PDF2D/LDR_2Dpdf_gridpoint_fullgrid.png', plot = p_ldr_full, width = 15, height = 10, units = "cm", dpi = 300)
#
#
#   # Save plots
#   ggsave('figure/PDF2D/full_2Dpdf_gridpoint.png', plot = p_full, width = 15, height = 10, units = "cm", dpi = 300)
#   ggsave('figure/PDF2D/HDR_2Dpdf_gridpoint.png', plot = p_hdr, width = 15, height = 10, units = "cm", dpi = 300)
#   ggsave('figure/PDF2D/LDR_2Dpdf_gridpoint.png', plot = p_ldr, width = 15, height = 10, units = "cm", dpi = 300)
#
#
#   # Plot HDR with grey for empty bins
#   p_hdr_grey <- ggplot(counts, aes(x = x, y = y, fill = hdr_value)) +
#     geom_tile() +
#     scale_fill_gradientn(colors = viridis(64), name = "Count",
#                          limits = c(0, 80),
#                          oob = scales::squish,
#                          na.value = "grey80") +   # grey for empty bins
#     labs(title = "HDR (Top 90%)", x = "Temperature [K]", y = "Precipitation [log2(mm/day+1)]") +
#     scale_x_continuous(limits = tas_limit, expand = c(0, 0)) +
#     scale_y_continuous(limits = pr_limit, expand = c(0, 0)) +
#     theme_bw() +
#     easy_center_title()
#   print(p_hdr_grey)
#
#   # Plot LDR with grey for empty bins
#   counts$ldr_full_value_grey <- NA
#   counts$ldr_full_value_grey[ldr_indices] <- counts$count[ldr_indices]
#   p_ldr_grey <- ggplot(counts, aes(x = x, y = y, fill = ldr_full_value_grey)) +
#     geom_tile() +
#     scale_fill_gradientn(colors = viridis(64), name = "Count",
#                          limits = c(0, 80),
#                          oob = scales::squish,
#                          na.value = "grey80") +   # grey for empty bins
#     labs(title = "LDR (Lowest 10%)", x = "Temperature [K]", y = "Precipitation [log2(mm/day+1)]") +
#     scale_x_continuous(limits = tas_limit, expand = c(0, 0)) +
#     scale_y_continuous(limits = pr_limit, expand = c(0, 0)) +
#     theme_bw() +
#     easy_center_title()
#   print(p_ldr_grey)
#
#   # Save these grey-empty plots
#   ggsave('figure/PDF2D/HDR_2Dpdf_gridpoint_greyempty.png', plot = p_hdr_grey, width = 15, height = 10, units = "cm", dpi = 300)
#   ggsave('figure/PDF2D/LDR_2Dpdf_gridpoint_greyempty.png', plot = p_ldr_grey, width = 15, height = 10, units = "cm", dpi = 300)
#
# }
#
#
# {
#
#   # Example grid point
#   lon <-  101
#   lat <- 90 + 4
#
#   # Fixed axis limits
#   tas_min <- 292
#   tas_max <- 305
#   tas_limit <- c(tas_min, tas_max)
#   pr_max <- max(log2(180))
#   pr_limit <- c(0, pr_max)
#
#   # Load time series data
#   nc <- nc_open('data/CMIP6_merged_all/EC-Earth3-CC/tas/tas_EC-Earth3-CC_19500101-21001230.nc')
#   var <- 'tas'
#   year_present <- 1970:2014
#   yyyy <- substr(as.character(nc.get.time.series(nc)), 1, 4)
#   iyyyy <- which(yyyy %in% year_present)
#   tas <- ncvar_get(nc, var, start = c(lon, lat, min(iyyyy)), count = c(1, 1, length(iyyyy)))
#
#   var <- 'pr'
#   nc <- nc_open('data/CMIP6_merged_all/EC-Earth3-CC/pr/pr_EC-Earth3-CC_19500101-21001230.nc')
#   pr <- ncvar_get(nc, var, start = c(lon, lat, min(iyyyy)), count = c(1, 1, length(iyyyy)))
#   pr <- log2(pr + 1)
#
#   # Create data frame
#   df <- data.frame(x = tas, y = pr)
#
#   # Bin definitions
#   tas_bins <- seq(tas_min, tas_max, length.out = 65)
#   pr_bins <- seq(0, pr_max, length.out = 65)
#
#   df$x_bin <- cut(df$x, breaks = tas_bins, include.lowest = TRUE)
#   df$y_bin <- cut(df$y, breaks = pr_bins, include.lowest = TRUE)
#
#   counts <- as.data.frame(table(df$x_bin, df$y_bin))
#   names(counts) <- c("x_bin", "y_bin", "count")
#
#   x_centers <- (tas_bins[-1] + tas_bins[-length(tas_bins)]) / 2
#   y_centers <- (pr_bins[-1] + pr_bins[-length(pr_bins)]) / 2
#   counts$x <- x_centers[as.numeric(counts$x_bin)]
#   counts$y <- y_centers[as.numeric(counts$y_bin)]
#
#   # Convert counts to densities
#   total_points <- sum(counts$count)
#   bin_area <- (tas_max - tas_min)/64 * (pr_max - 0)/64
#   counts$density <- counts$count / (total_points * bin_area)
#
#   # Replace zeros with NA for full PDF plot to show as white
#   counts$value <- ifelse(counts$density == 0, NA, counts$density)
#   max_density <- max(counts$density, na.rm = TRUE)
#
#   # Plot full PDF (density)
#   p_full <- ggplot(counts, aes(x = x, y = y, fill = value)) +
#     geom_tile() +
#     scale_fill_gradientn(colors = viridis(64), name = "Density",
#                          limits = c(0, max_density),
#                          oob = scales::squish,
#                          na.value = "white") +
#     labs(title = "Full PDF", x = "Temperature [K]", y = "Precipitation [log2(mm/day+1)]") +
#     scale_x_continuous(limits = tas_limit, expand = c(0, 0)) +
#     scale_y_continuous(limits = pr_limit, expand = c(0, 0)) +
#     theme_bw() +
#     easy_center_title()
#   print(p_full)
#
#   # Compute HDR and LDR
#   pdf_vector <- counts$density / sum(counts$density)
#   hdr_indices <- select_hdr_indices(pdf_vector, tau = 0.10)
#   all_indices <- seq_along(pdf_vector)
#   ldr_indices <- setdiff(all_indices, hdr_indices)
#
#   # Create HDR and LDR density columns
#   counts$hdr_density <- NA
#   counts$ldr_density <- NA
#   counts$hdr_density[hdr_indices] <- counts$density[hdr_indices]
#   counts$ldr_density[ldr_indices] <- counts$density[ldr_indices]
#
#   # Plot HDR
#   p_hdr <- ggplot(counts, aes(x = x, y = y, fill = hdr_density)) +
#     geom_tile() +
#     scale_fill_gradientn(colors = viridis(64), name = "Density",
#                          limits = c(0, max_density),
#                          oob = scales::squish,
#                          na.value = "white") +
#     labs(title = "HDR (Top 90%)", x = "Temperature [K]", y = "Precipitation [log2(mm/day+1)]") +
#     scale_x_continuous(limits = tas_limit, expand = c(0, 0)) +
#     scale_y_continuous(limits = pr_limit, expand = c(0, 0)) +
#     theme_bw() +
#     easy_center_title()
#   print(p_hdr)
#
#   # Filter to LDR bins with density > 0
#   ldr_df <- counts[ldr_indices, ]
#   ldr_df <- ldr_df[ldr_df$density > 0, ]
#
#   p_ldr <- ggplot(ldr_df, aes(x = x, y = y, fill = density)) +
#     geom_tile() +
#     scale_fill_gradientn(colors = viridis(64), name = 'Density', limits = c(0, max_density), oob = scales::squish) +
#     labs(title = 'LDR (Lowest 10%)', x = 'Temperature [K]', y = 'Precipitation [log2(mm/day+1)]') +
#     scale_x_continuous(limits = tas_limit, expand = c(0, 0)) +
#     scale_y_continuous(limits = pr_limit, expand = c(0, 0)) +
#     theme_bw() +
#     easy_center_title()
#   print(p_ldr)
#
#   # Plot LDR with full grid
#   counts$ldr_full_density <- NA
#   counts$ldr_full_density[ldr_indices] <- counts$density[ldr_indices]
#   p_ldr_full <- ggplot(counts, aes(x = x, y = y, fill = ldr_full_density)) +
#     geom_tile() +
#     scale_fill_gradientn(colors = viridis(64), name = "Density",
#                          limits = c(0, max_density),
#                          oob = scales::squish,
#                          na.value = "white") +
#     labs(title = "LDR (Lowest 10%) Full Grid", x = "Temperature [K]", y = "Precipitation [log2(mm/day+1)]") +
#     scale_x_continuous(limits = tas_limit, expand = c(0, 0)) +
#     scale_y_continuous(limits = pr_limit, expand = c(0, 0)) +
#     theme_bw() +
#     easy_center_title()
#   print(p_ldr_full)
#
#   # Plot HDR with grey for empty bins
#   p_hdr_grey <- ggplot(counts, aes(x = x, y = y, fill = hdr_density)) +
#     geom_tile() +
#     scale_fill_gradientn(colors = viridis(64), name = "Density",
#                          limits = c(0, max_density),
#                          oob = scales::squish,
#                          na.value = "grey80") +
#     labs(title = "HDR (Top 90%)", x = "Temperature [K]", y = "Precipitation [log2(mm/day+1)]") +
#     scale_x_continuous(limits = tas_limit, expand = c(0, 0)) +
#     scale_y_continuous(limits = pr_limit, expand = c(0, 0)) +
#     theme_bw() +
#     easy_center_title()
#   print(p_hdr_grey)
#
#   # Plot LDR with grey for empty bins
#   counts$ldr_full_density_grey <- NA
#   counts$ldr_full_density_grey[ldr_indices] <- counts$density[ldr_indices]
#   p_ldr_grey <- ggplot(counts, aes(x = x, y = y, fill = ldr_full_density_grey)) +
#     geom_tile() +
#     scale_fill_gradientn(colors = viridis(64), name = "Density",
#                          limits = c(0, max_density),
#                          oob = scales::squish,
#                          na.value = "grey80") +
#     labs(title = "LDR (Lowest 10%)", x = "Temperature [K]", y = "Precipitation [log2(mm/day+1)]") +
#     scale_x_continuous(limits = tas_limit, expand = c(0, 0)) +
#     scale_y_continuous(limits = pr_limit, expand = c(0, 0)) +
#     theme_bw() +
#     easy_center_title()
#   print(p_ldr_grey)
#
#   # Save plots
#   ggsave('figure/PDF2D2/full_2Dpdf_density.png', plot = p_full, width = 15, height = 10, units = "cm", dpi = 300)
#   ggsave('figure/PDF2D2/HDR_2Dpdf_density.png', plot = p_hdr, width = 15, height = 10, units = "cm", dpi = 300)
#   ggsave('figure/PDF2D2/LDR_2Dpdf_density.png', plot = p_ldr, width = 15, height = 10, units = "cm", dpi = 300)
#   ggsave('figure/PDF2D2/LDR_2Dpdf_density_fullgrid.png', plot = p_ldr_full, width = 15, height = 10, units = "cm", dpi = 300)
#   ggsave('figure/PDF2D2/HDR_2Dpdf_density_greyempty.png', plot = p_hdr_grey, width = 15, height = 10, units = "cm", dpi = 300)
#   ggsave('figure/PDF2D2/LDR_2Dpdf_density_greyempty.png', plot = p_ldr_grey, width = 15, height = 10, units = "cm", dpi = 300)
# }



