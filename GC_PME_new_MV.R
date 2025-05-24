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


range_var_final <- readRDS('ranges/range_var_final_allModelsPar_1950-2023_90deg_3v_PME.rds')

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
year_future <<- 2075:2100
# data directory
data_dir <<- 'data/CMIP6_merged_all/'

# List of the variable used
variables <- c('pr', 'tas', 'psl')

## 1.  bin-resolution choices
nbins3d <- 8     # 3-D joint
nbins2d <- 16    # 2-D pairs
nbins1d <- 32    # 1-D marginals

## 2.  model list
model_names <- scan("model_names_pr_tas_psl_perfect_model_without_duplicate.txt", what = "", quiet = TRUE)

## 3.  number of parallel workers
workers <- 1   # adapt to your machine

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
filename <- paste0(stamp, "_workspace_multiRes_3v.RData")
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
smooth_cost <- 1

GC_result <- tryCatch({
  GraphCutHellinger_nD_lat(
    pdf_models_future = pdf3_models_fut,
    h_dist            = h_dist_pres,
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


{
  # Extract the label attribution for the current smooth cost
  GC_labels <- GC_result$label_attribution

  # Convert the label matrix to a data frame for plotting
  label_df <- reshape2::melt(GC_labels, varnames = c("lon_idx", "lat_idx"), value.name = "label_attribution")

  # Explicitly assign correct longitude and latitude values
  label_df$lon <- lon[label_df$lon_idx]  # Map longitude indices to values
  label_df$lat <- lat[label_df$lat_idx]    # Map latitude indices to values

  # Force label_attribution to be a factor with levels in the exact order of model_names.
  label_df$label_attribution <- factor(label_df$label_attribution,
                                       levels = seq_along(model_names),
                                       labels = model_names)
  # Create a named color palette: each model name is explicitly mapped to its color.
  color_palette <- c(
    "#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd",
    "#8c564b", "#e377c2", "#7f7f7f", "#bcbd22", "#17becf",
    "#aec7e8", "#ffbb78", "#98df8a", "#ff9896", "#c5b0d5",
    "#c49c94", "#f7b6d2", "#c7c7c7", "#dbdb8d", "#9edae5",
    "#393b79", "#5254a3", "#6b6ecf"
  )
  # Name the palette vector with model_names (in the same order)
  names(color_palette) <- model_names

  p6 <- ggplot() +
    geom_tile(data = label_df, aes(x = lon, y = lat, fill = label_attribution)) +
    scale_fill_manual(
      values = color_palette,
      na.value = "white",
      guide = guide_legend(title = "Model Names", ncol = 1)
    ) +
    ggtitle(paste("Label GC Hellinger - Lambda:", 0.1)) +
    borders("world", colour = 'black', size = 0.12) +
    theme_bw() +
    coord_fixed(ratio = 1.3, xlim = c(-180, 180), ylim = c(-84, 90), expand = FALSE)+
    theme(
      legend.position = 'right',
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.background = element_blank(),
      legend.key.size = unit(0.5, 'cm'),
      legend.key.height = unit(0.5, 'cm'),
      legend.key.width = unit(0.5, 'cm'),
      legend.title = element_text(size = 10),
      legend.text = element_text(size = 8),
      plot.title = element_text(size = 16),
      plot.subtitle = element_text(size = 12, hjust = 0.5),
      axis.text = element_text(size = 10),
      axis.title = element_text(size = 12)
    ) +
    xlab('Longitude') +
    ylab('Latitude') +
    easy_center_title()

  p6
}

plot(pdf1_future$pr[188,90+45,,9])
summary(pdf1_future$psl[188,90+45,,9])


{
  library(plotly)

  # Example grid point indices
  lon_index <- 180+28 # Example longitude index
  lat_index <- 90+46  # Example latitude index

  # Extract the 512-bin PDF vector for the specific grid point
  pdf_vector <- pdf3_future[lon_index, lat_index, ,3]

  # Reshape the PDF vector into a 3D array of dimensions [8, 8, 8]
  nbins <- 8
  pdf_3d <- array(pdf_vector, dim = c(nbins, nbins, nbins))

  # Extract the ranges for the variables from range_var_final
  range_var <- range_var_final$ranges
  var1_min <- range_var$pr[lon_index, lat_index, 1]
  var1_max <- range_var$pr[lon_index, lat_index, 2]
  var2_min <- range_var$tas[lon_index, lat_index, 1]
  var2_max <- range_var$tas[lon_index, lat_index, 2]
  var3_min <- range_var$psl[lon_index, lat_index, 1]
  var3_max <- range_var$psl[lon_index, lat_index, 2]

  # Create bin edges for each variable
  x_bins <- seq(var1_min, var1_max, length.out = nbins + 1)
  y_bins <- seq(var2_min, var2_max, length.out = nbins + 1)
  z_bins <- seq(var3_min, var3_max, length.out = nbins + 1)

  # Create the coordinates for the centers of the bins
  x_centers <- (x_bins[-1] + x_bins[-length(x_bins)]) / 2
  y_centers <- (y_bins[-1] + y_bins[-length(y_bins)]) / 2
  z_centers <- (z_bins[-1] + z_bins[-length(z_bins)]) / 2

  # Expand the grid of coordinates
  grid <- expand.grid(x = x_centers, y = y_centers, z = z_centers)

  # Flatten the PDF array into a vector
  pdf_flat <- as.vector(pdf_3d)

  # Combine the coordinates with the PDF values
  plot_data <- data.frame(
    x = grid$x,
    y = grid$y,
    z = grid$z,
    value = pdf_flat
  )

  # Normalize PDF values for marker size
  normalized_pdf <- pdf_flat / max(pdf_flat, na.rm = TRUE)

  # Plot the 3D histogram
  fig <- plot_ly(
    data = plot_data,
    x = ~x,
    y = ~y,
    z = ~z,
    type = "scatter3d",
    mode = "markers",
    marker = list(
      size = ~normalized_pdf * 75,  # Adjust size scaling factor
      color = ~value,
      colorscale = "Viridis",
      showscale = TRUE
    ),
    text = ~paste("PDF Value:", round(value, 4))
  ) %>%
    layout(
      scene = list(
        xaxis = list(title = "Variable 1 (pr)"),
        yaxis = list(title = "Variable 2 (tas)"),
        zaxis = list(title = "Variable 3 (psl)")
      ),
      title = paste("3D PDF for Grid Point (Lon:", lon_index, ", Lat:", lat_index, ")")
    )

  # Show the plot
  fig
}
