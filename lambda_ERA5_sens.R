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
  seq(0.0, 0.150, by = 0.025),   # 0.025, 0.05, 0.075, 0.10, 0.125, 0.15
  seq(0.20,  0.50,  by = 0.05),    # 0.20, 0.25, 0.30, 0.35, 0.40, 0.45, 0.50
  seq(0.60,  2.00,  by = 0.10)     # 0.60 → 2.00 in steps of 0.10
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

  message("Producing H-dist3 map for lambda = ", lam)

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
    labs(subtitle = "Projection period: 1998-2023") +
    ggtitle(sprintf("GraphCut (lambda = %.3f) - Mean H = %.3f", lam, global_h)) +
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

  print(p_hdist$labels$title)
  charToRaw(p_hdist$labels$title)

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
# Hellinger metrics vs lambda (GC curves) + MMM reference lines
#   Metrics:
#     1) simple global mean H
#     2) area-weighted mean H
#     3) simple global mean gradient H
#     4) area-weighted mean gradient H
# ======================================================================

# gradient_hdist() must be defined
# GC_maps must be filled
# MMM_hdist_fut must be [lon x lat] matrix of MMM H-dist (future)
{

  # ---------------- GC part: compute metrics per lambda ------------------

  valid_idx <- which(!sapply(GC_maps, is.null))

  cos_weights <- cos(lat * pi/180)
  weight_matrix <- matrix(
    rep(cos_weights, each = length(lon)),
    nrow = length(lon),
    byrow = FALSE
  )

  lambda_vec             <- numeric(length(valid_idx))
  mean_plain             <- numeric(length(valid_idx))
  mean_weighted          <- numeric(length(valid_idx))
  mean_gradient_plain    <- numeric(length(valid_idx))
  mean_gradient_weighted <- numeric(length(valid_idx))

  for (k in seq_along(valid_idx)) {

    idx   <- valid_idx[k]
    entry <- GC_maps[[idx]]
    H     <- entry$GC_hdist_fut

    lambda_vec[k] <- entry$lambda

    # mask NA for H
    mask_H <- !is.na(H)
    Hv     <- H[mask_H]
    w      <- weight_matrix
    wv     <- w[mask_H]

    # 1) simple global mean H
    mean_plain[k] <- mean(Hv)

    # 2) area-weighted mean H
    mean_weighted[k] <- sum(Hv * wv) / sum(wv)

    # 3) gradient map and simple global mean gradient
    grad_map <- gradient_hdist(H)
    mean_gradient_plain[k] <- mean(grad_map, na.rm = TRUE)

    # 4) area-weighted mean gradient
    mask_G  <- !is.na(grad_map)
    Gv      <- grad_map[mask_G]
    wG      <- weight_matrix[mask_G]

    mean_gradient_weighted[k] <- sum(Gv * wG) / sum(wG)
  }

  df_lambda <- data.frame(
    lambda                 = lambda_vec,
    mean_plain_fut         = mean_plain,
    mean_weighted_fut      = mean_weighted,
    mean_grad_plain_fut    = mean_gradient_plain,
    mean_grad_weighted_fut = mean_gradient_weighted
  )

  df_lambda <- df_lambda[order(df_lambda$lambda), ]

  df_long <- data.frame(
    lambda        = df_lambda$lambda,
    mean_H_plain  = df_lambda$mean_plain_fut,
    mean_H_weighted  = df_lambda$mean_weighted_fut,
    grad_plain       = df_lambda$mean_grad_plain_fut,
    grad_weighted    = df_lambda$mean_grad_weighted_fut
  )

  # scaling factor to map gradient metrics onto H scale
  scale_factor <- diff(range(df_long$mean_H_plain)) / diff(range(df_long$grad_plain))

  # ---------------- MMM reference metrics (constants) --------------------

  # 1) H metrics
  mask_H_MMM <- !is.na(MMM_hdist_fut)
  Hv_MMM     <- MMM_hdist_fut[mask_H_MMM]
  wv_MMM     <- weight_matrix[mask_H_MMM]

  MMM_mean_H_plain    <- mean(Hv_MMM)
  MMM_mean_H_weighted <- sum(Hv_MMM * wv_MMM) / sum(wv_MMM)

  # 2) gradient metrics
  grad_MMM <- gradient_hdist(MMM_hdist_fut)

  mask_G_MMM  <- !is.na(grad_MMM)
  Gv_MMM      <- grad_MMM[mask_G_MMM]
  wG_MMM      <- weight_matrix[mask_G_MMM]

  MMM_mean_grad_plain    <- mean(grad_MMM, na.rm = TRUE)
  MMM_mean_grad_weighted <- sum(Gv_MMM * wG_MMM) / sum(wG_MMM)

  # values for right axis (scaled)
  MMM_mean_grad_plain_scaled    <- MMM_mean_grad_plain * scale_factor
  MMM_mean_grad_weighted_scaled <- MMM_mean_grad_weighted * scale_factor

  # ---------------------------- Plot -------------------------------------

  out_dir <- "figure/Labelling/lambda_sweep_seed1/comprarison/"
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  p_lambda <- ggplot(df_long, aes(x = lambda)) +
    # GC: left axis (H)
    geom_line(aes(y = mean_H_plain,     colour = "Mean H"), size = 1) +
    geom_point(aes(y = mean_H_plain,    colour = "Mean H"), size = 2) +
    geom_line(aes(y = mean_H_weighted,  colour = "Area-weighted mean H"), size = 1) +
    geom_point(aes(y = mean_H_weighted, colour = "Area-weighted mean H"), size = 2) +

    # GC: right axis (gradients, scaled)
    geom_line(aes(y = grad_plain    * scale_factor, colour = "Mean gradient H"), size = 1) +
    geom_point(aes(y = grad_plain   * scale_factor, colour = "Mean gradient H"), size = 2) +
    geom_line(aes(y = grad_weighted * scale_factor, colour = "Area-weighted mean gradient H"), size = 1) +
    geom_point(aes(y = grad_weighted* scale_factor, colour = "Area-weighted mean gradient H"), size = 2) +

    # MMM: reference horizontal lines (thinner, dashed, same colours)
    geom_hline(yintercept = MMM_mean_H_plain,
               colour = "black", linetype = "dashed", linewidth = 0.4, show.legend = FALSE) +
    geom_hline(yintercept = MMM_mean_H_weighted,
               colour = "red",   linetype = "dashed", linewidth = 0.4, show.legend = FALSE) +
    geom_hline(yintercept = MMM_mean_grad_plain_scaled,
               colour = "blue",  linetype = "dashed", linewidth = 0.4, show.legend = FALSE) +
    geom_hline(yintercept = MMM_mean_grad_weighted_scaled,
               colour = "darkgreen", linetype = "dashed", linewidth = 0.4, show.legend = FALSE) +

    # axes
    scale_y_continuous(
      name = "Hellinger distance (future)",
      sec.axis = sec_axis(~./scale_factor,
                          name = "Gradient Hellinger (future)")
    ) +
    scale_colour_manual(
      values = c(
        "Mean H"                          = "black",
        "Area-weighted mean H"            = "red",
        "Mean gradient H"                 = "blue",
        "Area-weighted mean gradient H"   = "darkgreen"
      ),
      name = NULL
    ) +
    scale_x_continuous(
      breaks = df_long$lambda,
      labels = format(df_long$lambda, digits = 3)
    ) +
    xlab("Lambda (smoothness weight)") +
    ggtitle("Hellinger metrics vs lambda (seed = 1)\nGC curves and MMM reference") +
    theme_bw() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      plot.title  = element_text(size = 14),
      legend.position = "right"
    )

  print(p_lambda)

  base_name <- file.path(out_dir, "Hdist_grad_metrics_vs_lambda_with_MMM_seed1")
  ggsave(paste0(base_name, ".pdf"),
         plot = p_lambda, width = 18, height = 10, units = "cm", dpi = 300)
  ggsave(paste0(base_name, ".png"),
         plot = p_lambda, width = 18, height = 10, units = "cm", dpi = 300)




  # ======================================================================
  # Assumes df_lambda already exists as:
  # df_lambda <- data.frame(
  #   lambda                 = lambda_vec,
  #   mean_plain_fut         = mean_plain,
  #   mean_weighted_fut      = mean_weighted,
  #   mean_grad_plain_fut    = mean_gradient_plain,
  #   mean_grad_weighted_fut = mean_gradient_weighted
  # )
  # and GC_maps, MMM_hdist_fut, gradient_hdist(), lon, lat exist.
  # ======================================================================

  # Sort by lambda (just to be safe)
  df_lambda <- df_lambda[order(df_lambda$lambda), ]

  # Precompute weights for MMM
  cos_weights <- cos(lat * pi/180)
  weight_matrix <- matrix(
    rep(cos_weights, each = length(lon)),
    nrow = length(lon),
    byrow = FALSE
  )

  # ---------------- MMM metrics (constants) --------------------

  # 1) H metrics
  mask_H_MMM <- !is.na(MMM_hdist_fut)
  Hv_MMM     <- MMM_hdist_fut[mask_H_MMM]
  wv_MMM     <- weight_matrix[mask_H_MMM]

  MMM_mean_H_plain    <- mean(Hv_MMM)
  MMM_mean_H_weighted <- sum(Hv_MMM * wv_MMM) / sum(wv_MMM)

  # 2) gradient metrics
  grad_MMM <- gradient_hdist(MMM_hdist_fut)

  mask_G_MMM  <- !is.na(grad_MMM)
  Gv_MMM      <- grad_MMM[mask_G_MMM]
  wG_MMM      <- weight_matrix[mask_G_MMM]

  MMM_mean_grad_plain    <- mean(grad_MMM, na.rm = TRUE)
  MMM_mean_grad_weighted <- sum(Gv_MMM * wG_MMM) / sum(wG_MMM)

  # ---------------- Data frames for plotting ------------------

  df_H <- data.frame(
    lambda            = df_lambda$lambda,
    mean_H_plain      = df_lambda$mean_plain_fut,
    mean_H_weighted   = df_lambda$mean_weighted_fut
  )

  df_grad <- data.frame(
    lambda              = df_lambda$lambda,
    mean_grad_plain     = df_lambda$mean_grad_plain_fut,
    mean_grad_weighted  = df_lambda$mean_grad_weighted_fut
  )

  out_dir <- "figure/Labelling/lambda_sweep_seed1/comprarison/"
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  # ======================================================================
  # 1) Plot for H (no gradient here)
  # ======================================================================

  p_H <- ggplot(df_H, aes(x = lambda)) +
    geom_line(aes(y = mean_H_plain,    colour = "Mean H"), size = 1) +
    geom_point(aes(y = mean_H_plain,   colour = "Mean H"), size = 2) +
    geom_line(aes(y = mean_H_weighted, colour = "Area-weighted mean H"), size = 1) +
    geom_point(aes(y = mean_H_weighted, colour = "Area-weighted mean H"), size = 2) +

    # MMM reference lines (same colours, thinner, dashed)
    geom_hline(yintercept = MMM_mean_H_plain,
               colour = "black", linetype = "dashed", linewidth = 0.4, show.legend = FALSE) +
    geom_hline(yintercept = MMM_mean_H_weighted,
               colour = "red",   linetype = "dashed", linewidth = 0.4, show.legend = FALSE) +

    scale_colour_manual(
      values = c(
        "Mean H"               = "black",
        "Area-weighted mean H" = "red"
      ),
      name = NULL
    ) +
    scale_x_continuous(
      breaks = df_H$lambda,
      labels = format(df_H$lambda, digits = 3)
    ) +
    xlab("Lambda (smoothness weight)") +
    ylab("Hellinger distance (future)") +
    ggtitle("Mean Hellinger distance vs lambda (GC) with MMM reference\n(seed = 1)") +
    theme_bw() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      plot.title  = element_text(size = 14),
      legend.position = "right"
    )

  print(p_H)

  base_name_H <- file.path(out_dir, "H_metrics_vs_lambda_with_MMM_seed1")
  ggsave(paste0(base_name_H, ".pdf"),
         plot = p_H, width = 18, height = 10, units = "cm", dpi = 300)
  ggsave(paste0(base_name_H, ".png"),
         plot = p_H, width = 18, height = 10, units = "cm", dpi = 300)

  # ======================================================================
  # 2) Plot for gradient(H) only
  # ======================================================================

  p_grad <- ggplot(df_grad, aes(x = lambda)) +
    geom_line(aes(y = mean_grad_plain,    colour = "Mean gradient H"), size = 1) +
    geom_point(aes(y = mean_grad_plain,   colour = "Mean gradient H"), size = 2) +
    geom_line(aes(y = mean_grad_weighted, colour = "Area-weighted mean gradient H"), size = 1) +
    geom_point(aes(y = mean_grad_weighted, colour = "Area-weighted mean gradient H"), size = 2) +

    # MMM reference lines for gradients
    geom_hline(yintercept = MMM_mean_grad_plain,
               colour = "blue",      linetype = "dashed", linewidth = 0.4, show.legend = FALSE) +
    geom_hline(yintercept = MMM_mean_grad_weighted,
               colour = "darkgreen", linetype = "dashed", linewidth = 0.4, show.legend = FALSE) +

    scale_colour_manual(
      values = c(
        "Mean gradient H"               = "blue",
        "Area-weighted mean gradient H" = "darkgreen"
      ),
      name = NULL
    ) +
    scale_x_continuous(
      breaks = df_grad$lambda,
      labels = format(df_grad$lambda, digits = 3)
    ) +
    xlab("Lambda (smoothness weight)") +
    ylab("Gradient of Hellinger distance (future)") +
    ggtitle("Gradient Hellinger metrics vs lambda (GC) with MMM reference\n(seed = 1)") +
    theme_bw() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      plot.title  = element_text(size = 14),
      legend.position = "right"
    )

  print(p_grad)

  base_name_G <- file.path(out_dir, "Grad_H_metrics_vs_lambda_with_MMM_seed1")
  ggsave(paste0(base_name_G, ".pdf"),
         plot = p_grad, width = 18, height = 10, units = "cm", dpi = 300)
  ggsave(paste0(base_name_G, ".png"),
         plot = p_grad, width = 18, height = 10, units = "cm", dpi = 300)

}


# ======================================================================
# Hellinger metrics vs lambda (PRESENT) for GC + MMM reference
#   Metrics:
#     1) simple global mean H
#     2) area-weighted mean H
#     3) simple global mean gradient H
#     4) area-weighted mean gradient H
# ======================================================================
{
  # Select only valid entries
  valid_idx <- which(!sapply(GC_maps, is.null))

  # Precompute latitude weights (area ~ cos(lat))
  cos_weights <- cos(lat * pi/180)
  weight_matrix <- matrix(
    rep(cos_weights, each = length(lon)),
    nrow = length(lon),
    byrow = FALSE
  )

  lambda_vec_pres             <- numeric(length(valid_idx))
  mean_plain_pres             <- numeric(length(valid_idx))
  mean_weighted_pres          <- numeric(length(valid_idx))
  mean_gradient_plain_pres    <- numeric(length(valid_idx))
  mean_gradient_weighted_pres <- numeric(length(valid_idx))

  for (k in seq_along(valid_idx)) {

    idx   <- valid_idx[k]
    entry <- GC_maps[[idx]]
    H     <- entry$GC_hdist_pres   # [lon x lat] PRESENT map

    lambda_vec_pres[k] <- entry$lambda

    # mask NA for H
    mask_H <- !is.na(H)
    Hv     <- H[mask_H]
    w      <- weight_matrix
    wv     <- w[mask_H]

    # 1) simple global mean H
    mean_plain_pres[k] <- mean(Hv)

    # 2) area-weighted mean H
    mean_weighted_pres[k] <- sum(Hv * wv) / sum(wv)

    # 3) gradient map and simple global mean gradient
    grad_map <- gradient_hdist(H)
    mean_gradient_plain_pres[k] <- mean(grad_map, na.rm = TRUE)

    # 4) area-weighted mean gradient
    mask_G  <- !is.na(grad_map)
    Gv      <- grad_map[mask_G]
    wG      <- weight_matrix[mask_G]

    mean_gradient_weighted_pres[k] <- sum(Gv * wG) / sum(wG)
  }

  df_lambda_pres <- data.frame(
    lambda                    = lambda_vec_pres,
    mean_plain_pres           = mean_plain_pres,
    mean_weighted_pres        = mean_weighted_pres,
    mean_grad_plain_pres      = mean_gradient_plain_pres,
    mean_grad_weighted_pres   = mean_gradient_weighted_pres
  )

  df_lambda_pres <- df_lambda_pres[order(df_lambda_pres$lambda), ]

  # ---------------- MMM metrics for PRESENT (constants) ------------------

  # 1) H metrics (present)
  mask_H_MMM_pres <- !is.na(MMM_hdist_pres)
  Hv_MMM_pres     <- MMM_hdist_pres[mask_H_MMM_pres]
  wv_MMM_pres     <- weight_matrix[mask_H_MMM_pres]

  MMM_mean_H_plain_pres    <- mean(Hv_MMM_pres)
  MMM_mean_H_weighted_pres <- sum(Hv_MMM_pres * wv_MMM_pres) / sum(wv_MMM_pres)

  # 2) gradient metrics (present)
  grad_MMM_pres <- gradient_hdist(MMM_hdist_pres)

  mask_G_MMM_pres  <- !is.na(grad_MMM_pres)
  Gv_MMM_pres      <- grad_MMM_pres[mask_G_MMM_pres]
  wG_MMM_pres      <- weight_matrix[mask_G_MMM_pres]

  MMM_mean_grad_plain_pres    <- mean(grad_MMM_pres, na.rm = TRUE)
  MMM_mean_grad_weighted_pres <- sum(Gv_MMM_pres * wG_MMM_pres) / sum(wG_MMM_pres)

  # ---------------- Data frames for plotting ----------------------------

  df_H_pres <- data.frame(
    lambda              = df_lambda_pres$lambda,
    mean_H_plain_pres   = df_lambda_pres$mean_plain_pres,
    mean_H_weighted_pres= df_lambda_pres$mean_weighted_pres
  )

  df_grad_pres <- data.frame(
    lambda                  = df_lambda_pres$lambda,
    mean_grad_plain_pres    = df_lambda_pres$mean_grad_plain_pres,
    mean_grad_weighted_pres = df_lambda_pres$mean_grad_weighted_pres
  )

  out_dir <- "figure/Labelling/lambda_sweep_seed1"
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  # ======================================================================
  # 1) Plot for H (present)
  # ======================================================================

  p_H_pres <- ggplot(df_H_pres, aes(x = lambda)) +
    geom_line(aes(y = mean_H_plain_pres,    colour = "Mean H (present)"), size = 1) +
    geom_point(aes(y = mean_H_plain_pres,   colour = "Mean H (present)"), size = 2) +
    geom_line(aes(y = mean_H_weighted_pres, colour = "Area-weighted mean H (present)"), size = 1) +
    geom_point(aes(y = mean_H_weighted_pres, colour = "Area-weighted mean H (present)"), size = 2) +

    # MMM reference lines (present)
    geom_hline(yintercept = MMM_mean_H_plain_pres,
               colour = "black", linetype = "dashed", linewidth = 0.4, show.legend = FALSE) +
    geom_hline(yintercept = MMM_mean_H_weighted_pres,
               colour = "red",   linetype = "dashed", linewidth = 0.4, show.legend = FALSE) +

    scale_colour_manual(
      values = c(
        "Mean H (present)"               = "black",
        "Area-weighted mean H (present)" = "red"
      ),
      name = NULL
    ) +
    scale_x_continuous(
      breaks = df_H_pres$lambda,
      labels = format(df_H_pres$lambda, digits = 3)
    ) +
    xlab("Lambda (smoothness weight)") +
    ylab("Hellinger distance (present)") +
    ggtitle("Mean Hellinger distance vs lambda (present, GC) with MMM reference\n(seed = 1)") +
    theme_bw() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      plot.title  = element_text(size = 14),
      legend.position = "right"
    )

  print(p_H_pres)

  base_name_H_pres <- file.path(out_dir, "H_metrics_vs_lambda_with_MMM_present_seed1")
  ggsave(paste0(base_name_H_pres, ".pdf"),
         plot = p_H_pres, width = 18, height = 10, units = "cm", dpi = 300)
  ggsave(paste0(base_name_H_pres, ".png"),
         plot = p_H_pres, width = 18, height = 10, units = "cm", dpi = 300)

  # ======================================================================
  # 2) Plot for gradient(H) (present)
  # ======================================================================

  p_grad_pres <- ggplot(df_grad_pres, aes(x = lambda)) +
    geom_line(aes(y = mean_grad_plain_pres,    colour = "Mean gradient H (present)"), size = 1) +
    geom_point(aes(y = mean_grad_plain_pres,   colour = "Mean gradient H (present)"), size = 2) +
    geom_line(aes(y = mean_grad_weighted_pres, colour = "Area-weighted mean gradient H (present)"), size = 1) +
    geom_point(aes(y = mean_grad_weighted_pres, colour = "Area-weighted mean gradient H (present)"), size = 2) +

    # MMM reference lines (present)
    geom_hline(yintercept = MMM_mean_grad_plain_pres,
               colour = "blue",      linetype = "dashed", linewidth = 0.4, show.legend = FALSE) +
    geom_hline(yintercept = MMM_mean_grad_weighted_pres,
               colour = "darkgreen", linetype = "dashed", linewidth = 0.4, show.legend = FALSE) +

    scale_colour_manual(
      values = c(
        "Mean gradient H (present)"               = "blue",
        "Area-weighted mean gradient H (present)" = "darkgreen"
      ),
      name = NULL
    ) +
    scale_x_continuous(
      breaks = df_grad_pres$lambda,
      labels = format(df_grad_pres$lambda, digits = 3)
    ) +
    xlab("Lambda (smoothness weight)") +
    ylab("Gradient of Hellinger distance (present)") +
    ggtitle("Gradient Hellinger metrics vs lambda (present, GC) with MMM reference\n(seed = 1)") +
    theme_bw() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      plot.title  = element_text(size = 14),
      legend.position = "right"
    )

  print(p_grad_pres)

  base_name_G_pres <- file.path(out_dir, "Grad_H_metrics_vs_lambda_with_MMM_present_seed1")
  ggsave(paste0(base_name_G_pres, ".pdf"),
         plot = p_grad_pres, width = 18, height = 10, units = "cm", dpi = 300)
  ggsave(paste0(base_name_G_pres, ".png"),
         plot = p_grad_pres, width = 18, height = 10, units = "cm", dpi = 300)
}


# ======================================================================
# Energy components vs lambda from GC_results_lambda
#   - Data cost
#   - Smooth cost
#   - Total cost = Data + Smooth
#   - Smooth cost / lambda  (on a separate plot)
# ======================================================================

# GC_results_lambda is a list like:
#   $lambda_0.025
#       $label_attribution
#       $`Data and smooth cost`$`Data cost`
#       $`Data and smooth cost`$`Smooth cost`

# -------- 1) Extract lambda and costs ---------------------------------
{
  # indices of non-NULL results
  valid_idx <- which(!sapply(GC_results_lambda, is.null))

  lambda_vec    <- numeric(length(valid_idx))
  data_cost_vec <- numeric(length(valid_idx))
  smooth_cost_vec <- numeric(length(valid_idx))

  for (k in seq_along(valid_idx)) {

    idx   <- valid_idx[k]
    res_k <- GC_results_lambda[[idx]]

    # lambda from list name, e.g. "lambda_0.025"
    lam_str <- sub("^lambda_", "", names(GC_results_lambda)[idx])
    lambda_vec[k] <- as.numeric(lam_str)

    dc <- res_k[["Data and smooth cost"]][["Data cost"]]
    sc <- res_k[["Data and smooth cost"]][["Smooth cost"]]

    data_cost_vec[k]   <- dc
    smooth_cost_vec[k] <- sc
  }

  # total energy
  total_cost_vec <- data_cost_vec + smooth_cost_vec

  df_energy <- data.frame(
    lambda      = lambda_vec,
    data_cost   = data_cost_vec,
    smooth_cost = smooth_cost_vec,
    total_cost  = total_cost_vec
  )

  # sort by lambda
  df_energy <- df_energy[order(df_energy$lambda), ]

  # ======================================================================
  # 2) Plot: Data, Smooth, Total cost vs lambda
  # ======================================================================

  out_dir <- "figure/Energy_lambda_sweep_seed1"
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  p_energy <- ggplot(df_energy, aes(x = lambda)) +
    geom_line(aes(y = data_cost,   colour = "Data cost"),   size = 1) +
    geom_point(aes(y = data_cost,  colour = "Data cost"),   size = 2) +
    geom_line(aes(y = smooth_cost, colour = "Smooth cost"), size = 1) +
    geom_point(aes(y = smooth_cost, colour = "Smooth cost"), size = 2) +
    geom_line(aes(y = total_cost,  colour = "Total cost"),  size = 1) +
    geom_point(aes(y = total_cost, colour = "Total cost"),  size = 2) +

    scale_colour_manual(
      values = c(
        "Data cost"   = "black",
        "Smooth cost" = "red",
        "Total cost"  = "blue"
      ),
      name = NULL
    ) +
    scale_x_continuous(
      breaks = df_energy$lambda,
      labels = format(df_energy$lambda, digits = 3)
    ) +
    xlab("Lambda (smoothness weight)") +
    ylab("Energy") +
    ggtitle("GraphCut energy components vs lambda (seed = 1)") +
    theme_bw() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      plot.title  = element_text(size = 14),
      legend.position = "right"
    )

  print(p_energy)

  base_energy <- file.path(out_dir, "Energy_components_vs_lambda_seed1")
  ggsave(paste0(base_energy, ".pdf"),
         plot = p_energy, width = 18, height = 10, units = "cm", dpi = 300)
  ggsave(paste0(base_energy, ".png"),
         plot = p_energy, width = 18, height = 10, units = "cm", dpi = 300)

  # ======================================================================
  # 3) Plot: Smooth cost normalized by lambda
  # ======================================================================

  # exclude lambda = 0 to avoid division by zero
  mask_pos_lambda <- df_energy$lambda > 0

  df_smooth_norm <- data.frame(
    lambda        = df_energy$lambda[mask_pos_lambda],
    smooth_over_l = df_energy$smooth_cost[mask_pos_lambda] /
      df_energy$lambda[mask_pos_lambda]
  )

  p_smooth_norm <- ggplot(df_smooth_norm, aes(x = lambda, y = smooth_over_l)) +
    geom_line(colour = "red", size = 1) +
    geom_point(colour = "red", size = 2) +
    scale_x_continuous(
      breaks = df_smooth_norm$lambda,
      labels = format(df_smooth_norm$lambda, digits = 3)
    ) +
    xlab("Lambda (smoothness weight)") +
    ylab("Smooth cost / lambda") +
    ggtitle("Smooth energy per unit lambda vs lambda (seed = 1)") +
    theme_bw() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      plot.title  = element_text(size = 14)
    )

  print(p_smooth_norm)

  base_smooth <- file.path(out_dir, "Smooth_cost_per_lambda_vs_lambda_seed1")
  ggsave(paste0(base_smooth, ".pdf"),
         plot = p_smooth_norm, width = 18, height = 10, units = "cm", dpi = 300)
  ggsave(paste0(base_smooth, ".png"),
         plot = p_smooth_norm, width = 18, height = 10, units = "cm", dpi = 300)
}

# ======================================================================
# Energy components vs lambda from GC_results_lambda
#   - Data cost
#   - Smooth cost
#   - Total cost = Data + Smooth
#   - Smooth cost / lambda
#   - Data cost + Smooth cost / lambda
# ======================================================================
{
  valid_idx <- which(!sapply(GC_results_lambda, is.null))

  lambda_vec      <- numeric(length(valid_idx))
  data_cost_vec   <- numeric(length(valid_idx))
  smooth_cost_vec <- numeric(length(valid_idx))

  for (k in seq_along(valid_idx)) {

    idx   <- valid_idx[k]
    res_k <- GC_results_lambda[[idx]]

    # lambda from list name, e.g. "lambda_0.025"
    lam_str <- sub("^lambda_", "", names(GC_results_lambda)[idx])
    lambda_vec[k] <- as.numeric(lam_str)

    dc <- res_k[["Data and smooth cost"]][["Data cost"]]
    sc <- res_k[["Data and smooth cost"]][["Smooth cost"]]

    data_cost_vec[k]   <- dc
    smooth_cost_vec[k] <- sc
  }

  total_cost_vec <- data_cost_vec + smooth_cost_vec

  df_energy <- data.frame(
    lambda      = lambda_vec,
    data_cost   = data_cost_vec,
    smooth_cost = smooth_cost_vec,
    total_cost  = total_cost_vec
  )

  df_energy <- df_energy[order(df_energy$lambda), ]

  # --- normalized smooth and normalized total --------------------------

  df_energy$smooth_over_lambda <- with(df_energy,
                                       ifelse(lambda > 0, smooth_cost / lambda, NA_real_)
  )

  df_energy$total_norm <- df_energy$data_cost + df_energy$smooth_over_lambda

  # ======================================================================
  # Plot: all five curves in one figure
  # ======================================================================

  out_dir <- "figure/Energy_lambda_sweep_seed1"
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  p_energy <- ggplot(df_energy, aes(x = lambda)) +
    geom_line(aes(y = data_cost,         colour = "Data cost"),   size = 1) +
    geom_point(aes(y = data_cost,        colour = "Data cost"),   size = 2) +

    geom_line(aes(y = smooth_cost,       colour = "Smooth cost"), size = 1) +
    geom_point(aes(y = smooth_cost,      colour = "Smooth cost"), size = 2) +

    geom_line(aes(y = total_cost,        colour = "Total cost"),  size = 1) +
    geom_point(aes(y = total_cost,       colour = "Total cost"),  size = 2) +

    geom_line(aes(y = smooth_over_lambda, colour = "Smooth cost / lambda"), size = 0.8) +
    geom_point(aes(y = smooth_over_lambda, colour = "Smooth cost / lambda"), size = 1.6) +

    geom_line(aes(y = total_norm,        colour = "Data + (Smooth / lambda)"), size = 0.8) +
    geom_point(aes(y = total_norm,       colour = "Data + (Smooth / lambda)"), size = 1.6) +

    scale_colour_manual(
      values = c(
        "Data cost"                    = "black",
        "Smooth cost"                  = "red",
        "Total cost"                   = "blue",
        "Smooth cost / lambda"         = "darkorange",
        "Data + (Smooth / lambda)"     = "darkgreen"
      ),
      name = NULL
    ) +
    scale_x_continuous(
      breaks = df_energy$lambda,
      labels = format(df_energy$lambda, digits = 3)
    ) +
    xlab("Lambda (smoothness weight)") +
    ylab("Energy") +
    ggtitle("GraphCut energy components vs lambda (seed = 1)") +
    theme_bw() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      plot.title  = element_text(size = 14),
      legend.position = "right"
    )

  print(p_energy)

  base_energy <- file.path(out_dir, "Energy_components_vs_lambda_seed1")
  ggsave(paste0(base_energy, ".pdf"),
         plot = p_energy, width = 18, height = 10, units = "cm", dpi = 300)
  ggsave(paste0(base_energy, ".png"),
         plot = p_energy, width = 18, height = 10, units = "cm", dpi = 300)
}



# ----------------------------------------------------------------------
# Compute local label-noise map and global means
# labels: [lon x lat] integer or factor matrix
# Uses 4-neighbourhood (N,S,E,W).
# Returns:
#   map             : per-cell mean mismatch with neighbours in [0,1]
#   mean_plain      : simple global mean of map
#   mean_weighted   : area-weighted mean (cos(lat))
# ----------------------------------------------------------------------
{
  label_noise_metric <- function(labels, lon, lat) {
    nlon <- nrow(labels)
    nlat <- ncol(labels)

    noise  <- matrix(0, nlon, nlat)
    denom  <- matrix(0, nlon, nlat)

    ileft  <- 1:(nlon - 1)
    iright <- 2:nlon
    itop   <- 1:(nlat - 1)
    ibot   <- 2:nlat

    # horizontal neighbours
    diff_lr <- labels[ileft, ] != labels[iright, ]
    noise[ileft, ]  <- noise[ileft, ]  + diff_lr
    noise[iright, ] <- noise[iright, ] + diff_lr
    denom[ileft, ]  <- denom[ileft, ]  + 1
    denom[iright, ] <- denom[iright, ] + 1

    # vertical neighbours
    diff_tb <- labels[, itop] != labels[, ibot]
    noise[, itop] <- noise[, itop] + diff_tb
    noise[, ibot] <- noise[, ibot] + diff_tb
    denom[, itop] <- denom[, itop] + 1
    denom[, ibot] <- denom[, ibot] + 1

    frac_map <- noise / denom  # fraction of disagreeing neighbours (0–1)

    # simple mean
    mean_plain <- mean(frac_map, na.rm = TRUE)

    # area-weighted mean
    cos_weights <- cos(lat * pi/180)
    weight_matrix <- matrix(
      rep(cos_weights, each = length(lon)),
      nrow = length(lon),
      byrow = FALSE
    )

    mean_weighted <- sum(frac_map * weight_matrix, na.rm = TRUE) /
      sum(weight_matrix, na.rm = TRUE)

    list(
      map           = frac_map,
      mean_plain    = mean_plain,
      mean_weighted = mean_weighted
    )
  }

  # ======================================================================
  # Label-noise vs lambda (purely structural diagnostic)
  # ======================================================================

  valid_idx <- which(!sapply(GC_results_lambda, is.null))

  lambda_vec        <- numeric(length(valid_idx))
  ln_mean_plain     <- numeric(length(valid_idx))
  ln_mean_weighted  <- numeric(length(valid_idx))

  for (k in seq_along(valid_idx)) {
    idx <- valid_idx[k]
    res_k <- GC_results_lambda[[idx]]

    # extract lambda from list name, e.g. "lambda_0.025"
    lam_str <- sub("^lambda_", "", names(GC_results_lambda)[idx])
    lambda_vec[k] <- as.numeric(lam_str)

    labels <- res_k$label_attribution  # [lon x lat]

    ln <- label_noise_metric(labels, lon = lon, lat = lat)

    ln_mean_plain[k]    <- ln$mean_plain
    ln_mean_weighted[k] <- ln$mean_weighted
  }

  df_label_noise <- data.frame(
    lambda              = lambda_vec,
    noise_plain         = ln_mean_plain,
    noise_weighted      = ln_mean_weighted
  )

  df_label_noise <- df_label_noise[order(df_label_noise$lambda), ]


  out_dir <- "figure/LabellingNoise_lambda_sweep_seed1"
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  p_lab <- ggplot(df_label_noise, aes(x = lambda)) +
    geom_line(aes(y = noise_plain,    colour = "Mean label-noise"), size = 1) +
    geom_point(aes(y = noise_plain,   colour = "Mean label-noise"), size = 2) +
    geom_line(aes(y = noise_weighted, colour = "Area-weighted mean label-noise"), size = 1) +
    geom_point(aes(y = noise_weighted, colour = "Area-weighted mean label-noise"), size = 2) +
    scale_colour_manual(
      values = c(
        "Mean label-noise"                = "black",
        "Area-weighted mean label-noise"  = "red"
      ),
      name = NULL
    ) +
    scale_x_continuous(
      breaks = df_label_noise$lambda,
      labels = format(df_label_noise$lambda, digits = 3)
    ) +
    xlab("Lambda (smoothness weight)") +
    ylab("Label-noise metric (fraction of disagreeing neighbours)") +
    ggtitle("Label-noise vs lambda (GC label maps)\n(seed = 1)") +
    theme_bw() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      plot.title  = element_text(size = 14),
      legend.position = "right"
    )

  print(p_lab)

  base_name <- file.path(out_dir, "Label_noise_vs_lambda_seed1")
  ggsave(paste0(base_name, ".pdf"),
         plot = p_lab, width = 18, height = 10, units = "cm", dpi = 300)
  ggsave(paste0(base_name, ".png"),
         plot = p_lab, width = 18, height = 10, units = "cm", dpi = 300)

}


# df_label_noise has: lambda, noise_weighted (or noise_plain)

N0 <- df_label_noise$noise_weighted[1]        # baseline at smallest lambda
r  <- df_label_noise$noise_weighted / N0
L_noise <- 10 * log10(r)

df_label_noise$rel_noise  <- r
df_label_noise$L_noise_dB <- L_noise

# Example target: -6 dB  (~75% reduction in noise power)
target_dB <- -6

cand_idx <- which(df_label_noise$L_noise_dB <= target_dB)
lambda_star <- df_label_noise$lambda[min(cand_idx)]
lambda_star


# avoid lambda = 0; use λ >= λ_min
fit_idx <- df_label_noise$lambda > 0
fit <- lm(log(noise_weighted) ~ lambda, data = df_label_noise[fit_idx, ])

lambda0 <- -1 / coef(fit)[["lambda"]]   # characteristic decay scale


df <- df_label_noise[order(df_label_noise$lambda), ]

lambda <- df$lambda
noise  <- df$noise_weighted

# ------------------------------------------------------------
# Rational decay model: N(lambda) = A / (1 + B*lambda) + C
# ------------------------------------------------------------

# Initial guesses
A0 <- max(noise) - min(noise)
B0 <- 5
C0 <- min(noise)

model <- nls(
  noise ~ A / (1 + B * lambda) + C,
  start = list(A = A0, B = B0, C = C0),
  control = list(maxiter = 500, warnOnly = TRUE)
)

summary(model)

# Extract fitted parameters
A <- coef(model)[["A"]]
B <- coef(model)[["B"]]
C <- coef(model)[["C"]]

# Fitted curve
lambda_dense <- seq(min(lambda), max(lambda), length.out = 300)
noise_fit <- A / (1 + B * lambda_dense) + C

# Plot
plot(lambda, noise, pch=19, col="black", cex=1.2,
     xlab="Lambda", ylab="Label-noise",
     main="Rational fit of label-noise vs lambda")
lines(lambda_dense, noise_fit, col="red", lwd=2)

legend("topright",
       legend=c("Empirical noise", "Rational fit"),
       col=c("black", "red"), pch=c(19, NA), lwd=c(NA,2))


lambda_from_frac <- function(frac, A, B, C) {
  # frac = N(λ*) / N0  (e.g. 0.25 → 75% reduction)
  (1/frac - 1) / B
}

lambda_50 <- lambda_from_frac(0.50, A, B, C)  # -3 dB (≈50% noise)
lambda_25 <- lambda_from_frac(0.25, A, B, C)  # -6 dB (25% noise)
lambda_10 <- lambda_from_frac(0.10, A, B, C)  # -10 dB (10% noise)

lambda_from_dB <- function(dB, A, B, C) {
  frac <- 10^(dB / 10)      # N(λ*) / N0
  (1/frac - 1) / B
}

lambda_m3  <- lambda_from_dB(-3,  A, B, C)  # ≈ 50% noise
lambda_m6  <- lambda_from_dB(-6,  A, B, C)  # ≈ 25% noise
lambda_m10 <- lambda_from_dB(-10, A, B, C)  # ≈ 10% noise


target_dB <- -6
lambda_star <- lambda_from_dB(target_dB, A, B, C)
noise_star  <- A / (1 + B * lambda_star) + C

# existing plot of points + fitted curve here...

abline(v = lambda_star, col = "blue", lty = 2)
abline(h = noise_star,  col = "blue", lty = 2)
text(lambda_star, par("usr")[4],
     labels = sprintf("λ* = %.3f (%.0f dB)", lambda_star, target_dB),
     pos = 3, col = "blue")


# ============================================================
# 0) Data
# ============================================================
df <- df_label_noise[order(df_label_noise$lambda), ]
lambda <- df$lambda
noise  <- df$noise_weighted

# ============================================================
# 1) Fit rational decay: N(λ) = A / (1 + B λ) + C
# ============================================================

A0 <- max(noise) - min(noise)
B0 <- 5
C0 <- min(noise)

model <- nls(
  noise ~ A / (1 + B * lambda) + C,
  start = list(A = A0, B = B0, C = C0),
  control = list(maxiter = 500, warnOnly = TRUE)
)

co <- coef(model)
A <- co["A"]; B <- co["B"]; C <- co["C"]

# baseline noise at λ = 0
N0 <- A + C

# fitted curve for plotting
lambda_dense <- seq(min(lambda), max(lambda), length.out = 300)
noise_fit    <- A / (1 + B * lambda_dense) + C

# ============================================================
# 2) Helper: λ* from target dB
# ============================================================

lambda_from_dB <- function(dB, A, B, C) {
  frac <- 10^(dB / 10)          # N(λ*) / N0
  (1/frac - 1) / B              # λ* = (1/f - 1) / B
}

targets_dB <- c(-3, -6, -10)

lambda_star <- sapply(targets_dB, lambda_from_dB, A = A, B = B, C = C)
names(lambda_star) <- paste0(targets_dB, " dB")

noise_star  <- A / (1 + B * lambda_star) + C   # N(λ*) on fitted curve

# ============================================================
# 3) Plot: empirical noise + fit + dB thresholds + λ values
# ============================================================

plot(lambda, noise,
     pch = 19, col = "black", cex = 1.2,
     xlab = "Lambda",
     ylab = "Label-noise",
     main = "Rational fit of label-noise vs lambda\nwith -3/-6/-10 dB thresholds")

lines(lambda_dense, noise_fit, col = "red", lwd = 2)

legend("topright",
       legend = c("Empirical noise", "Rational fit"),
       col    = c("black", "red"),
       pch    = c(19, NA),
       lwd    = c(NA, 2),
       bty    = "n")

cols <- c("-3 dB" = "blue", "-6 dB" = "darkgreen", "-10 dB" = "purple")

y_top <- par("usr")[4]
y_min <- par("usr")[3]

for (i in seq_along(targets_dB)) {

  dB       <- targets_dB[i]
  lab_dB   <- paste0(dB, " dB")
  lam_star <- lambda_star[i]
  N_star   <- noise_star[i]
  coli     <- cols[lab_dB]

  # horizontal line at noise level
  abline(h = N_star, col = coli, lty = 2)

  # vertical line at lambda*
  abline(v = lam_star, col = coli, lty = 2)

  # label at the top: "-6 dB"
  text(x = lam_star, y = y_top,
       labels = lab_dB,
       col = coli, pos = 3, cex = 0.8)

  # label with lambda value a bit above the horizontal line
  text(x = lam_star, y = N_star,
       labels = sprintf("lambda = %.3f", lam_star),
       col = coli, pos = 3, cex = 0.8)
}



# ----------------------------------------------------------------------
# Second-order label oscillation metric (step = 2)
# labels: [lon x lat] matrix of integer/factor labels
# lon, lat: numeric vectors
# Returns:
#   map           : [lon x lat] 0/1 oscillation flag
#   mean_plain    : simple mean over grid
#   mean_weighted : area-weighted mean (cos(lat))
# ----------------------------------------------------------------------
{
label_oscillation_metric <- function(labels, lon, lat, step = 2) {

  nlon <- nrow(labels)
  nlat <- ncol(labels)

  osc_flag <- matrix(0L, nlon, nlat)

  # ---- horizontal direction: (i-step, j), (i, j), (i+step, j) ----
  if (nlon > 2 * step) {
    for (i in (1 + step):(nlon - step)) {
      L1 <- labels[i - step, ]
      L2 <- labels[i,       ]
      L3 <- labels[i + step, ]

      change12 <- (L1 != L2)
      change23 <- (L2 != L3)

      bad <- change12 & change23 & !is.na(L1) & !is.na(L2) & !is.na(L3)

      # mark central cell as oscillatory where pattern has two changes
      osc_flag[i, bad] <- 1L
    }
  }

  # ---- vertical direction: (i, j-step), (i, j), (i, j+step) ----
  if (nlat > 2 * step) {
    for (j in (1 + step):(nlat - step)) {
      L1 <- labels[, j - step]
      L2 <- labels[, j       ]
      L3 <- labels[, j + step]

      change12 <- (L1 != L2)
      change23 <- (L2 != L3)

      bad <- change12 & change23 & !is.na(L1) & !is.na(L2) & !is.na(L3)

      # mark central cell as oscillatory where pattern has two changes
      osc_flag[bad, j] <- 1L
    }
  }

  # simple mean
  mean_plain <- mean(osc_flag, na.rm = TRUE)

  # area-weighted mean (cos(lat))
  cos_weights <- cos(lat * pi/180)
  weight_matrix <- matrix(
    rep(cos_weights, each = length(lon)),
    nrow = length(lon),
    byrow = FALSE
  )

  mask <- !is.na(osc_flag)
  ov   <- osc_flag[mask]
  wv   <- weight_matrix[mask]

  mean_weighted <- sum(ov * wv) / sum(wv)

  list(
    map           = osc_flag,
    mean_plain    = mean_plain,
    mean_weighted = mean_weighted
  )
}

  # ======================================================================
  # Second-order oscillation vs lambda
  # ======================================================================

  valid_idx <- which(!sapply(GC_results_lambda, is.null))

  lambda_vec_osc       <- numeric(length(valid_idx))
  osc_mean_plain       <- numeric(length(valid_idx))
  osc_mean_weighted    <- numeric(length(valid_idx))

  for (k in seq_along(valid_idx)) {
    idx   <- valid_idx[k]
    res_k <- GC_results_lambda[[idx]]

    # lambda from list name, e.g. "lambda_0.025"
    lam_str <- sub("^lambda_", "", names(GC_results_lambda)[idx])
    lambda_vec_osc[k] <- as.numeric(lam_str)

    labels <- res_k$label_attribution  # [lon x lat]

    osc <- label_oscillation_metric(labels, lon = lon, lat = lat, step = 2)

    osc_mean_plain[k]    <- osc$mean_plain
    osc_mean_weighted[k] <- osc$mean_weighted
  }

  df_label_osc <- data.frame(
    lambda          = lambda_vec_osc,
    osc_plain       = osc_mean_plain,
    osc_weighted    = osc_mean_weighted
  )

  df_label_osc <- df_label_osc[order(df_label_osc$lambda), ]

  out_dir <- "figure/LabellingOscillation_lambda_sweep_seed1"
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  p_osc <- ggplot(df_label_osc, aes(x = lambda)) +
    geom_line(aes(y = osc_plain,    colour = "Mean oscillation"), size = 1) +
    geom_point(aes(y = osc_plain,   colour = "Mean oscillation"), size = 2) +
    geom_line(aes(y = osc_weighted, colour = "Area-weighted mean oscillation"), size = 1) +
    geom_point(aes(y = osc_weighted, colour = "Area-weighted mean oscillation"), size = 2) +
    scale_colour_manual(
      values = c(
        "Mean oscillation"                = "black",
        "Area-weighted mean oscillation"  = "red"
      ),
      name = NULL
    ) +
    scale_x_continuous(
      breaks = df_label_osc$lambda,
      labels = format(df_label_osc$lambda, digits = 3)
    ) +
    xlab("Lambda (smoothness weight)") +
    ylab("Oscillation metric (2-step label changes)") +
    ggtitle("Second-order label oscillation vs lambda (step = 2)") +
    theme_bw() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      plot.title  = element_text(size = 14),
      legend.position = "right"
    )

  print(p_osc)

  base_name_osc <- file.path(out_dir, "Label_oscillation_vs_lambda_seed1")
  ggsave(paste0(base_name_osc, ".pdf"),
         plot = p_osc, width = 18, height = 10, units = "cm", dpi = 300)
  ggsave(paste0(base_name_osc, ".png"),
         plot = p_osc, width = 18, height = 10, units = "cm", dpi = 300)


  # ======================================================================
  # Compute -3 / -6 / -10 dB lambda thresholds from empirical curve
  # ======================================================================

  lambda_vals <- df_label_osc$lambda
  O_vals      <- df_label_osc$osc_weighted   # or osc_plain

  # baseline at smallest lambda
  O0 <- O_vals[1]

  targets_dB <- c(-3, -6, -10)

  compute_lambda_dB_empirical <- function(target_dB, lambda_vals, O_vals, O0) {
    # target oscillation
    T <- O0 * 10^(target_dB / 10)

    # find where the curve crosses T
    idx <- which(O_vals <= T)

    if (length(idx) == 0) {
      return(NA_real_)  # never crossed
    }

    i <- idx[1]   # first lambda where O <= T

    if (i == 1) {
      return(lambda_vals[1])
    }

    # linear interpolation between (i-1) and i
    l1 <- lambda_vals[i-1];  o1 <- O_vals[i-1]
    l2 <- lambda_vals[i];    o2 <- O_vals[i]

    # linear interpolation
    slope <- (o2 - o1) / (l2 - l1)
    lambda_star <- l1 + (T - o1) / slope

    lambda_star
  }

  lambda_star_emp <- sapply(
    targets_dB,
    compute_lambda_dB_empirical,
    lambda_vals = lambda_vals,
    O_vals      = O_vals,
    O0          = O0
  )

  names(lambda_star_emp) <- paste0(targets_dB, " dB")
  lambda_star_emp

}





# ----------------------------------------------------------------------
# patch_stats_small()
#   labels      : matrix [lon x lat] of integer/factor labels
#   min_size    : threshold size (in cells) for "small" patches
# Returns:
#   sizes       : vector of patch sizes (all patches)
#   n_small     : number of patches < min_size
#   frac_small  : n_small / total number of patches
# ----------------------------------------------------------------------
{
  # ----------------------------------------------------------------------
  # patch_stats_small_manual()
  #
  # labels      : matrix [lon x lat] of integer/factor labels
  #               (optionally with 0 or NA as "no land / no label")
  # min_size    : threshold size (in cells) for "small" patches
  # ignore_val  : value to be treated as background (e.g. 0), set to NULL if none
  #
  # Returns:
  #   sizes      : vector of patch sizes (all patches, in cells)
  #   n_small    : number of patches with size < min_size
  #   frac_small : n_small / total number of patches
  # ----------------------------------------------------------------------
  patch_stats_small_manual <- function(labels,
                                       min_size   = 4,
                                       ignore_val = 0) {

    nlon <- nrow(labels)
    nlat <- ncol(labels)

    # treat ignore_val as NA if specified
    if (!is.null(ignore_val)) {
      labels[labels == ignore_val] <- NA
    }

    visited <- matrix(FALSE, nlon, nlat)
    sizes   <- integer(0)

    is_valid <- function(i, j, lab) {
      i >= 1 && i <= nlon &&
        j >= 1 && j <= nlat &&
        !visited[i, j] &&
        !is.na(labels[i, j]) &&
        labels[i, j] == lab
    }

    for (i0 in 1:nlon) {
      for (j0 in 1:nlat) {

        if (visited[i0, j0] || is.na(labels[i0, j0])) next

        lab <- labels[i0, j0]

        qi <- i0
        qj <- j0
        visited[i0, j0] <- TRUE
        head <- 1L
        size <- 0L

        while (head <= length(qi)) {
          ci <- qi[head]
          cj <- qj[head]
          head <- head + 1L
          size <- size + 1L

          nbs <- list(
            c(ci - 1L, cj),
            c(ci + 1L, cj),
            c(ci, cj - 1L),
            c(ci, cj + 1L)
          )

          for (nb in nbs) {
            ni <- nb[1]
            nj <- nb[2]
            if (is_valid(ni, nj, lab)) {
              visited[ni, nj] <- TRUE
              qi <- c(qi, ni)
              qj <- c(qj, nj)
            }
          }
        }

        sizes <- c(sizes, size)
      }
    }

    if (length(sizes) == 0L) {
      n_small          <- 0L
      frac_small_patches <- 0
      frac_small_cells   <- 0
    } else {
      n_small            <- sum(sizes < min_size)
      frac_small_patches <- n_small / length(sizes)

      total_cells        <- sum(!is.na(labels))
      small_cells        <- sum(sizes[sizes < min_size])
      frac_small_cells   <- small_cells / total_cells
    }

    list(
      sizes              = sizes,
      n_small            = n_small,
      frac_small_patches = frac_small_patches,  # old metric (if you still want it)
      frac_small_cells   = frac_small_cells     # new, area-like metric
    )
  }



  # ======================================================================
  # Count small patches vs lambda  (manual connected components)
  # ======================================================================

  min_patch_size <- 9   # or 9, etc.

  valid_idx <- which(!sapply(GC_results_lambda, is.null))

  lambda_vec        <- numeric(length(valid_idx))
  frac_small_cells  <- numeric(length(valid_idx))

  for (k in seq_along(valid_idx)) {

    idx   <- valid_idx[k]
    res_k <- GC_results_lambda[[idx]]

    lam_str <- sub("^lambda_", "", names(GC_results_lambda)[idx])
    lambda_vec[k] <- as.numeric(lam_str)

    labels <- res_k$label_attribution

    ps <- patch_stats_small_manual(labels,
                                   min_size   = min_patch_size,
                                   ignore_val = 0)

    frac_small_cells[k] <- ps$frac_small_cells
  }

  df_patches <- data.frame(
    lambda           = lambda_vec,
    frac_small_cells = frac_small_cells
  )
  df_patches <- df_patches[order(df_patches$lambda), ]



  out_dir <- "figure/LabellingPatchSize_lambda_sweep_seed1"
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  p_patch <- ggplot(df_patches, aes(x = lambda, y = frac_small_cells)) +
    geom_line(colour = "black") +
    geom_point(colour = "black", size = 2) +
    scale_x_continuous(
      breaks = df_patches$lambda,
      labels = format(df_patches$lambda, digits = 3)
    ) +
    xlab("Lambda (smoothness weight)") +
    ylab(sprintf("Fraction of grid cells in patches < %d cells", min_patch_size)) +
    ggtitle(sprintf("Area fraction of very small patches vs lambda\n(min size = %d cells)", min_patch_size)) +
    theme_bw() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      plot.title  = element_text(size = 14)
    )


  print(p_patch)


  base_name <- file.path(out_dir, sprintf("Small_patches_vs_lambda_min%d", min_patch_size))
  ggsave(paste0(base_name, ".pdf"),
         plot = p_patch, width = 18, height = 10, units = "cm", dpi = 300)
  ggsave(paste0(base_name, ".png"),
         plot = p_patch, width = 18, height = 10, units = "cm", dpi = 300)

  frac0 <- df_patches$frac_small[1]
  target_ratio <- 0.10           # keep only 10% of original tiny patches

  target_frac <- frac0 * target_ratio

  idx <- which(df_patches$frac_small <= target_frac)
  lambda_star_small <- if (length(idx) == 0L) NA_real_ else df_patches$lambda[min(idx)]
  lambda_star_small
}
