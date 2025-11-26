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



