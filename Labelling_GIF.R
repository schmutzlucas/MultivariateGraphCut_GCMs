load('202502140016_my_workspace_ERA5_allModels_final_results.RData')

# ------------------------------------------------------------
# 0.  Packages ------------------------------------------------
# ------------------------------------------------------------
pkgs <- c("ggplot2", "reshape2", "gganimate", "viridisLite")
new  <- pkgs[!pkgs %in% installed.packages()[, "Package"]]
if(length(new)) install.packages(new, repos = "https://cloud.r-project.org")
lapply(pkgs, library, character.only = TRUE)

library(devtools)
lapply(list_of_packages, library, character.only = TRUE)
install_github("schmutzlucas/gcoWrapR")


# Loading local functions
source_code_dir <- 'functions/' #The directory where all functions are saved.
file_paths <- list.files(source_code_dir, full.names = T)
for(path in file_paths){source(path)}

# ------------------------------------------------------------
# 1. Run the graph-cut one swap at a time, keep every frame ---
# ------------------------------------------------------------
n_frames     <- 8                     # 0 … 7
smooth_cost  <- 0.3
GC_results   <- vector("list", n_frames)

lab_prev     <- NULL                  # first run → random map
cum_iter     <- 0                     # cumulative α-β swaps so far

for (k in seq_len(n_frames)) {

  iterations <- if (k == 1) 0 else 1  # 0 for the raw random map, then 1

  res <- GraphCutHellinger_nD_lat_gif(
    pdf_models_future = pdf_models_future,
    h_dist            = h_dist,
    weight_data       = 1,
    weight_smooth     = smooth_cost,
    nBins             = nbins1d^3,
    iterations        = iterations,
    labelling         = lab_prev,
    lat               = lat,
    seed              = 1,
    verbose           = TRUE,
    rebuild           = (k == 1)        # build XPtrs only once
  )

  cum_iter             <- cum_iter + iterations   # update total passes
  GC_results[[k]]      <- res                     # store in order
  names(GC_results)[k] <- cum_iter                # "0","1","2",…

  lab_prev <- res$label_attribution               # next initialise
  gc()
}



# ------------------------------------------------------------
# 2.  Make one frame at a time  -------------------------------
# ------------------------------------------------------------
library(ggplot2)
library(viridisLite)

# your colour palette ---------------------------------------------------
color_palette <- c(
  "#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd",
  "#8c564b", "#e377c2", "#7f7f7f", "#bcbd22", "#17becf",
  "#aec7e8", "#ffbb78", "#98df8a", "#ff9896", "#c5b0d5",
  "#c49c94", "#f7b6d2", "#c7c7c7", "#dbdb8d", "#9edae5",
  "#393b79", "#5254a3", "#6b6ecf"
)
names(color_palette) <- model_names     # 22 names ←→ 22 colours

# helper that wraps longitudes 0…359 into –180…180 ----------------------
wrap_lon <- function(x) ifelse(x > 180, x - 360, x)

# function that builds **your** static map ------------------------------
make_frame <- function(lab_matrix, title_text = NULL) {

  df <- reshape2::melt(lab_matrix,
                       varnames = c("lon_idx","lat_idx"),
                       value.name = "label_attribution")

  df$lon  <- wrap_lon(lon[df$lon_idx])
  df$lat  <- lat[df$lat_idx]
  df$label_attribution <- factor(df$label_attribution,
                                 levels = seq_along(model_names),
                                 labels = model_names)

  ggplot() +
    geom_tile(data = df,
              aes(x = lon, y = lat, fill = label_attribution)) +
    scale_fill_manual(values = color_palette,
                      guide = "none", drop = FALSE) +
    borders("world", colour = "black", size = 0.12) +
    coord_fixed(ratio = 1.3,
                xlim = c(-180, 180),
                ylim = c(-84,  90),
                expand = FALSE) +
    theme_bw() +
    theme(panel.grid = element_blank(),
          axis.title  = element_blank(),
          axis.text   = element_blank(),
          axis.ticks  = element_blank(),
          plot.title  = element_text(hjust = .5, size = 12)) +
    labs(title = if (is.null(title_text)) "" else title_text)
}

# # ------------------------------------------------------------
# # 3.  Render PNGs and assemble the GIF ------------------------
# # ------------------------------------------------------------
# png_dir <- file.path(tempdir(), "frames_gc")
# dir.create(png_dir, showWarnings = FALSE)
#
# # Full-HD size expressed in inches at 300 dpi
# dpi_png <- 300
# w_in    <- 1920 / dpi_png          # 6.4″
# h_in    <- 1080 / dpi_png          # 3.6″
#
# for (k in seq_along(GC_results)) {
#   lab  <- GC_results[[k]]$label_attribution
#   ggsave(file.path(png_dir, sprintf("frame_%03d.png", k - 1)),
#          plot   = make_frame(lab, sprintf("Graphcut iteration %d", k - 1)),
#          width  = w_in, height = h_in,
#          units  = "in", dpi = dpi_png)
# }
#
#
#
# # turn the stack into a GIF (or MP4) ------------------------------------
# png_files <- list.files(png_dir, pattern = "png$", full.names = TRUE)
#
# w_px <- 1920
# h_px <- 1080            # or 960 if you want a perfect 2:1 world map
#
# gifski(png_files,
#        gif_file = "figure/graphcut_evolution.gif",
#        width    = w_px,          # <-- keep these two lines
#        height   = h_px,
#        delay    = 1)
#
#
# # For MP4 instead of GIF:
# # av::av_encode_video(png_files, "figure/graphcut_evolution.mp4",
# #                     framerate = 1)
#
# cat("GIF written to", normalizePath("figure/graphcut_evolution030.gif"), "\n")



# ------------------------------------------------------------
# 3.  Render PNGs – one file per frame  -----------------------
# ------------------------------------------------------------

# create the output folder once
png_dir <- file.path("figure", "GraphCutEvolution")
dir.create(png_dir, recursive = TRUE, showWarnings = FALSE)

# Full-HD raster size: 1920 × 1080 px @ 300 dpi
dpi_png <- 300
w_in    <- 1920 / dpi_png    # 6.4 in
h_in    <- 1080 / dpi_png    # 3.6 in

for (k in seq_along(GC_results)) {
  lab  <- GC_results[[k]]$label_attribution
  file <- file.path(png_dir, sprintf("frame_%03d.png", k - 1))

  ggsave(filename = file,
         plot     = make_frame(lab, sprintf("Graph-cut iteration %d", k - 1)),
         width    = w_in,
         height   = h_in,
         units    = "in",
         dpi      = dpi_png)
}

cat("✅", length(GC_results), "PNG frames written to",
    normalizePath(png_dir), "\n")
