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
# 1.  Run the graph-cut for the chosen iteration counts -------
# ------------------------------------------------------------
iters        <- c(0, 1, 5, 10)   # frames you want
smooth_cost  <- 0.5
GC_results   <- list()

for(i in seq_along(iters)){
  it <- iters[i]
  GC_results[[as.character(it)]] <- GraphCutHellinger_nD_lat_gif(
    pdf_models_future = pdf_models_future,
    h_dist            = h_dist,
    weight_data       = 1,
    weight_smooth     = smooth_cost,
    nBins             = nbins1d^3,
    iterations        = it,
    lat               = lat,
    seed              = 1,
    verbose           = TRUE,
    rebuild           = (i == 1)     # rebuild only on the first call
  )
  gc()                              # tidy up between runs
}

# ------------------------------------------------------------
# 2.  Convert label matrices to a single long data frame ------
# ------------------------------------------------------------
# fixed colour palette (same as your example, truncated to n models)
color_palette <- c(
  "#1f77b4","#ff7f0e","#2ca02c","#d62728","#9467bd",
  "#8c564b","#e377c2","#7f7f7f","#bcbd22","#17becf",
  "#aec7e8","#ffbb78","#98df8a","#ff9896","#c5b0d5",
  "#c49c94","#f7b6d2","#c7c7c7","#dbdb8d","#9edae5",
  "#393b79","#5254a3","#6b6ecf"
)[seq_along(model_names)]
names(color_palette) <- model_names

lon_vec <- lon   # for clarity
lat_vec <- lat

frame_df <- do.call(
  rbind,
  lapply(names(GC_results), function(step){
    lab <- GC_results[[step]]$label_attribution
    df  <- melt(lab, varnames = c("lon_idx","lat_idx"), value.name = "label")
    df$lon   <- lon_vec[df$lon_idx]
    df$lat   <- lat_vec[df$lat_idx]
    df$label <- factor(df$label,
                       levels = seq_along(model_names),
                       labels = model_names)
    df$step  <- sprintf("%05d", as.integer(step))   # nice ordering
    df
  })
)

# ------------------------------------------------------------
# 3.  Build the animation ------------------------------------
# ------------------------------------------------------------
plot_base <- ggplot(frame_df, aes(lon, lat, fill = label)) +
  geom_tile() +
  scale_fill_manual(values = color_palette, guide = "none") +  # <- no legend
  coord_fixed(ratio = 1.3,
              xlim  = c(-180, 180),
              ylim  = c(-84,  90),
              expand = FALSE) +
  theme_void()   # <- removes axes, background, titles, etc.

anim <- plot_base +
  transition_states(step, transition_length = 0, state_length = 1) +
  enter_fade() + exit_fade()

animate(anim,
        fps      = 2,                           # 2 frames / s  = 2 s total
        nframes  = length(iters),
        width    = 1000,
        height   = 450,
        renderer = gifski_renderer("graphcut_evolution.gif"))

# ------------------------------------------------------------
# 4.  Result --------------------------------------------------
# ------------------------------------------------------------
cat("GIF written to", normalizePath("graphcut_evolution.gif"), "\n")
