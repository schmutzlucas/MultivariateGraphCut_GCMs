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

range_var_final <- readRDS('ranges/range_var_final_allModelsPar_1950-2023_90deg_3v.rds')

# Setting global variables
lon <- -180:179
lat <- -90:90

# # Setting global variables
# # 2° grid from your my_grid_2deg.txt:
# lon <- seq(-180, 178, by = 2)  # length = 180
# lat <- seq(-90,   90,  by = 2)  # length =  91


# Temporal ranges
year_present <<- 1950:1975
year_future <<- 2075:2100
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
  tmp <- compute_nd_pdf_optimized_0centered(variables, model_names, data_dir, year_present, year_future,
                                            lon, lat, aperm(abind(range_var_final$ranges, along = 4), c(1, 2, 4, 3)), nbins1d, workers = 4)
})
cat("Time taken for compute_nd_pdf_optimized: ", format_time(time_optimized["elapsed"]), "\n")


# Get the current date and time
current_time <- Sys.time()

# Format the date and time as a string in the format 'yyyymmddhhmm'
formatted_time <- format(current_time, "%Y%m%d%H%M")

# Concatenate the formatted time string with your desired filename
filename <- paste0(formatted_time, "_my_workspace_PME_allModels_beforeOptim_3v.RData")

# Save the workspace using the generated filename
save.image(file = filename, compress = FALSE)


# — subset tmp$present / tmp$future to only the "no-dup" models ——————

# 1. read in both lists
full_models  <- read.table('model_names_pr_tas_psl_perfect_model.txt',              stringsAsFactors = FALSE)[[1]]
short_models <- read.table('model_names_pr_tas_psl_perfect_model_without_duplicate.txt', stringsAsFactors = FALSE)[[1]]

# 2. find their positions in the full list
keep_idx <- match(short_models, full_models)
if (any(is.na(keep_idx))) {
  stop("Some models in the 'without_duplicate' list were not found in the full list:\n",
       paste(short_models[is.na(keep_idx)], collapse = ", "))
}

# 3. subset your 4-D PDF arrays along the model dimension (4th dim)
pdf_present_all <- tmp$present[,,, keep_idx]
pdf_future_all  <- tmp$future[,,, keep_idx]

# 4. update all_models to the shorter list
all_models <- as.list(short_models)


# --- initialize structures ---------------------------------------------------

# list of all models
all_models <- read.table('model_names_pr_tas_psl_perfect_model_without_duplicate.txt')
all_models <- as.list(all_models[['V1']])

# dimension sizes
lon_size <- length(lon)
lat_size <- length(lat)

# total number of joint‐pdf bins
nbins_total <- nbins1d^length(variables)

# for full H‐distance: include every bin
all_bins <- seq_len(nbins_total)
selected_indices_all <- lapply(seq_len(lon_size), function(i)
  lapply(seq_len(lat_size), function(j) all_bins)
)

# smoothness weights to sweep
smooth_costs <- c(0.05, 0.1, 0.6, 0.7, 0.8, 0.9, 1, 1.1, 1.2, 1.3)

# pre-allocate result lists
GC_result_list               <- setNames(vector("list", length(smooth_costs)), as.character(smooth_costs))
GC_hdist_present_list        <- GC_hdist_future_list        <- GC_partial_hdist_future_list <- GC_result_list
MMM_hdist_present_list       <- list()
MMM_hdist_future_list        <- list()
MMM_partial_hdist_future_list<- list()
MMM_present_list             <- list()
MMM_future_list              <- list()

# parallel leave-one-out with future.apply ------------------------------

library(future.apply)
options(future.globals.maxSize = 16 * 1024^3)


# 1) choose as many workers as you can afford in RAM
plan(multisession, workers = 4)

# 2) wrap one reference’s entire workflow in a function
process_ref <- function(m) {
  reference_name <- all_models[[m]]
  other_models   <- all_models[-m]

  # slice PDFs
  pdf_ref_present    <- pdf_present_all[,,,m]
  pdf_models_present <- pdf_present_all[,,,-m]
  pdf_ref_future     <- pdf_future_all[,,,m]
  pdf_models_future  <- pdf_future_all[,,,-m]

  # full H-dist maps
  h_dist_pres <- compute_partial_hdist(pdf_ref_present, pdf_models_present, selected_indices_all)
  h_dist_fut  <- compute_partial_hdist(pdf_ref_future,  pdf_models_future,  selected_indices_all)

  # MMM pdfs + H-dist
  MMM_pres <- apply(pdf_models_present, c(1,2,3), mean)
  MMM_fut  <- apply(pdf_models_future,  c(1,2,3), mean)

  lon_size <- dim(h_dist_pres)[1]
  lat_size <- dim(h_dist_pres)[2]

  MMM_hdist_pres <- MMM_hdist_fut <- matrix(NA, lon_size, lat_size)
  for(i in seq_len(lon_size)) for(j in seq_len(lat_size)) {
    MMM_hdist_pres[i,j] <- sqrt(sum((sqrt(MMM_pres[i,j,]) - sqrt(pdf_ref_present[i,j,]))^2)) / sqrt(2)
    MMM_hdist_fut[i,j]  <- sqrt(sum((sqrt(MMM_fut[i,j, ]) - sqrt(pdf_ref_future[i,j, ]))^2)) / sqrt(2)
  }
  MMM_hdist_pres[is.nan(MMM_hdist_pres)] <- 0
  MMM_hdist_fut[is.nan(MMM_hdist_fut)]   <- 0

  # partial H-dist on outer 10%
  ldr_idx <- lapply(seq_len(lon_size), function(i)
    lapply(seq_len(lat_size), function(j)
      select_hdr_indices(pdf_ref_future[i,j,], tau = 0.10)
    )
  )
  partial_fut <- compute_partial_hdist(pdf_ref_future, pdf_models_future, ldr_idx)

  # GraphCut sweep
  GC_res   <- list()
  GC_hp    <- list()
  GC_hf    <- list()
  GC_hpf   <- list()

  for(lambda in smooth_costs) {
    key <- as.character(lambda)
    gc_out <- tryCatch({
      GraphCutHellinger_nD_lat(
        pdf_models_future = pdf_models_future,
        h_dist            = h_dist_pres,
        weight_data       = 1,
        weight_smooth     = lambda,
        nBins             = nbins_total,
        lat               = lat,
        seed              = 1,
        verbose           = FALSE,
        rebuild           = FALSE
      )
    }, error = function(e) NULL)

    if (!is.null(gc_out)) {
      # unpack
      hgp <- hgf <- hpp <- matrix(NA, lon_size, lat_size)
      for(l in seq_along(other_models)) {
        mask         <- (gc_out$label_attribution == l)
        slice_pres   <- h_dist_pres[  ,  , l]
        slice_fut    <- h_dist_fut[   ,  , l]
        slice_pfut   <- partial_fut[  ,  , l]
        hgp[mask]   <- slice_pres[mask]
        hgf[mask]   <- slice_fut[mask]
        hpp[mask]   <- slice_pfut[mask]
      }
      GC_res[[key]] <- gc_out
      GC_hp[[key]]  <- hgp
      GC_hf[[key]]  <- hgf
      GC_hpf[[key]] <- hpp
    }
  }

  # return a named list of everything for this ref
  list(
    ref                = reference_name,
    MMM_hdist_pres     = MMM_hdist_pres,
    MMM_hdist_fut      = MMM_hdist_fut,
    MMM_partial_fut    = partial_fut,
    MMM_pdf_pres       = MMM_pres,
    MMM_pdf_fut        = MMM_fut,
    GC_results         = GC_res,
    GC_hdist_pres_maps = GC_hp,
    GC_hdist_fut_maps  = GC_hf,
    GC_part_fut_maps   = GC_hpf
  )
}

# 3) launch them in parallel
all_outputs <- future_lapply(seq_along(all_models), process_ref, future.seed = TRUE)

# name them by model
names(all_outputs) <- all_models

# 4) unpack into your *_list objects
for(out in all_outputs) {
  ref <- out$ref

  MMM_hdist_present_list[[ref]]         <- out$MMM_hdist_pres
  MMM_hdist_future_list[[ref]]          <- out$MMM_hdist_fut
  MMM_partial_hdist_future_list[[ref]]  <- out$MMM_partial_fut
  MMM_present_list[[ref]]               <- out$MMM_pdf_pres
  MMM_future_list[[ref]]                <- out$MMM_pdf_fut

  for(lambda in names(out$GC_results)) {
    GC_result_list[[lambda]][[ref]]               <- out$GC_results[[lambda]]
    GC_hdist_present_list[[lambda]][[ref]]        <- out$GC_hdist_pres_maps[[lambda]]
    GC_hdist_future_list[[lambda]][[ref]]         <- out$GC_hdist_fut_maps[[lambda]]
    GC_partial_hdist_future_list[[lambda]][[ref]] <- out$GC_part_fut_maps[[lambda]]
  }
}

plan(sequential)

# done – now all *_list objects filled exactly as in the serial loop.

# Get the current date and time
current_time <- Sys.time()

# Format the date and time as a string in the format 'yyyymmddhhmm'
formatted_time <- format(current_time, "%Y%m%d%H%M")

# Concatenate the formatted time string with your desired filename
filename <- paste0(formatted_time, "_my_workspace_PME_allModels_results.RData")

# Save the workspace using the generated filename
save.image(file = filename, compress = FALSE)


# the two reference names you want to drop
drop_models <- c("CMCC-CM2-SR5", "CMCC-ESM2")

# tmp1 will now be a copy of GC_hdist_future_list
# but with those two references removed from every cost
tmp1 <- lapply(GC_hdist_future_list, function(cost_list) {
  cost_list[ ! names(cost_list) %in% drop_models ]
})

tmp2 <- MMM_hdist_future_list[ ! names(MMM_hdist_future_list) %in% drop_models ]



# Extract the label attribution for the current smooth cost
GC_labels <- GC_result_list$`0.1`$`MIROC-ES2L`$label_attribution

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

#Violin partial H by ref
{
  library(ggplot2)
  library(dplyr)

  # Crée le dossier 'figure' s'il n'existe pas
  if (!dir.exists("figure")) dir.create("figure")

  # Liste des références disponibles
  references_done <- names(MMM_partial_hdist_future_list)

  for (ref_name in references_done) {

    # MMM
    mmm_vals <- as.vector(MMM_partial_hdist_future_list[[ref_name]])
    df_ref <- data.frame(
      Hellinger = mmm_vals,
      Method = "MMM",
      SmoothCost = "MMM",
      Reference = ref_name
    )

    # Tous les GraphCut
    for (cost in c("0.05", "0.1", "0.6", "1", "2")) {
      if (!is.null(GC_partial_hdist_future_list[[cost]][[ref_name]])) {
        gc_vals <- as.vector(GC_partial_hdist_future_list[[cost]][[ref_name]])
        df_ref <- rbind(df_ref, data.frame(
          Hellinger = gc_vals,
          Method = "GraphCut",
          SmoothCost = cost,
          Reference = ref_name
        ))
      }
    }

    df_ref <- na.omit(df_ref)
    df_ref$MethodLabel <- ifelse(df_ref$Method == "MMM", "MMM", paste0("GC-", df_ref$SmoothCost))
    df_ref$MethodLabel <- factor(df_ref$MethodLabel, levels = c("GC-0.05", "GC-0.1", "GC-0.6", "GC-1", "GC-2", "MMM"))

    # Plot
    p_ref <- ggplot(df_ref, aes(x = MethodLabel, y = Hellinger, fill = MethodLabel)) +
      geom_violin(scale = "area", adjust = 1.2, width = 0.7, alpha = 0.85) +
      scale_fill_manual(values = c(
        "GC-0.05" = "#1b9e77",
        "GC-0.1"  = "#d95f02",
        "GC-0.6"  = "#7570b3",
        "GC-1"    = "#e7298a",
        "GC-2"    = "#66a61e",
        "MMM"     = "#e6ab02"
      )) +
      labs(
        title = paste("Partial Hellinger Distance (Future) - Reference:", ref_name),
        x = "Method",
        y = "Partial Hellinger Distance",
        fill = "Method"
      ) +
      theme_minimal(base_size = 14) +
      theme(
        axis.text.x = element_text(size = 10),
        plot.title = element_text(size = 16, face = "bold"),
        legend.position = "none"
      )

    # Nettoyer le nom du fichier
    clean_name <- gsub("[^A-Za-z0-9]", "", ref_name)
    file_base <- paste0("figure/Violin_Hdist_BC_22model_", clean_name)

    # Sauvegarde
    ggsave(paste0(file_base, "_partial.pdf"), plot = p_ref, width = 20, height = 15, units = "cm", dpi = 300)
    ggsave(paste0(file_base, "_partial.png"), plot = p_ref, width = 20, height = 15, units = "cm", dpi = 300)

    cat("Saved partial Hellinger violin plot for reference:", ref_name, "\n")
  }

}


# Violins aggrégés par méthode ----------------------------------------------
{

  library(ggplot2)
  library(dplyr)
  library(RColorBrewer)

  if (!dir.exists("figure")) dir.create("figure")

  # 1) Gather the smooth costs that actually exist
  smooth_costs <- names(GC_hdist_future_list)
  smooth_costs <- as.character(sort(as.numeric(smooth_costs)))

  # 2) Build the big data.frame across MMM + all GC-cost methods
  df_all <- do.call(rbind, lapply(names(MMM_hdist_future_list), function(ref_name) {
    # MMM
    mmm_vals <- as.vector(MMM_hdist_future_list[[ref_name]])
    df <- data.frame(
      Hellinger  = mmm_vals,
      Method     = "MMM",
      SmoothCost = "MMM",
      Reference  = ref_name,
      stringsAsFactors = FALSE
    )
    # each GraphCut cost
    for (cost in smooth_costs) {
      gc_map <- GC_hdist_future_list[[cost]][[ref_name]]
      if (!is.null(gc_map)) {
        df <- rbind(df, data.frame(
          Hellinger  = as.vector(gc_map),
          Method     = "GraphCut",
          SmoothCost = cost,
          Reference  = ref_name,
          stringsAsFactors = FALSE
        ))
      }
    }
    return(df)
  }))

  df_all <- na.omit(df_all)

  # 3) Build a single label factor that orders first the GC methods, then MMM
  method_labels <- c(paste0("GC-", smooth_costs), "MMM")
  df_all$MethodLabel <- with(df_all,
                             ifelse(Method=="MMM", "MMM", paste0("GC-", SmoothCost))
  )
  df_all$MethodLabel <- factor(df_all$MethodLabel, levels = method_labels)

  # 4) Plot
  # generate a GC palette of the correct length
  library(RColorBrewer)
  n_gc <- length(smooth_costs)
  gc_colors <- colorRampPalette(brewer.pal(8, "Set2"))(n_gc)
  names(gc_colors) <- paste0("GC-", smooth_costs)

  all_colors <- c(gc_colors, MMM = "#e6ab02")

  p_agg <- ggplot(df_all, aes(x = MethodLabel, y = Hellinger, fill = MethodLabel)) +
    geom_violin(scale = "area", trim = TRUE, adjust = 1.5, alpha = 0.85, width = 0.7) +
    stat_summary(fun = mean, geom = "point", shape = 20, size = 2.2,
                 color = "black", position = position_dodge(width = 0.7)) +
    stat_summary(fun = median, geom = "crossbar", width = 0.4,
                 color = "red", fatten = 1, position = position_dodge(width = 0.7)) +
    scale_fill_manual(values = all_colors) +
    labs(
      title    = "Distribution of Hellinger Distance (Future)",
      subtitle = "Aggregated across all reference models",
      x        = "Method",
      y        = "Hellinger Distance",
      fill     = "Method"
    ) +
    theme_minimal(base_size = 14) +
    theme(
      axis.text.x     = element_text(size = 11, angle = 45, hjust = 1),
      plot.title      = element_text(size = 16, face = "bold"),
      legend.position = "none"
    )


  file_base <- "figure/Aggregated_Hdist_BC_22models"
  ggsave(paste0(file_base, ".pdf"), p_agg, width = 20, height = 15, units = "cm", dpi = 300)
  ggsave(paste0(file_base, ".png"), p_agg, width = 20, height = 15, units = "cm", dpi = 300)

  cat("Saved aggregated Hellinger distance plot to", file_base, "\n")

}


# Violins aggrégés par méthode sans jumeaux
{

  library(ggplot2)
  library(dplyr)
  library(RColorBrewer)

  if (!dir.exists("figure")) dir.create("figure")

  # 1) Gather the smooth costs that actually exist
  smooth_costs <- names(GC_hdist_future_list)
  smooth_costs <- as.character(sort(as.numeric(smooth_costs)))

  # 2) Build the big data.frame across MMM + all GC-cost methods
  df_all <- do.call(rbind, lapply(names(tmp2), function(ref_name) {
    # MMM
    mmm_vals <- as.vector(tmp2[[ref_name]])
    df <- data.frame(
      Hellinger  = mmm_vals,
      Method     = "MMM",
      SmoothCost = "MMM",
      Reference  = ref_name,
      stringsAsFactors = FALSE
    )
    # each GraphCut cost
    for (cost in smooth_costs) {
      gc_map <- tmp1[[cost]][[ref_name]]
      if (!is.null(gc_map)) {
        df <- rbind(df, data.frame(
          Hellinger  = as.vector(gc_map),
          Method     = "GraphCut",
          SmoothCost = cost,
          Reference  = ref_name,
          stringsAsFactors = FALSE
        ))
      }
    }
    return(df)
  }))

  df_all <- na.omit(df_all)

  # 3) Build a single label factor that orders first the GC methods, then MMM
  method_labels <- c(paste0("GC-", smooth_costs), "MMM")
  df_all$MethodLabel <- with(df_all,
                             ifelse(Method=="MMM", "MMM", paste0("GC-", SmoothCost))
  )
  df_all$MethodLabel <- factor(df_all$MethodLabel, levels = method_labels)

  # 4) Plot
  # generate a GC palette of the correct length
  library(RColorBrewer)
  n_gc <- length(smooth_costs)
  gc_colors <- colorRampPalette(brewer.pal(8, "Set2"))(n_gc)
  names(gc_colors) <- paste0("GC-", smooth_costs)

  all_colors <- c(gc_colors, MMM = "#e6ab02")

  p_agg <- ggplot(df_all, aes(x = MethodLabel, y = Hellinger, fill = MethodLabel)) +
    geom_violin(scale = "area", trim = TRUE, adjust = 1.5, alpha = 0.85, width = 0.7) +
    stat_summary(fun = mean, geom = "point", shape = 20, size = 2.2,
                 color = "black", position = position_dodge(width = 0.7)) +
    stat_summary(fun = median, geom = "crossbar", width = 0.4,
                 color = "red", fatten = 1, position = position_dodge(width = 0.7)) +
    scale_fill_manual(values = all_colors) +
    labs(
      title    = "Distribution of Hellinger Distance (Future)",
      subtitle = "Aggregated across all reference models",
      x        = "Method",
      y        = "Hellinger Distance",
      fill     = "Method"
    ) +
    theme_minimal(base_size = 14) +
    theme(
      axis.text.x     = element_text(size = 11, angle = 45, hjust = 1),
      plot.title      = element_text(size = 16, face = "bold"),
      legend.position = "none"
    )


  file_base <- "figure/Aggregated_Hdist_BC_22models"
  ggsave(paste0(file_base, ".pdf"), p_agg, width = 20, height = 15, units = "cm", dpi = 300)
  ggsave(paste0(file_base, ".png"), p_agg, width = 20, height = 15, units = "cm", dpi = 300)

  cat("Saved aggregated Hellinger distance plot to", file_base, "\n")

}

{
  # --------------------------------------------------------------------------
  # Build Global / Land / Ocean violins by masking, without spatial joins
  # --------------------------------------------------------------------------

  library(ggplot2)
  library(dplyr)
  library(RColorBrewer)
  library(maps)

  # ensure figure directory
  if (!dir.exists("figure")) dir.create("figure")

  # 0) correct mask
  lon_rep <- rep(lon, times = length(lat))
  lat_rep <- rep(lat, each   = length(lon))
  reg     <- map.where("world", lon_rep, lat_rep)
  ocean_mask <- matrix(is.na(reg), nrow=length(lon), ncol=length(lat))
  land_mask  <- !ocean_mask

  # 1) build data.frame by slicing each mat:
  df_list <- list()
  for(ref in names(tmp2)) {
    # MMM global
    mmm <- tmp2[[ref]]
    df_list[[paste(ref,"MMM","Global")]] <- data.frame(
      Hellinger=as.vector(mmm),
      MethodLabel="MMM",
      Region="Global",
      Reference=ref
    )
    # MMM land
    df_list[[paste(ref,"MMM","Land")]] <- data.frame(
      Hellinger=mmm[land_mask],
      MethodLabel="MMM",
      Region="Land",
      Reference=ref
    )
    # MMM ocean
    df_list[[paste(ref,"MMM","Ocean")]] <- data.frame(
      Hellinger=mmm[ocean_mask],
      MethodLabel="MMM",
      Region="Ocean",
      Reference=ref
    )

    # GC for each λ
    for(cost in smooth_costs) {
      label <- paste0("GC-",cost)
      gc  <- tmp1[[cost]][[ref]]  # a lon×lat matrix
      if(is.null(gc)) next

      df_list[[paste(ref,label,"Global")]] <- data.frame(
        Hellinger=as.vector(gc),
        MethodLabel=label,
        Region="Global",
        Reference=ref
      )
      df_list[[paste(ref,label,"Land")]] <- data.frame(
        Hellinger=gc[land_mask],
        MethodLabel=label,
        Region="Land",
        Reference=ref
      )
      df_list[[paste(ref,label,"Ocean")]] <- data.frame(
        Hellinger=gc[ocean_mask],
        MethodLabel=label,
        Region="Ocean",
        Reference=ref
      )
    }
  }

  df_all <- bind_rows(df_list)

  # now make sure factors are right:
  df_all$MethodLabel <- factor(df_all$MethodLabel,
                               levels = c(paste0("GC-",smooth_costs),"MMM"))
  df_all$Region      <- factor(df_all$Region,
                               levels = c("Global","Land","Ocean"))

  # 2) prepare colours
  gc_cols <- colorRampPalette(brewer.pal(8, "Set2"))(length(smooth_costs))
  names(gc_cols) <- paste0("GC-", smooth_costs)
  all_cols <- c(gc_cols, MMM = "#e6ab02")

  # 3) plot
  p_region <- ggplot(df_all, aes(x = MethodLabel, y = Hellinger, fill = MethodLabel)) +
    geom_violin(scale = "area", trim = TRUE, adjust = 1.5, alpha = 0.8, width = 0.7) +
    stat_summary(fun = mean, geom = "point", shape = 20, size = 2, color = "black",
                 position = position_dodge(width = 0.7)) +
    stat_summary(fun = median, geom = "crossbar", width = 0.4, color = "red", fatten = 1,
                 position = position_dodge(width = 0.7)) +
    scale_fill_manual(values = all_cols) +
    facet_grid(. ~ Region) +
    labs(
      title    = "Hellinger Distance by Method & Region",
      subtitle = "Global vs. Land vs. Ocean",
      x        = "Method",
      y        = "Hellinger Distance",
      fill     = "Method"
    ) +
    theme_minimal(base_size = 14) +
    theme(
      axis.text.x     = element_text(angle = 45, hjust = 1),
      plot.title      = element_text(size = 16, face = "bold"),
      legend.position = "none"
    )

  # 4) save
  ggsave("figure/Aggregated_Hdist_By_Region.pdf", p_region,
         width = 30, height = 10, units = "cm", dpi = 300)
  ggsave("figure/Aggregated_Hdist_By_Region.png", p_region,
         width = 30, height = 10, units = "cm", dpi = 300)

  cat("Saved region-separated H-distance violins to figure/Aggregated_Hdist_By_Region.*\n")


}



# violin panel with 22 ref : H dist
{
  library(ggplot2)
  library(dplyr)
  library(patchwork)

  # Liste des références disponibles
  references_done <- names(MMM_hdist_future_list)

  # Créer tous les plots dans une liste
  plot_list <- list()

  for (ref_name in references_done) {

    mmm_vals <- as.vector(MMM_hdist_future_list[[ref_name]])
    df_ref <- data.frame(
      Hellinger = mmm_vals,
      Method = "MMM",
      SmoothCost = "MMM",
      Reference = ref_name
    )

    for (cost in c("0.05", "0.1", "0.6", "1", "2")) {
      if (!is.null(GC_hdist_future_list[[cost]][[ref_name]])) {
        gc_vals <- as.vector(GC_hdist_future_list[[cost]][[ref_name]])
        df_ref <- rbind(df_ref, data.frame(
          Hellinger = gc_vals,
          Method = "GraphCut",
          SmoothCost = cost,
          Reference = ref_name
        ))
      }
    }

    df_ref <- na.omit(df_ref)
    df_ref$MethodLabel <- ifelse(df_ref$Method == "MMM", "MMM", paste0("GC-", df_ref$SmoothCost))
    df_ref$MethodLabel <- factor(df_ref$MethodLabel, levels = c("GC-0.05", "GC-0.1", "GC-0.6", "GC-1", "GC-2", "MMM"))

    # Create the individual plot
    p <- ggplot(df_ref, aes(x = MethodLabel, y = Hellinger, fill = MethodLabel)) +
      geom_violin(scale = "area", adjust = 1.2, width = 0.7, alpha = 0.85) +
      stat_summary(fun = mean, geom = "point", shape = 20, size = 2.2, color = "black", position = position_dodge(width = 0.7)) +
      stat_summary(fun = median, geom = "crossbar", width = 0.4, color = "red", fatten = 1, position = position_dodge(width = 0.7)) +
      scale_fill_manual(values = c(
        "GC-0.05" = "#1b9e77",
        "GC-0.1"  = "#d95f02",
        "GC-0.6"  = "#7570b3",
        "GC-1"    = "#e7298a",
        "GC-2"    = "#66a61e",
        "MMM"     = "#e6ab02"
      )) +
      labs(
        title = ref_name,
        x = NULL,
        y = NULL
      ) +
      theme_minimal(base_size = 10) +
      theme(
        axis.text.x = element_text(size = 7, angle = 90, hjust = 1, vjust = 0.5),
        axis.text.y = element_text(size = 7),
        plot.title = element_text(size = 10, face = "bold"),
        legend.position = "none"
      )

    plot_list[[ref_name]] <- p
  }

  # Combine all plots with patchwork
  big_plot <- wrap_plots(plot_list, ncol = 8) +
    plot_annotation(
      title = "Hellinger Distance (Future) by Reference and Method",
      theme = theme(
        plot.title = element_text(size = 18, face = "bold", hjust = 0.5)
      )
    )

  # Save in high-resolution (4K scale, roughly)
  ggsave("figure/AllReferences_Hellinger_Grid_22models.png", plot = big_plot,
         width = 1920/96, height = 1080/96, dpi = 300, units = "in")  # 4K: 3840×2160 pixels
  ggsave("figure/AllReferences_Hellinger_Grid_22models.pdf", plot = big_plot,
         width = 1920/96, height = 1080/96, dpi = 300, units = "in")

  cat("✅ Multi-panel violin plot saved as 4K image and PDF.\n")

}



# Summary table
{
  # Résumé des statistiques Hellinger par méthode et par référence
  summary_stats <- df_all %>%
    group_by(Reference, MethodLabel) %>%
    summarise(
      Mean = mean(Hellinger, na.rm = TRUE),
      Median = median(Hellinger, na.rm = TRUE),
      .groups = "drop"
    )

  # Ajouter la ligne globale (toutes références confondues)
  global_summary <- df_all %>%
    group_by(MethodLabel) %>%
    summarise(
      Mean = mean(Hellinger, na.rm = TRUE),
      Median = median(Hellinger, na.rm = TRUE)
    ) %>%
    mutate(Reference = "ALL")

  # Combiner
  summary_stats <- bind_rows(summary_stats, global_summary)

  # Afficher le tableau
  print(summary_stats)

  # Optionnel : export CSV
  # write.csv(summary_stats, "summary_Hellinger_stats.csv", row.names = FALSE)

}

# Summary table
{
  # Summary stats from tmp2 (MMM) and tmp1 (GraphCut) -----------------------

  library(dplyr)

  # 1) per-reference, per-method summary
  stats_list <- list()

  # MMM entries
  for(ref in names(tmp2)) {
    vals <- as.vector(tmp2[[ref]])
    stats_list[[length(stats_list)+1]] <- data.frame(
      Reference   = ref,
      MethodLabel = "MMM",
      Mean        = mean(vals,   na.rm = TRUE),
      Median      = median(vals, na.rm = TRUE),
      stringsAsFactors = FALSE
    )
  }

  # GraphCut entries
  for(cost in smooth_costs) {
    label <- paste0("GC-", cost)
    for(ref in names(tmp1[[cost]])) {
      mat   <- tmp1[[cost]][[ref]]
      if (is.null(mat)) next
      vals  <- as.vector(mat)
      stats_list[[length(stats_list)+1]] <- data.frame(
        Reference   = ref,
        MethodLabel = label,
        Mean        = mean(vals,   na.rm = TRUE),
        Median      = median(vals, na.rm = TRUE),
        stringsAsFactors = FALSE
      )
    }
  }

  # bind per-ref results
  summary_stats <- bind_rows(stats_list)

  # 2) global summary (across *all* references) for each method
  global_list <- list()

  # MMM global
  all_mmm <- unlist(lapply(tmp2, as.vector))
  global_list[[length(global_list)+1]] <- data.frame(
    Reference   = "ALL",
    MethodLabel = "MMM",
    Mean        = mean(all_mmm,   na.rm = TRUE),
    Median      = median(all_mmm, na.rm = TRUE),
    stringsAsFactors = FALSE
  )

  # GraphCut global
  for(cost in smooth_costs) {
    label <- paste0("GC-", cost)
    mats  <- tmp1[[cost]]
    all_gc <- unlist(lapply(mats, as.vector))
    global_list[[length(global_list)+1]] <- data.frame(
      Reference   = "ALL",
      MethodLabel = label,
      Mean        = mean(all_gc,   na.rm = TRUE),
      Median      = median(all_gc, na.rm = TRUE),
      stringsAsFactors = FALSE
    )
  }

  global_summary <- bind_rows(global_list)

  # 3) combine
  summary_stats <- bind_rows(summary_stats, global_summary)

  # 4) display (and optionally save)
  print(summary_stats)
  # write.csv(summary_stats, "summary_Hellinger_stats.csv", row.names = FALSE)

}


# Summary table by median
{
  library(ggplot2)
  library(dplyr)
  library(tidyr)

  # Identifier la méthode gagnante par référence (plus petite médiane)
  summary_stats_highlight <- summary_stats %>%
    group_by(Reference) %>%
    mutate(
      IsBest = Median == min(Median, na.rm = TRUE),
      FontFace = ifelse(IsBest, "bold", "plain")
    ) %>%
    ungroup()

  # Heatmap avec texte noir et valeur gagnante en gras
  p_median <- ggplot(summary_stats_highlight, aes(x = MethodLabel, y = Reference, fill = Median)) +
    geom_tile(color = "white") +
    geom_text(
      aes(label = sprintf("%.3f", Median), fontface = FontFace),
      size = 3, color = "black"
    ) +
    scale_fill_viridis_c(option = "C", name = "Median\nHellinger") +
    labs(
      title = "Median Hellinger Distance by Method and Reference",
      subtitle = "Best per row is shown in bold",
      x = "Method",
      y = "Reference"
    ) +
    theme_minimal(base_size = 13) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      plot.title = element_text(face = "bold")
    )

  # Affichage + sauvegarde
  print(p_median)
  ggsave("figure/PME_NOBC_Summary_Hellinger_Heatmap_Median_Highlighted.pdf", plot = p_median,
         width = 20, height = 15, units = "cm", dpi = 300)
  ggsave("figure/PME_NOBC_Summary_Hellinger_Heatmap_Median_Highlighted.png", plot = p_median,
         width = 20, height = 15, units = "cm", dpi = 300)
}

# Summary table by mean
{
  library(ggplot2)
  library(dplyr)
  library(tidyr)

  # Reordonner les références avec "ALL" à la fin
  summary_stats <- summary_stats %>%
    mutate(
      Reference = factor(Reference, levels = c(setdiff(unique(Reference), "ALL"), "ALL")),
      MethodLabel = factor(MethodLabel, levels = c("GC-0.05", "GC-0.1", "GC-0.6", "GC-1", "GC-2", "MMM"))
    )

  # Identifier la meilleure méthode par référence (plus petite moyenne)
  summary_stats_highlight <- summary_stats %>%
    group_by(Reference) %>%
    mutate(
      IsBest = Mean == min(Mean, na.rm = TRUE),
      FontFace = ifelse(IsBest, "bold", "plain")
    ) %>%
    ungroup()

  # Heatmap avec mise en évidence
  p_mean <- ggplot(summary_stats_highlight, aes(x = MethodLabel, y = Reference, fill = Mean)) +
    geom_tile(color = "white") +
    geom_text(aes(label = sprintf("%.3f", Mean), fontface = FontFace), size = 3, color = "black") +
    scale_fill_viridis_c(option = "C", name = "Mean\nHellinger") +
    labs(
      title = "Mean Hellinger Distance by Method and Reference",
      subtitle = "Best per row is shown in bold",
      x = "Method",
      y = "Reference"
    ) +
    theme_minimal(base_size = 13) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      plot.title = element_text(face = "bold")
    )

  # Affichage
  print(p_mean)

  # Sauvegarde
  ggsave("figure/Summary_Hellinger_Heatmap_Mean_Highlighted.pdf", plot = p_mean, width = 20, height = 15, units = "cm", dpi = 300)
  ggsave("figure/Summary_Hellinger_Heatmap_Mean_Highlighted.png", plot = p_mean, width = 20, height = 15, units = "cm", dpi = 300)
}


# Aggregated partial Hellinger for smooth 1 only
{
  library(ggplot2)
  library(dplyr)

  if (!dir.exists("figure")) dir.create("figure")

  # Aggregate only for GC-1 and MMM
  df_all <- do.call(rbind, lapply(names(MMM_partial_hdist_future_list), function(ref_name) {
    mmm_vals <- as.vector(MMM_partial_hdist_future_list[[ref_name]])
    df <- data.frame(
      Hellinger = mmm_vals,
      Method = "MMM",
      SmoothCost = "MMM",
      Reference = ref_name
    )

    if (!is.null(GC_partial_hdist_future_list[["1"]][[ref_name]])) {
      gc_vals <- as.vector(GC_partial_hdist_future_list[["1"]][[ref_name]])
      df <- rbind(df, data.frame(
        Hellinger = gc_vals,
        Method = "GraphCut",
        SmoothCost = "1",
        Reference = ref_name
      ))
    }

    return(df)
  }))

  df_all <- na.omit(df_all)
  df_all$MethodLabel <- ifelse(df_all$Method == "MMM", "MMM", "GC-1")
  df_all$MethodLabel <- factor(df_all$MethodLabel, levels = c("GC-1", "MMM"))

  p_agg <- ggplot(df_all, aes(x = MethodLabel, y = Hellinger, fill = MethodLabel)) +
    geom_violin(scale = "area", trim = TRUE, adjust = 1.5, alpha = 0.85, width = 0.7) +
    stat_summary(fun = mean, geom = "point", shape = 20, size = 2.2, color = "black", position = position_dodge(width = 0.7)) +
    stat_summary(fun = median, geom = "crossbar", width = 0.4, color = "red", fatten = 1, position = position_dodge(width = 0.7)) +
    scale_fill_manual(values = c(
      "GC-1" = "#d95f02",
      "MMM"    = "#e6ab02"
    )) +
    labs(
      title = "Partial Hellinger Distance (Future)",
      subtitle = "Aggregated across all references - Smooth = 1",
      x = "Method",
      y = "Partial Hellinger Distance",
      fill = "Method"
    ) +
    theme_minimal(base_size = 14) +
    theme(
      axis.text.x = element_text(size = 11),
      plot.title = element_text(size = 16, face = "bold"),
      legend.position = "none"
    )

  file_base <- "figure/Aggregated_HdistPartial_BC_22models_Smooth1_noBC"
  ggsave(paste0(file_base, ".pdf"), plot = p_agg, width = 20, height = 15, units = "cm", dpi = 300)
  ggsave(paste0(file_base, ".png"), plot = p_agg, width = 20, height = 15, units = "cm", dpi = 300)

  cat("✅ Aggregated plot for Smooth = 1 saved at", file_base, "\n")
}


{
  library(ggplot2)
  library(dplyr)

  # Filtered references
  references_done <- names(MMM_hdist_future_list)

  # Initialize dataframe to collect everything
  df_all <- data.frame()

  for (ref_name in references_done) {

    # Get MMM values
    mmm_vals <- as.vector(MMM_hdist_future_list[[ref_name]])
    df_tmp <- data.frame(
      Hellinger = mmm_vals,
      Method = "MMM",
      SmoothCost = "MMM",
      Reference = ref_name
    )

    # Get GC-1 values
    if (!is.null(GC_hdist_future_list[["1"]][[ref_name]])) {
      gc_vals <- as.vector(GC_hdist_future_list[["1"]][[ref_name]])
      df_tmp <- rbind(df_tmp, data.frame(
        Hellinger = gc_vals,
        Method = "GraphCut",
        SmoothCost = "1",
        Reference = ref_name
      ))
    }

    df_all <- rbind(df_all, df_tmp)
  }

  df_all <- na.omit(df_all)
  df_all$MethodLabel <- ifelse(df_all$Method == "MMM", "MMM", "GC-1")
  df_all$MethodLabel <- factor(df_all$MethodLabel, levels = c("GC-1", "MMM"))

  # One violin plot per method, faceted by reference
  p <- ggplot(df_all, aes(x = MethodLabel, y = Hellinger, fill = MethodLabel)) +
    geom_violin(scale = "area", adjust = 1.2, width = 0.7, alpha = 0.85) +
    stat_summary(fun = mean, geom = "point", shape = 20, size = 2.2, color = "black", position = position_dodge(width = 0.7)) +
    stat_summary(fun = median, geom = "crossbar", width = 0.4, color = "red", fatten = 1, position = position_dodge(width = 0.7)) +
    scale_fill_manual(values = c(
      "GC-1" = "#d95f02",
      "MMM"    = "#e6ab02"
    )) +
    labs(
      title = "Hellinger Distance (Future) by Reference – Smooth = 1",
      x = "Method",
      y = "Hellinger Distance"
    ) +
    theme_minimal(base_size = 12) +
    theme(
      axis.text.x = element_text(size = 10, angle = 0, hjust = 0.5),
      axis.text.y = element_text(size = 10),
      plot.title = element_text(size = 14, face = "bold")
    ) +
    facet_wrap(~Reference, ncol = 9)

  # Save plot
  ggsave("figure/Hellinger_Comparison_Smooth_1_AllRefs_nobc.png", plot = p,
         width = 20, height = 10, units = "in", dpi = 300)
  # Save plot
  ggsave("figure/Hellinger_Comparison_Smooth_1_AllRefs_nobc.pdf", plot = p,
         width = 20, height = 10, units = "in", dpi = 300)

  cat("✅ Violin plot for Smooth = 0.1 and MMM saved.\n")

}