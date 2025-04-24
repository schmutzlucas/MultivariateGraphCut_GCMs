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


pdf_present_all <- tmp$present
pdf_future_all  <- tmp$future


# --- initialize structures ---------------------------------------------------

# list of all models
all_models <- model_names

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
smooth_costs <- c(0.05, 0.1, 0.6, 1, 2)

# pre-allocate result lists
GC_result_list               <- setNames(vector("list", length(smooth_costs)), as.character(smooth_costs))
GC_hdist_present_list        <- GC_hdist_future_list        <- GC_partial_hdist_future_list <- GC_result_list
MMM_hdist_present_list       <- list()
MMM_hdist_future_list        <- list()
MMM_partial_hdist_future_list<- list()
MMM_present_list             <- list()
MMM_future_list              <- list()

# --- main leave-one-out loop ------------------------------------------------
for (m in seq_along(all_models)) {
  reference_name <- all_models[[m]]
  other_models   <- all_models[-m]

  cat("=== Processing reference:", reference_name, "===\n")

  # slice out ref vs. others
  pdf_ref_present    <- pdf_present_all[,,,m]
  pdf_models_present <- pdf_present_all[,,,-m]
  pdf_ref_future     <- pdf_future_all[,,,m]
  pdf_models_future  <- pdf_future_all[,,,-m]

  # compute full H-distance maps
  h_dist_present <- compute_partial_hdist(pdf_ref_present, pdf_models_present, selected_indices_all)
  h_dist_future  <- compute_partial_hdist(pdf_ref_future,  pdf_models_future,  selected_indices_all)

  # compute MMM PDFs (cellwise mean over the other_models)
  MMM_present <- apply(pdf_models_present,    c(1,2,3), mean)
  MMM_future  <- apply(pdf_models_future,     c(1,2,3), mean)

  # compute H-distance of MMM vs. ref
  MMM_hdist_present <- matrix(NA, lon_size, lat_size)
  MMM_hdist_future  <- matrix(NA, lon_size, lat_size)
  for (i in seq_len(lon_size)) {
    for (j in seq_len(lat_size)) {
      MMM_hdist_present[i,j] <- sqrt(sum((sqrt(MMM_present[i,j,]) - sqrt(pdf_ref_present[i,j,]))^2)) / sqrt(2)
      MMM_hdist_future[i,j]  <- sqrt(sum((sqrt(MMM_future[i,j, ]) - sqrt(pdf_ref_future[i,j, ]))^2)) / sqrt(2)
    }
  }
  MMM_hdist_present[is.nan(MMM_hdist_present)] <- 0
  MMM_hdist_future[is.nan(MMM_hdist_future)]   <- 0

  # select “low-density” bins for partial HDist (tau = 0.10)
  ldr_indices <- lapply(seq_len(lon_size), function(i)
    lapply(seq_len(lat_size), function(j)
      select_hdr_indices(pdf_ref_future[i,j,], tau = 0.10)
    )
  )
  partial_hdist_future <- compute_partial_hdist(pdf_ref_future, pdf_models_future, ldr_indices)

  # store MMM results
  MMM_hdist_present_list[[reference_name]]        <- MMM_hdist_present
  MMM_hdist_future_list[[reference_name]]         <- MMM_hdist_future
  MMM_partial_hdist_future_list[[reference_name]] <- partial_hdist_future
  MMM_present_list[[reference_name]]              <- MMM_present
  MMM_future_list[[reference_name]]               <- MMM_future

  # --- GraphCut sweep over smoothness weights -------------------------------
  for (λ in smooth_costs) {
    λ_key <- as.character(λ)
    result <- tryCatch({
      GraphCutHellinger_nD_lat(
        pdf_models_future = pdf_models_future,
        h_dist            = h_dist_present,
        weight_data       = 1,
        weight_smooth     = λ,
        nBins             = nbins_total,
        lat               = lat,
        seed              = 1,
        verbose           = TRUE,
        rebuild           = TRUE
      )
    }, error = function(e) {
      cat("  [Error] λ =", λ, ":", e$message, "\n")
      NULL
    })

    if (!is.null(result)) {
      # unpack label‐map into H-distance maps
      hdist_gc_pres <- matrix(NA, lon_size, lat_size)
      hdist_gc_fut  <- matrix(NA, lon_size, lat_size)
      partial_gc_fut<- matrix(NA, lon_size, lat_size)

      for (lab in seq_along(other_models)) {
        mask <- (result$label_attribution == lab)
        hdist_gc_pres[mask]   <- h_dist_present[mask, lab]
        hdist_gc_fut[mask]    <- h_dist_future[mask, lab]
        partial_gc_fut[mask]  <- partial_hdist_future[mask, lab]
      }

      # store per-λ, per-ref
      GC_result_list[[λ_key]][[reference_name]]                <- result
      GC_hdist_present_list[[λ_key]][[reference_name]]         <- hdist_gc_pres
      GC_hdist_future_list[[λ_key]][[reference_name]]          <- hdist_gc_fut
      GC_partial_hdist_future_list[[λ_key]][[reference_name]]  <- partial_gc_fut
    }

    gc()  # free memory
  }

  # --- checkpoint after each reference -------------------------------------
  cp_name <- gsub("[^A-Za-z0-9]", "", reference_name)
  if (!dir.exists("checkpoints_2deg")) dir.create("checkpoints_2deg")
  saveRDS(
    list(
      GC_result_list,
      GC_hdist_present_list, GC_hdist_future_list, GC_partial_hdist_future_list,
      MMM_hdist_present_list, MMM_hdist_future_list, MMM_partial_hdist_future_list,
      MMM_present_list, MMM_future_list
    ),
    file = sprintf("checkpoints_2deg/PME_%s_%s.rds",
                   cp_name, format(Sys.time(), "%Y%m%d_%H%M%S")),
    compress = FALSE
  )
  cat("  → checkpoint saved for", reference_name, "\n\n")
}

cat("All references processed. Ready for aggregation or plotting.\n")
