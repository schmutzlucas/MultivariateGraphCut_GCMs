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


range_var_final <- readRDS('ranges/range_var_summer_ERA5_1950-2023_3v.rds')

# ------------------------------------------------------------------
# A. build permutation that converts 0?359 ? -180?+179 order
lon_file        <- 0:359
lon_adj         <- ifelse(lon_file >= 180, lon_file - 360, lon_file)
lon_order_adj   <- order(lon_adj)                # length 360

# ------------------------------------------------------------------
# B. re-order every variable?s range matrix
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
year_future <<- 1999:2024

# Seasonal window (user-settable): month/day
# Example summer: April 15 to October 14
season_start_md <- c(4, 15)
season_end_md   <- c(10, 14)

# data directory
data_dir <<- 'data/CMIP6_merged_all'


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
cat("-> building PDFs and means ...\n")
t0 <- Sys.time()
tmp <- compute_nd_pdf_multi_seasonal(
  variables, model_names,
  data_dir,
  year_present, year_future,
  lon, lat,
  range_var_gc,
  season_start_md = season_start_md,
  season_end_md   = season_end_md,
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
cat("? workspace saved to", filename, "\n")


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
#  2) build ?all?bins? index list for full Hellinger
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
smooth_cost <- 0.01

GC_result <- tryCatch({
  GraphCutHellinger_nD_lat(
    pdf_models_future = pdf3_models_fut,   # using ?present? PDFs for labeling
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
  cat("??  GraphCut failed at smooth_cost =", smooth_cost, ":\n", e$message, "\n")
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
}
