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

range_var_final <- readRDS('ranges/range_var_final_GreenwichCentered_1950-2023_90deg_3v.rds')

# Setting global variables
lon <- 0:10
lat <- 0:10
lon_size <- length(lon)
lat_size <- length(lat)
# Temporal ranges
year_present <<- 1950:1975
year_future <<- 1998:2023
# data directory
data_dir <<- 'data/CMIP6_merged_all/'

# List of the variable used
variables <- c('pr', 'tas', 'psl')

## 1.  bin-resolution choices
nbins3d <- 8     # 3-D joint
nbins2d <- 16    # 2-D pairs
nbins1d <- 32    # 1-D marginals

## 2.  model list
model_names <- scan("model_names_pr_tas_psl.txt", what = "", quiet = TRUE)

## 3.  number of parallel workers
workers <- 32   # adapt to your machine

## 4.  call the multi-resolution histogram builder
cat("→ building PDFs and means …\n")
t0 <- Sys.time()
tmp <- compute_nd_pdf_multi(
  variables, model_names,
  data_dir,
  year_present, year_future,
  lon, lat,
  aperm(abind(range_var_final$ranges, along = 4), c(1, 2, 4, 3)),
  nbins3d = nbins3d,
  nbins2d = nbins2d,
  nbins1d = nbins1d,
  workers  = workers)
cat("   done in", round(difftime(Sys.time(), t0, units = "mins"), 2), "minutes\n")

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

nmodels2 <- length(model_names)  # after dropping ref

# --------------------------------------------------------------------
#  2) build “all‐bins” index list for full Hellinger
# --------------------------------------------------------------------
nbins_total3d <- nbins3d^3
selected_indices_all <- lapply(seq_len(nlon), function(i) {
  lapply(seq_len(nlat), function(j) {
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
  GraphCutHellinger_nD(
    pdf_models_future = pdf3_models_pres,   # using “present” PDFs for labeling
    h_dist            = h_dist_pres,        # datacost = Hellinger(pres)
    weight_data       = 1,
    weight_smooth     = smooth_cost,
    nBins             = nbins_total3d,
    seed              = 1,
    verbose           = TRUE,
    rebuild           = TRUE
  )
}, error = function(e) {
  cat("⚠️  GraphCut failed at smooth_cost =", smooth_cost, ":\n", e$message, "\n")
  NULL
})

# result in GC_result$label_attribution, etc.
