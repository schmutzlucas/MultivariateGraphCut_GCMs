# Install and load necessary libraries
list_of_packages <- read.table("package_list.txt", sep="\n")$V1
new.packages <- list_of_packages[!(list_of_packages %in% installed.packages()[,"Package"])]
if(length(new.packages))
  install.packages(new.packages, repos = "https://cloud.r-project.org")

library(devtools)
lapply(list_of_packages, library, character.only = TRUE)
install_github("thaos/gcoWrapR")


# Loading local functions
source_code_dir <- 'functions/' #The directory where all functions are saved.
file_paths <- list.files(source_code_dir, full.names = T)
for(path in file_paths){source(path)}



range_var_final <- readRDS('ranges/range_var_final_allModelsPar_1950-2022_new.rds')

# Setting global variables
lon <- -180:179
lat <- -90:90
# Temporal ranges
year_present <<- 1977:1999
year_future <<- 2000:2022
# data directory
data_dir <<- 'data/CMIP6_merged_all/'

# Bins for the kde
nbins1d <<- 32



# List of the variable used
# variables <- c('tas', 'tasmax',
#                'tasmin', 'pr',
#                'psl', 'hur',
#                'huss')
# variables <- c('pr', 'tas', 'tasmin', 'tasmax')
variables <- c('pr', 'tas')

# Obtains the list of models from the model names or from a file
# # Method 1
# dir_path <- paste0('data/CMIP6_merged_all/')
# model_names <- list.dirs(dir_path, recursive = FALSE)
# model_names <- basename(model_names)

# Method 2
# model_names <- c(  'ERA5',
#                                   'MIROC6',
#                                   'IPSL-CM6A-LR',
#                                   'NorESM2-MM',
#                                   'UKESM1-0-LL')

# Method 3
model_names <- read.table('model_names_long.txt')
model_names <- as.list(model_names[['V1']])
# Index of the reference
ref_index <<- 1


tmp <- OpenAndHist2D_range(model_names, variables, year_present, range_var_final)

pdf_matrix <- tmp[[1]]
kde_models <- pdf_matrix[ , , , -ref_index]
kde_ref <- pdf_matrix[ , , , ref_index]
rm(pdf_matrix)
range_matrix <- tmp[[2]]
x_breaks <- tmp[[3]]
y_breaks <- tmp[[4]]


tmp <- OpenAndHist2D_range(model_names, variables, year_future , range_var_final)

pdf_matrix <- tmp[[1]]
kde_models_future <- pdf_matrix[ , , , -ref_index]
kde_ref_future <- pdf_matrix[ , , , ref_index]
range_matrix_future <- tmp[[2]]
x_breaks_future <- tmp[[3]]
y_breaks_future <- tmp[[4]]
rm(pdf_matrix)
rm(tmp)

# Choose the reference in the models
reference_name <<- model_names[ref_index]
model_names <<- model_names[-ref_index]


# Open and average the models for the selected time periods
tmp <- OpenAndAverageCMIP6(
  model_names, variables, year_present, year_future
)
models_list <- tmp[[1]]
models_matrix <- tmp[[2]]
rm(tmp)


tmp <- OpenAndAverageCMIP6(
  reference_name, variables, year_present, year_future
)
reference_list <- tmp[[1]]
reference_matrix <- tmp[[2]]
rm(tmp)

models_matrix_nrm <- list()
models_matrix_nrm <- NormalizeVariables(models_matrix, variables, 'StdSc')

reference_matrix_nrm <- list()
reference_matrix_nrm <- NormalizeVariables(reference_matrix, variables, 'StdSc')


# Get the current date and time

current_time <- Sys.time()

# Format the date and time as a string in the format 'yyyymmddhhmm'
formatted_time <- format(current_time, "%Y%m%d%H%M")

# Concatenate the formatted time string with your desired filename
filename <- paste0(formatted_time, "_my_workspace_ERA5_allModels_beforeOptim.RData")

# Save the workspace using the generated filename
save.image(file = filename, compress = FALSE)


# Computing the sum of hellinger distances between models and reference --> used as datacost
h_dist <- array(data = 0, dim = c(length(lon), length(lat),
                                  length(model_names)))
h_dist_unchecked <- array(data = 0, dim = c(length(lon), length(lat),
                                            length(model_names)))

# Loop through variables and models
m <- 1
for (model_name in model_names) {
  for (i in seq_along(lon)) {
    for (j in seq_along(lat)) {
      # Compute Hellinger distance
      h_dist_unchecked[i, j, m] <- sqrt(sum((sqrt(kde_models[i, j, , m]) - sqrt(kde_ref[i, j, ]))^2)) / sqrt(2)
    }
  }
  m <- m + 1
}

hist(h_dist_unchecked)
# Replace NaN with 0
h_dist[,,] <- replace(h_dist_unchecked[,,], is.nan(h_dist_unchecked), 0)
hist(h_dist)
rm(h_dist_unchecked)


# Computing the sum of hellinger distances between models and reference --> used as datacost
h_dist_future <- array(data = 0, dim = c(length(lon), length(lat),
                                         length(model_names)))

h_dist_unchecked <- array(data = 0, dim = c(length(lon), length(lat),
                                            length(model_names)))

# Loop through variables and models
m <- 1
for (model_name in model_names) {
  for (i in seq_along(lon)) {
    for (j in seq_along(lat)) {
      # Compute Hellinger distance
      h_dist_unchecked[i, j, m] <- sqrt(sum((sqrt(kde_models_future[i, j, , m]) - sqrt(kde_ref_future[i, j, ]))^2)) / sqrt(2)
    }
  }
  m <- m + 1
}
hist(h_dist_unchecked)
# Replace NaN with 0
h_dist_future[,,] <- replace(h_dist_unchecked[,,], is.nan(h_dist_unchecked), 0)
hist(h_dist_future)
rm(h_dist_unchecked)




# Graphcut hellinger labelling
GC_result_hellinger_mixte <- list()
GC_result_hellinger_mixte <- GraphCutHellinger2D(kde_ref = kde_ref,
                                                           models_smoothcost = models_matrix_nrm$future,
                                                           h_dist = h_dist,
                                                           weight_data = 1,
                                                           weight_smooth = 1,
                                                           verbose = TRUE)


# Graphcut hellinger labelling
GC_result_hellinger_new2 <- list()
GC_result_hellinger_new2 <- GraphCutHellinger2D_new2(pdf_ref = kde_ref,
                                                     kde_models = kde_models,
                                                     pdf_models_future = kde_models_future,
                                                     h_dist = h_dist,
                                                     weight_data = 1,
                                                     weight_smooth = 1,
                                                     nBins = nbins1d^2,
                                                     verbose = TRUE,
                                                     rebuild = TRUE)


# Get the current date and time

current_time <- Sys.time()

# Format the date and time as a string in the format 'yyyymmddhhmm'
formatted_time <- format(current_time, "%Y%m%d%H%M")

# Concatenate the formatted time string with your desired filename
filename <- paste0(formatted_time, "_my_workspace_ERA5_allModels_mixte.RData")

# Save the workspace using the generated filename
save.image(file = filename, compress = FALSE)




GC_hellinger_projections <- list()
j <- 1
for(var in variables){
  for(l in 0:(length(model_names))){
    islabel <- which(GC_result_hellinger_new2$label_attribution == l)
    GC_hellinger_projections[[var]][islabel] <- models_matrix$future[,,(l),j][islabel]
  }
  j <- j + 1
}
GC_hellinger_projections$tas <- matrix(GC_hellinger_projections$tas, nrow = 360)
GC_hellinger_projections$pr <- matrix(GC_hellinger_projections$pr, nrow = 360)





h_dist_map <- array(NA, dim = dim(GC_result_hellinger_new2$label_attribution))

for(j in seq_along(variables)){
  for(l in 0:(length(model_names))){
    islabel <- which(GC_result_hellinger_new2$label_attribution == l)
    h_dist_map[islabel] <- h_dist_future[,,(l)][islabel]
  }
}

test_df <- melt(h_dist_map, c("lon", "lat"), value.name = "Bias")

test_df$lon[test_df$lon > 180] <- test_df$lon[test_df$lon > 180] - 360


p5 <- ggplot() +
  geom_tile(data=test_df, aes(x=lon, y=lat-90, fill=Bias))+
  labs(subtitle = 'Projection period : 1999 - 2022')+
  ggtitle(paste0('GraphCut Hellinger', ': Mean Hellinger distance = ', round(mean(h_dist_map), 2)))+
  scale_fill_gradient(low = "white", high = "#015a8c", limits = c(0.1, 0.70), oob = scales::squish)+
  borders("world", colour = 'black', lwd = 0.12) +
  scale_x_continuous(, expand = c(0, 0)) +
  scale_y_continuous(, expand = c(0,0))+
  theme(legend.position = 'bottom')+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(panel.background = element_blank())+
  xlab('Longitude')+
  ylab('Latitude') +
  labs(fill='Hellinger \nDistance')+
  theme_bw()+
  theme(legend.key.size = unit(1, 'cm'), #change legend key size
        legend.key.height = unit(1.4, 'cm'), #change legend key height
        legend.key.width = unit(0.4, 'cm'), #change legend key width
        legend.title = element_text(size=16), #change legend title font sizen
        legend.text = element_text(size=12))+ #change legend text font size
  theme(plot.title = element_text(size=24),
        plot.subtitle = element_text(size = 20,hjust=0.5),
        axis.text=element_text(size=14),
        axis.title=element_text(size=16),)+
  easy_center_title()
p5


name <- paste0('figure/H_dist_future_GC_hellinger_new_26models8')
ggsave(paste0(name, '.pdf'), plot = p5, width = 35, height = 25, units = "cm", dpi = 300)
ggsave(paste0(name, '.png'), plot = p5, width = 35, height = 25, units = "cm", dpi = 300)



# Map of the bias that works
# Creating the dataframe
# For tas
bias_tmp <- GC_hellinger_projections$tas -
  reference_list$future$tas[[1]]

test_df <- melt(bias_tmp, c("lon", "lat"), value.name = "Bias")

# Settings the limits
# limit <- max(abs(c(min(test_df$Bias), max(test_df$Bias))))
# limit <- round_any(limit, 10, floor)

limit <- 8
limits <- c(-8, 8)
v_limits <- as.numeric(format(seq(-limit, limit, len=6), digits = 3))

test_df$Bias[test_df$Bias < -limit] <- -limit
test_df$Bias[test_df$Bias > limit] <- limit

p <- ggplot() +
  geom_tile(data=test_df, aes(x=lon, y=lat-90, fill=Bias),)+
  ggtitle(paste0('GraphCut Hellinger'))+
  labs(subtitle = paste0('Average Bias : ', round(mean(abs(bias_tmp)), 2)))+
  scale_fill_gradientn(colours = rev(RColorBrewer::brewer.pal(11, "RdBu")),
                       breaks = v_limits,
                       limits = limits)+
  # guide = guide_colourbar(
  #   barwidth = 0.4,
  #   ticks.colour = 'black',
  #   ticks.linewidth = 1,
  #   frame.colour = 'black',
  #   draw.ulim = TRUE,
  #   draw.llim = TRUE),)
  borders("world2", colour = 'black', lwd = 0.12) +
  scale_x_continuous(, expand = c(0, 0)) +
  scale_y_continuous(, expand = c(0,0))+
  theme(legend.position = 'bottom')+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(panel.background = element_blank())+
  xlab('Longitude')+
  ylab('Latitude') +
  labs(fill='Bias [K]')+
  theme_bw()+
  theme(legend.key.size = unit(1, 'cm'), #change legend key size
        legend.key.height = unit(1.4, 'cm'), #change legend key height
        legend.key.width = unit(0.4, 'cm'), #change legend key width
        legend.title = element_text(size=16), #change legend title font size
        legend.text = element_text(size=12))+ #change legend text font size
  theme(plot.title = element_text(size=24),
        plot.subtitle = element_text(size = 20,hjust=0.5),
        axis.text=element_text(size=14),
        axis.title=element_text(size=16),)+
  easy_center_title()
p
mean(abs(bias_tmp))

name <- paste0('figure/GC_Hellinger_bias_tas_26models_new7')
ggsave(paste0(name, '.pdf'), plot = p, width = 35, height = 25, units = "cm", dpi = 300)
ggsave(paste0(name, '.png'), plot = p, width = 35, height = 25, units = "cm", dpi = 300)



# Map of the bias that works
# Creating the dataframe
# For pr
bias_tmp <- GC_hellinger_projections$pr - reference_list$future$pr[[1]]

# Ensure you are removing NA values in the computation
rm(avg_bias)
avg_bias <- round(mean(abs(bias_tmp), na.rm = TRUE), 2)

# Settings the limits
limit <- c(-10, 10)

# Settings the limits

v_limits <- as.numeric(format(seq(min(limit), max(limit), len=5), digits = 3))
bias_tmp[bias_tmp < min(limit)] <- min(limit)
bias_tmp[bias_tmp > max(limit)] <- max(limit)

test_df <- melt(bias_tmp, c("lon", "lat"), value.name = "Bias")

p <- ggplot() +
  geom_tile(data=test_df, aes(x=lon, y=lat-90, fill=Bias),)+
  ggtitle(paste0('GraphCut Hellinger'))+
  labs(subtitle = paste0('Average Bias : ', avg_bias))+
  scale_fill_gradientn(colours = rev(RColorBrewer::brewer.pal(11, "RdBu")),
                       breaks = v_limits,
                       limits = limit)+
  # guide = guide_colourbar(
  #   barwidth = 0.4,
  #   ticks.colour = 'black',
  #   ticks.linewidth = 1,
  #   frame.colour = 'black',
  #   draw.ulim = TRUE,
  #   draw.llim = TRUE),)
  borders("world2", colour = 'black', lwd = 0.12) +
  scale_x_continuous(, expand = c(0, 0)) +
  scale_y_continuous(, expand = c(0,0))+
  theme(legend.position = 'bottom')+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(panel.background = element_blank())+
  xlab('Longitude')+
  ylab('Latitude') +
  labs(fill='[mm/day]')+
  theme_bw()+
  theme(legend.key.size = unit(1, 'cm'), #change legend key size
        legend.key.height = unit(1.4, 'cm'), #change legend key height
        legend.key.width = unit(0.4, 'cm'), #change legend key width
        legend.title = element_text(size=16), #change legend title font size
        legend.text = element_text(size=12))+ #change legend text font size
  theme(plot.title = element_text(size=24),
        plot.subtitle = element_text(size = 20,hjust=0.5),
        axis.text=element_text(size=14),
        axis.title=element_text(size=16),)+
  easy_center_title()
p
mean(abs(bias_tmp))

name <- paste0('figure/GC_Hellinger_new_bias_pr_26models_new7')
ggsave(paste0(name, '.pdf'), plot = p, width = 35, height = 25, units = "cm", dpi = 300)
ggsave(paste0(name, '.png'), plot = p, width = 35, height = 25, units = "cm", dpi = 300)





h_dist_map <- array(NA, dim = dim(GC_result_hellinger_new2$label_attribution))

for(j in seq_along(variables)){
  for(l in 0:(length(model_names))){
    islabel <- which(GC_result_hellinger_new2$label_attribution == l)
    h_dist_map[islabel] <- h_dist[,,(l)][islabel]
  }
}

test_df <- melt(h_dist_map, c("lon", "lat"), value.name = "Bias")

p5 <- ggplot() +
  geom_tile(data=test_df, aes(x=lon, y=lat-90, fill=Bias))+
  labs(subtitle = 'Projection period : 1977 - 1999')+
  ggtitle(paste0('GraphCut Hellinger', ': Mean Hellinger distance = ', round(mean(h_dist_map), 2)))+
  scale_fill_gradient(low = "white", high = "#015a8c", limits = c(0.1, 0.70), oob = scales::squish)+
  borders("world2", colour = 'black', lwd = 0.12) +
  scale_x_continuous(, expand = c(0, 0)) +
  scale_y_continuous(, expand = c(0,0))+
  theme(legend.position = 'bottom')+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(panel.background = element_blank())+
  xlab('Longitude')+
  ylab('Latitude') +
  labs(fill='Hellinger \nDistance')+
  theme_bw()+
  theme(legend.key.size = unit(1, 'cm'), #change legend key size
        legend.key.height = unit(1.4, 'cm'), #change legend key height
        legend.key.width = unit(0.4, 'cm'), #change legend key width
        legend.title = element_text(size=16), #change legend title font sizen
        legend.text = element_text(size=12))+ #change legend text font size
  theme(plot.title = element_text(size=24),
        plot.subtitle = element_text(size = 20,hjust=0.5),
        axis.text=element_text(size=14),
        axis.title=element_text(size=16),)+
  easy_center_title()
p5

name <- paste0('figure/H_dist_future_GC_hellinger_new_26models_present7')
ggsave(paste0(name, '.pdf'), plot = p5, width = 35, height = 25, units = "cm", dpi = 300)
ggsave(paste0(name, '.png'), plot = p5, width = 35, height = 25, units = "cm", dpi = 300)



n1 <- dim(kde_models)[1]
n2 <- dim(kde_models)[2]
n3 <- dim(kde_models)[3]
n4 <- dim(kde_models)[4]

MMM_KDE_future <- array(0, dim = c(n1, n2, n3))

for (i in 1:n1) {
  for (j in 1:n2) {
    for (k in 1:n3) {
      MMM_KDE_future[i,j,k] <- mean(kde_models_future[i,j,k,])
    }
  }
}



# Computing the sum of hellinger distances between models and reference --> used as datacost
MMM_h_dist_future <- array(data = 0, dim = dim(h_dist_map))
h_dist_unchecked <- array(data = 0, dim = dim(h_dist_map))

# Loop through grip-points
for (i in 1:360) {
  for (j in 1:181) {
    # Compute Hellinger distance
    h_dist_unchecked[i,j] <- sqrt(sum((sqrt(MMM_KDE_future[i,j, ]) - sqrt(kde_ref[i, j, ]))^2)) / sqrt(2)

  }
}

# Replace NaN with 0
MMM_h_dist_future[, ] <- replace(h_dist_unchecked[,], is.nan(h_dist_unchecked), 0)

MMM <- list()
MMM$tas <- apply(models_matrix$future[,,,2], c(1, 2), mean)

MMM$pr <- apply(models_matrix$future[,,,1], c(1, 2), mean)