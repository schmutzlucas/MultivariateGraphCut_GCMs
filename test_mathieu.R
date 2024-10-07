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

# Initialize an empty list to store the results
GC_result_hellinger_test <- list()

# Compute the weights based on the current lambdas
weight_data <- 1
weight_smooth <- 1

# Use tryCatch to handle errors
GC_result_hellinger_test <- tryCatch({
  # Call the GraphCutHellinger2D_new2 function and store the result
  GraphCutHellinger2D_new2(pdf_ref = kde_ref,
                           kde_models = kde_models,
                           pdf_models_future = kde_models_future,
                           h_dist = h_dist,
                           weight_data = weight_data,
                           weight_smooth = weight_smooth,
                           nBins = nbins1d^2,
                           verbose = TRUE,
                           rebuild = FALSE)
}, error = function(e) {
  # Handle the error, e.g., by printing a message and returning NULL or a specific error value
  cat("An error occurred in iteration", i, "with lambda =", lambdas[i], "\nError message:", e$message, "\n")
  # Optionally, return a specific value indicating the error
  NULL  # Or any other value that makes sense in your context
})




h_dist_future_lambda11 <- list()
avg_hdist_future_lambda11 <- list()

tmp <- array(NA, dim = dim(GC_result_hellinger_test$label_attribution))
for(j in seq_along(variables)){
  for(l in 0:(length(model_names))){
    islabel <- which(GC_result_hellinger_test$label_attribution == l)
    tmp[islabel] <- h_dist_future[,,(l)][islabel]

  }
}

h_dist_future_lambda11 <- tmp
avg_hdist_future_lambda11 <- mean(h_dist_future_lambda01[[i]])



test_df <- melt(tmp, c("lon", "lat"), value.name = "Bias")

p5 <- ggplot() +
  geom_tile(data=test_df, aes(x=lon, y=lat-90, fill=Bias))+
  labs(subtitle = 'Projection period : 2000-20022')+
  ggtitle(paste0('GraphCut Hellinger new', ': Mean Hellinger distance = ', round(mean(h_dist_map), 2)))+
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

name <- paste0('figure/H_dist_future_GC_hellinger_new_26models')
ggsave(paste0(name, '.pdf'), plot = p5, width = 35, height = 25, units = "cm", dpi = 300)
ggsave(paste0(name, '.png'), plot = p5, width = 35, height = 25, units = "cm", dpi = 300)
