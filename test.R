# Generate two time series of random values
x <- rnorm(1e5)
y <- x + 2

# Estimate the probability density functions of the time series using kernel density estimation
fx <- density(x, from = -4, to = 11, n = 100)
fy <- density(y, from = -4, to = 11, n = 100)

# Compute the Hellinger distance between the probability density functions
h_dist <- sqrt(sum((sqrt(fx$y) - sqrt(fy$y))^2)) / sqrt(2)


library(distr)
ts1 <- ts(c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10))
ts2 <- ts(c(1, 4, 9, 16, 25, 36, 49, 64, 81, 100))
# Define two probability distributions
p1 <- density(ts1)
p2 <- density(ts2)

# Compute the Hellinger distance
hellinger_distance <- sqrt(sum((sqrt(p1$y) - sqrt(p2$y))^2)) / sqrt(2)



library(MASS)

# Generate two multivariate normal distributions with different means and covariance matrices
mu1 <- c(1, 2)
mu2 <- c(12, 22)
sigma1 <- matrix(c(1, 0.5, 0.5, 2), ncol = 2)
sigma2 <- matrix(c(1, 0.5, 0.5, 2), ncol = 2)
dist1 <- mvrnorm(1000, mu1, sigma1)
dist2 <- mvrnorm(1000, mu2, sigma2)

# Compute the Mahalanobis distance
mahalanobis_dist <- mahalanobis(dist1, dist2, cov = (sigma1 + sigma2) / 2)

# Compute the Euclidean distance
euclidean_dist <- dist(rbind(mu1, mu2))

# Print the distances
cat("Mahalanobis distance:", mahalanobis_dist, "\n")
cat("Euclidean distance:", euclidean_dist, "\n")

library(ggplot2)

# Plot the first variable against the second variable of dist1
ggplot(data.frame(dist2), aes(x = dist2[,1], y = dist2[,2])) +
  geom_point() +
  xlab("Variable 1") +
  ylab("Variable 2")


# Load the `transport` package
library(transport)

# Generate two univariate distributions with different means and variances
x <- rnorm(100, mean = 1, sd = 1)
y <- rnorm(100, mean = 1, sd = 1)

# Calculate the Wasserstein distance between the two distributions
dist <- wasserstein1d(x, y)

# Print the results
print(dist)


# Load required packages
library(ggplot2)
library(dplyr)
library(transport)

# Set correlation coefficient and generate sample data
set.seed(123)
n <- 1000
correlation <- 0.5
temp_model1 <- rnorm(n, mean = 20, sd = 2)
precip_model1 <- rnorm(n, mean = 50 + correlation * (temp_model1 - 20), sd = 10)
temp_model2 <- rnorm(n, mean = 18, sd = 3)
precip_model2 <- rnorm(n, mean = 55 + correlation * (temp_model2 - 18), sd = 15)

# Combine the data into a single data frame
data <- data.frame(
  model = rep(c("Model 1", "Model 2"), each = n),
  temperature = c(temp_model1, temp_model2),
  precipitation = c(precip_model1, precip_model2)
)

# Compute the joint distribution for each model
joint_dist_model1 <- data %>% filter(model == "Model 1") %>%
  ggplot(aes(x = temperature, y = precipitation)) +
  geom_point(alpha = 0.5) +
  labs(title = "Joint Distribution of Temperature and Precipitation - Model 1")

joint_dist_model2 <- data %>% filter(model == "Model 2") %>%
  ggplot(aes(x = temperature, y = precipitation)) +
  geom_point(alpha = 0.5) +
  labs(title = "Joint Distribution of Temperature and Precipitation - Model 2")



{
# Load required packages
library(diverse)

# Generate sample data for two climate models
set.seed(123)
n <- 1000
correlation <- 0.5

# Climate model 1
temp_model1 <- rnorm(n, mean = 20, sd = 2)
precip_model1 <- rnorm(n, mean = 50 + correlation * (temp_model1 - 20), sd = 10)
joint_dist_model1 <- cbind(temp_model1, precip_model1)

# Climate model 2
temp_model2 <- rnorm(n, mean = 18, sd = 3)
precip_model2 <- rnorm(n, mean = 55 + correlation * (temp_model2 - 18), sd = 15)
joint_dist_model2 <- cbind(temp_model2, precip_model2)

library(MASS) # Load MASS package for kde2d function

# Combine the two time series into a single matrix
joint_dist <- rbind(joint_dist_model1, joint_dist_model2)

# Compute the kernel density estimate of the joint PDF
pdf_est <- kde2d(joint_dist[,1], joint_dist[,2])

# Plot the joint PDF
contour(pdf_est, xlab = "Temperature", ylab = "Precipitation", main = "Joint PDF")

  library(diverse) # Load diverse package for JSD function

# Compute the kernel density estimate of the joint PDF for model 1
pdf_est_model1 <- kde2d(joint_dist_model1[,1], joint_dist_model1[,2])

# Compute the kernel density estimate of the joint PDF for model 2
pdf_est_model2 <- kde2d(joint_dist_model2[,1], joint_dist_model2[,2])

# Compute the Jensen-Shannon Divergence between the two models
jsd <- JSD(cbind(c(pdf_est_model1$z), c(pdf_est_model2$z)))

# Print the JSD value
cat("Jensen-Shannon Divergence between the two models:", jsd, "\n")




}


{
  library(ks) # Load ks package for kde3d function

# Climate model 1
n <- 1000
correlation <- 0.5
temp_model1 <- rnorm(n, mean = 20, sd = 2)
precip_model1 <- rnorm(n, mean = 50 + correlation * (temp_model1 - 20), sd = 10)
humidity_model1 <- rnorm(n, mean = 70 + correlation * (temp_model1 - 20), sd = 5)
joint_dist_model1 <- cbind(temp_model1, precip_model1, humidity_model1)

# Climate model 2
temp_model2 <- rnorm(n, mean = 18, sd = 3)
precip_model2 <- rnorm(n, mean = 55 + correlation * (temp_model2 - 18), sd = 15)
humidity_model2 <- rnorm(n, mean = 75 + correlation * (temp_model2 - 18), sd = 10)
joint_dist_model2 <- cbind(temp_model2, precip_model2, humidity_model2)

# Combine the two time series into a single matrix
joint_dist <- rbind(joint_dist_model1, joint_dist_model2)

# Compute the kernel density estimate of the joint PDF
pdf_est <- kde3d(joint_dist[,1], joint_dist[,2], joint_dist[,3], n = 50)

# Plot the joint PDF
persp(pdf_est, xlab = "Temperature", ylab = "Precipitation", zlab = "Humidity", main = "Joint PDF")

}

{library(ks) # Load ks package for kde3d function

# Climate model 1
n <- 10000
correlation <- 0.5
temp_model1 <- rnorm(n, mean = 20, sd = 2)
precip_model1 <- rnorm(n, mean = 50 + correlation * (temp_model1 - 20), sd = 10)
humidity_model1 <- rnorm(n, mean = 70 + correlation * (temp_model1 - 20), sd = 5)
joint_dist_model1 <- cbind(temp_model1, precip_model1, humidity_model1)

# Climate model 2
temp_model2 <- rnorm(n, mean = 18, sd = 3)
precip_model2 <- rnorm(n, mean = 55 + correlation * (temp_model2 - 18), sd = 15)
humidity_model2 <- rnorm(n, mean = 75 + correlation * (temp_model2 - 18), sd = 10)
joint_dist_model2 <- cbind(temp_model2, precip_model2, humidity_model2)

# Combine the two time series into a single matrix
joint_dist <- rbind(joint_dist_model1, joint_dist_model2)

# Compute the kernel density estimate of the joint PDF
pdf_est <- kde3d(joint_dist[,1], joint_dist[,2], joint_dist[,3], n = 80)


}

{library(ks) # Load ks package for kde function

# Climate model 1
n <- 100
correlation <- 0.5
temp_model1 <- rnorm(n, mean = 20, sd = 2)
precip_model1 <- rnorm(n, mean = 50 + correlation * (temp_model1 - 20), sd = 10)
humidity_model1 <- rnorm(n, mean = 70 + correlation * (temp_model1 - 20), sd = 5)
wind_model1 <- rnorm(n, mean = 10 + correlation * (temp_model1 - 20), sd = 2)
joint_dist_model1 <- cbind(temp_model1, precip_model1, humidity_model1, wind_model1)

# Climate model 2
temp_model2 <- rnorm(n, mean = 18, sd = 3)
precip_model2 <- rnorm(n, mean = 55 + correlation * (temp_model2 - 18), sd = 15)
humidity_model2 <- rnorm(n, mean = 75 + correlation * (temp_model2 - 18), sd = 10)
wind_model2 <- rnorm(n, mean = 12 + correlation * (temp_model2 - 18), sd = 3)
joint_dist_model2 <- cbind(temp_model2, precip_model2, humidity_model2, wind_model2)

# Combine the two time series into a single matrix
joint_dist <- rbind(joint_dist_model1, joint_dist_model2)

# Compute the kernel density estimate of the joint PDF
pdf_est <- kde(joint_dist)

# Plot the joint PDF

}

