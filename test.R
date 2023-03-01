# Generate two time series of random values
x <- rnorm(1e8)
y <- rnorm(1e8, mean = 1000)

# Estimate the probability density functions of the time series using kernel density estimation
fx <- density(x)
fy <- density(y)

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

