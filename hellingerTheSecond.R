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