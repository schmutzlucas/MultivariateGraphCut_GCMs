

n1 <- dim(kde_models)[1]
n2 <- dim(kde_models)[2]
n3 <- dim(kde_models)[3]
n4 <- dim(kde_models)[4]

MMM_KDE <- array(0, dim = c(n1, n2, n3))

for (i in 1:n1) {
  for (j in 1:n2) {
    for (k in 1:n3) {
      MMM_KDE[i,j,k] <- mean(kde_models[i,j,k,])
    }
  }
}
