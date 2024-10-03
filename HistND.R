# code for computing multivar histograms

library("SBCK")

# needs SBCK
get_histNd <- function(X, bw=NULL, bo=NULL, limvar=NULL){
  # dim(X) = times x n var = times x d
  # length(bw) = d
  # length(bo) = d
  # dim(limvar) = 2 x d where limvar[,v] is the desired range for variable v

  d <- dim(X)[2]
  sh <- SparseHist(X, bin_width=bw, bin_origin = bo)
  shbw <- sh$bin_width
  shbo <- sh$bin_origin

  limshc <- apply(sh$c, 2,FUN=range)
  N_in_shc <- c()
  for(v in 1:d){
    N_in_shc <- c(N_in_shc, round( ((limshc[2,v]+(shbw[v]/2)) - (limshc[1,v]-(shbw[v]/2)) ) / shbw[v]))
  }

  coordH <- list()
  for(v in 1:d){
    coordH[[v]] <- seq(limshc[1,v],limshc[2,v], by=shbw[v])
  }

  indcell <- array(NaN, dim=c(length(sh$c[,1]), d))
  for(i in 1:length(sh$c[,1])){
    for(v in 1:d){
      indcell[i,v] <- which(abs((coordH[[v]])-sh$c[i,v])<1e-9)
    }
  }

  matH <- array(0, dim = N_in_shc) # N_in_shc is a vector
  matH[indcell] <- sh$p

  totsize <- c() ; leftadd <- c() ; rightadd <- c()

  if(!is.null(limvar)){

    if(prod(limvar[1,]<limshc[1,])==TRUE && prod(limvar[2,]>limshc[2,])){
      # False <- 0 ; True <- 1
      # => if prod == FALSE (=0), at least one component is False
      # => if prod == TRUE (=0), all elements are TRUE

      for(v in 1:d){
        limleft <- coordH[[v]][1] - (shbw[v]/2) ; naddleft <- 0
        while(abs(limleft-limvar[1,v])>1e-6){
          naddleft <- naddleft + 1 ; limleft <- limleft - shbw[v]
        }

        limright <- coordH[[v]][(N_in_shc[v])] + (shbw[v]/2)
        naddright <- 0
        while(abs(limright-limvar[2,v])>1e-6){
          naddright <- naddright + 1
          limright <- limright + shbw[v]
        }

        if(naddleft>0){
          newcoordleft <- coordH[[v]][1] - seq(from = naddleft, to=1, by=-1) * shbw[v]
        }else{newcoordleft <- coordH[[v]][1]}
        if(naddright>0){
          newcoordright <- coordH[[v]][(N_in_shc[v])] + seq(from = 1, to=naddright, by=1) * shbw[v]
        }else{newcoordright <- coordH[[v]][(N_in_shc[v])]}
        if(naddleft>0 & naddright>0){coordH[[v]] <- c(newcoordleft, coordH[[v]], newcoordright)}
        if(naddleft>0 & naddright==0){coordH[[v]] <- c(newcoordleft, coordH[[v]])}
        if(naddleft==0 & naddright>0){coordH[[v]] <- c(coordH[[v]], newcoordright)}

        totsize <- c(totsize, (naddleft+naddright))
        leftadd <- c(leftadd, naddleft) ; rightadd <- c(rightadd, naddright)
      } # end for var

      begin <- end <- c()
      for(v in 1:d){
        begin <- c(begin, (leftadd[v]+1))
        end <- c(end, (length(coordH[[v]])-rightadd[v]))
      }
      matHinc0 <- array(0, dim=(dim(matH)+totsize))

      L <- list(matHinc0)
      for(v in 1:d){
        L[[v+1]] <- begin[v]:end[v]
      }
      L[[d+2]] <- matH
      matHinc <- do.call("[<-", L)
      matH <- matHinc
      rm(matHinc0, matHinc) ; gc()

    } else{
      cat("limvar [1,] is NOT < limshc ==> limit of the histogram is not changed\n")
    }

  }

  return(list(p=matH, coord = coordH, leftadd=leftadd, rightadd=rightadd))
}

#############################
# Test with data
#############################
library(MASS)

######### parameters for generating the data
n <- 1000
nd <- 2
m <- c(0,0) ; S <- matrix(c(1,0.3,0.3,1),2,2) ; S

bw <- c(0.1, 0.1)
bo <- c(0,0)

# Generating artificial data
X <- mvrnorm(n, mu=m, Sigma=S) ; cov(X)

limvar <- array(NaN, dim=c(2, length(m)))
for(var in 1:nd){
  limvar[1,var] <- floor(range(X[,var])[1])
  limvar[2,var] <- ceiling(range(X[,var])[2])
}
limvar

H <- get_histNd(X, bw, bo, limvar=NULL)
sum(H$p)

library(fields)
image.plot(H$coord[[1]], H$coord[[2]], H$p, xlab="V1", ylab="V2")

library(ggplot2)
library(plotly)

n <- 10000
nd <- 3
mean_vector <- c(0, 0, 0)
cov_matrix <- matrix(c(1, 0.3, 0.2, 0.3, 1, 0.5, 0.2, 0.5, 1), nrow=3)

data_3d <- mvrnorm(n, mu=mean_vector, Sigma=cov_matrix)
bin_width <- c(0.1, 0.1, 0.1)
bin_origin <- c(0, 0, 0)

hist_3d <- get_histNd(data_3d, bw=bin_width, bo=bin_origin)

x_coords <- hist_3d$coord[[1]]
y_coords <- hist_3d$coord[[2]]
z_coords <- hist_3d$coord[[3]]
values <- hist_3d$p

plot_data <- expand.grid(X = x_coords, Y = y_coords, Z = z_coords)
plot_data$Value <- as.vector(values)

plot_data <- subset(plot_data, Value > 0)

plot_ly(data = plot_data, x = ~X, y = ~Y, z = ~Z, size = ~Value, color = ~Value, colors = c("blue", "red")) %>%
  add_markers() %>%
  layout(scene = list(xaxis = list(title = 'X'),
                      yaxis = list(title = 'Y'),
                      zaxis = list(title = 'Z')),
         title <- "3D Histogram Visualization")
