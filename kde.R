
# Kernel density estimation -----------------------------------------------

library(tidyverse)
library(Rcpp)
library(microbenchmark)
library(tictoc)
sourceCpp("kde_estimation.cpp")

#Generate random data
sample <-  rweibull(10000, 3)

#Give values to estimation function
x_grid <- seq(min(sample), max(sample), length.out = 512)

#Maximal smoothing principle window parameter
h <- 3*sd(sample)*(0.28/(35*length(sample)))^(1/5)
dens_sample<- kde_cpp(sample, x_grid, 0.01)

par(mfrow=c(1,2))
#Plot results
hist(sample, probability = T, breaks = 100, main = "KDE with h=0.01")
lines(x_grid,dens_sample, col = "red")
