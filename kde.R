
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
h <- 1.06*sd(sample)*length(sample)^(-1/5)
dens_sample<- kde_cpp(sample, x_grid, h)

#Plot results
hist(sample, probability = T, breaks = 200)
lines(x_grid,dens_sample)

