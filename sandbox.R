
# Libraries ---------------------------------------------------------------
library(data.table)
library(tidyverse)
library(Rcpp)
library(tictoc)
sourceCpp("kde_estimation.cpp")
sourceCpp("logconctest.cpp")

# Read data ---------------------------------------------------------------

test0  <- data.table::fread("C:\\Users\\bscuser\\Projects\\wcet-test0\\long-test0.txt",header = F,
                                   col.names = "X")

tic()
l_test <- function(h,s){logconctest_cpp(kde_cpp(s, x_grid, h))}


# Prepare tail sample -----------------------------------------------------

sample_data <- sample(test0$X, size = 1000, replace = T)
tail_size <- c(50, 100, 200, 500)
threshold <- sort(sample_data)[1000-tail_size[4]]
sample_tail <- sample_data[sample_data > threshold]/threshold

#Grids for KDE 
x_grid <- seq(min(sample_tail), max(sample_tail), length.out = 500)
h_grid <- seq(0.0001, 0.1, length.out = 1000)

hist(sample_tail, probability = T, breaks = 70)
lines(x_grid, kde_cpp(sample_tail, x_grid, 0.009))

ggplot(data = data.frame(sample_tail), aes(x = sample_tail))+
  geom_histogram(aes(y=..density..), colour = "black", fill ="#E69F00")+
  geom_line(aes(x = x_grid), y = kde_cpp(sample_tail, x_grid, 0.009), colour = "purple4")+
  ggtitle("Histogram for sample of the tail with KDE N=500")+
  theme_bw()
# Binary Search -----------------------------------------------------------

#Variable for binary search
max_index <- length(h_grid)
low_index <- 0
index <- 0
while(abs(max_index - low_index) >1){
  index <- floor((max_index-low_index)/2) +low_index
  if(l_test(h_grid[index], sample_tail)==0){
    low_index <- index
  }
  else {
    max_index <- index
  }
}

if(max_index<low_index){index <- index+1}
#Bandwidth value to compare with other samples
h_0 <- h_grid[index]


# Bootstrap and p-avlue ---------------------------------------------------


h_star <- c()

#Perform again the test for 1000 bootstrap samples of tail
# and store the value of h_star
for(i in 1:1000){
    
  sample_tail_2 <- sample(sample_tail, length(sample_tail), replace = T)
  x_grid <- seq(min(sample_tail_2), max(sample_tail_2), length.out = 500)
  h_grid <- seq(0.0001, 0.1, length.out = 1000)
  max_index <- length(h_grid)
  low_index <- 0
  index <- 0
    while(abs(max_index - low_index) >1){
      index <- floor((max_index-low_index)/2) +low_index
      if(l_test(h_grid[index], sample_tail_2)==0){
        low_index <- index
      }
      else {
        max_index <- index
      }
    }
  
  if(max_index<low_index){index <- index+1}
  h_star[i]<- h_grid[index]
}

#The p value is simply the proportion of h_star values above h_0
p <- length(h_star[h_star > h_0])/length(h_star)
toc()

ggplot(data = data.frame(h_star), aes(x = h_star))+
  geom_histogram(aes(y=..density..), colour = "black", fill ="#E69F00")+
  geom_vline(xintercept=h_0)+
  scale_x_continuous(breaks = c(pretty(h_star), h_0), labels = c(pretty(h_star), "h_0"))+
  ggtitle("Histogram of h star for sample of the tail N=500")+
  theme_bw()




