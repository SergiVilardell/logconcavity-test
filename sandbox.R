
# Libraries ---------------------------------------------------------------
library(data.table)
library(Rcpp)
library(tictoc)
sourceCpp("kde_estimation.cpp")
sourceCpp("logconctest.cpp")

# Read data ---------------------------------------------------------------

test0  <- data.table::fread("C:\\Users\\bscuser\\Projects\\wcet-test0\\long-test0.txt",header = F,
                                   col.names = "X")

l_test <- function(h){logconctest_cpp(kde_cpp(sample_tail, x_grid, h))}


sample_data <- sample(test0$X, size = 1000, replace = T)
tail_size <- c(50, 100, 200, 500)
threshold <- sort(sample_data)[1000-tail_size[4]]
sample_tail <- sample_data[sample_data > threshold]/threshold
x_grid <- seq(min(sample_tail), max(sample_tail), length.out = 512)
h_grid <- seq(0.001, 0.1, length.out = 500)
max_index <- length(h_grid)
low_index <- 0
index <- 0

while(abs(max_index - low_index) >1){
  index <- floor((max_index-low_index)/2) +low_index
  if(l_test(h_grid[index])==0){
    low_index <- index
  }
  else {
    max_index <- index
  }
}

if(max_index<low_index){index <- index+1}
h_0 <- h_grid[index]
h_star <- c()

for(i in 1:1000){
    
  sample_tail_2 <- sample(sample_tail, length(sample_tail_2), replace = T)
  x_grid <- seq(min(sample_tail_2), max(sample_tail_2), length.out = 512)
  h_grid <- seq(0.001, 0.1, length.out = 500)
  max_index <- length(h_grid)
  low_index <- 0
  index <- 0
    while(abs(max_index - low_index) >1){
      index <- floor((max_index-low_index)/2) +low_index
      if(l_test(h_grid[index])==0){
        low_index <- index
      }
      else {
        max_index <- index
      }
    }
  
  if(max_index<low_index){index <- index+1}
  h_star[i]<- h_grid[index]
}

p <- length(h_star[h_star > h_0])/length(h_star)



