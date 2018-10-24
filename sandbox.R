
# Libraries ---------------------------------------------------------------
library(data.table)
library(viridis)
library(tidyverse)
library(Rcpp)
library(tictoc)
sourceCpp("kde_estimation.cpp")
sourceCpp("logconctest.cpp")

# Read data ---------------------------------------------------------------

test0  <- data.table::fread("C:\\Users\\bscuser\\Projects\\wcet-test0\\long-test0.txt",header = F,
                                   col.names = "X")

#Compute empirical cumulative density function
a <-ecdf(test0$X)

#Get the value at quantule 1-10^-6 for latter comparison 
emp_quantile <- as.numeric(quantile(test0$X, 1-10^-6))

lweibullt<-function(x,u,a,b) length(x)*log(a)+length(x)*log(b)+(b-1)*sum(log(x))-a*sum(x^b-u^b)
qweibullt<-function(p,u,a,b) ((-log(1-p)/a)+u^b)^(1/b)

tic()

total_p_values <- list()
total_quantiles <- list()
tail_size <- c(50, 100, 200, 500)


# Prepare tail sample -----------------------------------------------------
for (k in 1:1000) {
  print(k)
  sample_data <- sample(test0$X, size = 1000, replace = T)
  ecdf_test_sample <- ecdf(sample_data)
  p_values <- c()
  quantiles <- c()
  for(j in 1:length(tail_size)){

    threshold <- sort(sample_data)[1000-tail_size[j]]
    sample_tail <- sample_data[sample_data > threshold]/threshold
    
    #Fit
    aux<-function(p){-lweibullt(sample_tail,u=1,a=1, b=p)}
    solfit <- nlm(aux, 10)
    shape0 <- solfit$estimate
    
    aux<-function(p){-lweibullt(sample_tail,u=1,a=p, b=1)}
    solfit <- nlm(aux, 10)
    scale0 <- solfit$estimate
    
    aux<-function(p){-lweibullt(sample_tail,u=1,a=p[1], b=p[2])}
    solfit <- nlm(aux, c(scale0,shape0))
    scale <- solfit$estimate[1]
    shape <- solfit$estimate[2]
    
    #Compute quantile
    lambda <- ecdf_test_sample(threshold)
    p <- (1-10^-6 - lambda)/(1 - lambda)
    quantiles[j] <-  qweibullt(p, u = 1, b = shape, a = scale)*threshold/emp_quantile
    
    #Grids for KDE 
    h_grid <- seq(0.0001, 0.1, length.out = 1000)
    x_grid <- seq(min(sample_tail), max(sample_tail), length.out = length(sample_tail))
    # Binary Search -----------------------------------------------------------
    
    #Variable for binary search
    max_index <- length(h_grid)
    low_index <- 0
    index <- 0
    while(abs(max_index - low_index) >1){
      index <- floor((max_index-low_index)/2) +low_index
      if(logconctest_cpp(kde_cpp(sample_tail, x_grid, h_grid[index]))==0){
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
    
    sum <- 0
    
    #Perform again the test for 1000 bootstrap samples of tail
    # and store the value of h_star
    for(i in 1:1000){
      
      sample_tail_2 <- sample(sample_tail, length(sample_tail), replace = T)
      x_grid <- seq(min(sample_tail_2), max(sample_tail_2), length.out = length(sample_tail_2))
      sum <- sum + logconctest_cpp(kde_cpp(sample_tail_2, x_grid,h_0))
    }
    
    #The p value is simply the proportion of h_star values above h_0
    p_values[j] <- sum/1000
  }
  
  total_p_values[[k]] <- p_values
  total_quantiles[[k]] <- quantiles
}

toc()

sim_data <- data.frame(p_values = unlist(total_p_values), quantiles = unlist(total_quantiles))
sim_data$tail <- as.factor(rep(tail_size, 1000))
sim_data$sim <- as.factor(rep(seq(1,1000), each=4))
sim_data$test <- as.factor(rep(0,4000))

ggplot(data = sim_data, aes( y = p_values, group = tail, x = tail, fill = tail))+
  geom_boxplot()+
  scale_fill_viridis(discrete = T)+
  theme_bw()+
  ggtitle("Distribution of p-values for log-concavity test with different tail sizes")


write.csv(sim_data, "sim_data_t.csv", row.names = F)
ggplot(data = data.frame(sample_tail), aes(x = sample_tail))+
  geom_histogram(aes(y=..density..), colour = "black", fill ="#E69F00")+
  geom_line(aes(x = x_grid), y = kde_cpp(sample_tail, x_grid, 0.009), colour = "purple4")+
  ggtitle("Histogram for sample of the tail with KDE N=500")+
  theme_bw()


