
# Libraries ---------------------------------------------------------------
library(data.table)
library(viridis)
library(tidyverse)
library(Rcpp)
library(gridExtra)
library(tictoc)
sourceCpp("kde_estimation.cpp")
sourceCpp("logconctest.cpp")

# Read data ---------------------------------------------------------------

test_file <- c("/home/bill/BSC/data/long-test8.txt", "/home/bill/BSC/data/long-test9.txt")


# Functions ---------------------------------------------------------------

lweibullt<-function(x,u,a,b) length(x)*log(a)+length(x)*log(b)+(b-1)*sum(log(x))-a*sum(x^b-u^b)
qweibullt<-function(p,u,a,b) ((-log(1-p)/a)+u^b)^(1/b)

#Tailweibull fitting
quantile_tailweib <- function(sample_tail){
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
  return(qweibullt(p, u = 1, b = shape, a = scale))
}

#Binary Search
binary_search_h <- function(h_grid, x_grid, sample_tail){
  max_index <- length(h_grid)
  low_index <- 0
  index <- 0
  while(abs(max_index - low_index) >1){
    index <- floor((max_index-low_index)/2) +low_index
    if(logconctest_cpp(kde_cpp(sample_tail, x_grid, h_grid[index]), x_grid)==0){
      low_index <- index
    }
    else {
      max_index <- index
    }
  }
  if(max_index<low_index){index <- index+1}
  return(index)
}

# Test --------------------------------------------------------------------

for(m in 1:2){ 
  
  #Read data
  test0  <- data.table::fread(test_file[m],header = F,
                                     col.names = "X")
  
  #Compute empirical cumulative density function
  a <-ecdf(test0$X)
  
  #Get the value at quantule 1-10^-6 for latter comparison 
  emp_quantile <- as.numeric(quantile(test0$X, 1-10^-6))

  total_p_values <- list()
  total_quantiles <- list()
  tail_size <- c(50, 100, 200, 500)
  tail_5 <- c(3, 5, 10, 25)
  sample_list <- list()
  
  # Prepare tail sample -----------------------------------------------------
  
  for (k in 1:1000) {
    print(k)
    sample_data <- sample(test0$X, size = 1000, replace = T)
    sample_list[[k]] <- sample_data
    ecdf_test_sample <- ecdf(sample_data)
    p_values <- c()
    quantiles <- c()

    for(j in 1:length(tail_size)){
      threshold <- sort(sample_data)[1000-tail_size[j]]
      sample_tail <- sort(sample_data[sample_data > threshold]/threshold)
      stdev <- sd(sample_tail)
      quantiles[j] <-  quantile_tailweib(sample_tail)*threshold/emp_quantile
      
      #Grids for KDE 
      h_grid <- seq(0.0001, 0.1, length.out = 1000)
      x_grid <- seq(min(sample_tail), max(sample_tail), length.out = length(sample_tail))
      
      # Binary Search -----------------------------------------------------------
      
      #Bandwidth value to compare with other samples
      index <- binary_search_h(h_grid, sample_tail, sample_tail)
      h_0 <- h_grid[index]
      
      # Bootstrap and p-avlue ---------------------------------------------------
      
      sum <- 0
  
      #Perform again the test for 1000 bootstrap samples of tail
      # and store the value of h_star
      nsim <- 1000
      for(i in 1:nsim){
        s <- head(sort(sample(sample_tail, length(sample_tail), replace = T)),-tail_5[j])
        soroll <- h_0*rnorm(length(s))
        s_boot <-sort(sqrt(stdev^2/(stdev^2 + h_0^2))*(s +soroll)) 
        x_grid <- seq(min(s_boot), max(s_boot), length.out = length(s_boot))
        sum <- sum + logconctest_cpp(kde_cpp(s_boot, s_boot,h_0), s_boot)
      }
      
      #The p value is simply the proportion of h_star values above h_0
      p_values[j] <- 1-sum/nsim
    }
    
    total_p_values[[k]] <- p_values
    total_quantiles[[k]] <- quantiles
  }

  sim_data <- data.frame(p_values = unlist(total_p_values), quantiles = unlist(total_quantiles))
  sim_data$tail <- as.factor(rep(tail_size, 1000))
  sim_data$sim <- as.factor(rep(seq(1,1000), each=4))
  sim_data$test <- as.factor(rep(m-1,4000))
  sample_df <- data.frame(sample = unlist(sample_list))
  write.csv(sample_df, paste(paste("samples_", m+8, sep = ""),".csv", sep = ""), row.names = F)
  write.csv(sim_data, paste(paste("sim_data_t_", m+8, sep = ""),".csv", sep = ""), row.names = F)
}


# Plots -------------------------------------------------------------------

filt_sim_data %>% filter(tail == 200)


p1 <- ggplot(data = sim_data, aes( y = quantiles, group = tail, x = tail, fill = tail))+
  geom_boxplot(show.legend = F)+
  geom_hline(yintercept =  1, col = "red") +
  scale_fill_viridis(discrete = T)+
  theme_bw()+
  ylim(c(0.95,1.2))+
  facet_wrap(~test, nrow = 1)+
  ggtitle("Distribution of quantiles for log-concavity test with different tail sizes")

filt_sim_data <- sim_data %>% 
  filter(p_values <0.01)

p2 <- ggplot(data = filt_sim_data, aes( y = quantiles, group = tail, x = tail, fill = tail))+
  geom_boxplot(show.legend = F)+
  geom_hline(yintercept =  1, col = "red") +
  scale_fill_viridis(discrete = T)+
  theme_bw()+
  ylim(c(0.95,1.2))+
  facet_wrap(~test, nrow = 1)+
  ggtitle("Distribution of quantiles for log-concavity test with different tail sizes")

grid.arrange(p1, p2, nrow=1)


ggplot(data = sim_data, aes( y = p_values, group = tail, x = tail, fill = tail))+
  geom_boxplot()+
  geom_hline(yintercept =  0.05, col = "red") +
  scale_fill_viridis(discrete = T)+
  theme_bw()+
  ggtitle("Distribution of p-values for log-concavity test with different tail sizes")

ggplot(data = data.frame(sample_tail), aes(x = sample_tail))+
  geom_histogram(aes(y=..density..), colour = "black", fill ="#E69F00")+
  geom_line(aes(x = x_grid), y = kde_cpp(sample_tail, x_grid, 0.009), colour = "purple4")+
  ggtitle("Histogram for sample of the tail with KDE N=500")+
  theme_bw()





