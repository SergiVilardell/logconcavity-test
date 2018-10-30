library(logcondens)
library(data.table)
library(viridis)
library(tidyverse)

# Functions ---------------------------------------------------------------


lweibullt<-function(x,u,a,b) length(x)*log(a)+length(x)*log(b)+(b-1)*sum(log(x))-a*sum(x^b-u^b)
tweibull_mle <- function(x){
	aux<-function(p){-lweibullt(x,u=1,a=1, b=p)}
	solfit <- nlm(aux, 10)
	shape0 <- solfit$estimate
	
	aux<-function(p){-lweibullt(x,u=1,a=p, b=1)}
	solfit <- nlm(aux, 10)
	scale0 <- solfit$estimate
	
	aux<-function(p){-lweibullt(x,u=1,a=p[1], b=p[2])}
	solfit <- nlm(aux, c(scale0,shape0))
	scale <- solfit$estimate[1]
	shape <- solfit$estimate[2]
	return(lweibullt(x, u = 1, a = scale, b = shape))
}


# Data --------------------------------------------------------------------

test_file <- c("/home/bill/BSC/data/long-test0.txt",
			   "/home/bill/BSC/data/long-test1.txt",
			   "/home/bill/BSC/data/long-test2.txt",
			   "/home/bill/BSC/data/long-test3.txt",
			   "/home/bill/BSC/data/long-test4.txt",
			   "/home/bill/BSC/data/long-test5.txt",
			   "/home/bill/BSC/data/long-test6.txt",
			   "/home/bill/BSC/data/long-test7.txt",
			   "/home/bill/BSC/data/long-test8.txt", 
			   "/home/bill/BSC/data/long-test9.txt")

test_data <- data.table::fread(test_file[1],header = F,
							   col.names = "X")


tail_size <- c(50,100,200,500)
sample_data <- sort(sample(test_data$X, 1000, replace = T))
threshold <- sample_data[1000-tail_size[4]]
sample_tail <- sample_data[sample_data > threshold]/threshold

res <- logConDens(sample_tail, smoothed = F)
l <- tweibull_mle(sample_tail)
l_logcon <- res$L

f <- exp(res$phi)
hist(sample_tail, probability = T)
lines(res$x, f, col = "red")
