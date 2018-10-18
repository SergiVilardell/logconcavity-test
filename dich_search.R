a <- rep(0,50)
b <- rep(1, 5550)

data <- c(a,b)
max_index <- length(data)
low_index <- 0
index <- 0
while(abs(max_index - low_index) >1){
  index <- floor((max_index-low_index)/2) +low_index
  if(data[index]==0){
    low_index <- index
  }
  else {
    max_index <- index
  }
}

print(index)
