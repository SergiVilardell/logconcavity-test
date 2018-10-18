# Generate log-concave data -----------------------------------------------

data <- dweibull(seq(0,1, length.out = 10000), 3, 1)

plot(data)
first_derivative <- diff(data)
log_derivative <- diff(diff(log(data)))
p <- first_derivative/data[1:length(data)-1]


plot(p[-c(1:100)], main = "g'(x)/g(x) is monotone decreasing")
plot(log_derivative[-c(1:100)],  main = "second derivative of ln(g(x)) is negative")
abline(h = 0)
