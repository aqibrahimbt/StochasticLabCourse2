## Loading libraries and dependencies
options(warn=-1)
options(repr.plot.width=6, repr.plot.height=4)
packages = c("tidyverse", "EnvStats")

## Check to see if package is available and load else install the package and its dependencies
package.check <- lapply(packages, FUN = function(x) {
  if (!require(x, character.only = TRUE)) {
    install.packages(x, dependencies = TRUE)
    library(x, character.only = TRUE)
  }
})


## Inversion Method
# Computations
set.seed(250)

# Change default to Wichmann-Hill
RNGkind(kind = "Wichmann-Hill", normal.kind = NULL)

#Inversion method to simulate a binomial random variable
N = 1000
n = 10
p = 0.4
u <- runif(N)

bins <- .bincode(u, breaks = c(0, pbinom(0:10, 10, 0.4)), right = F, include.lowest = T)

## Simulate binomial random variables using the Inversion Method
binomial_rv_inversion <- numeric()
for(i in 1:N){
  binomial_rv_inversion[i] <- bins[i]-1
}


#Simulation of a binomial random variable by simulating corresponding Bernoulli random variables by inversion method
binomial_bernoulli_rv_inversion <- numeric()
for (i in 1:N){
  v <- runif(n)
  binomial_bernoulli_rv_inversion[i] <- sum(v < p)
}


#Simulation of a binomial random variable with the rbinom function
binomial_rv_rbinom <- rbinom(N, n, p)

## Visualization of Generated Samples
#Plot the histograms of all three samples on one panel
binomial_rv_inversion <- data.frame(binomial_rv_inversion)
binomial_rv_inversion$method <- rep("Binomial Random Vaariables", 1000)
colnames(binomial_rv_inversion) <- c("rand_num", "method")

binomial_bernoulli_rv_inversion <- data.frame(binomial_bernoulli_rv_inversion)
binomial_bernoulli_rv_inversion$method <- rep("Bernoulli Random Variables", 1000)
colnames(binomial_bernoulli_rv_inversion) <- c("rand_num", "method")

binomial_rv_rbinom <- data.frame(binomial_rv_rbinom)
binomial_rv_rbinom$method <- rep("Samples from 'rbinom'", 1000)
colnames(binomial_rv_rbinom) <- c("rand_num", "method")

samples <- rbind(binomial_rv_inversion, binomial_bernoulli_rv_inversion, binomial_rv_rbinom)

ggplot(samples, aes(x = rand_num, fill = method)) +
  geom_density() + labs(title = "Empirical PDF of Generated Samples",caption = "Simulated data using Inversion Methods and Wichmann-Hill generator",
                                                        x = "Generated Sample",y = "Density", colour = "Method") + theme_classic()

## Accept Reject Method
#Switch the random number generator back to its default
RNGkind(kind = "default", normal.kind = NULL)

## density of the normally distributed random variables
f <- function(x){
  ((2*pi)^(-1/2))*exp(-(x^2)/2)
}

## density of the standard Cauchy distribution
g <- function(x){
  (pi*(1 + x^2))^(-1)
}

# Determine the best value of the constant c, such that f(x) <= c * g(x)
x <- -2000:2000
(c <- max(f(x)/g(x)))

# Generating Cauchy distributed random variables using inversion method
N <- 10000
j <- 0

rand_num <- numeric()
while(length(rand_num) != N){
  w <- runif(1)
  cauchy <- tan((w-(1/2)) * pi)
  U <- runif(1)
  if(U * c * g(cauchy) <= f(cauchy)){
    rand_num[j] <- cauchy
    j <- j + 1
  }
}


## Histogram of the obtained sample with the standard normal density
k <- rnorm(N)
samples <- data.frame(rand_num, k)

ggplot(samples) +
  geom_histogram(aes( x = rand_num, y = ..density.., colour = rand_num), colour ="white") +
  geom_density(aes(x = k), colour = "blue") + 
  labs(title = "Empirical PDF of Generated Samples from Cauchy distribution",caption = "Simulated data using Accept-Reject Method and default generator",
       x = "Generated Sample",y = "Relative Frequency", colour = "Method") + theme_classic()

#QQ-plot
ggplot(data = samples, mapping = aes(sample = rand_num)) +
  stat_qq() + 
  labs(title = "QQ-Plot of Generated Samples",caption = "QQ-plot of the simulated standard normal distributed random variables ",
       x = "Theoretical",y = "Sample") + theme_classic()