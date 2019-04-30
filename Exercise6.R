#specify the packages of interest
options(warn=-1)
options(repr.plot.width=6, repr.plot.height=4)

packages = c("readr", "ggplot2", "tidyverse")

## Check to see if package is available and load else install the package and its dependencies
package.check <- lapply(packages, FUN = function(x) {
  if (!require(x, character.only = TRUE)) {
    install.packages(x, dependencies = TRUE)
    library(x, character.only = TRUE)
  }
})

StudentsPerformance <- read_csv("data/StudentsPerformance.csv")

math_score <-StudentsPerformance$`math score`
reading_score <-StudentsPerformance$`reading score`
writing_score <-StudentsPerformance$`writing score`

## Kernel functions
uniform <- function(x) (abs(x) < 1) * 0.5
triangular <- function(x) (abs(x) < 1) * (1 - abs(x))
gaussian <- function(x) 1/sqrt(2*pi) * exp(-(x^2)/2)
epanechnikov <- function(x) (abs(x) < 1) * (0.75 * (1 - x^2))


k.names <- list("uniform", "triangular", "epanechnikov", "gaussian")
k.funtions <- list(uniform, triangular, epanechnikov, gaussian)

#KDE = Kernel density estimator
KDE <- function(sample, bandwith, kernel){
  K <- k.funtions[[which(k.names == kernel)]]
  a <- function(x){
    t <- (x - sample)/bandwith
    (1/(bandwith*length(sample)))*(sum(K(t)))
  }
  return(Vectorize(a))
}

ggplot(StudentsPerformance, aes(x = `math score`)) + 
  stat_function(fun = KDE(math_score, 2, "epanechnikov"), aes(colour = "2")) +
  stat_function(fun = KDE(math_score, 5, "epanechnikov"), aes(colour = "5")) +
  stat_function(fun = KDE(math_score, 8, "epanechnikov"), aes(colour = "8")) +
  stat_function(fun = KDE(math_score, 11, "epanechnikov"), aes(colour = "11")) +
  theme(legend.position = c(0.96, 0.95),legend.justification = c("right", "top"), legend.title = element_text()) +
  labs(color = "Bandwidth") + labs(title = "Epanechnikov Kernel function with different bandwidths",
                                   x = "Math Score",y = "Frequency", colour = "Bandwith") + theme_classic()


#plot: bandwidth = 11, with different kernel
ggplot(StudentsPerformance, aes(x = `math score`)) + 
  stat_function(fun = KDE(math_score, 11, "epanechnikov"), aes(colour = "epanechnikov")) +
  stat_function(fun = KDE(math_score,11, "uniform"), aes(colour = "uniform")) +
  stat_function(fun = KDE(math_score, 11, "triangular"), aes(colour = "triangular")) +
  stat_function(fun = KDE(math_score, 11, "gaussian"), aes(colour = "gaussian")) +
  theme(legend.position = c(0.96, 0.95),legend.justification = c("right", "top"), legend.title = element_text()) +
  labs(color = "Kernel") + labs(title = "Kernel Density plot with different kernel functions",
                                x = "Math Score",y = "Relative Frequency", colour = "Kernel") + theme_classic()


#Implement the cross-validation criterion to find the optimal bandwidth

#cross-validation function from the lecture 6
CV = function(sample, bandwidth){
  n = length(sample)
  
  f.hat = KDE(sample, bandwidth, "gaussian")
  
  f.hat.sqr = function(x){f.hat(x)**2} 
  
  A = integrate(f.hat.sqr, lower = Inf, upper = Inf) 
  
  B = (sum(gaussian(outer(sample, sample, FUN = "-")/bandwidth)) - n * gaussian (0)) * 2/(n *(n - 1)* bandwidth)
  
  C = A$value - B
  
  return(C)
}

CV.math = Vectorize(function(bandwidth){CV(math_score, bandwidth)})
CV.reading = Vectorize(function(bandwidth){CV(reading_score, bandwidth)})
CV.writing = Vectorize(function(bandwidth){CV(writing_score, bandwidth)})

# # # Optimize objective function
CV.math.opt = optimize(CV.math, interval = c(1, 20))$minimum
CV.reading.opt = optimize(CV.reading, interval = c(1,15))$minimum
CV.writing.opt = optimize(CV.writing, interval = c(1, 20))$minimum


#bw.ucv
bw.ucv.math.opt = bw.ucv(math_score, lower = 1, upper = 20)
bw.ucv.reading.opt = bw.ucv(reading_score, lower = 1, upper = 15)
bw.ucv.writing.opt = bw.ucv(writing_score, lower = 1, upper = 20)
c(bw.ucv.math.opt, bw.ucv.reading.opt, bw.ucv.writing.opt)


#bw.bcv
bw.bcv.math.opt = bw.bcv(math_score, lower = 1, upper = 20)
bw.bcv.reading.opt = bw.bcv(reading_score, lower = 1, upper = 15)
bw.bcv.writing.opt = bw.bcv(writing_score, lower = 1, upper = 20)
c(bw.bcv.math.opt, bw.bcv.reading.opt, bw.bcv.writing.opt)


## Completed preparatory Course
students_participated <- StudentsPerformance %>%
  filter(`test preparation course` == "completed")

math_score_c = students_participated$`math score`
reading_score_c = students_participated$`reading score`
writing_score_c = students_participated$`writing score`

CV.math_score_c = Vectorize(function(bandwidth){CV(math_score_c, bandwidth)})
CV.reading_score_c = Vectorize(function(bandwidth){CV(reading_score_c, bandwidth)})
CV.writing_score_c = Vectorize(function(bandwidth){CV(writing_score_c, bandwidth)})

#optimal bandwidth for the group "completed"
CV.math.opt.c = optimize(CV.math_score_c, interval = c(1, 20))$minimum
CV.reading.opt.c = optimize(CV.reading_score_c, interval = c(1, 15))$minimum
CV.writing.opt.c = optimize(CV.writing_score_c, interval = c(1, 20))$minimum
# #KDEs

math_score_kde_c <- KDE(students_participated$`math score`, CV.math.opt.c, "gaussian") 
reading_score_kde_c <- KDE(students_participated$`reading score`, CV.reading.opt.c, "gaussian")
writing_score_kde_c <- KDE(students_participated$`writing score`, CV.writing.opt.c, "gaussian")

## Didn't complete preparasatory course
students_notparticipated <- StudentsPerformance %>%
  filter(`test preparation course` == "none")

math_score_p = students_notparticipated$`math score`
reading_score_p = students_notparticipated$`reading score`
writing_score_p = students_notparticipated$`writing score`

CV.math_score_p = Vectorize(function(bandwidth){CV(math_score_p, bandwidth)})
CV.reading_score_p = Vectorize(function(bandwidth){CV(reading_score_p, bandwidth)})
CV.writing_score_p = Vectorize(function(bandwidth){CV(writing_score_p, bandwidth)})

#optimal bandwidth for the group "not participated"
CV.math.opt.p = optimize(CV.math_score_p, interval = c(1, 20))$minimum
CV.reading.opt.p = optimize(CV.reading_score_p, interval = c(1, 15))$minimum
CV.writing.opt.p = optimize(CV.writing_score_p, interval = c(1, 20))$minimum

# #KDEs
math_score_kde_p <- KDE(students_notparticipated$`math score`, CV.math.opt.p, "gaussian") 
reading_score_kde_p <- KDE(students_notparticipated$`reading score`, CV.reading.opt.p, "gaussian")
writing_score_kde_p <- KDE(students_notparticipated$`writing score`, CV.writing.opt.p, "gaussian")


#plots
ggplot(StudentsPerformance, aes(x = `math score`)) + 
  stat_function(fun = math_score_kde_c, aes(color = "completed")) +
  stat_function(fun = math_score_kde_p, aes(color = "none")) +
  scale_colour_manual(name="Test preparation course", values = c("red", "purple")) +
  theme(legend.position = c(0.96, 0.95),legend.justification = c("right", "top"), legend.title = element_text()) + labs(title = "",
                                                                                                                        x = "Math Score",y = "y", colour = "Preparatory Course") + theme_classic()
ggplot(StudentsPerformance, aes(x = `reading score`)) + 
  stat_function(fun = reading_score_kde_c, aes(color = "completed")) +
  stat_function(fun = reading_score_kde_p, aes(color = "none")) +
  scale_colour_manual(name="Test preparation course", values = c("red", "purple")) +
  theme(legend.position = c(0.96, 0.95),legend.justification = c("right", "top"), legend.title = element_text()) + labs(title = "",
                                                                                                                        x = "Reading Score",y = "y", colour = "Preparatory Course") + theme_classic()

ggplot(StudentsPerformance, aes(x = `writing score`)) + 
  stat_function(fun = writing_score_kde_c, aes(color = "completed")) +
  stat_function(fun = writing_score_kde_p, aes(color = "none")) +
  scale_colour_manual(name="Test preparation course", values = c("red", "purple")) +
  theme(legend.position = c(0.96, 0.95),legend.justification = c("right", "top"), legend.title = element_text()) + labs(title = "",
                                                                                                                        x = "Writing Score",y = "y", colour = "Preparatory Course") + theme_classic()


