## Loading library and dependencies
options(warn=-1)
options(repr.plot.width=6, repr.plot.height=4)

## Normality and Poisonness
student_mat <- read.csv('data/student-mat.csv')

group1 <- data.frame(student_mat$G1, rep('G1', nrow(student_mat)))
colnames(group1) <- c("grades", "types")
group2 <- data.frame(student_mat$G2, rep('G2', nrow(student_mat)))
colnames(group2) <- c("grades", "types")
group3 <- data.frame(student_mat$G3, rep('G3', nrow(student_mat)))
colnames(group3) <- c("grades", "types")
all_groups <- rbind(group1, group2, group3)

#Normal distributed?
ggplot(data = all_groups, mapping = aes(sample = grades)) + 
  geom_density(aes(x = grades), fill = "chartreuse") +
  facet_wrap(. ~types) + 
  labs(title = "Density plot of student grades",caption = "Grades data of students in math and portuguese language courses in secondary school"
       ,x = "Grades",y = "Density") + theme_classic()


ggplot(data = all_groups, mapping = aes(sample = grades)) + 
  stat_qq(distribution = stats::qnorm, dparams = list(mean = mean(all_groups$grades), sd = sd(all_groups$grades))) +
  geom_abline(alpha = 0.25) +
  facet_wrap(. ~types) + 
  labs(title = "QQ-Plot plot of student grades",caption = "Grades data of students in math and portuguese language courses in secondary school"
       ,x = "Theoretical",y = "Sample") + theme_classic()

#Poisson distributed?
ggplot(data = all_groups, mapping = aes(sample = grades)) + 
  stat_qq(distribution = stats::qpois, dparams = list(lambda = mean(all_groups$grades))) +
  geom_abline(alpha = 0.25) +
  facet_wrap(. ~types) + labs(title = "QQ-Plot plot of student grades",caption = "Grades data of students in math and portuguese language courses in secondary school"
                              ,x = "Theoretical",y = "Sample") + theme_classic()

par(mfrow=c(1,1))
n=length(all_groups$grades)
(x=table(all_groups$grades))
k=as.numeric(names(x))
plot(k,log(x)+lfactorial(k))


## generalized linear model
model_1 <- glm(formula = G1 ~. -G2 -G3, family = poisson, data = student_mat) 
summary(model_1)

#Pearson residuals for model_1
pearson_residual <- residuals(model_1, "pearson")
pearson_residual <- data.frame(pearson_residual)

## QQ-Plot visualization
ggplot(data = pearson_residual, mapping = aes(sample = pearson_residual)) + 
  stat_qq(distribution = stats::qnorm, dparams = list(mean = mean(pearson_residual$pearson_residual),
                                                      sd = sd(pearson_residual$pearson_residual))) +
  geom_abline(alpha = 0.25) + labs(title = "QQ-Plot plot of Pearson Residuals"
                                   ,x = "Theoretical",y = "Sample") + theme_classic()


## Function to compute the Anscombe Residuals
anscombe <- function(y, mu){
  (3 *(y ** (2/3)- mu ** (2/3)))/2*(mu ** (1/6))
}


#Anscombe residuals for model_1
anscombe_residual <- anscombe(student_mat$G1, model_1$fitted.values)
anscombe_residual <- data.frame(anscombe_residual)

## QQ-Plot visualization
ggplot(data = anscombe_residual, mapping = aes(sample = anscombe_residual)) + 
  stat_qq(distribution = stats::qnorm, dparams = list(mean = mean(anscombe_residual$anscombe_residual), 
                                                      sd = sd(anscombe_residual$anscombe_residual))) +
  geom_abline(alpha = 0.25) + labs(title = "QQ-Plot plot of Anscombe Residuals"
                                   ,x = "Theoretical",y = "Sample") + theme_classic()

# Residual analysis for model 1
par(mfrow = c(2, 2))
plot(model_1)

model_2 <- glm(formula = G1 ~ sex + Fedu + studytime + failures + schoolsup + famsup + goout -G2 -G3, 
               family = poisson, data = student_mat) 
summary(model_2)


#Analysis of deviance
anova(model_2, model_1, test = "Chisq")


#model 3 and comparison with model 2
model_3 <- glm(formula = G1 ~ sex + Fedu + studytime + failures + schoolsup + famsup + Walc , family = poisson, data = student_mat) 
summary(model_3)

# Residual analysis for model 1
par(mfrow = c(2, 2))
plot(model_3)