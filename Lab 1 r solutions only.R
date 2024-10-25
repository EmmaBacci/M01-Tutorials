############################
#  Get data and package    #
############################

rm(list=ls( ))
#remove old data in memory


######## Exercise 1 #############
#################################


library(lmtest)
library(sandwich)
library(mfx) # Marginal Effects, Odds Ratios and Incidence Rate Ratios for GLMs
library(censReg) # for Tobit
#load libraries


load("/Users/bacciemma/Documents/Unifern/Period 1/JobCorps_PC1.RData")
load("JobCorps_PC1.RData")
#load data

View(JC) #spreadsheet-style data viewer
str(JC)
#short description of the content of the dataset

attach(JC)
# every variable in loaded into work space
# Assuming df is your dataframe
set.seed(123)  # For reproducibility
noise <- rnorm(n = nrow(JC), mean = 0, sd = 0.1)  # Adjust sd for the amount of noise
JC <- JC + noise

#1
summary(assignment)

#2
ols1 <- lm(pworky3 ~ assignment)
coeftest(ols1, vcov = vcovHC)

plot(assignment, pworky3, pch=19)
abline(ols1$coefficients[1],ols1$coefficients[2])

#3.
ols2 <- lm(pworky3 ~ assignment + female + hsdegree + geddegree)
summary(ols2)


#4.
ols3 <- lm(earny3 ~ assignment + female + hsdegree + geddegree)
summary(ols3)


#5.
coeftest(ols3, vcov = vcovHC)

#6.
ols4 <- lm(earny3 ~ assignment + female) 
waldtest(ols4, ols3, vcov = vcovHC)

#7
worky3 <- (pworky3 > 0)
age2 <- age^2
femalehaschild <- female*haschild
ols5 <- lm(worky3 ~ assignment + haschild + femalehaschild + female + age + age2)
coeftest(ols5, vcov = vcovHC)

#8.
print(ols5$coefficients["age"] + ols5$coefficients["age2"]*2*20) #marginal effect of age


#9.
print(ols5$coefficients["female"] + ols5$coefficients["haschild"] + ols5$coefficients["femalehaschild"])


#10.
hist(ols5$fitted.values, col = "blue") #histogram
summary(ols5$fitted.values)


#11.
probit1 <- glm(worky3 ~ assignment + haschild + femalehaschild + female + age + age2, family=binomial(probit))
summary(probit1)


#12.
probiteffect <- pnorm(cbind(1,assignment,1,1,1,age,age2)%*%probit1$coefficients) - pnorm(cbind(1,assignment,0,0,0,age,age2)%*%probit1$coefficients)
mean(probiteffect)


#13.
probitmfx(probit1, data = JC, atmean = TRUE) # average marginal ("treatment") effects
probitmfx(probit1, data = JC, atmean = FALSE)


#15.
model <- earny3 ~ assignment + haschild + femalehaschild + female + age + age2
tobit <- censReg(model,data=JC)
summary(tobit)

ols6 <- lm(model)
summary(ols6)
coeftest(ols6, vcov = vcovHC)


#16.
xb <- cbind(1, assignment, haschild, femalehaschild, female, age, age2)%*%tobit$estimate[1:7]
sigma <- exp(tobit$estimate[8])
predictions <- (pnorm((xb)/sigma))*xb+sigma*(dnorm((xb)/sigma))

summary(fitted(ols6))

# create a new plotting window and set the plotting area into a 1*2 array
par(mfrow = c(1, 2))
hist(predictions, col = "red", main = "Tobit", xlab = "predictions") #plot histograms
hist(fitted(ols6), col = "blue", main = "LPM", xlab = "predictions")


#17.
xb1 <- cbind(1,assignment, 1, 1, 1, age, age2)%*%tobit$estimate[1:7]
xb2 <- cbind(1,assignment, 0, 0, 0, age, age2)%*%tobit$estimate[1:7]
sigma <- exp(tobit$estimate[8])
predictions1 <- (pnorm( (xb1) / sigma))*xb1+sigma*(dnorm( (xb1) / sigma))
predictions2 <- (pnorm( (xb2) / sigma))*xb2+sigma*(dnorm( (xb2) / sigma))
tobiteffect <- mean(predictions1-predictions2)
print(tobiteffect)
print(ols6$coefficients["female"] + ols6$coefficients["haschild"] + ols6$coefficients["femalehaschild"])





######## Exercise 2 #############
#################################

#Remove previous data from the memory
rm(list=ls())


#load libraries
library(quantreg)
library(lmtest)
library(sandwich)

#Simple example for a quantile regression
x <- seq(0,10, by = 0.01); n <- length(x)
beta_0 <- 1; beta_1 <- 0.5
y <- beta_0 + beta_1*x + rnorm(n) #homoscedastic
y <- beta_0 + beta_1*x + 0.1*x*rnorm(n) #heteroscedastic
plot(x, y, pch = 3)
abline(beta_0, beta_1, col="red", lty = 2)

model <- (y ~ x)
ols <- lm(model)
abline(ols$coefficients[1], ols$coefficients[2], col="red", lty = 1)
summary(ols)

qr0.5 <- rq(model, tau = 0.5)
qr0.2 <- rq(model, tau = 0.2)
qr0.8 <- rq(model, tau = 0.8)

abline(qr0.5$coefficients[1], qr0.5$coefficients[2], col="blue")
abline(qr0.2$coefficients[1], qr0.2$coefficients[2], col="blue")
abline(qr0.8$coefficients[1], qr0.8$coefficients[2], col="blue")




load("growth_PC1.RData")
attach(growth)
str(growth)
View(growth)


#1.
model <- (y_net ~ lgdp2 + fse2 + fhe2 
          + mse2 + mhe2 + lexp2 + lintr2 + 
            gedy2 + Iy2 + gcony2 + lblakp2 + pol2 + ttrad2)
ols <- lm(model)
summary(ols)
coeftest(ols, vcov = vcovHC) 

#3.
growth <- rq(model, tau = 0.5)
summary.rq(growth, se = "boot", R = 1000)
#R = 1000 number of bootstrap replications

#4
model2 <- (y_net ~ lgdp2 + mse2 + mhe2 + lblakp2)
growth0.5 <- rq(model2, tau = 0.5)
summary.rq(growth0.5, se = "boot", R = 1000)

growth0.1<-rq(model2, tau = 0.1)
summary.rq(growth0.1, se = "boot", R = 1000)

#5.
plot(summary(rq(model2,tau = 1:9/10)))

#8.
n <- 99
x <- rnorm(n)
x <- rcauchy(n, location = 0, scale = 1) #Cauchy distribution (heavy-tailed)
hist(x)

mean(x)
median(x)
sort(x)[50]

