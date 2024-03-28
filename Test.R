library(tidyverse)
library(Rcpp)
library(RcppArmadillo)
library(foreach)
library(doParallel)
library(Deriv)
library(nimble)


path <- "/Users/kevin-imac/Desktop/Github - Repo/"
if(! file.exists(path)){
  path <- "/Users/kevinkvp/Desktop/Github Repo/"
}

dat <- read.table(paste0(path, "HW5MCMC/coal.dat"), header = TRUE)

sourceCpp(paste0(path, "HW5MCMC/src/main.cpp"))
test <- update_theta(lambda1 = 10, lambda2 = 5, dat = dat$disasters)
sum(test)
testBase <- rep(NA, 111)
for(i in 1:111){
  testBase[i] <- (1/1) ^ sum(dat$disasters[1:i])
}

testTheta <- rep(NA, 100000)
for(i in 1:100000){
  testTheta[i] <- update_theta(lambda1 = 1, lambda2 = 1, dat = dat$disasters)
}
hist(testTheta)
table(testTheta)

set.seed(10)
result <- gibbsGamma(iter = 100000, dat = dat$disasters) 
result <- gibbsGamma(iter = 100000, dat = dat$disasters) 

plot(result[, 1], type = "l")

result <- gibbsHalfN(iter = 50000, s2_1 = 1, s2_2 = 1, dat = dat$disasters)

for(i in 2:111){
  c(var(dat$disasters[1:i]), var(dat$disasters[-(1:i)])) %>% print()
}



plot(result[, 1], type = "l")
plot(result[, 2], type = "l")
plot(result[, 3], type = "l")

dat


update_theta(theta_old = 10, lambda1 = 0.074, lambda2 = 0.075, alpha = 10, dat = dat$disasters)


sampResult <- rep(NA, 100000)
for(i in 1:100000){
  sampResult[i] <- update_theta(lambda1 = 0.074, lambda2 = 0.075, alpha = 10, dat = dat$disasters)
}

table(sampResult)

testMat <- update_theta(lambda1 = 0.074, lambda2 = 0.075, alpha = 10, dat = dat$disasters)
testMat
sum(dat$disasters)
5^191
log(rgamma(1, 3, 1)^191)

20^140 * 40^51/sum(rep(20^140 * 40^51, 10))
