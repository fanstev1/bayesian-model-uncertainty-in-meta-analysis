library(MCMCpack)
library(mvtnorm)
library(truncnorm)
library(metafor)
#setwd("C:/Users/user/Downloads/")
setwd("Z:/Steve tmp/R sessions simulation")
load("sim_reduced.RData")
load("sim_reduced_2.RData")

#for (i in 599999:400000) {
#for (i in 400001:547000) {
#for (i in 435001:547000) {
543000
for (i in 543001:547000) {
  y    <- yi.lst[[i]]
  var_y<- vi.lst[[i]]
  nobs <- length(y)
  
  beta<- rnorm(1, mean= sum(y*(1/var_y)/sum(1/var_y)), sd= sqrt(sum(1/var_y)))
  lambda<- rnorm(1)
  b<- rnorm(nobs)
  var_b<- rinvgamma(1, shape= .5, scale= .5*nobs)
  
  parms.1<- list(beta= beta,
                 lambda= lambda,
                 b= b,
                 var_b= var_b,
                 ## additional parameters to keep track 
                 s= lambda * b,
                 sd.random= abs(lambda) * sqrt(var_b),
                 or= exp(beta),
                 fixed.effect= (lambda==0)
  )
  parms.ini<- list(parms.1)
  tmp<- bayesian.mdl.uncertainty.meta.analysis(y= y, 
                                              var_y= var_y, 
                                              parms.ini= parms.ini, 
                                              mu_beta= 0, 
                                              var_beta= 100, 
                                              alpha_0= .5, 
                                              theta_0= .5 * nobs, 
                                              var_lambda_0= var.prior[i], 
                                              prob.random= .5,
                                              nthin= 10, 
                                              nburn= 5000, 
                                              niter= 10000
  ) 
  mcmc.res.lst[[i]]<- summary(tmp)
  HPDinterval.lst[[i]]<- HPDinterval(tmp)
  
  print(i)
  
  if (i %% 1000 == 0) {
    rm(list= c("tmp", "beta", "lambda", "b", "var_b", "parms.1", "parms.ini", "y", "var_y", "nobs"))
    save(list = c("i", "mcmc.res.lst", "HPDinterval.lst"), file = "sim_reduced_2.RData")
  }
}

