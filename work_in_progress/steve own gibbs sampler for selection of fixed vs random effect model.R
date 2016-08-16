library(MCMCpack)
library(mvtnorm)
library(truncnorm)
setwd("I:/Users/sfan/Google Drive/MSc HSR thesis work/program/bayesian meta analysis")
source("gibbs sampler for bayesian model uncertainty in meta analysis.R")

# save.image("steve work in progress.RData")


## model Y_i = mu + lambda * b_i + epislon_i
## prior density: 
##    beta ~ N(mu_beta, var_beta), interpreted as overall
##    lambda ~ ZI-N(0, var_lambda_0, prob.random)
##    b ~ MVN(mu_b_0, var_b), var_b ~ inverse gamma shape (alpha) and scale (theta)

## data
y<- as.numeric(eff.summary$yi)
var_y<- as.numeric(eff.summary$vi)  
nobs<- length(y)

# initial 
set.seed(45)
beta<- rnorm(1, mean= sum(y*(1/var_y)/sum(1/var_y)), sd= sqrt(sum(1/var_y)))
lambda<- rnorm(1)
b<- rnorm(nobs)
var_b<- rinvgamma(1, shape= .5, scale= .5)

parms.1<- list(beta= beta,
             lambda= lambda,
             b= b,
             var_b= var_b,
             ## additional parameters to keep track 
             sd.random= abs(lambda) * sqrt(var_b),
             or= exp(beta),
             random.absence= (lambda==0),
             # model checking
             dev= rep(NA, 1), 
             p.inv= rep(NA, nobs),
             y.rep= rep(NA, nobs),
             p.smaller= rep(NA, nobs),
             dev.rep= rep(NA)
             )


set.seed(546)
beta<- rnorm(1, mean= sum(y*(1/var_y)/sum(1/var_y)), sd= sqrt(sum(1/var_y)))
lambda<- rnorm(1)
b<- rnorm(nobs)
var_b<- rinvgamma(1, shape= .5, scale= .5)

parms.2<- list(beta= beta,
               lambda= lambda,
               b= b,
               var_b= var_b,
               ## additional parameters to keep track 
               sd.random= abs(lambda) * sqrt(var_b),
               or= exp(beta),
               random.absence= (lambda==0),
               # model checking
               dev= rep(NA, 1), 
               p.inv= rep(NA, nobs),
               y.rep= rep(NA, nobs),
               p.smaller= rep(NA, nobs),
               dev.rep= rep(NA)
)



set.seed(132)
beta<- rnorm(1, mean= sum(y*(1/var_y)/sum(1/var_y)), sd= sqrt(sum(1/var_y)))
lambda<- rnorm(1)
b<- rnorm(nobs)
var_b<- rinvgamma(1, shape= .5, scale= .5)

parms.3<- list(beta= beta,
               lambda= lambda,
               b= b,
               var_b= var_b,
               ## additional parameters to keep track 
               sd.random= abs(lambda) * sqrt(var_b),
               or= exp(beta),
               random.absence= (lambda==0),
               # model checking
               dev= rep(NA, 1), 
               p.inv= rep(NA, nobs),
               y.rep= rep(NA, nobs),
               p.smaller= rep(NA, nobs),
               dev.rep= rep(NA)
)


parms.ini<- list(parms.1, parms.2, parms.3)

#nthin<- 10
#nburn<- 5000
#niter<- 5000
#nchain<- length(parms.ini)

#mcmc.out<- NULL
#mcmc.out.lst<- list(mcmc.out= NULL, mcmc.out= NULL, mcmc.out= NULL)
mcmc.out.lst<- bayesian.mdl.uncertainty.meta.analysis(y= as.numeric(eff.summary$yi), 
                                                     var_y= as.numeric(eff.summary$vi), 
                                                     parms.ini= parms.ini, 
                                                     mu_beta= 0, 
                                                     var_beta= 10000, 
                                                     alpha_0= .5, 
                                                     theta_0= .5 * nobs, 
                                                     var_lambda_0= 25, 
                                                     prob.random= .5,
                                                     nthin= 10, 
                                                     nburn= 500, 
                                                     niter= 500, 
                                                     #nchain= 3, 
                                                     mcmc.out.lst= list(mcmc.out= NULL, mcmc.out= NULL, mcmc.out= NULL)
                                                     ) 

#mcmc.out<- mcmc(mcmc.out)
summary(mcmc.out.lst)
plot(mcmc.out)


psedoBF.meta<- function(mcmc.out, pos= grep("p.inv", colnames(mcmc.out))){
  cpo<- 1/apply(mcmc.out[, pos], 2, mean) # conditional predictive ordinate
  return(-2* sum(log(predictive.ordinate)) )
}


