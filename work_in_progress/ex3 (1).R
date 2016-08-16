library(metafor)
library(coda)
#library(R2WinBUGS)
#library(R2OpenBUGS)
setwd("I:/Users/sfan/Google Drive/MSc HSR thesis work/program/bayesian meta analysis")
source("gibbs sampler for bayesian model uncertainty in meta analysis.R")

# DATA 3: Meta analysis EXAMPLE from Bayesian Data Analysis by Andrew Gelman
meta.data<- data.frame(study= 1:22,
                       studyName= paste(1:22), 
                       nc= c(39, 116, 93, 1520, 365, 52, 939, 471, 282, 1921, 583, 266, 293, 883, 147, 213, 122, 154, 134, 218, 364, 674), 
                       rc= c(3, 14, 11, 127, 27, 6, 152, 48, 37, 188, 52, 47, 16, 45, 31, 38, 12, 6, 3, 40, 43, 39),
                       nt= c(38, 114, 69, 1533, 355, 59, 945, 632, 278, 1916, 873, 263, 291, 858, 154, 207, 251, 151, 174, 209, 391, 680),
                       rt= c(3, 7, 5, 102, 28, 4, 98, 60, 25, 138, 64, 45, 9, 57, 25, 33, 28, 8, 6, 32, 27, 22))

attach(meta.data)

eff.summary<- escalc(measure= "OR",
                     ai= rt,
                     n1i= nt,
                     ci= rc,
                     n2i= nc, 
                     data= meta.data)

eff.summary<- cbind(eff.summary, logor.ctl=log(rc/(nc-rc))) # adding logit(event rate) in the control group

detach(meta.data)
save.image("ex3.RData")


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
               fixed.effect= (lambda==0),
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
               fixed.effect= (lambda==0),
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
               fixed.effect= (lambda==0),
               # model checking
               dev= rep(NA, 1), 
               p.inv= rep(NA, nobs),
               y.rep= rep(NA, nobs),
               p.smaller= rep(NA, nobs),
               dev.rep= rep(NA)
)


parms.ini<- list(parms.1, parms.2, parms.3)


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
                                                      nburn= 5000, 
                                                      niter= 10000, 
                                                      print.progress= TRUE, 
                                                      mcmc.out.lst= list(mcmc.out= NULL, mcmc.out= NULL, mcmc.out= NULL)
) 
save(list = ls(all.names = TRUE), file = "ex3.RData")

summary(mcmc.out.lst[, c("beta", "or", "sd.random", "fixed.effect")])
plot(mcmc.out.lst[, c("beta", "or", "sd.random", "fixed.effect")])
gelman.diag(mcmc.out.lst[, c("beta", "or", "sd.random", "fixed.effect")]) 
HPDinterval(mcmc.out.lst[, c("beta", "or", "sd.random", "fixed.effect")])

