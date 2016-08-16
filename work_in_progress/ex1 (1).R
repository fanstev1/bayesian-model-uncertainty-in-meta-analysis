library(metafor)
library(coda)
library(R2WinBUGS)
#library(R2OpenBUGS)
setwd("I:/Users/sfan/Google Drive/MSc HSR thesis work/program/bayesian meta analysis")
source("gibbs sampler for bayesian model uncertainty in meta analysis.R")

# DATA 1: EXAMPLE 1 FROM COMPREHENSIVE DECISION ANALYTICAL MODELLING IN ECONOMIC EVALUATION
meta.data<- data.frame(study= 1:6,
                       studyName= c("NAIA3005", "NAI30010", "NAIA/B2009", "WV15673", "WV15697", "WV15799"), 
                       nc= c(554, 423, 144, 268, 251, 462), 
                       rc= c(34, 40, 9, 19, 6, 34),
                       nt= c(553, 414, 144, 268, 252, 493),
                       rt= c(11, 7, 3, 3, 3, 4))

attach(meta.data)

eff.summary<- escalc(measure= "OR",
                     ai= rt,
                     n1i= nt,
                     ci= rc,
                     n2i= nc, 
                     data= meta.data)

eff.summary<- cbind(eff.summary, logor.ctl=log(rc/(nc-rc))) # adding logit(event rate) in the control group

detach(meta.data)

save(list = ls(all.names = TRUE), file = "ex1.RData")


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
                                                      #nchain= 3, 
                                                      mcmc.out.lst= list(mcmc.out= NULL, mcmc.out= NULL, mcmc.out= NULL)
) 
save(list = ls(all.names = TRUE), file = "ex1.RData")
     
summary(mcmc.out.lst)
plot(mcmc.out.lst)
gelman.diag(mcmc.out.lst) 
HPDinterval(mcmc.out.lst) 
