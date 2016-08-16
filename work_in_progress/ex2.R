library(metafor)
library(coda)
library(R2WinBUGS)
#library(R2OpenBUGS)
setwd("I:/Users/sfan/Google Drive/MSc HSR thesis work/program/bayesian meta analysis")
source("gibbs sampler for bayesian model uncertainty in meta analysis.R")

# DATA 2: SYSTEMATIC REVIEWS IN HEALTH CARE CASE STUDY 1 pg 303
meta.data<- data.frame(study= 1:6,
                       studyName= c("Breart 1992", "Breart 1992", "Gagnon 1997", "Hodnett 1989", "Kennell 1991", "Langer 1998"), 
                       nc= c(131, 664, 204, 73, 200, 363), 
                       rc= c(62, 319, 142, 43, 55, 303),
                       nt= c(133, 656, 209, 72, 212, 361),
                       rt= c(55, 281, 139, 30, 24, 205))

attach(meta.data)

eff.summary<- escalc(measure= "OR",
                     ai= rt,
                     n1i= nt,
                     ci= rc,
                     n2i= nc, 
                     data= meta.data)

eff.summary<- cbind(eff.summary, logor.ctl=log(rc/(nc-rc))) # adding logit(event rate) in the control group

detach(meta.data)

save(list = ls(all.names = TRUE), file = "ex2.RData")


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
save(list = ls(all.names = TRUE), file = "ex2.RData")

summary(mcmc.out.lst)
plot(mcmc.out.lst)
gelman.diag(mcmc.out.lst) 
HPDinterval(mcmc.out.lst)
