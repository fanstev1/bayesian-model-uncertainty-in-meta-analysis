# objective: conduct a simulation to examine the estimated posterior model probability using unit-information prior for the parameter lambda.

# the number of studies= 3, 5, 10, 15, 20
# within study sigma^2 ~ inverse gamma(alpha= 2, theta= 1)
# condition on sigma^2, beta ~ N(.5, sigma^2)
library(MCMCpack)
library(mvtnorm)
library(truncnorm)
library(metafor)
#setwd("I:/Users/sfan/Google Drive/MSc HSR thesis work/program/bayesian meta analysis")
#setwd("I:/Users/sfan/Google Drive/MSc HSR thesis work/program/R sessions simulation")
setwd("T:/Steve tmp/R sessions simulation")
source("gibbs sampler for bayesian model uncertainty in meta analysis.R")
source("functions to calculate unit information.R")


#i= 3
#j= 2
#k= 2


nstudy<- c(3, 5, 10, 15, 20)
var_s<- c(# 0, no between study hetergeneity
          .2, .5, 1, 2)
random.eff.prob<- c(0, .25, .5, 1)
nsim<- 10000

# generate data
sim.data<- NULL
for (k in 1:length(random.eff.prob)) {
  for (j in 1:length(var_s)) {
    for (i in 1:length(nstudy)) {
      
      p    <- rbinom(n= nstudy[i]*nsim, size= 1, prob= random.eff.prob[k])
      var_y<- rinvgamma(n= nstudy[i]*nsim, shape= 2, scale= 1)
      s    <- rnorm(n= nstudy[i]*nsim, sd= sqrt(var_s[j]), mean= 0)
      y    <- rnorm(n= nstudy[i]*nsim, sd= sqrt(var_y)   , mean= 1 + p*s)
      
      tmp<- data.frame(nstudy= nstudy[i],
                       var.random= var_s[j],
                       rdm.prob= random.eff.prob[k],
                       sim.idx= rep(1:nsim, each= nstudy[i]),
                       yi= y,
                       vi= var_y,
                       true.s= s,
                       p= p
      )
      sim.data<- rbind(sim.data, tmp)
      rm(list= c("y", "s", "var_y", "tmp"))
      #save(list = ls(all.names = TRUE), file = "sim.RData")
      print(c(random.eff.prob[k], var_s[j], nstudy[i]))
  
    }
  }
}
save(list = ls(all.names = TRUE), file = "sim.RData")

# freqentist 

dd<- unique(sim.data[, c("nstudy", "var.random", "rdm.prob", "sim.idx")])
dd<- cbind(sid= 1:nrow(dd), dd)
dd<- merge(x= dd, y= sim.data, by= c("nstudy", "var.random", "rdm.prob", "sim.idx"))
sim.data<- dd[order(dd$sid, dd$nstudy, dd$var.random, dd$rdm.prob, dd$sim.idx),]
save(list = ls(all.names = TRUE), file = "sim.RData")

res.re<- with(sim.data, 
              by(sim.data, INDICES = list(sid),
                 function(tmp) rma(yi= tmp$yi, vi= tmp$vi, method= "REML", control= list(maxiter= 10000, threshold= 1E-7, stepadj=0.5)) )
)
res.re.bk<- res.re
#dd<- sapply(res.re.bk, function(dd) dd$tau2)
#idx.zero<- which(dd==0)
save(list = ls(all.names = TRUE), file = "sim.RData")

# theoretically equal to zero
llk<- sapply(res.re.bk, function(rma.obj) llk.prime.var_s(y= rma.obj$yi, var_y= rma.obj$vi, var_s= rma.obj$tau2))
idx.zero<- which(abs(llk)> 0.005)
for (i in idx.zero) {

  sign.zero<- sign(llk.prime.var_s(y= res.re[[i]]$yi, var_y= res.re[[i]]$vi, var_s= res.re[[i]]$tau2))
  ll<- res.re[[i]]$tau2 - 0.001
  repeat {
    sign.ll<- sign(llk.prime.var_s(y= res.re[[i]]$yi, var_y= res.re[[i]]$vi, var_s= ll))
    if (sign.zero!= sign.ll) {break}
    else {
      ll<- ll-.001
    }
  }
  res.re[[i]]$tau2<- uniroot(llk.prime.var_s, interval= c(ll, 0), y= res.re[[i]]$yi, var_y= res.re[[i]]$vi, tol= 1E-10)$root
  rm(list= c("sign.zero", "ll", "sign.ll"))  

  #print(i)
}
# check
# tmp<- sapply(res.re, function(rma.obj) llk.prime.var_s(y= rma.obj$yi, var_y= rma.obj$vi, var_s= rma.obj$tau2))
# which(abs(tmp)> 0.005)
# rm(list= c("xx", "tmp", "i", "j"))
var.prior<- sapply(res.re, function(rma.obj) var.unit.info(y= rma.obj$yi, var_y= rma.obj$vi, var_s= rma.obj$tau2))
summary(var.prior)

var.prior.bk<- var.prior
var.prior[which(var.prior.bk<= 0)]<- median(var.prior.bk)

## check ##
#dd<- unique(sim.data[, c("nstudy", "var.random", "sim.idx", "sid")])
#table(dd$var.random)
#table(dd$nstudy[which(var.prior<0)], dd$var.random[which(var.prior<0)])
#llk<- sapply(res.re, function(rma.obj) llk.prime.var_s(y= rma.obj$yi, var_y= rma.obj$vi, var_s= rma.obj$tau2))


#res.re<- rma(yi= y, vi=var_y, method= "REML", control= list(maxiter= 500, threshold= 1E-10))
#summary(res.re)
#forest(res.re)

### memory reducing strategy ##
yi.lst<- lapply(res.re, function(rma.obj) as.numeric(rma.obj$yi) )
vi.lst<- lapply(res.re, function(rma.obj) as.numeric(rma.obj$vi) )
rm("res.re")

# mcmc.res.lst<-  vector("list", length(yi.lst))
# HPDinterval.lst<- vector("list", length(yi.lst))
for (i in 8808:length(yi.lst)) {
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
    #save(list = c("mcmc.res.lst", "HPDinterval.lst", "var.prior", "yi.lst", "vi.lst"), 
    #     file = "C:/Users/user/Downloads/sim_reduced.RData")
    save(list = ls(all.names = TRUE), file = "C:/Users/user/Downloads/sim_reduced.RData")
  }
}
# save(list = ls(all.names = TRUE), file = "C:/Users/user/Downloads/sim.RData")



HPDinterval(ex4.subgrp)




summary(mcmc.out.lst[, c("fixed.effect", "sd.random", "beta", "lambda", "s1", "s2", "s3")])
plot(mcmc.out.lst[, c("sd.random", "beta", "lambda", "s1", "s2", "s3", "b1", "b2", "b3")])
table(unlist(mcmc.out.lst[, c("fixed.effect")]))
# when var_lambda= 25
#   0     1 
#1251 28749 

#gelman.diag(mcmc.out.lst) 
#HPDinterval(mcmc.out.lst)

var.prior<- var.unit.info(y= y, var_y= var_y, var_s= res.re$tau2)
mcmc.out.unit<- bayesian.mdl.uncertainty.meta.analysis.unit(y= y, 
                                                      var_y= var_y, 
                                                      parms.ini= parms.ini, 
                                                      mu_beta= 0, 
                                                      var_beta= 10000, 
                                                      var_lambda_0= var.prior, 
                                                      prob.random= .5,
                                                      nthin= 10, 
                                                      nburn= 5000, 
                                                      niter= 10000
) 
plot(mcmc.out.unit[, c("beta", "lambda", "s1", "s2", "s3")])
summary(mcmc.out.unit[, c("fixed.effect", "beta", "lambda", "s1", "s2", "s3")])
table(unlist(mcmc.out.unit[, c("fixed.effect")]))
#    0     1 
#10674 19326

mcmc.out.unit2<- bayesian.mdl.uncertainty.meta.analysis.unit(y= y, 
                                                            var_y= var_y, 
                                                            parms.ini= parms.ini, 
                                                            mu_beta= 0, 
                                                            var_beta= 10000, 
                                                            var_lambda_0= var.prior*3, 
                                                            prob.random= .5,
                                                            nthin= 10, 
                                                            nburn= 5000, 
                                                            niter= 10000
) 
plot(mcmc.out.unit2[, c("beta", "lambda", "s1", "s2", "s3")])
summary(mcmc.out.unit2[, c("fixed.effect", "beta", "lambda", "s1", "s2", "s3")])
table(unlist(mcmc.out.unit2[, c("fixed.effect")]))
#   0     1 
#7722 22278 
save(list = ls(all.names = TRUE), file = "sim.RData")



