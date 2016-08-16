library(R2OpenBUGS)
library(coda)
setwd("I:/Users/sfan/Google Drive/MSc HSR thesis work/program/bayesian meta analysis")

meta.data<- data.frame(study= 1:6,
                       studyName= c("NAIA3005", "NAI30010", "NAIA/B2009", "WV15673", "WV15697", "WV15799"), 
                       nc= c(554, 423, 144, 268, 251, 462), 
                       rc= c(34, 40, 9, 19, 6, 34),
                       nt= c(553, 414, 144, 268, 252, 493),
                       rt= c(11, 7, 3, 3, 3, 4))

meta.analysis.fe<- function() {
  for (i in 1:nStudy){
    rc[i] ~ dbin(pc[i], nc[i])
    rt[i] ~ dbin(pt[i], nt[i])
    
    logit(pc[i])<- mu[i] # study effect
    #mu[i] ~ dnorm(0, 1E-6)
    mu[i] ~ dflat()
    
    logit(pt[i])<- mu[i] + tx.effect
  }    
  #tx.effect ~ dnorm(0, 1E-6)
  tx.effect ~ dflat()
  
  # OR
  OR<- exp(tx.effect)  
}
write.model(meta.analysis.fe, "bayes_meta_analysis_fe.txt")

input.param.fe<- c("OR", "tx.effect", "mu")

input.init.1<- list(mu= rep(1, nrow(meta.data)), 
                    tx.effect= 2)
input.init.2<- list(mu= rep(-1, nrow(meta.data)), 
                    tx.effect= -2)
input.init.3<- list(mu= rep(0, nrow(meta.data)), 
                    tx.effect= 0)
input.init.fe<- list(input.init.1, input.init.2, input.init.3)




# Random effect on the effect measure (log_OR)
meta.analysis.re<- function() {
  for (i in 1:nStudy){
    rc[i] ~ dbin(pc[i], nc[i])
    rt[i] ~ dbin(pt[i], nt[i])
    
    # study effect
    logit(pc[i])<- mu[i] 
    #mu[i] ~ dnorm(pbo.effect, tau.pbo.effect)
    mu[i] ~ dnorm(0, 1E-1)
    
    logit(pt[i])<- mu[i] + delta[i]
    delta[i] ~ dnorm(tx.effect, tau.tx.effect)
  }
  
  # prior distributions
  tx.effect ~ dnorm(0, 1E-1)
  tau.tx.effect ~ dgamma(.001, .001)
  #pbo.effect ~ dnorm(0, 1E-6)
  #tau.pbo.effect ~ dgamma(.001, .001)

  # OR
  OR<- exp(tx.effect)

  # probability of event
  #odds.c<- exp(pbo.effect)
  #odds.t<- exp(pbo.effect + tx.effect)
  #p.c<- odds.c/(1+odds.c)
  #p.t<- odds.t/(1+odds.t)
}
write.model(meta.analysis.re, "bayes_meta_analysis_re.txt")
#file.show("bayes_meta_analysis_re.txt")

input.param.re<- c("OR", "tx.effect", "mu", "delta")

input.init.1<- list(mu= rep(1, nrow(meta.data)), 
                  delta= rep(-1, nrow(meta.data)),
                  tx.effect= 2, 
                  tau.tx.effect= 10)
input.init.2<- list(mu= rep(-1, nrow(meta.data)), 
                    delta= rep(1, nrow(meta.data)),
                    tx.effect= -2, 
                    tau.tx.effect= 1)
input.init.3<- list(mu= rep(0, nrow(meta.data)), 
                    delta= rep(0, nrow(meta.data)),
                    tx.effect= 0, 
                    tau.tx.effect= 5)


#input.init.1<- list(mu= rep(1, nrow(meta.data)), 
#                  delta= rep(-1, nrow(meta.data)),
#                  pbo.effect= -1, 
#                  tau.pbo.effect=10, 
#                  tx.effect= 2, 
#                  tau.tx.effect= 10)
#input.init.2<- list(mu= rep(-1, nrow(meta.data)), 
#                    delta= rep(1, nrow(meta.data)),
#                    pbo.effect= 1, 
#                    tau.pbo.effect=1, 
#                    tx.effect= -1, 
#                    tau.tx.effect= 1)
#input.init.3<- list(mu= rep(0, nrow(meta.data)), 
#                    delta= rep(0, nrow(meta.data)),
#                    pbo.effect= 0, 
#                    tau.pbo.effect=5, 
#                    tx.effect= 0, 
#                    tau.tx.effect= 5)

input.init.re<- list(input.init.1, input.init.2, input.init.3)






input.data<- list(nStudy= nrow(meta.data),
                  nc= meta.data$nc,
                  nt= meta.data$nt,
                  rc= meta.data$rc,
                  rt= meta.data$rt)

res.fe<- bugs(data= input.data,
              inits= input.init.fe,
              parameters.to.save= input.param.fe,
              n.iter= 20000,
              model.file= "bayes_meta_analysis_fe.txt",
              n.chains= 3,
              n.thin= 10,
              codaPkg = TRUE)
  codaobject.fe <- read.bugs(res.fe)
  summary(codaobject.fe)
  plot(codaobject.fe)
  # check for convergence
  gelman.diag(codaobject.fe)  

res.re<- bugs(data= input.data,
              inits= input.init.re,
              parameters.to.save= input.param.re,
              n.iter= 20000,
              model.file= "bayes_meta_analysis_re.txt",
              n.chains= 3,
              n.thin= 10,
              codaPkg = TRUE)
  codaobject.re <- read.bugs(res.re)
  summary(codaobject.re)
  plot(codaobject.re)
  HPDinterval(codaobject.re)
  # check for convergence
  gelman.diag(codaobject.re)
