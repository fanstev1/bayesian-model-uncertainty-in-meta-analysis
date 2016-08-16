library(R2OpenBUGS)
library(coda)
setwd("I:/Users/sfan/Google Drive/MSc HSR thesis work/program/bayesian meta analysis")

meta.data<- data.frame(study= 1:6,
                       studyName= c("NAIA3005", "NAI30010", "NAIA/B2009", "WV15673", "WV15697", "WV15799"), 
                       nc= c(554, 423, 144, 268, 251, 462), 
                       rc= c(34, 40, 9, 19, 6, 34),
                       nt= c(553, 414, 144, 268, 252, 493),
                       rt= c(11, 7, 3, 3, 3, 4))

randomEfft.meta.analysis<- function() {
  for (i in 1:nStudy){
    rc[i] ~ dbin(pc[i], nc[i])
    rt[i] ~ dbin(pt[i], nt[i])
    
    logit(pc[i])<- mu[i]
    logit(pt[i])<- mu[i] + delta[i]
    
    mu[i] ~ dnorm(pbo.effect, tau.pbo.effect)
    delta[i] ~ dnorm(tx.effect, tau.tx.effect)
  }
  
  # prior distributions
  pbo.effect ~ dnorm(0, 1E-6)
  tau.pbo.effect ~ dgamma(.001, .001)
  
  tx.effect ~ dnorm(0, 1E-6)
  tau.tx.effect ~ dgamma(.001, .001)
  
  # OR
  odds.c<- exp(pbo.effect)
  odds.t<- exp(pbo.effect + tx.effect)
  OR<- exp(tx.effect)
  
  # probability of event
  p.c<- odds.c/(1+odds.c)
  p.t<- odds.t/(1+odds.t)
}
write.model(randomEfft.meta.analysis, "randomEfft.meta.analysis.txt")
file.show("randomEfft.meta.analysis.txt")



input.data<- list(nStudy= nrow(meta.data),
                  nc= meta.data$nc,
                  nt= meta.data$nt,
                  rc= meta.data$rc,
                  rt= meta.data$rt)

input.init.1<- list(mu= rep(1, nrow(meta.data)), 
                  delta= rep(-1, nrow(meta.data)),
                  pbo.effect= -1, 
                  tau.pbo.effect=10, 
                  tx.effect= 2, 
                  tau.tx.effect= 10)
input.init.2<- list(mu= rep(-1, nrow(meta.data)), 
                    delta= rep(1, nrow(meta.data)),
                    pbo.effect= 1, 
                    tau.pbo.effect=1, 
                    tx.effect= -1, 
                    tau.tx.effect= 1)
input.init.3<- list(mu= rep(0, nrow(meta.data)), 
                    delta= rep(0, nrow(meta.data)),
                    pbo.effect= 0, 
                    tau.pbo.effect=5, 
                    tx.effect= 0, 
                    tau.tx.effect= 5)
input.init<- list(input.init.1, input.init.2, input.init.3)

input.param<- c("odds.c", "odds.t", "OR", "p.c", "p.t", "pbo.effect", "tx.effect", "tau.tx.effect")

output<- bugs(data= input.data,
              inits= input.init,
              parameters.to.save= input.param,
              n.iter= 20000,
              model.file= "randomEfft.meta.analysis.txt",
              n.chains= 3,
              n.thin= 10,
              codaPkg = TRUE)
codaobject <- read.bugs(output)
plot(codaobject)
summary(codaobject)

# check for convergence
gelman.diag(codaobject)
