library(metafor)
library(coda)
library(R2WinBUGS)
#library(R2OpenBUGS)
setwd("I:/Users/sfan/Google Drive/MSc HSR thesis work/program/bayesian meta analysis")

# DATA 1: EXAMPLE 1 FROM COMPREHENSIVE DECISION ANALYTICAL MODELLING IN ECONOMIC EVALUATION
meta.data<- data.frame(study= 1:6,
                       studyName= c("NAIA3005", "NAI30010", "NAIA/B2009", "WV15673", "WV15697", "WV15799"), 
                       nc= c(554, 423, 144, 268, 251, 462), 
                       rc= c(34, 40, 9, 19, 6, 34),
                       nt= c(553, 414, 144, 268, 252, 493),
                       rt= c(11, 7, 3, 3, 3, 4))

# DATA 2: SYSTEMATIC REVIEWS IN HEALTH CARE CASE STUDY 1 pg 303
#meta.data<- data.frame(study= 1:6,
#                       studyName= c("Breart 1992", "Breart 1992", "Gagnon 1997", "Hodnett 1989", "Kennell 1991", "Langer 1998"), 
#                       nc= c(131, 664, 204, 73, 200, 363), 
#                       rc= c(62, 319, 142, 43, 55, 303),
#                       nt= c(133, 656, 209, 72, 212, 361),
#                       rt= c(55, 281, 139, 30, 24, 205))
attach(meta.data)

eff.summary<- escalc(measure= "OR",
                     ai= rt,
                     n1i= nt,
                     ci= rc,
                     n2i= nc, 
                     data= meta.data)

eff.summary<- cbind(eff.summary, logor.ctl=log(rc/(nc-rc))) # adding logit(event rate) in the control group

detach(meta.data)



meta.analysis.re.ZI_half<- function() {
  for (i in 1:nStudy){
    yi[i] ~ dnorm(theta[i], tau.y[i])
    theta[i]<- mu.theta + lambda * b[i]
    tau.y[i]<- 1/vi[i]
  }
  # prior on the overall effect
  mu.theta ~ dnorm(0, 1E-6)
  
  # prior distributions on lamda and b[i] -> induce half t distribution on the variance component of the random study effect
    lambda ~ dflat()  
    PI<- 3.141592653589793238462643383279502884197169399375105820974944592307816406286
    #loglik<-  (equals(lambda, 0) * log(1-p0) + 
    #                            (1-equals(lambda, 0))* ( 
    #                              log(p0) +  0.5*log(2/PI) - log(prior.scale) -0.5 * pow(lambda/prior.scale, 2) ) 
    #                          )
    loglik<- log( equals(lambda, 0) * (1-p0) +
                    (1-equals(lambda, 0))* ( p0 * pow(2*PI, -.5) * pow(prior.scale, -1) * exp(-.5 * pow(lambda/prior.scale, 2)) )
      
    )
    C<- 100000
    loglik.c<-  -1*loglik + C       
    dummy <- 0
    dummy ~ dpois(loglik.c)
  #  p0 ~ dunif(0, 1)
  #lambda  ~ dnorm(0, tau.lambda)
  #tau.lambda<- pow(prior.scale, -2)
  #p0 ~ dbern(prior.re.prob)
  
  for (i in 1:nStudy){
    b[i] ~ dnorm(0, tau.b)
  }
  prior.scale.tau.b<- .5 * nStudy
  tau.b ~ dgamma(.5, prior.scale.tau.b)  # non-central chi squared distribution with 1 df
  sigma.theta<- abs(lambda)/sqrt(tau.b)
  
  
  # OR
  or.overall<- exp(mu.theta)
}
write.model(meta.analysis.re.ZI_half, "meta_analysis_re_ZI_half.txt")
# file.show("meta_analysis_re_ZI_half.txt")


input.param.re<- c("theta", "mu.theta", "sigma.theta", "or.overall", "lambda")

input.init.1<- list(mu.theta= rnorm(1), 
                    b= rnorm(nrow(meta.data)),
                    lambda= rnorm(1), 
                    #theta= rnorm(nrow(meta.data)),
                    #p0= as.integer(0),
                    tau.b= runif(1)                   
)
input.init.2<- list(mu.theta= rnorm(1), 
                    b= rnorm(nrow(meta.data)),
                    lambda= rnorm(1), 
                    #theta= rnorm(nrow(meta.data)),
                    #p0= as.integer(1),
                    tau.b= runif(1)                   
)
input.init.3<- list(mu.theta= rnorm(1), 
                    b= rnorm(nrow(meta.data)),
                    lambda= rnorm(1), 
                    #theta= rnorm(nrow(meta.data)),
                    #p0= as.integer(1),
                    tau.b= runif(1)                   
)
input.init.re<- list(input.init.1, input.init.2, input.init.3)






input.data<- list(nStudy= as.integer(nrow(meta.data)),
                  yi= as.numeric(eff.summary$yi),
                  vi= as.numeric(eff.summary$vi),
                  prior.scale= 25,
                  p0= .5
                  )


res.re<- bugs(data= input.data,
              inits= input.init.re,
              parameters.to.save= input.param.re,
              n.iter= 6000,
              model.file= "meta_analysis_re_ZI_half.txt",
              n.chains= 3,
              n.thin= 10,
              n.burnin= 1000,
              bugs.directory= "I:/Users/sfan/winbugs14",
              #OpenBUGS.pgm= "C:/Program Files (x86)/OpenBUGS/OpenBUGS323/OpenBUGS.exe",
              #debug= TRUE,
              over.relax= TRUE,
              codaPkg = TRUE)

codaobject.re <- read.bugs(res.re)
  # check for convergence
  gelman.diag(codaobject.re)  

summary(codaobject.re)
plot(codaobject.re)
