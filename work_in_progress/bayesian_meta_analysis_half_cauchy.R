library(metafor)
library(coda)
library(R2WinBUGS)
#library(R2OpenBUGS)
setwd("I:/Users/sfan/Google Drive/MSc HSR thesis work/program/bayesian meta analysis")

# DATA 1: EXAMPLE 1 FROM COMPREHENSIVE DECISION ANALYTICAL MODELLING IN ECONOMIC EVALUATION
#meta.data<- data.frame(study= 1:6,
#                       studyName= c("NAIA3005", "NAI30010", "NAIA/B2009", "WV15673", "WV15697", "WV15799"), 
#                       nc= c(554, 423, 144, 268, 251, 462), 
#                       rc= c(34, 40, 9, 19, 6, 34),
#                       nt= c(553, 414, 144, 268, 252, 493),
#                       rt= c(11, 7, 3, 3, 3, 4))

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

meta.analysis.re.half.cauchy<- function() {
  for (i in 1:nStudy){
    yi[i] ~ dnorm(theta[i], tau.y[i])
    theta[i]<- mu.theta + lambda * b[i]
    tau.y[i]<- 1/vi[i]
  }
  # prior on the overall effect
  mu.theta ~ dnorm(0, 1E-6)
  
  # prior distributions on lamda and b[i] -> induce half t distribution on the variance component of the random study effect
  lambda  ~ dnorm(0, tau.lambda)
  tau.lambda<- pow(prior.scale, -2)
  
  for (i in 1:nStudy){
    b[i] ~ dnorm(0, tau.b)
  }  
  tau.b ~ dgamma(.5, .5)  # chi squared distribution with 1 df
  sigma.theta<- abs(lambda)/sqrt(tau.b) # sd for study random effect
  
  
  # OR
  or.overall<- exp(mu.theta)
}

write.model(meta.analysis.re.half.cauchy, "meta.analysis.re.half.cauchy.txt")
#file.show("meta.analysis.re.half.cauchy.txt")

input.param.re<- c("theta", "mu.theta", "sigma.theta", "or.overall")

input.init.1<- list(mu.theta= rnorm(1), 
                    b= rnorm(nrow(meta.data)),
                    lambda= rnorm(1), 
                    tau.b= runif(1), 
                    theta= rnorm(nrow(meta.data))
                    )
input.init.2<- list(mu.theta= rnorm(1), 
                    b= rnorm(nrow(meta.data)),
                    lambda= rnorm(1), 
                    tau.b= runif(1), 
                    theta= rnorm(nrow(meta.data))
                    )
input.init.3<- list(mu.theta= rnorm(1), 
                    b= rnorm(nrow(meta.data)),
                    lambda= rnorm(1), 
                    tau.b= runif(1), 
                    theta= rnorm(nrow(meta.data))
                    )
input.init.re<- list(input.init.1, input.init.2, input.init.3)






input.data<- list(nStudy= nrow(meta.data),
                  yi= as.numeric(eff.summary$yi),
                  vi= as.numeric(eff.summary$vi),
                  prior.scale= 25)

res.re<- bugs(data= input.data,
              inits= input.init.re,
              parameters.to.save= input.param.re,
              n.iter= 60000,
              model.file= "meta.analysis.re.half.cauchy.txt",
              n.chains= 3,
              n.thin= 10,
              n.burnin= 10000,
              bugs.directory= "I:/Users/sfan/winbugs14",
              #debug= TRUE,
              codaPkg = TRUE)
codaobject.re <- read.bugs(res.re)
  # check for convergence
  gelman.diag(codaobject.re)  

summary(codaobject.re)
plot(codaobject.re)








# adjust for underlying risk of influenza


meta.analysis.re.half.cauchy<- function() {
  for (i in 1:nStudy){
    yi[i] ~ dnorm(theta[i], tau.y[i])
    theta[i]<- mu.theta + lambda * b.adj[i] + gamma * (logor.ctl[i] - mean(logor.ctl[]))
    tau.y[i]<- 1/vi[i]
  }
  # prior on the overall effect
  mu.theta ~ dnorm(0, 1E-6)
  gamma ~ dnorm(0, 1E-6)
  
  # prior distributions on lamda and b[i] -> induce half t distribution on the variance component of the random study effect
  lambda  ~ dnorm(0, tau.lambda)
  tau.lambda<- pow(prior.scale, -2)
  
  for (i in 1:nStudy){
    b.adj[i] ~ dnorm(0, tau.b)
  }  
  tau.b ~ dgamma(.5, .5)  # chi squared distribution with 1 df
  sigma.theta<- abs(lambda)/sqrt(tau.b)
  
  
  # OR
  or.overall<- exp(mu.theta)
}

write.model(meta.analysis.re.half.cauchy, "meta.analysis.re.half.cauchy.adj.txt")

input.param.re<- c("theta", "mu.theta", "sigma.theta", "or.overall", "gamma")

input.init.1<- list(mu.theta= rnorm(1), 
                    b.adj= rnorm(nrow(meta.data)),
                    lambda= rnorm(1), 
                    tau.b= runif(1), 
                    gamma= rnorm(1),
                    theta= rnorm(nrow(meta.data))
)
input.init.2<- list(mu.theta= rnorm(1), 
                    b.adj= rnorm(nrow(meta.data)),
                    lambda= rnorm(1), 
                    tau.b= runif(1), 
                    gamma= rnorm(1),
                    theta= rnorm(nrow(meta.data))
)
input.init.3<- list(mu.theta= rnorm(1), 
                    b.adj= rnorm(nrow(meta.data)),
                    lambda= rnorm(1), 
                    tau.b= runif(1),
                    gamma= rnorm(1),
                    theta= rnorm(nrow(meta.data))
)
input.init.re<- list(input.init.1, input.init.2, input.init.3)


input.data<- list(nStudy= nrow(meta.data),
                  yi= as.numeric(eff.summary$yi),
                  vi= as.numeric(eff.summary$vi),
                  logor.ctl= as.numeric(eff.summary$logor.ctl),
                  prior.scale= 25)

res.re.adj<- bugs(data= input.data,
              inits= input.init.re,
              parameters.to.save= input.param.re,
              n.iter= 1000000,
              model.file= "meta.analysis.re.half.cauchy.adj.txt",
              n.chains= 3,
              n.thin= 200,
              n.burnin= 200000,
              bugs.directory= "I:/Users/sfan/winbugs14",
              #debug= TRUE,
              codaPkg = TRUE)
codaobject.re.adj <- read.bugs(res.re.adj)
# check for convergence
gelman.diag(codaobject.re.adj)  

summary(codaobject.re.adj)
plot(codaobject.re.adj)
