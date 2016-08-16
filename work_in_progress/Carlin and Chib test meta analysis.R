library(metafor)
library(coda)
library(R2WinBUGS)
#library(R2OpenBUGS)
setwd("I:/Users/sfan/Google Drive/MSc HSR thesis work/program/bayesian meta analysis")
source("Carlin and Chib gibbs sampler for bayesian model uncertainty in meta analysis.R")

meta.analysis.fe<- function() {
  for (i in 1:nStudy){
    y[i] ~ dnorm( delta, tau_y[i])
  }
  delta ~ dnorm(0, 1E-6)
  #tx.effect ~ dflat()
  
  # OR
  OR<- exp(delta)  
}
write.model(meta.analysis.fe, "bayes_meta_analysis_normal_approx_fe.txt")
meta.analysis.re<- function() {
  for (i in 1:nStudy){
    y[i] ~ dnorm( delta[i], tau_y[i])
    delta[i]<- beta + lambda*b[i]
    b[i] ~ dnorm(0, tau.b)
  }
  
  # prior distributions
  lambda ~ dnorm(0, tau.lambda)
  tau.lambda<- pow(prior.scale, -2)
  beta ~ dnorm(0, 1E-6)
  tau.b ~ dgamma(.5, .5)
  
  # OR
  OR<- exp(beta)
  sigma.study.effect<- abs(lambda)/sqrt(tau.b)
  for (i in 1:nStudy){
    s[i]<- lambda*b[i]
  }
}
write.model(meta.analysis.re, "bayes_meta_analysis_normal_approx_re.txt")

#meta.data<- data.frame(study= 1:22,
#                       studyName= paste(1:22), 
#                       nc= c(39, 116, 93, 1520, 365, 52, 939, 471, 282, 1921, 583, 266, 293, 883, 147, 213, 122, 154, 134, 218, 364, 674), 
#                       rc= c(3, 14, 11, 127, 27, 6, 152, 48, 37, 188, 52, 47, 16, 45, 31, 38, 12, 6, 3, 40, 43, 39),
#                       nt= c(38, 114, 69, 1533, 355, 59, 945, 632, 278, 1916, 873, 263, 291, 858, 154, 207, 251, 151, 174, 209, 391, 680),
#                       rt= c(3, 7, 5, 102, 28, 4, 98, 60, 25, 138, 64, 45, 9, 57, 25, 33, 28, 8, 6, 32, 27, 22))

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

## data
y<- as.numeric(eff.summary$yi)
var_y<- as.numeric(eff.summary$vi)  
nobs<- length(y)

input.data<- list(nStudy= nobs, y= y, tau_y= 1/var_y)


input.param.fe<- c("delta")

input.init.1<- list(delta= 2)
input.init.2<- list(delta= -2)
input.init.3<- list(delta= 0)
input.init.fe<- list(input.init.1, input.init.2, input.init.3)

res.fe<- bugs(data= input.data,
              inits= input.init.fe,
              parameters.to.save= input.param.fe,
              n.iter= 20000,
              model.file= "bayes_meta_analysis_normal_approx_fe.txt",
              n.chains= 3,
              n.thin= 10,
              bugs.directory = "I:/Users/sfan/WinBUGS14/",
              #debug= TRUE,
              codaPkg = TRUE)
plot(codaobject.fe <- read.bugs(res.fe))
summary(codaobject.fe)
summary(codaobject.fe[, "delta"])$statistics



# Random effect on the effect measure (log_OR)

input.param.re<- c("tau.b", "beta", "lambda","b", "sigma.study.effect", "s")

input.init.1<- list(b= rep(1, nrow(meta.data)), 
                    #delta= rep(-1, nrow(meta.data)),
                    beta= 2, 
                    lambda= rnorm(n= 1, sd= 5),
                    tau.b= 10)
input.init.2<- list(b= rep(-1, nrow(meta.data)), 
                    #delta= rep(1, nrow(meta.data)),
                    beta= -2, 
                    lambda= rnorm(n= 1, sd= 5),
                    tau.b= 1)
input.init.3<- list(b= rep(0, nrow(meta.data)), 
                    #delta= rep(0, nrow(meta.data)),
                    beta= 0, 
                    lambda=rnorm(n= 1, sd= 5),
                    tau.b= 5)
input.init.re<- list(input.init.1, input.init.2, input.init.3)


input.data<- list(nStudy= nobs, y= y, tau_y= 1/var_y, prior.scale= 25)

res.re<- bugs(data= input.data,
              inits= input.init.re,
              parameters.to.save= input.param.re,
              n.iter= 20000,
              model.file= "bayes_meta_analysis_normal_approx_re.txt",
              n.chains= 3,
              n.thin= 10,
              bugs.directory = "I:/Users/sfan/WinBUGS14/",
              #debug= TRUE,
              codaPkg = TRUE)


plot( codaobject.re <- read.bugs(res.re))
summary(codaobject.re)$statistics
(summary(codaobject.re[, "tau.b"])$statistics["Mean"]/summary(codaobject.re[, "tau.b"])$statistics["SD"])^2
summary(codaobject.re[, "tau.b"])$statistics["SD"]^2/summary(codaobject.re[, "tau.b"])$statistics["Mean"]

as.numeric(summary(codaobject.re[, paste("b[", 1:nobs, "]", sep= "")])$statistics[,"Mean"])
as.numeric(summary(codaobject.re[, paste("b[", 1:nobs, "]", sep= "")])$statistics[,"SD"])



### testing 


prior.parm.mdl1<- list(delta= list(mu= 0, var= 100))
# pseduo prior for Model 1 parameters when the model 2 is chosen
pseudo.prior.parm.mdl1<- list(delta= list(mu= summary(codaobject.fe)$statistics["delta", "Mean"], 
                                          var= (5*summary(codaobject.fe)$statistics["delta", "SD"])^2
)
)



prior.parm.mdl2<- list(beta   = list(mu= 0, var= 100),
                       lambda = list(mu= 0, var= 625),
                       var_b  = list(alpha= .5, theta= .5*length(y))
)
# pseudo prior for hyperparameters in model 2 given model 1
pseudo.prior.parm.mdl2<- list(beta  = list(mu= summary(codaobject.re)$statistics["beta", "Mean"], 
                                           var= (5*summary(codaobject.re)$statistics["beta", "SD"])^2
                                           ),

                              lambda= list(mu= as.numeric(summary(codaobject.re[,"lambda"])$statistics["Mean"]), 
                                           var= (5*as.numeric(summary(codaobject.re[,"lambda"])$statistics["SD"]))^2
                                           ),
                              
                              b     = list(mu  = as.numeric(summary(codaobject.re)$statistics[paste("b[", 1:nobs, "]", sep= ""),"Mean"]),
                                           var = (5*as.numeric(summary(codaobject.re)$statistics[paste("b[", 1:nobs, "]", sep= ""),"SD"]))^2
                                           ),
                              
                              var_b = list(alpha= 1, theta= 2)
                              )

# initial value
parm1<- list(mdl.idx= 1, 
             parm1= list(delta= 0), 
             parm2= list(beta= 0, lambda= rnorm(1, sd= 25), b= rnorm(length(y), sd= 10), var_b= 100)
            )
parm2<- list(mdl.idx= 1, 
             parm1= list(delta= 0), 
             parm2= list(beta= 0, lambda= rnorm(1, sd= 25), b= rnorm(length(y), sd= 10), var_b= 100)
             )
parm3<- list(mdl.idx= 1, 
             parm1= list(delta= 0), 
             parm2= list(beta= 0, lambda= rnorm(1, sd= 25), b= rnorm(length(y), sd= 10), var_b= 100)
             )

parm.ini<- list(parm= parm1, parm= parm2, parm= parm3)

mcmc.out.lst<- bayes.mdl.uncertainty.meta.analysis(y= y, var_y, parm.ini= parm.ini,
                                                   niter= 10000, nburn= 5000, nthin= 10,
                                                   prior.prob.fixed= .5,
                                                   prior.parm.mdl1= prior.parm.mdl1, 
                                                   pseudo.prior.parm.mdl1= pseudo.prior.parm.mdl1,
                                                   prior.parm.mdl2= prior.parm.mdl2, 
                                                   pseudo.prior.parm.mdl2= pseudo.prior.parm.mdl2) 

#save(list = ls(all.names = TRUE), file = "Carlin_Chib_test.RData")
table(unlist(mcmc.out.lst[,"mdl.idx"]))
plot(mcmc.out.lst)


dd<- data.frame(chain.idx= rep( 1:3, each=10000),
                mdl.idx= unlist(mcmc.out.lst[,"mdl.idx"]),
                s1= unlist(mcmc.out.lst[,"parm2.s1"]),
                s2= unlist(mcmc.out.lst[,"parm2.s2"]),
                s3= unlist(mcmc.out.lst[,"parm2.s3"]),
                s4= unlist(mcmc.out.lst[,"parm2.s4"]),
                s5= unlist(mcmc.out.lst[,"parm2.s5"]),
                s6= unlist(mcmc.out.lst[,"parm2.s5"]),
                var_s= unlist(mcmc.out.lst[,"parm2.var_s"]),
                lambda= unlist(mcmc.out.lst[,"parm2.lambda"])
                )

plot(dd$s1[which(dd$mdl.idx==2 & dd$chain.idx==1)], col= "red", type= "l")
lines(dd$s1[which(dd$mdl.idx==2 & dd$chain.idx==2)], col= "blue", type= "l")
lines(dd$s1[which(dd$mdl.idx==2 & dd$chain.idx==3)], col= "green", type= "l")

plot(dd$s2[which(dd$mdl.idx==2 & dd$chain.idx==1)], col= "red", type= "l")
lines(dd$s2[which(dd$mdl.idx==2 & dd$chain.idx==2)], col= "blue", type= "l")
lines(dd$s2[which(dd$mdl.idx==2 & dd$chain.idx==3)], col= "green", type= "l")

plot( dd$s3[which(dd$mdl.idx==2 & dd$chain.idx==1)], col= "red", type= "l")
lines(dd$s3[which(dd$mdl.idx==2 & dd$chain.idx==2)], col= "blue", type= "l")
lines(dd$s3[which(dd$mdl.idx==2 & dd$chain.idx==3)], col= "green", type= "l")

plot( dd$var_s[which(dd$chain.idx==1)], col= "red", type= "l", ylim= c(0, 140))
lines(dd$var_s[which(dd$chain.idx==2)], col= "blue", type= "l")
lines(dd$var_s[which(dd$chain.idx==3)], col= "green", type= "l")

tapply(dd$s1, dd$mdl.idx, summary)
apply(dd, 2, tapply, dd$mdl, hist)

rand.effect.mdl.idx<- which(unlist(mcmc.out.lst[,"mdl.idx"])==2)
s1<- unlist(mcmc.out.lst[,"parm2.s1"])
s2<- unlist(mcmc.out.lst[,"parm2.s2"])
s3<- unlist(mcmc.out.lst[,"parm2.s3"])
s4<- unlist(mcmc.out.lst[,"parm2.s4"])
s5<- unlist(mcmc.out.lst[,"parm2.s5"])
s6<- unlist(mcmc.out.lst[,"parm2.s5"])
var_s<- unlist(mcmc.out.lst[,"parm2.var_s"])

summary(s1[rand.effect.mdl.idx])
summary(s2[rand.effect.mdl.idx])
summary(s3[rand.effect.mdl.idx])
summary(s4[rand.effect.mdl.idx])
summary(s5[rand.effect.mdl.idx])
summary(s6[rand.effect.mdl.idx])

plot(density(s1[rand.effect.mdl.idx]))
plot(density(s2[rand.effect.mdl.idx]))
plot(density(s3[rand.effect.mdl.idx]))
plot(density(s4[rand.effect.mdl.idx]))
plot(density(s5[rand.effect.mdl.idx]))
plot(density(s6[rand.effect.mdl.idx]))
plot(density(var_s[rand.effect.mdl.idx]))
plot(density(abs(dd$lambda[rand.effect.mdl.idx])))
### end of testing
hist(abs(unlist(mcmc.out.lst[, c("parm2.lambda")])) * unlist(mcmc.out.lst[, "parm2.var_b"]))

length(dd$s1)

table(mcmc.out.lst[,"mdl.idx"])
table(unlist(mcmc.out.lst[,"mdl.idx"]))
