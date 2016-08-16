library(metafor)
library(coda)
library(R2WinBUGS)
#library(R2OpenBUGS)
setwd("I:/Users/sfan/Google Drive/MSc HSR thesis work/program/bayesian meta analysis")
source("gibbs sampler for bayesian model uncertainty in meta analysis.R")
source("functions to calculate unit information.R")
# load("ex4.RData")

# DATA 4: Meta analysis EXAMPLE from Bayesian Data Analysis by Andrew Gelman
meta.data.af<- data.frame(study= 1:4,
                          studyName= c("CIBIS-II", "MERIT-HF", "SENIORS", "US-Carvedilol"), 
                          subgrp= c("AF"),
                          nc= c(264, 282, 237, 52), 
                          rc= c(43, 31, 51, 6),
                          nt= c(257, 274, 227, 84),
                          rt= c(42, 30, 38, 4))

meta.data.sr<- data.frame(study= 1:4,
                          studyName= c("CIBIS-II", "MERIT-HF", "SENIORS", "US-Carvedilol"), 
                          subgrp= c("SR"),
                          nc= c(1004, 1569, 444, 346), 
                          rc= c(170, 160, 84, 25),
                          nt= c(1014, 1563, 451, 612),
                          rt= c(105, 102, 77, 18))

#attach(meta.data.af)
eff.summary.af<- escalc(measure= "OR",
                        ai= rt,
                        n1i= nt,
                        ci= rc,
                        n2i= nc, 
                        data= meta.data.af)
eff.summary.sr<- escalc(measure= "OR",
                        ai= rt,
                        n1i= nt,
                        ci= rc,
                        n2i= nc, 
                        data= meta.data.sr)

eff.summary<- rbind(eff.summary.af, eff.summary.sr)

eff.summary$af <- ifelse(eff.summary$subgrp == "AF", 1, 0)
eff.summary$sr <- ifelse(eff.summary$subgrp == "SR", 1, 0)
#detach(meta.data.af)


#eff.summary<- cbind(eff.summary, logor.ctl=log(rc/(nc-rc))) # adding logit(event rate) in the control group

save.image("ex4.RData")

fe.af<- rma(method= "FE",
            measure= "OR",
            ai= rt,
            n1i= nt,
            ci= rc,
            n2i= nc,
            subset = (subgrp== "AF"),
            data= eff.summary)


re.af<- rma(method= "REML",
            measure= "OR",
            ai= rt,
            n1i= nt,
            ci= rc,
            n2i= nc,
            subset = (subgrp== "AF"),
            data= eff.summary)


fe.sr<- rma(method= "FE",
            measure= "OR",
            ai= rt,
            n1i= nt,
            ci= rc,
            n2i= nc,
            subset = (subgrp== "SR"),
            data= eff.summary)

re.sr<- rma(method= "REML",
            measure= "OR",
            ai= rt,
            n1i= nt,
            ci= rc,
            n2i= nc,
            subset = (subgrp== "SR"),
            data= eff.summary)

par(mfrow= c(2,2))
radial(fe.af, main = "Fixed-Effects Model")
radial(re.af, main = "Random-Effects Model")
radial(fe.sr, main = "Fixed-Effects Model")
radial(re.sr, main = "Random-Effects Model")


fe.overall<- rma(yi, vi, method= "FE", data= eff.summary)
fe.subgrp<- rma(yi, vi, method= "FE", intercept = FALSE, mods= cbind(af, sr), data= eff.summary)

re.subgrp<- rma(yi, vi, method= "REML", intercept = FALSE, mods= cbind(af, sr), data= eff.summary)

preds.overall<- predict(fe.overall)
preds.subgrp<- predict(fe.subgrp, newmods=  as.matrix(data.frame(af= c(1, 0), sr= c(0, 1))))

forest(eff.summary$yi, eff.summary$vi, atransf = exp, 
       ylim = c(-3.5, 11), xlim = c(-7, 7),
       slab = paste(eff.summary$studyName, " (",eff.summary$subgrp, ")", sep= ""))
abline(h = 0, lwd= 2)
abline(h= 4.5, lty= 3)
addpoly(c(preds.subgrp$pred, preds.overall$pred), 
        sei = c(preds.subgrp$se, preds.overall$se), 
        atransf = exp, mlab = c("AF", "SR", "Overall"))
text(-7, 10, "Study (Subgroup)", pos = 4, font = 2)
text(7, 10, "Odds-Ratio [95% CI]", pos = 2, font = 2)



## data

y<- as.numeric(eff.summary$yi)
var_y<- as.numeric(eff.summary$vi)
subgrp<- ifelse(eff.summary$subgrp=="AF", 1, 2)
nobs<- length(unique(eff.summary$study))

# initial 
set.seed(45)
beta<- rnorm(2, mean= sum(y*(1/var_y)/sum(1/var_y)), sd= sqrt(sum(1/var_y)))
lambda<- rnorm(2)
b<- rnorm(length(as.numeric(eff.summary$yi)))
var_b<- rinvgamma(2, shape= .5, scale= .5)

parms.1<- list(beta= beta,
               lambda= lambda,
               b= b,
               var_b= var_b,
               ## additional parameters to keep track 
               sd.random= abs(lambda) * sqrt(var_b),
               or= exp(beta),
               fixed.effect= (lambda==0),
               diff.subgrp= NA,
               or.subgrp= NA,
               effsmaller.subgrp= NA
               # model checking
               #dev= rep(NA, 1), 
               #p.inv= rep(NA, nobs),
               #y.rep= rep(NA, nobs),
               #p.smaller= rep(NA, nobs),
               #dev.rep= rep(NA)
)
set.seed(145)
beta<- rnorm(2, mean= sum(y*(1/var_y)/sum(1/var_y)), sd= sqrt(sum(1/var_y)))
lambda<- rnorm(2)
b<- rnorm(length(as.numeric(eff.summary$yi)))
var_b<- rinvgamma(2, shape= .5, scale= .5)

parms.2<- list(beta= beta,
               lambda= lambda,
               b= b,
               var_b= var_b,
               ## additional parameters to keep track 
               sd.random= abs(lambda) * sqrt(var_b),
               or= exp(beta),
               fixed.effect= (lambda==0),
               diff.subgrp= NA,
               or.subgrp= NA,
               effsmaller.subgrp= NA
               # model checking
               #dev= rep(NA, 1), 
               #p.inv= rep(NA, nobs),
               #y.rep= rep(NA, nobs),
               #p.smaller= rep(NA, nobs),
               #dev.rep= rep(NA)
)

set.seed(245)
beta<- rnorm(2, mean= sum(y*(1/var_y)/sum(1/var_y)), sd= sqrt(sum(1/var_y)))
lambda<- rnorm(2)
b<- rnorm(length(as.numeric(eff.summary$yi)))
var_b<- rinvgamma(2, shape= .5, scale= .5)

parms.3<- list(beta= beta,
               lambda= lambda,
               b= b,
               var_b= var_b,
               ## additional parameters to keep track 
               sd.random= abs(lambda) * sqrt(var_b),
               or= exp(beta),
               fixed.effect= (lambda==0),
               diff.subgrp= NA,
               or.subgrp= NA,
               effsmaller.subgrp= NA
               # model checking
               #dev= rep(NA, 1), 
               #p.inv= rep(NA, nobs),
               #y.rep= rep(NA, nobs),
               #p.smaller= rep(NA, nobs),
               #dev.rep= rep(NA)
)

parms.ini<- list(parms.1, parms.2, parms.3)

# AF subgroup (= 1)
llk<- llk.prime.var_s(y= re.af$yi, var_y= re.af$vi, var_s= re.af$tau2)
  sign.zero<- sign(llk.prime.var_s(y= re.af$yi, var_y= re.af$vi, var_s= re.af$tau2))
  ll<- re.af$tau2 - 0.001
  repeat {
    sign.ll<- sign(llk.prime.var_s(y= re.af$yi, var_y= re.af$vi, var_s= ll))
    if (sign.zero!= sign.ll) {break}
    else {
      ll<- ll-.001
    }
  }
  re.af$tau2<- uniroot(llk.prime.var_s, interval= c(ll, 0), y= re.af$yi, var_y= re.af$vi, tol= 1E-10)$root
  rm(list= c("sign.zero", "ll", "sign.ll"))  
var.unit.info(y= re.af$yi, var_y= re.af$vi, var_s= re.af$tau2)

# SR subgroup (= 2)
var.unit.info(y= re.sr$yi, var_y= re.sr$vi, var_s= re.sr$tau2)

var.unit.info(y= re.subgrp$yi, var_y= re.subgrp$vi, var_s= re.subgrp$tau2)

ex4.subgrp<- bayesian.mdl.uncertainty.meta.analysis.subgrp(y= as.numeric(eff.summary$yi), 
                                                           var_y= as.numeric(eff.summary$vi), 
                                                           subgrp= ifelse(eff.summary$subgrp=="AF", 1, 2),
                                                           parms.ini= parms.ini, 
                                                           mu_beta= rep(0, 2), 
                                                           var_beta= rep(10000, 2), 
                                                           alpha_0= rep(.5, 2), 
                                                           theta_0= rep(.5*nobs, 2), 
                                                           var_lambda_0= rep(0.1556011, 2), 
                                                           prob.random= rep(.5, 2),
                                                           nthin= 10, 
                                                           nburn= 5000, 
                                                           niter= 10000
)



gelman.diag(ex4.subgrp) 
autocorr.plot(ex4.subgrp)
HPDinterval(ex4.subgrp)
summary(ex4.subgrp[, c("or1", "or2", "or.subgrp", "fixed.effect1", "fixed.effect2", "sd.random1", "sd.random2")])
plot(ex4.subgrp[, c("or1", "or2", "or.subgrp", "fixed.effect1", "fixed.effect2", "sd.random1", "sd.random2")])
autocorr.plot(ex4.subgrp[, c("or1", "or2", "or.subgrp", "fixed.effect1", "fixed.effect2", "sd.random1", "sd.random2")])
HPDinterval(ex4.subgrp[, c("or1", "or2", "or.subgrp", "fixed.effect1", "fixed.effect2", "sd.random1", "sd.random2")])


par(mfcol= c(3, 1))
plot(density(unlist(ex4.subgrp[, c("or2")])), main= "Treatment effect in HF patients with SR", xlim= c(0, 1.5), xlab= "Odds ratio")
#      abline(v= mean(unlist(ex4.subgrp[, c("or2")])), lwd= 2)
plot(density(unlist(ex4.subgrp[, c("or1")])), main= "Treatment effect in HF patients with AF", xlim= c(0, 1.5), xlab= "Odds ratio")
#      abline(v= mean(unlist(ex4.subgrp[, c("or1")])), lwd= 2)
plot(density(unlist(ex4.subgrp[, c("or.subgrp")])), main= "Differential treatment effect (SR vs. AF)", xlim= c(0, 1.5), xlab= "Odds ratio")
#      abline(v= mean(unlist(ex4.subgrp[, c("or.subgrp")])), lwd= 2)

# SR subgroup (= 2)
mean(unlist(ex4.subgrp[, c("fixed.effect2")])>0)/(1-mean(unlist(ex4.subgrp[, c("fixed.effect2")])>0))

mean(unlist(ex4.subgrp[, c("fixed.effect1")])>0)/(1-mean(unlist(ex4.subgrp[, c("fixed.effect1")])>0))
save.image("ex4.RData")







# common tau^2
re<- rma(yi, vi, method= "REML",
         mods= ~factor(subgrp)-1,
         data= eff.summary)



rma(yi, vi, mods = cbind(ablat), data = dat)