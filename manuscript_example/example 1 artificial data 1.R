library(metafor)
library(coda)
#library(R2WinBUGS)
#library(R2OpenBUGS)
#setwd("I:/users/sfan/Google Drive/MSc HSR thesis work/document/")
#source("I:/Users/sfan/Google Drive/MSc HSR thesis work/program/bayesian meta analysis/gibbs sampler for bayesian model uncertainty in meta analysis.R")
#source("I:/Users/sfan/Google Drive/MSc HSR thesis work/program/bayesian meta analysis/functions to calculate unit information.R")
setwd("L:/Steve Fan/Documents/Research/model uncertainty manuscript example")
source("gibbs sampler for bayesian model uncertainty in meta analysis.R")
source("functions to calculate unit information.R")
op<- par()
load("artificial_data_ex.RData")

hypo.data<- data.frame(study= 1:5,
                       yi= c(0.44591340, -0.88931133, -1.58538866, -1.34807315, -0.46941765),
                       vi= c(0.532505845, 0.325584765, 0.194581121, 0.415367965,0.056434210))


hypo.data.2<- hypo.data.1<- hypo.data
hypo.data.1$vi<- hypo.data$vi*2
hypo.data.2$vi<- hypo.data$vi


# conduct random effect model analysis
re.1<- rma(method= "REML", yi= yi, vi= vi, data= hypo.data.1)
fe.1<- rma(method= "FE", yi= yi, vi= vi, data= hypo.data.1)

re.2<- rma(method= "REML", yi= yi, vi= vi, data= hypo.data.2)
fe.2<- rma(method= "FE", yi= yi, vi= vi, data= hypo.data.2)

#par(mfrow= c(2, 2))
pdf("artificial_data_ex.pdf", onefile = TRUE, paper= "letter", width = 7, height= 7)
par(mar=c(4,4,1,2))
par(mfrow= c(2, 1))
forest(hypo.data.1$yi, hypo.data.1$vi, xlab= expression(paste("Heterogeneity: Q= 4.3, df= 4 (P= 0.37); ", I^2, "= 12.3%")), digits= 3, 
       rows= 8:4, ylim= c(1.5, 11), xlim= c(-8, 8))
addpoly(fe.1, row= 3, mlab= "FE Model", digits = 3)
addpoly(re.1, row= 2, mlab= "RE Model", digits = 3)
text(x= -8, y= 9.2, "Artificial data 1", cex= 2, adj= c(0, 0))
#text(x= -5.4, y= 0.25, labels= expression(paste("Heterogeneity: Q= 4.3, df= 4 (P= 0.37); ", I^2, "= 12.3%")), pos= 3, adj= c(0, 0))


forest(hypo.data.2$yi, hypo.data.2$vi, xlab= expression(paste("Heterogeneity: Q= 8.6, df= 4 (P= 0.07); ", I^2, "= 51.7%")), digits= 3, 
       rows= 8:4, ylim= c(1.5, 11), xlim= c(-8, 8), at= seq(from= -4, to= 4, by= 2))
addpoly(fe.2, row= 3, mlab= "FE Model", digits = 3)
addpoly(re.2, row= 2, mlab= "RE Model", digits = 3)
text(x= -8, y= 9.2, "Artificial data 2", cex= 2, adj= c(0, 0))
#text(x= -5.15, y= 0.25, labels= expression(paste("Heterogeneity: Q= 8.6, df= 4 (P= 0.07); ", I^2, "= 51.7%")), pos= 3, adj= c(0, 0))
dev.off()




# analyze the first hypothetical data
y    <- as.numeric(hypo.data.1$yi)
var_y<- as.numeric(hypo.data.1$vi)
nobs <- length(unique(hypo.data.1$study))

# initializaton
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
               fixed.effect= (lambda==0)
               # model checking
               #dev= rep(NA, 1), 
               #p.inv= rep(NA, nobs),
               #y.rep= rep(NA, nobs),
               #p.smaller= rep(NA, nobs),
               #dev.rep= rep(NA)
)

set.seed(145)
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
               fixed.effect= (lambda==0)
               # model checking
               #dev= rep(NA, 1), 
               #p.inv= rep(NA, nobs),
               #y.rep= rep(NA, nobs),
               #p.smaller= rep(NA, nobs),
               #dev.rep= rep(NA)
)

set.seed(345)
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
               fixed.effect= (lambda==0)
               # model checking
               #dev= rep(NA, 1), 
               #p.inv= rep(NA, nobs),
               #y.rep= rep(NA, nobs),
               #p.smaller= rep(NA, nobs),
               #dev.rep= rep(NA)
)
parms.ini<- list(parms.1, parms.2, parms.3)


set.seed(43234)
mcmc.out.lst.1<- bayesian.mdl.uncertainty.meta.analysis(y= y, 
                                                      var_y= var_y, 
                                                      parms.ini= parms.ini, 
                                                      mu_beta= 0, 
                                                      var_beta= 100, 
                                                      alpha_0= .5, 
                                                      theta_0= .5 * nobs, 
                                                      var_lambda_0= var.unit.info(y= re.1$yi, var_y= re.1$vi, var_s= re.1$tau2), 
                                                      prob.random= .5,
                                                      nthin= 10, 
                                                      nburn= 5000, 
                                                      niter= 10000, 
                                                      print.progress= FALSE) 
# Bayesian fixed model specification for dataset 1
set.seed(43234)
mcmc.out.fixed.1<- bayesian.mdl.uncertainty.meta.analysis(y= y, 
                                                          var_y= var_y, 
                                                          parms.ini= parms.ini, 
                                                          mu_beta= 0, 
                                                          var_beta= 100, 
                                                          alpha_0= .5, 
                                                          theta_0= .5 * nobs, 
                                                          var_lambda_0= var.unit.info(y= re.1$yi, var_y= re.1$vi, var_s= re.1$tau2), 
                                                          prob.random= 0,
                                                          nthin= 10, 
                                                          nburn= 5000, 
                                                          niter= 10000, 
                                                          print.progress= FALSE) 

# Bayesian random model specification for dataset 1
set.seed(2134)
mcmc.out.randm.1<- bayesian.mdl.uncertainty.meta.analysis(y= y, 
                                                          var_y= var_y, 
                                                          parms.ini= parms.ini, 
                                                          mu_beta= 0, 
                                                          var_beta= 100, 
                                                          alpha_0= .5, 
                                                          theta_0= .5 * nobs, 
                                                          var_lambda_0= var.unit.info(y= re.1$yi, var_y= re.1$vi, var_s= re.1$tau2), 
                                                          prob.random= 1,
                                                          nthin= 10, 
                                                          nburn= 5000, 
                                                          niter= 10000, 
                                                          print.progress= FALSE) 
save.image("artificial_data_ex.RData")


plot(mcmc.out.lst.1[, c("beta", "sd.random", "fixed.effect")])
gelman.diag(mcmc.out.lst.1[, c("beta", "sd.random", "fixed.effect")]) 

summary(mcmc.out.lst.1[, c("beta", "sd.random", "fixed.effect")])
HPDinterval(mcmc.out.lst.1[, c("beta", "sd.random", "fixed.effect")])
summary(fe.1); plot(fe.1)
summary(re.1)

summary(mcmc.out.fixed.1[, c("beta", "sd.random", "fixed.effect")])
HPDinterval(mcmc.out.fixed.1[, c("beta", "sd.random", "fixed.effect")])


summary(mcmc.out.randm.1[, c("beta", "sd.random", "fixed.effect")])
HPDinterval(mcmc.out.randm.1[, c("beta", "sd.random", "fixed.effect")])











# analyze the second hypothetical data
y    <- as.numeric(hypo.data.2$yi)
var_y<- as.numeric(hypo.data.2$vi)
nobs <- length(unique(hypo.data.2$study))

# initializaton
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
               fixed.effect= (lambda==0)
               # model checking
               #dev= rep(NA, 1), 
               #p.inv= rep(NA, nobs),
               #y.rep= rep(NA, nobs),
               #p.smaller= rep(NA, nobs),
               #dev.rep= rep(NA)
)

set.seed(145)
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
               fixed.effect= (lambda==0)
               # model checking
               #dev= rep(NA, 1), 
               #p.inv= rep(NA, nobs),
               #y.rep= rep(NA, nobs),
               #p.smaller= rep(NA, nobs),
               #dev.rep= rep(NA)
)

set.seed(345)
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
               fixed.effect= (lambda==0)
               # model checking
               #dev= rep(NA, 1), 
               #p.inv= rep(NA, nobs),
               #y.rep= rep(NA, nobs),
               #p.smaller= rep(NA, nobs),
               #dev.rep= rep(NA)
)
parms.ini<- list(parms.1, parms.2, parms.3)


set.seed(43234)
mcmc.out.lst.2<- bayesian.mdl.uncertainty.meta.analysis(y= y, 
                                                      var_y= var_y, 
                                                      parms.ini= parms.ini, 
                                                      mu_beta= 0, 
                                                      var_beta= 100, 
                                                      alpha_0= .5, 
                                                      theta_0= .5 * nobs, 
                                                      var_lambda_0= var.unit.info(y= re.2$yi, var_y= re.2$vi, var_s= re.2$tau2), 
                                                      prob.random= .5,
                                                      nthin= 10, 
                                                      nburn= 5000, 
                                                      niter= 10000, 
                                                      print.progress= FALSE) 
set.seed(43234)
mcmc.out.fixed.2<- bayesian.mdl.uncertainty.meta.analysis(y= y, 
                                                        var_y= var_y, 
                                                        parms.ini= parms.ini, 
                                                        mu_beta= 0, 
                                                        var_beta= 100, 
                                                        alpha_0= .5, 
                                                        theta_0= .5 * nobs, 
                                                        var_lambda_0= var.unit.info(y= re.2$yi, var_y= re.2$vi, var_s= re.2$tau2), 
                                                        prob.random= 0,
                                                        nthin= 10, 
                                                        nburn= 5000, 
                                                        niter= 10000, 
                                                        print.progress= FALSE)
set.seed(43234)
mcmc.out.randm.2<- bayesian.mdl.uncertainty.meta.analysis(y= y, 
                                                        var_y= var_y, 
                                                        parms.ini= parms.ini, 
                                                        mu_beta= 0, 
                                                        var_beta= 100, 
                                                        alpha_0= .5, 
                                                        theta_0= .5 * nobs, 
                                                        var_lambda_0= var.unit.info(y= re.2$yi, var_y= re.2$vi, var_s= re.2$tau2), 
                                                        prob.random= 1,
                                                        nthin= 10, 
                                                        nburn= 5000, 
                                                        niter= 10000, 
                                                        print.progress= FALSE) 
save.image("artificial_data_ex.RData")


plot(mcmc.out.lst.2[, c("beta", "sd.random", "fixed.effect")])
gelman.diag(mcmc.out.lst.2[, c("beta", "sd.random", "fixed.effect")]) 

summary(mcmc.out.lst.2[, c("beta", "sd.random", "fixed.effect")])
HPDinterval(mcmc.out.lst.2[, c("beta", "sd.random", "fixed.effect")])
summary(re.2); forest(re.2, digits=3)


post.fixed.mdl.prob<- summary(mcmc.out.lst.2[, "fixed.effect"])$statistics[1]
post.fixed.mdl.prob/(1-post.fixed.mdl.prob)

save.image("artificial_data_ex.RData")
post.fixed.mdl.prob<- summary(mcmc.out.lst.1[, "fixed.effect"])$statistics[1]
