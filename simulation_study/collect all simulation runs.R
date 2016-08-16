library(MCMCpack)
library(mvtnorm)
library(truncnorm)
library(metafor)
#setwd("C:/Users/sfan/Downloads/")
setwd("Z:/Steve tmp/R sessions simulation")
load("sim_reduced.RData")


mcmc.res.lst.original<- mcmc.res.lst
HPDinterval.lst.original<- HPDinterval.lst


load("sim_reduced_5.RData")
for (j in 100001:200000) {
  mcmc.res.lst.original[[j]]<- mcmc.res.lst[[j]]
  HPDinterval.lst.original[[j]]<- HPDinterval.lst[[j]]
}


load("sim_reduced_4.RData")
for (j in 200001:300000) {
  mcmc.res.lst.original[[j]]<- mcmc.res.lst[[j]]
  HPDinterval.lst.original[[j]]<- HPDinterval.lst[[j]]
}
for (j in 347001:399999) {
  mcmc.res.lst.original[[j]]<- mcmc.res.lst[[j]]
  HPDinterval.lst.original[[j]]<- HPDinterval.lst[[j]]
}
# 400000 is missing

load("sim_reduced_6.RData")
for (j in 300001:347000) {
  mcmc.res.lst.original[[j]]<- mcmc.res.lst[[j]]
  HPDinterval.lst.original[[j]]<- HPDinterval.lst[[j]]
}

load("sim_reduced_2.RData")
for (j in 400000:500000) {
  mcmc.res.lst.original[[j]]<- mcmc.res.lst[[j]]
  HPDinterval.lst.original[[j]]<- HPDinterval.lst[[j]]
}
for (j in 547001:599999) {
  mcmc.res.lst.original[[j]]<- mcmc.res.lst[[j]]
  HPDinterval.lst.original[[j]]<- HPDinterval.lst[[j]]
}

load("sim_reduced_7.RData")
for (j in 500001:547000) {
  mcmc.res.lst.original[[j]]<- mcmc.res.lst[[j]]
  HPDinterval.lst.original[[j]]<- HPDinterval.lst[[j]]
}

load("sim_reduced_3.RData")
for (j in 600000:800000) {
  mcmc.res.lst.original[[j]]<- mcmc.res.lst[[j]]
  HPDinterval.lst.original[[j]]<- HPDinterval.lst[[j]]
}
