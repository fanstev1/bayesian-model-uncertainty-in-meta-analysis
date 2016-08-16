library(metafor)
setwd("I:/Users/sfan/Google Drive/MSc HSR thesis work/program/bayesian meta analysis")

meta.data<- data.frame(study= 1:6,
                       studyName= c("NAIA3005", "NAI30010", "NAIA/B2009", "WV15673", "WV15697", "WV15799"), 
                       nc= c(554, 423, 144, 268, 251, 462), 
                       rc= c(34, 40, 9, 19, 6, 34),
                       nt= c(553, 414, 144, 268, 252, 493),
                       rt= c(11, 7, 3, 3, 3, 4))

# assess the effect measure 
# ai - number of events in treated group
# n1i- number of subjects in treated group
# ci - number of events in control group
# n2i- number of subjects in control gorup
eff.summary<- escalc(measure= "OR",
                     ai= meta.data$rt,
                     n1i= meta.data$nt,
                     ci= meta.data$rc,
                     n2i= meta.data$nc,
                     data= meta.data)

pooled.res.fe<- rma(method= "FE",
                    measure= "OR",
                    ai= rt,
                    n1i= nt,
                    ci= rc,
                    n2i= nc,
                    data= eff.summary)

pooled.res.re<- rma(method= "REML",
                    measure= "OR",
                    ai= rt,
                    n1i= nt,
                    ci= rc,
                    n2i= nc, 
                    data= eff.summary)

summary(pooled.res.fe)
summary(pooled.res.re)

forest(pooled.res.fe, atransf = exp)
forest(pooled.res.re, atransf = exp)
funnel(pooled.res.fe)
