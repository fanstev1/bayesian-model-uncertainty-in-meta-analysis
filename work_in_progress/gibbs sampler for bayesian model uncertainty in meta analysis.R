library(MCMCpack)
library(mvtnorm)
library(truncnorm)

bayesian.mdl.uncertainty.meta.analysis.subgrp<- 
  function(y, var_y, subgrp, parms.ini, 
           mu_beta, var_beta, 
           alpha_0, theta_0, var_lambda_0, prob.random,
           nthin= 1, nburn, niter, nchain= length(parms.ini), print.progress= FALSE) {
    ## model Y_i = mu + lambda * b_i + epislon_i
    ## prior density: 
    ##    beta ~ N(mu_beta, var_beta), interpreted as overall
    ##    lambda ~ ZI-N(0, var_lambda_0, prob.random)
    ##    b ~ MVN(mu_b_0, var_b), var_b ~ inverse gamma shape (alpha) and scale (theta)
    
    mcmc.out.lst<- vector("list", nchain)
    grpID<- unique(subgrp)
    nsubgrp<- length(grpID) 
    iter.idx<- 0
    while (iter.idx <= (niter+nburn)*nthin){
      for (k in 1:nchain) {
        parms<- parms.ini[[k]]
        
        for (j in grpID) {
          parms$beta[j]   <- gs_cd_beta(b= parms$b[which(subgrp==j)], 
                                        lambda= parms$lambda[j], 
                                        mu_beta= mu_beta[j], 
                                        var_beta= var_beta[j], 
                                        y= y[which(subgrp==j)], 
                                        var_y= var_y[which(subgrp==j)])
          parms$b[which(subgrp==j)]<- gs_cd_bi(beta= parms$beta[j], lambda= parms$lambda[j], var_b= parms$var_b[j], y= y[which(subgrp==j)], var_y= var_y[which(subgrp==j)])
          parms$var_b[j]     <- gs_cd_var_b(b= parms$b[which(subgrp==j)], alpha_0= alpha_0[j], theta_0= theta_0[j])
          parms$lambda[j]    <- gs_cd_lambda.norm(mu= parms$beta[j], b=parms$b[which(subgrp==j)], var_lambda_0= var_lambda_0[j], prob.random= prob.random[j], y= y[which(subgrp==j)], var_y= var_y[which(subgrp==j)])
          #parms$lambda<- gs_cd_lambda.truncnorm(mu= parms$beta, b=parms$b, var_lambda_0= var_lambda_0, prob.random= prob.random, y= y, var_y= var_y)
          
          
          ## additional parameters to keep track 
          parms$sd.random[j]<- abs(parms$lambda[j]) * sqrt(parms$var_b[j])
          parms$or[j]<- exp(parms$beta[j])
          parms$fixed.effect[j]<- (parms$lambda[j]==0)
          
        }
        parms$diff.subgrp<- parms$beta[2] - parms$beta[1]
        parms$or.subgrp<- exp(parms$diff.subgrp)
        parms$effsmaller.subgrp<- (parms$or.subgrp > 1)
        ## deviance 
        #parms$dev       <- -2*dmvnorm(y, mean= parms$beta + parms$lambda * parms$b, sigma= diag(var_y), log= TRUE )
        #parms$p.inv     <- 1/dnorm(y, mean= parms$beta + parms$lambda * parms$b, sd= sqrt(var_y))
        
        # for model checking
        #parms$y.rep     <- rmvnorm(1, mean= parms$beta + parms$lambda * parms$b, sigma= diag(var_y))
        #parms$dev.rep   <- -2*dmvnorm(parms$y.rep, mean= parms$beta + parms$lambda * parms$b, sigma= diag(var_y), log= TRUE )
        #parms$p.smaller <- (y< parms$y.rep)  
        
        parms.ini[[k]]<- parms
        
        if (iter.idx> nburn * nthin & iter.idx %% nthin == 0) {
          mcmc.out.lst[[k]]<- rbind(mcmc.out.lst[[k]], unlist(parms))
          
        }
      }
      if (print.progress) print(iter.idx)
      iter.idx<- 1+iter.idx
    }
    for (k in 1:nchain) { mcmc.out.lst[[k]]<- mcmc(mcmc.out.lst[[k]], start= 1+nburn, thin= nthin) }
    return(mcmc.out.lst<- mcmc.list(mcmc.out.lst))
  }



bayesian.mdl.uncertainty.meta.analysis<- 
  function(y, var_y, parms.ini, 
           mu_beta, var_beta, alpha_0= .5, theta_0= .5*length(y), var_lambda_0= 25, prob.random= .5,
           nthin= 1, nburn, niter, nchain= length(parms.ini), print.progress= FALSE) {
    ## model Y_i = mu + lambda * b_i + epislon_i
    ## prior density: 
    ##    beta ~ N(mu_beta, var_beta), interpreted as overall
    ##    lambda ~ ZI-N(0, var_lambda_0, prob.random)
    ##    b ~ MVN(mu_b_0, var_b), var_b ~ inverse gamma shape (alpha) and scale (theta)
    
    # prior hyperparameter for beta 
    #mu_beta<- 0        # prior mean for beta
    #var_beta<- 10000   # prior variance for beta
    
    # prior hyperparameter for var_b
    # b_i is normally distributed with zero mean and var_b
    # var_b has inverse gamma distribution with shape (alpha_0) and scale (theta_0)
    #alpha_0<- .5
    #theta_0<- .5*nobs   # scale
    
    # prior hyperparameter for lambda
    #var_lambda_0<- 25   # scale
    #prob.random<- .5         # probabily of presence of random effect (i.e. )
    mcmc.out.lst<- vector("list", nchain)
    
    iter.idx<- 0
    while (iter.idx <= (niter+nburn)*nthin){
      for (k in 1:nchain) {
        parms<- parms.ini[[k]]
        
        ######## core of gibbs sampling #############
        parms$beta  <- gs_cd_beta(b= parms$b, lambda= parms$lambda, mu_beta= mu_beta, var_beta= var_beta, y= y, var_y= var_y)
        parms$b     <- gs_cd_bi(beta= parms$beta, lambda= parms$lambda, var_b= parms$var_b, y= y, var_y= var_y)
        parms$var_b <- gs_cd_var_b(b= parms$b, alpha_0= alpha_0, theta_0= theta_0)
        parms$lambda<- gs_cd_lambda.norm(mu= parms$beta, b=parms$b, var_lambda_0= var_lambda_0, prob.random= prob.random, y= y, var_y= var_y)
        #############################################
        
        
        
        
        ## additional parameters to keep track 
        parms$s         <- parms$b * parms$lambda
        parms$sd.random <- abs(parms$lambda) * sqrt(parms$var_b)
        parms$or        <- exp(parms$beta)
        parms$fixed.effect<- (parms$lambda==0)
        
        ## deviance 
       # parms$dev       <- -2*dmvnorm(y, mean= parms$beta + parms$lambda * parms$b, sigma= diag(var_y), log= TRUE )
       # parms$p.inv     <- 1/dnorm(y, mean= parms$beta + parms$lambda * parms$b, sd= sqrt(var_y))
        
        # for model checking
        #parms$y.rep     <- rmvnorm(1, mean= parms$beta + parms$lambda * parms$b, sigma= diag(var_y))
        #parms$dev.rep   <- -2*dmvnorm(parms$y.rep, mean= parms$beta + parms$lambda * parms$b, sigma= diag(var_y), log= TRUE )
        #parms$p.smaller <- (y< parms$y.rep)  
        
        parms.ini[[k]]<- parms
        
        if (iter.idx> nburn * nthin & iter.idx %% nthin == 0) {
          mcmc.out.lst[[k]]<- rbind(mcmc.out.lst[[k]], unlist(parms))
          
        }
      }
      if (print.progress) print(iter.idx)
      iter.idx<- 1+iter.idx
    }
    
    
    for (k in 1:nchain) { 
      mcmc.out.lst[[k]]<- mcmc(mcmc.out.lst[[k]], start= 1+nburn, thin= nthin) 
    }
    return(mcmc.out.lst<- mcmc.list(mcmc.out.lst))
}



bayesian.mdl.uncertainty.meta.analysis.unit<- 
  function(y, var_y, parms.ini, 
           mu_beta, var_beta, var_lambda_0= 25, prob.random= .5,
           nthin= 1, nburn, niter, nchain= length(parms.ini), print.progress= FALSE) {
    ## model Y_i = mu + lambda * b_i + epislon_i
    ## prior density: 
    ##    beta ~ N(mu_beta, var_beta), interpreted as overall
    ##    lambda ~ ZI-N(0, var_lambda_0, prob.random)
    ##    b ~ MVN(0, 1), var_b ~ inverse gamma shape (alpha) and scale (theta)
    
    # prior hyperparameter for beta 
    #mu_beta<- 0        # prior mean for beta
    #var_beta<- 10000   # prior variance for beta
    
    # prior hyperparameter for var_b
    # b_i is normally distributed with zero mean and var_b
    # var_b has inverse gamma distribution with shape (alpha_0) and scale (theta_0)
    #alpha_0<- .5
    #theta_0<- .5*nobs   # scale
    
    # prior hyperparameter for lambda
    #var_lambda_0<- 25   # scale
    #prob.random<- .5         # probabily of presence of random effect (i.e. )
    mcmc.out.lst<- vector("list", nchain)
    
    iter.idx<- 0
    while (iter.idx <= (niter+nburn)*nthin){
      for (k in 1:nchain) {
        parms<- parms.ini[[k]]
        
        parms$beta  <- gs_cd_beta(b= parms$b, lambda= parms$lambda, mu_beta= mu_beta, var_beta= var_beta, y= y, var_y= var_y)
        parms$b     <- gs_cd_bi(beta= parms$beta, lambda= parms$lambda, var_b= 1, y= y, var_y= var_y)
        #parms$var_b <- gs_cd_var_b(b= parms$b, alpha_0= alpha_0, theta_0= theta_0)
        #parms$lambda<- gs_cd_lambda.norm(mu= parms$beta, b=parms$b, var_lambda_0= var_lambda_0, prob.random= prob.random, y= y, var_y= var_y)
        parms$lambda<- gs_cd_lambda.truncnorm(mu= parms$beta, b=parms$b, var_lambda_0= var_lambda_0, prob.random= prob.random, y= y, var_y= var_y)
        
        
        ## additional parameters to keep track 
        parms$s         <- parms$b * parms$lambda
        parms$sd.random <- abs(parms$lambda) * sqrt(parms$var_b)
        parms$or        <- exp(parms$beta)
        parms$fixed.effect<- (parms$lambda==0)
        ## deviance 
        #parms$dev       <- -2*dmvnorm(y, mean= parms$beta + parms$lambda * parms$b, sigma= diag(var_y), log= TRUE )
        #parms$p.inv     <- 1/dnorm(y, mean= parms$beta + parms$lambda * parms$b, sd= sqrt(var_y))
        
        # for model checking
        #parms$y.rep     <- rmvnorm(1, mean= parms$beta + parms$lambda * parms$b, sigma= diag(var_y))
        #parms$dev.rep   <- -2*dmvnorm(parms$y.rep, mean= parms$beta + parms$lambda * parms$b, sigma= diag(var_y), log= TRUE )
        #parms$p.smaller <- (y< parms$y.rep)  
        
        parms.ini[[k]]<- parms
        
        if (iter.idx> nburn * nthin & iter.idx %% nthin == 0) {
          mcmc.out.lst[[k]]<- rbind(mcmc.out.lst[[k]], unlist(parms))
          
        }
      }
      if (print.progress) print(iter.idx)
      iter.idx<- 1+iter.idx
    }
    for (k in 1:nchain) { mcmc.out.lst[[k]]<- mcmc(mcmc.out.lst[[k]], start= 1+nburn, thin= nthin) }
    return(mcmc.out.lst<- mcmc.list(mcmc.out.lst))
  }

# conditional density for mu given Y, lambda, bi
gs_cd_beta<- function(mu_beta, var_beta, y, var_y, b, lambda) {
  #nobs<- length(y)
  # update hyperparameter
  var_beta.update<- (sum(var_y^(-1)) + var_beta^(-1))^(-1)
  mu_beta.update<- var_beta.update * sum((y- lambda*b)/var_y)
  mu.update<- rnorm(1, mean= mu_beta.update, sd= sqrt(var_beta.update))
  
  return(as.numeric(mu.update))
}
# parms$mu<- gs_cd_mu(b= parms$b, lambda= parms$lambda, mu_beta, var_beta, y, var_y)


# conditional density for b given Y, mu, lambda, bi, var_b
gs_cd_bi<- function(y, var_y, beta, lambda, var_b) {
  #nobs<- length(y)
  # update hyperparameter
  var_b.update<- ((lambda/sqrt(var_y))^2 + var_b^(-1))^(-1)
  mu_b.update<- var_b.update * (lambda * (y-beta)/var_y)
  b.update<- rmvnorm(1, mean= mu_b.update, sigma= diag(var_b.update) )  

  return(as.numeric(b.update))
}
# parms$b<- gs_cd_bi(y, var_y, mu= parms$mu, lambda= parms$lambda, var_b= parms$var_b)

# conditional density for var_b given Y, mu, lambda, bi
gs_cd_var_b<- function(b, alpha_0, theta_0) {
  return( rinvgamma(1, shape= alpha_0+.5*length(b), scale= theta_0 + .5*sum(b^2)) )
}
# parms$var_b<- gs_cd_var_b(b= parms$b)

# conditional density for lambda given Y, mu, lambda, bi
# prior density for lambda is normal
gs_cd_lambda.norm<- function(y, var_y, mu, b, var_lambda_0, prob.random){
  
  var_lambda<- (sum((b/sqrt(var_y))^2) + var_lambda_0^-1)^-1
  mu_lambda<- var_lambda * sum( b * (y - mu)/var_y)
  prob.random.update<- prob.random * dnorm(0, mean= 0, sd= sqrt(var_lambda_0)) /
    ( prob.random * dnorm(0, mean= 0, sd= sqrt(var_lambda_0)) + 
        (1-prob.random) * dnorm(0, mean= mu_lambda, sd= sqrt(var_lambda)) )
  
  #m<- rbinom(1, size= 1, prob= prob.random.update)
  lambda.update<- rbinom(1, size= 1, prob= prob.random.update) * 
    rnorm(1, mean= mu_lambda, sd= sqrt(var_lambda) )
  
  return(as.numeric(lambda.update))
}
# gs_cd_lambda(y, var_y, mu= parms$mu, b=parms$b)



# conditional density for lambda given Y, mu, lambda, bi
# prior density is a zero-inflated truncated normal with zero mean and variance var_lambda_0
gs_cd_lambda.truncnorm<- function(y, var_y, mu, b, var_lambda_0, prob.random){
  
  var_lambda<- (sum((b/sqrt(var_y))^2) + var_lambda_0^-1)^-1
  mu_lambda<- var_lambda * sum( b * (y - mu)/var_y)
  prob.random.update<- prob.random * dtruncnorm(0, a= 0, mean= 0, sd= sqrt(var_lambda_0)) /
    ( prob.random * dtruncnorm(0, a= 0, mean= 0, sd= sqrt(var_lambda_0)) + 
        (1-prob.random) * dtruncnorm(0, a= 0, mean= mu_lambda, sd= sqrt(var_lambda)) )
  
  #  m<- rbinom(1, size= 1, prob= prob.random.update)
  lambda.update<- rbinom(1, size= 1, prob= prob.random.update) * 
    rtruncnorm(n= 1, a=0, mean = mu_lambda, sd = sqrt(var_lambda))
  
  return(as.numeric(lambda.update))
}
