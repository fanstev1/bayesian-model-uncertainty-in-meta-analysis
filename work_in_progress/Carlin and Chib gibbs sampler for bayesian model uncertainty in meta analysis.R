library(MCMCpack)
library(mvtnorm)
#library(R2WinBUGS)


# Gibbs sampler for Carlin and Chib's method
# Model 1: Y_i ~ N(delta, sigma_i^2)
# Model 2: Y_i ~ N(beta + S_i, sigma_i^2) and S_i ~ N(0, sigma_s^2)

bayes.mdl.uncertainty.meta.analysis<- function(y, var_y, parm.ini, 
                                               niter, nburn, nthin, 
                                               nchain= length(parm.ini),
                                               prior.prob.fixed= .5,
                                               prior.parm.mdl1, pseudo.prior.parm.mdl1,
                                               prior.parm.mdl2, pseudo.prior.parm.mdl2,
                                               print.progress= FALSE) {
  # this is the main program implementing the Gibbs sampler of Carlin and Chib's method
  # parm is list, containing Model index, parm1 (list) and parm2 (list)
  # parm1 is a list variable containing all parameters of model 1
  # parm2 is a list variable containing all parameters of model 2
  
  mcmc.out.lst<- vector("list", nchain)
  
  iter.idx<- 0
  while (iter.idx <= (niter+nburn)*nthin){
    for (k in 1:nchain) {
      parm<- parm.ini[[k]]
    
      # update Model 1 parmeters
      parm$parm1$delta<- mdl1_post_delta(parm1= parm$parm1, 
                                         mdl= parm$mdl.idx, 
                                         y= y, var_y= var_y, 
                                         prior.parm= prior.parm.mdl1, 
                                         pseudo.prior.parm= pseudo.prior.parm.mdl1)
      
      
      # update Model 2 parmeters
      parm$parm2$beta <- mdl2_post_beta(parm2= parm$parm2, 
                                        mdl= parm$mdl.idx, 
                                        y= y, var_y= var_y, 
                                        prior.parm= prior.parm.mdl2, 
                                        pseudo.prior.parm= pseudo.prior.parm.mdl2)  
      
      
      parm$parm2$lambda<- mdl2_post_lambda( parm2= parm$parm2, 
                                            mdl= parm$mdl.idx, 
                                            y= y, var_y= var_y, 
                                            prior.parm= prior.parm.mdl2, 
                                            pseudo.prior.parm= pseudo.prior.parm.mdl2)
      
      parm$parm2$b    <- mdl2_post_b( parm2= parm$parm2, 
                                      mdl= parm$mdl.idx, 
                                      y= y, var_y= var_y, 
                                      prior.parm= prior.parm.mdl2, 
                                      pseudo.prior.parm= pseudo.prior.parm.mdl2)
  
      
      parm$parm2$var_b<- mdl2_post_var_b( parm2= parm$parm2, 
                                          mdl= parm$mdl.idx, 
                                          y= y, var_y= var_y, 
                                          prior.parm= prior.parm.mdl2, 
                                          pseudo.prior.parm= pseudo.prior.parm.mdl2)
      
      
      # update model probability
      parm$mdl.idx  <- gs_post_mdl(parm= parm, y= y, var_y= var_y, prior.prob.mdl1= prior.prob.fixed, 
                                   prior.parm.mdl1= prior.parm.mdl1, 
                                   prior.parm.mdl2= prior.parm.mdl2, 
                                   pseudo.prior.parm.mdl1= pseudo.prior.parm.mdl1, 
                                   pseudo.prior.parm.mdl2= pseudo.prior.parm.mdl2)
      
      parm.ini[[k]]<- parm
      
      if (iter.idx> nburn * nthin & iter.idx %% nthin == 0) {
        parm<- unlist(parm)
        parm2.s<- as.numeric(parm[which(names(parm) %in% paste("parm2.b", 1:nobs, sep= ""))] * parm["parm2.lambda"])
        names(parm2.s)<- paste("parm2.s", 1:nobs, sep= "")
        parm2.var_s<- as.numeric(parm[which(names(parm)=="parm2.var_b")] * abs(parm["parm2.lambda"]))
        names(parm2.var_s)<- "parm2.var_s"
        parm<- c(parm, parm2.s, parm2.var_s)
        
        mcmc.out.lst[[k]]<- rbind(mcmc.out.lst[[k]], parm)
      }
    }
    
    if (print.progress) print(iter.idx)
    iter.idx<- 1+iter.idx
  }
  
  for (k in 1:nchain) { mcmc.out.lst[[k]]<- mcmc(mcmc.out.lst[[k]], start= 1+nburn, thin= nthin) }
  return(mcmc.out.lst<- mcmc.list(mcmc.out.lst))
  #plot( mcmc.out.lst<- mcmc(mcmc.out.lst)  )
  #table(mcmc.out.lst[,"mdl.idx"])
}


gs_post_mdl<- function(parm, y, var_y, prior.prob.mdl1= .5, 
                       prior.parm.mdl1, pseudo.prior.parm.mdl1, 
                       prior.parm.mdl2, pseudo.prior.parm.mdl2) {
  
  log.marginal.lik.ratio<- 
    # given model 2
    dmvnorm(y, mean= parm$parm2$beta + parm$parm2$lambda * parm$parm2$b, sigma= diag(var_y), log= TRUE) + 
              # pseudo prior for Model 1 parameters when model 2 is assigned
              dnorm(parm$parm1$delta, mean= pseudo.prior.parm.mdl1$delta$mu, sd= sqrt(pseudo.prior.parm.mdl1$delta$var), log= TRUE) + 
    dnorm(parm$parm2$beta, mean= prior.parm.mdl2$beta$mu, sd= sqrt(prior.parm.mdl2$beta$var), log= TRUE) +
    dnorm(parm$parm2$lambda, mean= prior.parm.mdl2$lambda$mu, sd= sqrt(prior.parm.mdl2$lambda$var), log= TRUE) +
    sum(dnorm(parm$parm2$b, mean= 0, sd= sqrt(parm$parm2$var_b), log= TRUE) ) +
    log(dinvgamma(parm$parm2$var_b, shape= prior.parm.mdl2$var_b$alpha , scale= prior.parm.mdl2$var_b$theta)) +
    log(1-prior.prob.mdl1) -
    # given model 1
    ( dmvnorm(y, mean= rep(parm$parm1$delta, length(y)), sigma= diag(var_y), log= TRUE) +
        dnorm(parm$parm1$delta, mean= prior.parm.mdl1$delta$mu, sd= sqrt(prior.parm.mdl1$delta$var), log= TRUE) +
        # pseudo prior for Model 1 parameters when model 2 is assigned
        dnorm(parm$parm2$beta, mean= pseudo.prior.parm.mdl2$beta$mu, sd= sqrt(pseudo.prior.parm.mdl2$beta$var), log= TRUE) +
        dnorm(parm$parm2$lambda, mean= pseudo.prior.parm.mdl2$lambda$mu, sd= sqrt(pseudo.prior.parm.mdl2$lambda$var), log= TRUE) +
        sum(dnorm(parm$parm2$b, mean= pseudo.prior.parm.mdl2$b$mu, sd= sqrt(pseudo.prior.parm.mdl2$b$var), log= TRUE)) +
        log(dinvgamma(parm$parm2$var_b, shape= pseudo.prior.parm.mdl2$var_b$alpha, scale= pseudo.prior.parm.mdl2$var_b$theta)) +
        log(prior.prob.mdl1))
  
  prob.mdl1.update<- 1/(1+exp(log.marginal.lik.ratio)) # probably of model 1 (fixed-effect)
  
  return( 1+rbinom(n= 1, size= 1, prob= 1-prob.mdl1.update) )
}


mdl2_post_var_b<- function(parm2, mdl, y, var_y, prior.parm, pseudo.prior.parm){
  # parm.mdl2 is a list variable containing all parameters of model 2
  # prior.parm is a list containing all hyperparameters of beta when model= 1
  # pseudo.prior.parm is a list containing all hyperparameters of beta when model= 2
  
  if (mdl==2) {
    alpha.update<- prior.parm$var_b$alpha + 0.5 * length(parm2$b)
    theta.update<- prior.parm$var_b$theta + 0.5 * sum(parm2$b^2)
  }
  else {
    alpha.update<- pseudo.prior.parm$var_b$alpha
    theta.update<- pseudo.prior.parm$var_b$theta
  }
  return( rinvgamma(1, shape= alpha.update, scale= theta.update) )
}

mdl2_post_b<- function(parm2, mdl, y, var_y, prior.parm, pseudo.prior.parm){
  # parm.mdl2 is a list variable containing all parameters of model 2
  # prior.parm is a list containing all hyperparameters of beta when model= 1
  # pseudo.prior.parm is a list containing all hyperparameters of beta when model= 2
  
  if (mdl==2) {
    var.b.update  <- ( 1/parm2$var_b + parm2$lambda^2/var_y)^(-1)
    mu.b.update   <- var.b.update * parm2$lambda * (y - parm2$beta)/var_y 
  }
  else {
    var.b.update  <- pseudo.prior.parm$b$var
    mu.b.update   <- pseudo.prior.parm$b$mu
  }
  return( as.numeric( rmvnorm(n= 1, mean= mu.b.update, sigma= diag(var.b.update)) ) )
}

mdl2_post_lambda<- function(parm2, mdl, y, var_y, prior.parm, pseudo.prior.parm){
  # parm.mdl2 is a list variable containing all parameters of model 2
  # prior.parm is a list containing all hyperparameters of beta when model= 1
  # pseudo.prior.parm is a list containing all hyperparameters of beta when model= 2
  
  if (mdl==2) {
    sigma.lambda.update<- ( 1/prior.parm$lambda$var + sum(parm2$b^2/var_y) )^(-1/2)
    mu.lambda.update   <- sigma.lambda.update^2 * sum(parm2$b * (y - parm2$beta)/var_y)
  }
  else {
    sigma.lambda.update<- sqrt(pseudo.prior.parm$b$var)
    mu.lambda.update   <- pseudo.prior.parm$b$mu
  }
  return( rnorm(n= 1, mean= mu.lambda.update, sd= sigma.lambda.update) )
}

mdl2_post_beta<- function(parm2, mdl, y, var_y, prior.parm, pseudo.prior.parm){
  # parm.mdl2 is a list variable containing all parameters of model 2
  # prior.parm is a list containing all hyperparameters of beta when model= 1
  # pseudo.prior.parm is a list containing all hyperparameters of beta when model= 2
  
  if (mdl==2) {
    sigma.beta.update<- ( 1/prior.parm$beta$var + sum(1/var_y) )^(-1/2)
    mu.beta.update   <- sigma.beta.update^2 * sum((y - parm2$lambda* parm2$b)/var_y) 
  }
  else {
    sigma.beta.update<- sqrt(pseudo.prior.parm$beta$var)
    mu.beta.update   <- pseudo.prior.parm$beta$mu
  }
  return( rnorm(n= 1, mean= mu.beta.update, sd= sigma.beta.update) )
}



mdl1_post_delta<- function(parm1, mdl, y, var_y, prior.parm, pseudo.prior.parm){
  # parm1 is a list variable containing all parameters of model 1
  # prior.parm.mdl1 is a list containing all hyperparameters of delta when model= 1
  # pseudo.prior.parm.mdl1 is a list containing all hyperparameters of delta when model= 2
  
  if (mdl==1) {
    sigma.delta.update<- ( 1/prior.parm$delta$var + sum(1/var_y) )^(-1/2)
    mu.delta.update   <- sigma.delta.update^2 * sum(y/var_y) 
  }
  else {
    sigma.delta.update<- sqrt(pseudo.prior.parm$delta$var)
    mu.delta.update   <- pseudo.prior.parm$delta$mu
  }
  return( rnorm(n= 1, mean= mu.delta.update, sd= sigma.delta.update) )
}

















gs_post_s.mdl2<- function(parm.mdl2, mdl, y, var_y, prior.parm, pseudo.prior.parm){
  # parm.mdl2 is a list variable containing all parameters of model 2
  # prior.parm is a list containing all hyperparameters of beta when model= 1
  # pseudo.prior.parm is a list containing all hyperparameters of beta when model= 2
  
  if (mdl==2) {
    sigma.s.update<- ( (1/parm.mdl2$sigma_s)^2 + 1/var_y)^(-1/2)
    mu.s.update<- sigma.s.update^2 * (y-parm.mdl2$beta)/var_y 
  }
  else {
    sigma.s.update<- pseudo.prior.parm$s$sigma
    mu.s.update<- pseudo.prior.parm$s$mu
  }
  output<- as.numeric(rmvnorm(n= 1, mean= mu.s.update, sigma= diag(sigma.s.update^2)))
  return(output)
}

gs_post_sigma_s.mdl2<- function(parm.mdl2, mdl, y, var_y, prior.parm, pseudo.prior.parm){
  # parm.mdl2 is a list variable containing all parameters of model 2
  # prior.parm is a list containing all hyperparameters of beta when model= 1
  # pseudo.prior.parm is a list containing all hyperparameters of beta when model= 2
  
  if (mdl==2) {
    alpha.update<- prior.parm$sigma_s$alpha + 0.5*length(parm.mdl2$s)
    theta.update<- prior.parm$sigma_s$theta + 0.5*sum(parm.mdl2$s^2)
  }
  else {
    alpha.update<- pseudo.prior.parm$sigma_s$alpha
    theta.update<- pseudo.prior.parm$sigma_s$theta
  }
  output<- sqrt(rinvgamma(1, shape= alpha.update, scale= theta.update ))
  return( output )
}






