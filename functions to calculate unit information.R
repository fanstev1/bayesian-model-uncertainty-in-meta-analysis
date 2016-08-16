beta.fun<- function(y, var_y, var_s){
  wgt<- 1/(var_s + var_y)
  wgt<- wgt/sum(wgt)
  return( sum(y*wgt) )
}

beta.prime<- function(y, var_y, var_s){
  out<- sum( y * (var_s+var_y)^(-1) * sum( (var_s+var_y)^(-2) )/sum( 1/(var_s + var_y) )^2 ) -
                              sum( y*(var_s+var_y)^(-2)/sum( (var_s+var_y)^(-1) ) )
  return( out )
}
beta.2prime<- function(y, var_y, var_s){
  #sum(1/(var_s+var_y))
  out<- 2*sum(y*(var_s+var_y)^(-3))/sum(1/(var_s+var_y)) - 
    2*sum(y*(var_s+var_y)^(-2))*sum((var_s+var_y)^(-2))/sum(1/(var_s+var_y)) -
    2*sum(y/(var_s+var_y))*sum((var_s+var_y)^(-3))/sum(1/(var_s+var_y))^2 +
    2*sum(y/(var_s+var_y))*sum((var_s+var_y)^(-2))^2/sum((var_s+var_y)^(-1))^3
    
  return(out)
}

# first derivative of log-likelihood w.r.t variance
llk.prime.var_s<- function(y, var_y, var_s){
  out<- .5*sum((var_s+var_y)^(-2))/ sum((var_s+var_y)^(-1)) -
          .5*sum((var_s+var_y)^(-1)) -
          .5*sum( 2*(y-beta.fun(y, var_y, var_s))*-beta.prime(y, var_y, var_s)/(var_s+var_y) -
                    ( (y-beta.fun(y, var_y, var_s))/(var_s+var_y) )^2 )
  return(out)
}

# first derivative of log-likelihood w.r.t variance
fisher.info.var_s<- function(y, var_y, var_s){
  
  beta.hat<- beta.fun(y, var_y, var_s)
  beta.prime.hat<- beta.prime(y, var_y, var_s)
  beta.2prime.hat<- beta.2prime(y, var_y, var_s)
  
  tmp<- -4*(y- beta.hat)*(-beta.prime.hat)/(var_s + var_y)^2 +
    2*(y- beta.hat)*(-beta.2prime.hat)/(var_s + var_y) +
    2*beta.prime.hat^2/(var_s+var_y) +
    2*(y- beta.hat)^2/(var_s+var_y)^3
  
  out<- -.5*sum(tmp) - 
    sum((var_s+var_y)^(-3))/sum(1/(var_s+var_y)) +
    0.5*sum((var_s+var_y)^(-2))^2/ sum(1/(var_s+var_y))^2 +
    0.5*sum((var_s+var_y)^(-2))
  
  return( out )
  
}
var.unit.info<- function(y, var_y, var_s){
  wgt<- sum(1/var_y)
  
  #if (var_s==0) var_s<- uniroot(llk.prime.sigma_s, interval= c(-5, 0), y= y, var_y= var_y)
  out<- (-fisher.info.var_s(y, var_y, var_s)/wgt)^(-1)
  return( out )
  
}
#beta.fun(y= y, var_y= var_y, var_s= res.re$tau2)
#beta.prime(y= y, var_y= var_y, var_s= res.re$tau2)
#beta.2prime(y= y, var_y= var_y, var_s= res.re$tau2)
#llk.prime.sigma_s(y= y, var_y= var_y, var_s= res.re$tau2)
#fisher.info.sigma_s(y= y, var_y= var_y, var_s= res.re$tau2)
#var.unit.info(y= y, var_y= var_y, var_s= res.re$tau2)
#var.unit.info(y= c(-0.86, -0.33, -0.47,-0.50, 0.28, -0.04, -0.80, -0.19, -0.49), 
#              var_y= c(0.325732899, 0.307692308, 0.123609394, 0.062972292, 0.29154519, 0.075700227, 0.613496933, 0.017866714, 0.077399381),
#              var_s= -.0078)






####### full likelihood #############
#full.llk.prime<- function(y, var_y, beta, var_s) {
#  mat<- matrix(0, nrow= 2, ncol= 1)
#  mat[1, 1]<- sum( (y- beta)/(var_s + var_y) )
#  mat[2, 1]<- .5 * ( sum( ((y-beta)/(var_s+var_y))^2 ) - sum(1/(var_s+var_y)) )
#  return( mat )
#}
#obs.fisher.info.full<- function(y, var_y, beta, var_s) {
#  mat<- matrix(0, nrow= 2, ncol= 2)
#  mat[1, 1]<- -sum(1/(var_s + var_y))
#  mat[1, 2]<- -sum( (y- beta)/(var_s + var_y)^2)
#  mat[2, 1]<- mat[1, 2]
#  mat[2, 2]<- 0.5*sum((var_s+var_y)^(-2)) - sum((y-beta)^2/(var_s+var_y)^3)
#    
#  return( mat )
#}

#old<- matrix(c(1, .1), nrow= 2, ncol= 1)
#i= 1
#repeat{
#  new<- old - solve(obs.fisher.info.full(y= y, var_y= var_y, beta= old[1, 1], var_s= old[2, 1] )) %*%
#                full.llk.prime(y= y, var_y= var_y, beta= old[1, 1], var_s= old[2, 1])
#  #old<- new
#  if (sum(as.vector(new-old)^2)< 1E-6){ break}
#  else {old<- new; i<- 1+i}
#}







#sigma.unres.mle<- function(y, var_y, var_s.init, threshold= 1E-10){
#  old<- var_s.init
#  
#  repeat {
#    new<- old - llk.prime.sigma_s(y= y, var_y= var_y, var_s= old)/fisher.info.sigma_s(y= y, var_y= var_y, var_s= old)
#    if (abs(new-old)< threshold) {break}
#    else {
#      old<- new
#    }
#  }
#  return(new)
#}
#sigma.unres.mle(y= y, var_y= var_y, var_s.init= .09, threshold= 1E-10)