#=======================================================
#Sgen functions adapted from samSim
#
#
#=======================================================



Sgencompute <- function(S, a,b, Smsy ) {
	#modified from samsim sGenOptimSmsyum
  
  prt <- S * exp(a - b * S)
  epsilon <- log(Smsy) - log(prt)
  nLogLike <- sum(dnorm(epsilon, 0, 0.01, log = T))
  return(list(nSS = sum(nLogLike)))
  
}



sGenSolver <- function(a,b, Smsy) {
  #gives the min Ricker log-likelihood
  if(a>0){
    fnSGen <- function(S, a, b, Smsy) -1.0 * Sgencompute(S, a, b, Smsy)$nSS
    fit <- optimize(f = fnSGen, interval = c(0, ((a / b) * (0.5 - 0.07 * a))),
                 a=a,b=b, Smsy = Smsy)
  }else{
    fit <-list(minimum=NA)
  }

  return(list(fit = fit$minimum))
}