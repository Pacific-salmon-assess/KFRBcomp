#==============================================================
#Simple simulation evaluation of Ricker model with time-varying parameters
#Catarina Wor
#May 2022
#==============================================================


library(KFfuncs)
compile("TMBmodels/Ricker_tva_Smax_ratiovar.cpp")
dyn.load(dynlib("TMBmodels/Ricker_tva_Smax_ratiovar"))

#ao<-3
#b<-1/10000
#ER<-0.40
#fec= c(0,.1,.3,.5,.1)
#sig <- .5
#siga <-.3
#
#ps<-0:35000
#pr<-ao*ps*exp(-b*ps)
#sr <- simulateSRrandom()
#plot(sr$S,log(sr$S/sr$R))
#
#par(mfrow=c(2,1))
#plot(sr$S,sr$R, xlim=c(0,35000), ylim=c(0,35000))
#lines(ps,pr,col="red")
#plot(log(sr$a),type="b")




simulateSRrandom <- function(ao=3, b=1/10000, ER=0.4, fec=c(0,0,0,1,0), sig=.5, siga=.2, nobs=40,CapScalar=5 ){

    yrs <- 1:100
    S <- NULL
    R <- NULL
    Robs <- NULL
    a <- NULL
    Smsy <- NULL
    Umsy <- NULL
    Sgen <- NULL
    
    a[1:5] <- ao
    Seq <- ao/b
    S[1] <- Seq
    R[1] <- Seq*exp(ao-b*Seq)
  
    for(y in 2:length(fec)){
      S[y] <- 0    
      for(j in seq_along(fec)){
        if((y-j) > 0){
          S[y] <- S[y] + R[y-j]*(1-ER)*fec[j]
        }else{
          S[y] <- S[y] + R[1]*(1-ER)*fec[j]
        }
      }
      R[y]<-S[y]*exp(a[y]-b*S[y])        
    }
  
    
    for(y in (length(fec)+1):max(yrs)){  
      S[y]<-0
      for(j in seq_along(fec)){      	
        S[y] <- S[y] + R[y-j]*(1-ER)*fec[j]   	
      } 
      #truncated lognrmal error
      #a[y] <- a[y-1] + (qnorm(runif(1, 0.1, 0.9),0,siga))        
      a[y] <- a[y-1] + rnorm(1,0,siga) -(.5*siga^2)
      # avoid alpha that is super negative
      if(a[y]<0.1){
        a[y]<- 0.1
      }

      #R[y] <- S[y]*exp(a[y]-b*S[y]) *exp(qnorm(runif(1, 0.1, 0.95),0,sig))
      R[y] <- S[y]*exp(a[y]-b*S[y]) *exp(rnorm(1,0,sig)-(.5*sig^2))#-(.5*sig^2)
      if(!is.na(R[y]) &  R[y] > Seq*CapScalar){
        R[y] <- Seq*CapScalar
      }
      if(!is.na(R[y]) & R[y]<5){
      	R[y] <- NA
      }
      
    }

    for(n in seq_along(a)){

      if(a[n]>0){
        Smsy[n] <- (1 - gsl::lambert_W0(exp(1 - a[n]))) /b
        Umsy[n] <- .5 * a[n] - 0.07 * a[n]^2
        Sgen[n] <- unlist(mapply(sGenSolver,a=a[n],Smsy=Smsy[n], b=b))
      }else{
        Smsy[n] <- NA
        Umsy[n] <- NA
        Sgen[n] <- NA
      }
    } 
    


  outdf <- data.frame(R=R,
	S=S,
	a=a,
  Umsy=Umsy,
  Smsy=Smsy,
  Sgen=Sgen,
	logR_S=log(R/S))

  return(outdf[(100-nobs):100,])
}


simulateSRtrend <- function(ao=3, b=1/10000, ER=0.4, fec=c(0,0,0,1,0), sig=.5, 
	siga=.3, nobs=40, CapScalar=5, trend="decline",lowsca=.5,hisca=2, ampsc=.5){

    yrs <- 1:100
    S <- NULL
    R <- NULL
    a <- NULL
    Smsy <- NULL
    Umsy <- NULL
    Sgen <- NULL
    
    a[1:5] <- ao
    Seq <- ao/b
    S[1] <- Seq
    R[1] <- Seq*exp(ao-b*Seq)
    
    #trends in a
    #start trend 10 yrs before data collection started
    trlen <- (length(yrs)-(nobs+9)):length(yrs)  

    if(trend=="decline"){
    
      amin <- ao*lowsca
      atrend <- seq(ao,amin, length=nobs+10)
      a[6:(length(yrs)-(nobs+10))] <- ao
      a[trlen] <- atrend
    
    }else if(trend=="increase"){
    
      amax<-ao*hisca
      atrend<-seq(ao,amax, length=nobs+10)
    
    }else if(trend=="sine"){
       
      atrend<- ao + (ao*ampsc)* sin(2*pi/((nobs+10)/2)*trlen)

    }else if(trend=="regime"){
        
      atrend<- ao *rep(rep(c(.5,2),each=10),length.out=length(trlen))
        
    }
    a[6:(length(yrs)-(nobs+10))]<-ao
    a[trlen] <- atrend
    

    Smsy <- NULL
    Umsy <- NULL
    Sgen <- NULL
    for(n in seq_along(a)){

      if(a[n]>0){
        Smsy[n] <- (1 - gsl::lambert_W0(exp(1 - a[n]))) /b
        Umsy[n] <- .5 * a[n] - 0.07 * a[n]^2
        Sgen[n] <- unlist(mapply(sGenSolver,a=a[n],Smsy=Smsy[n], b=b))
      }else{
        Smsy[n] <- NA
        Umsy[n] <- NA
        Sgen[n] <- NA
      }
    } 
    
      
  
    for(y in 2:length(fec)){
      S[y]<-0    
      for(j in seq_along(fec)){
        if((y-j)>0){
          S[y] <- S[y] + R[y-j]*(1-ER)*fec[j]
        }else{
          S[y] <- S[y] + R[1]*(1-ER)*fec[j]
        }
      }
      R[y]<-S[y]*exp(a[y]-b*S[y])        
    }
  
    
    for(y in 6:max(yrs)){ 
      S[y]<-0 
      for(j in seq_along(fec)){      	
        S[y] <- S[y] + R[y-j]*(1-ER)*fec[j]   	
      }        
      
      R[y] <- S[y]*exp(a[y]-b*S[y])*exp(rnorm(1,0,sig))
      
      if(!is.na(R[y]) & R[y] > Seq*CapScalar){
        R[y] <- Seq*CapScalar
      }
      
      if(!is.na(R[y]) & R[y]<5){
      	R[y] <- NA
      }
    }

  outdf <- data.frame(R=R,
	S=S,
	a=a, 
  Umsy=Umsy,
  Smsy=Smsy,
  Sgen=Sgen,
	logR_S=log(R/S))

  return(outdf[(100-nobs+1):100,])
}





runrandomsims <- function(nsim=100, ao=2.5, b=1/30000, ER=0.0, plot_progress=TRUE, trend="random walk",
	fec= c(0,.1,.3,.5,.1), sig=.5, siga=.2, nobs=40, CapScalar=5,lowsca=.5,hisca=2, ampsc=.5, seed=sample.int(100000, 1)){


  #create empty list
  Sim <- list()
  cta <- list()
  ctsmax <- list()
  ctsig <- list()
  dlmKF <- list()
  RB <- list()
  dlmKFalpha <- list()
  RBalpha <- list()
  RBalphamc <- list()
  #holtKF <- list()
  #holtKFalpha <- list()
  tmbholtKF <- list()
  tmbholtKFalpha <- list()
  stanRB <-list()
  stanRBa <-list()
  stanGP <-list()
  stanGPa <- list()

  for(i in seq_len(nsim)){

    set.seed(seed+i)

    if(trend=="random walk"){
      s <- simulateSRrandom(ao=ao, b=b, ER=ER, fec= fec, sig=sig, siga=siga, nobs=nobs, CapScalar=CapScalar )

    }else if(trend=="decline"){
      s <- simulateSRtrend(ao=ao, b=b, ER=ER, fec= fec, sig=sig, siga=NA, nobs=nobs,
           CapScalar=CapScalar, trend="decline",lowsca=lowsca,hisca=hisca, ampsc=ampsc )
   
    }else if(trend=="increase"){
      s <- simulateSRtrend(ao=ao, b=b, ER=ER, fec= fec, sig=sig, siga=NA, nobs=nobs,
           CapScalar=CapScalar, trend="increase",lowsca=lowsca,hisca=hisca, ampsc=ampsc )

    }else if(trend=="sine"){
      s <- simulateSRtrend(ao=ao, b=b, ER=ER, fec= fec, sig=sig, siga=NA, nobs=nobs,
           CapScalar=CapScalar, trend="sine",lowsca=lowsca,hisca=hisca, ampsc=ampsc )
      
    }else if(trend=="regime"){
       s <- simulateSRtrend(ao=ao, b=b, ER=ER, fec= fec, sig=sig, siga=NA, nobs=nobs,
           CapScalar=CapScalar, trend="regime",lowsca=lowsca,hisca=hisca, ampsc=ampsc )
    }else{
      stop(paste("trend patten",trend,"not defined"))
    }

  
    #s <- simulateSRrandom(ao=ao, b=b, ER=ER, fec= fec, sig=sig, siga=siga, nobs=nobs, CapScalar=CapScalar )
    s$extinct<-is.na(s$R)
    Sim[[i]]<-s
  
    if(sum(is.na(s$R))>0){
      cta[[i]]<-NA
      ctsmax[[i]]<-NA
      ctsig[[i]]<-NA
      dlmKF[[i]]<-NA
      RB[[i]]<-NA
      dlmKFalpha[[i]]<-NA
      RBalpha[[i]]<-NA
      RBalphamc[[i]]<-NA
      #holtKF[[i]]<-NA
      #holtKFalpha[[i]]<-NA
      tmbholtKF[[i]]<-NA
      tmbholtKFalpha[[i]]<-NA
      stanRB[[i]]<-NA
      stanRBa[[i]]<-NA
      stanGP[[i]]<-NA
      stanGPa[[i]]<-NA
      next;
    }
 
    srm <- lm(s$logR_S ~ s$S)

    if(is.na(srm$coefficients[[2]])){
      next;
    }

    if(plot_progress){
      plot(s$logR_S~ s$S, main=paste(i))
      abline(srm)
    }
   
    cta[[i]]<-rep(srm$coefficients[[1]],length(s$R))
    ctsmax[[i]]<-1/-srm$coefficients[[2]]
    ctsig[[i]]<-summary(srm)$sigma

    SRdata<-list(obs_logRS=s$logR_S,obs_S=s$S, prbeta1=1.5,
    prbeta2=1.5)
  
    #Model 1 - TMB
    parameters<- list(
      alphao=srm$coefficients[[1]],
      logSmax = log(1/ifelse(-srm$coefficients[[2]]<0,1e-08,-srm$coefficients[2])),
      rho=.5,
      logvarphi=0,
      alpha=rep(srm$coefficients[1],length(s$R))
      )    

    obj <- MakeADFun(SRdata,parameters,DLL="Ricker_tva_Smax_ratiovar",random="alpha")#,lower = -Inf, upper = Inf)
     newtonOption(obj, smartsearch=FALSE)

    opt <- nlminb(obj$par,obj$fn,obj$gr)
  
    sdrep <- summary(sdreport(obj))


    Smsyrb <- (1 - gsl::lambert_W0(exp(1 - sdrep[rownames(sdrep)=="alpha",1]))) /sdrep[rownames(sdrep)=="beta",1]
    #print("Sgen Recursive Bayes")
    #Sgenrb <- unlist(mapply(sGenSolver,a=sdrep[rownames(sdrep)=="alpha",1],
    #  Smsy=Smsyrb, b=sdrep[rownames(sdrep)=="beta",1]))
  
    #MCMC
    fitmcmc1 <- tmbstan(obj, chains = 4, iter = 2000,
                init="random",
                lower=c(-5,4,0,-6),
                 upper=c(10,16,1,6),
                 control = list(adapt_delta = 0.95,max_treedepth = 15), warmup = 500,  thin = 1)
    
    mc <- extract(fitmcmc1, 
                inc_warmup=FALSE, permuted=FALSE)
    
    mcmc<-reshape::melt(mc, as.is=TRUE)
    mcmca<-mcmc[grep("alpha\\[",mcmc$parameters),]
    mcmca$b <- 1/exp(mcmc$value[mcmc$parameters=="logSmax"])
    
    fit_summary <- summary(fitmcmc1)

    umsyrbmcmc <- .5 * mcmca$value - 0.07 * (mcmca$value)^2
    Smsyrbmcmc <- (1 - gsl::lambert_W0(exp(1 - mcmca$value)))/mcmca$b
    #print("Sgen Recursive Bayes MCMC")
    #Sgenrbmcmc <- unlist(mapply(sGenSolver,a=mcmca$value,
    #  Smsy=Smsyrbmcmc, b=mcmca$b))
    mcmcrefs<-data.frame(umsy=umsyrbmcmc,
      Smsy=Smsyrbmcmc )
   
    RB[[i]] <- list(sdrep=sdrep, convergence=opt$convergence, message=opt$message, 
     mcmc= fitmcmc1, mcmcsummary=  fit_summary, mcmcRP=mcmcrefs )
    RBalpha[[i]] <- sdrep[which(rownames(sdrep)=="alpha"),1]
    RBalphamc[[i]] <- fit_summary$summary[grep("alpha\\[",rownames(fit_summary$summary)),"50%"] 
  
    #Model 2 - tv a and static b
    SRdata2 <- data.frame(byr=seq_along(s$S),
      spwn=s$S,
      rec=s$R)
    avarydlm <-fitDLM(data=SRdata2, alpha_vary = TRUE, beta_vary = FALSE)
    umsydlm <- .5 * avarydlm$results$alpha - 0.07 * (avarydlm$results$alpha)^2;
    Smsydlm <- (1 - gsl::lambert_W0(exp(1 - avarydlm$results$alpha))) /-avarydlm$results$beta    
    
    #print("Sgen KF dlm")    
    #Sgendlm <- unlist(mapply(sGenSolver,a=avarydlm$results$alpha,
    #  Smsy=Smsydlm, b=-avarydlm$results$beta))
   

    dlmKF[[i]] <- list(results=avarydlm$results, alpha=avarydlm$results$alpha, sigobs=avarydlm$sd.est[1], 
      siga=avarydlm$sd.est[2], beta=-avarydlm$results$beta[1], smax=-1/avarydlm$results$beta[1],
       message=avarydlm$message, convergence=avarydlm$convergence, Smsy=Smsydlm, umsy=umsydlm)
    dlmKFalpha[[i]] <- avarydlm$results$alpha

    #Model 3 Carrie's KF - this seem to be broken
    
    #initial <- list()
    #initial$mean.a <- srm$coefficients[1]
    #initial$var.a <- .5
    #initial$b <- srm$coefficients[2]
    #initial$ln.sig.e <- log(.5)
    #initial$ln.sig.w <- log(.5)
    #initial$Ts <- 0
    #initial$EstB <- TRUE
    #
    #holtKFfit <- kf.rw(initial=initial,x=s$S,y=s$logR_S)
    #
    #holtKF[[i]]<-list(alpha=holtKFfit$smoothe.mean.a, alpha.var=holtKFfit$smoothe.var.a,
    # sigobs=holtKFfit$sig.e, siga=holtKFfit$sig.w, beta=holtKFfit$b, 
    #  smax=1/holtKFfit$b, convergence=holtKFfit$Report$convergence, message= holtKFfit$Report$message,
    #  filter.alpha=holtKFfit$post.mean.a, filter.vara=holtKFfit$post.var.a) 
    #holtKFalpha[[i]]<-holtKFfit$smoothe.mean.a

    #Model 3.1 Carrie's KF TMB
    
    rekf <- kfTMB(data=s, silent = FALSE, control = TMBcontrol())
    kfrep <- summary(sdreport(rekf$tmb_obj))
    alphakftmb <- kfrep[which(rownames(kfrep)=="smoothemeana"),1]
    

    umsykftmb <- .5 * alphakftmb - 0.07 * (alphakftmb)^2;
    Smsykftmb <- (1 - gsl::lambert_W0(exp(1 - alphakftmb)))/-kfrep[which(rownames(kfrep)=="b"),1]
    #print("Sgen KF TMB") 
    #Sgenkftmb <- unlist(mapply(sGenSolver,a=alphakftmb,
    #  Smsy=Smsykftmb, b=-kfrep[which(rownames(kfrep)=="b"),1]))
  

    tmbholtKF[[i]]<-list(obj=rekf$tmb_obj, sdrep=kfrep, rep=rekf$tmb_obj$report(),message=rekf$model$message,
    convergence=rekf$model$convergence, umsy=umsykftmb, Smsy= Smsykftmb)
    tmbholtKFalpha[[i]]<-kfrep[which(rownames(kfrep)=="smoothemeana"),1]

    #Model 4 Stan RB
    stan_rb<-rstan::stan(file=here('code','stancode','ricker_linear_varying_a.stan'), data=list(R_S = s$logR_S,
                                                 N=nrow(s),
                                                 TT=as.numeric(factor(seq_len(nrow(s)))),
                                                 S=c((s$S))),                                                                                                   
            pars = c('log_a','b','sigma_e','sigma_a'),
            control = list(adapt_delta = 0.95,max_treedepth = 15), warmup = 500, chains = 4, iter = 2000, thin = 1)

    params_stan_rb<-rstan::extract(stan_rb)
    
    umsystanrb <- .5 * params_stan_rb$log_a - 0.07 * (params_stan_rb$log_a)^2
    Smsystanrb <- matrix(NA, nrow=nrow(umsystanrb), ncol=ncol(umsystanrb))
    #print("Sgen stan Recursive Bayes")
    #Sgenstanrb <- matrix(NA, nrow=nrow(umsystanrb), ncol=ncol(umsystanrb))

    for(n in 1:ncol(Smsystanrb)){
      Smsystanrb[,n] <- (1 - gsl::lambert_W0(exp(1 - params_stan_rb$log_a[,n])))/params_stan_rb$b
      #Sgenstanrb[,n] <- unlist(mapply(sGenSolver,a=params_stan_rb$log_a[,n],
      #Smsy=Smsystanrb[,n], b=params_stan_rb$b))
    }
    
    stanRB[[i]]<-list(stanfit=stan_rb,mcmcsummary=summary(stan_rb)$summary, umsy=umsystanrb, 
      Smsy=Smsystanrb)
    stanRBa[[i]]<-summary(stan_rb)$summary[grep("log_a\\[",rownames(summary(stan_rb)$summary)), "50%"] 



    #Model 4 Stan Gaussian Process    
    stan_GP=rstan::stan(file=here('code','stancode','ricker_linear_varying_a_GP.stan'),data=list(R_S = s$logR_S,
                                                                             N=nrow(s),
                                                                             TT=as.numeric(factor(seq_len(nrow(s)))),
                                                                             S=c(s$S)),
            pars = c('log_a','b','sigma_e'),
            control = list(adapt_delta = 0.95,max_treedepth = 15), warmup = 500, chains = 4, iter = 2000, thin = 1)
    
    params_stan_gp <- rstan::extract(stan_GP)
    

    umsyGP <- .5 * params_stan_gp$log_a - 0.07 * (params_stan_gp$log_a)^2
    SmsyGP <- matrix(NA, nrow=nrow(umsyGP), ncol=ncol(umsyGP))
    #SgenGP <- matrix(NA, nrow=nrow(umsyGP), ncol=ncol(umsyGP))
    #print("Sgen stan Gaussian Process")
    for(n in 1:ncol(SmsyGP)){
      SmsyGP[,n] <- (1 - gsl::lambert_W0(exp(1 - params_stan_gp$log_a[,n])))/params_stan_gp$b
    #  SgenGP[,n] <- unlist(mapply(sGenSolver,a=params_stan_gp$log_a[,n],
    #  Smsy=SmsyGP[,n], b=params_stan_gp$b))
    }
    

    stanGP[[i]]<-list(stanfit=stan_GP,mcmcsummary=summary(stan_GP)$summary, umsy=umsyGP, Smsy=SmsyGP)
    stanGPa[[i]]<-summary(stan_GP)$summary[grep("log_a\\[",rownames(summary(stan_GP)$summary)), "50%"] 


	}

  
	return(list(
    seed=seed,
    #s=s, #for testing
    Simtrend = Sim,
    ctatrend = cta,
    ctsmaxtrend = ctsmax,
    ctsigtrend=ctsig,
    dlmKFtrend = dlmKF,
    RBtrend = RB,
    dlmKFalphatrend = dlmKFalpha,
    RBalphatrend = RBalpha,
    RBalphamctrend = RBalphamc,
    #holtKFtrend = holtKF,
    #holtKFalphatrend = holtKFalpha, 
    tmbholtKFtrend = tmbholtKF,
    tmbholtKFalphatrend = tmbholtKFalpha,
    stanRBtrend = stanRB,
    stanRBatrend = stanRBa,
    stanGPtrend = stanGP,
    stanGPatrend = stanGPa
    ))



}


  

calculatepbias<- function(simresult, Smax, sig, siga, Bayesstat="median"){

  nsim <-length(simresult$Simtrend)
  sima <- lapply(simresult$Simtrend, function(x)x$a)
  valid <- unlist(lapply(simresult$Simtrend, function(x)sum(x$extinct)))<1

  
  #bias of estimators
  lmabias <- list()
  dlmabias <- list()
  rbabias <- list()
  rbmcabias <- list()
  #holtabias <- list()
  tmbholtabias <- list()
  stanrbabias <- list()
  stangpabias <- list()

  vs <- 0

  convlm <- sapply(simresult$ctsmaxtrend,  function(x)ifelse(anyNA(x),NA,rep(0,)))
  convdlm <- sapply(simresult$dlmKFtrend,  function(x)ifelse(anyNA(x),NA,x$convergence))
  convRB <- sapply(simresult$RBtrend,  function(x)ifelse(anyNA(x),NA,x$convergence))
  #convholt <- sapply(simresult$holtKFtrend,  function(x)ifelse(anyNA(x),NA,x$convergence))
  convholttmb <- sapply(simresult$tmbholtKFtrend,  function(x)ifelse(anyNA(x),NA,x$convergence))
  


  for(n in 1:nsim){
    if(valid[n]){
      vs<-vs+1
      lmabias[[n]] <-((simresult$ctatrend[[n]]-sima[[n]])/sima[[n]]*100)
      dlmabias[[n]]<-((simresult$dlmKFalphatrend[[n]]-sima[[n]])/sima[[n]]*100)
      rbabias[[n]]<-((simresult$RBalphatrend[[n]]-sima[[n]])/sima[[n]]*100)
      rbmcabias[[n]]<-((simresult$RBalphamctrend[[n]]-sima[[n]])/sima[[n]]*100)
      #holtabias[[n]]<-((simresult$holtKFalphatrend[[n]]-sima[[n]])/sima[[n]]*100)
      tmbholtabias[[n]]<-((simresult$tmbholtKFalphatrend[[n]]-sima[[n]])/sima[[n]]*100)
      arbsatn<-as.vector(simresult$stanRBtrend[[n]]$mcmcsummary[grep("log_a\\[",rownames(simresult$stanRBtrend[[n]]$mcmcsummary)), "50%"]) 
      stanrbabias[[n]]<-(arbsatn-sima[[n]])/sima[[n]]*100
      agpsatn<-as.vector(simresult$stanGPtrend[[n]]$mcmcsummary[grep("log_a\\[",rownames(simresult$stanGPtrend[[n]]$mcmcsummary)), "50%"]) 
      stangpabias[[n]]<-(agpsatn-sima[[n]])/sima[[n]]*100

      if(Bayesstat=="median"){
        
        rbmcabias[[n]]<-((simresult$RBalphamctrend[[n]]-sima[[n]])/sima[[n]]*100)
        arbsatn<-as.vector(simresult$stanRBtrend[[n]]$mcmcsummary[grep("log_a\\[",rownames(simresult$stanRBtrend[[n]]$mcmcsummary)), "50%"]) 
        stanrbabias[[n]]<-(arbsatn-sima[[n]])/sima[[n]]*100
        agpsatn<-as.vector(simresult$stanGPtrend[[n]]$mcmcsummary[grep("log_a\\[",rownames(simresult$stanGPtrend[[n]]$mcmcsummary)), "50%"])   
        stangpabias[[n]]<-(agpsatn-sima[[n]])/sima[[n]]*100
      
      }else if(Bayesstat=="mean"){
        
        arbsatn<-as.vector(simresult$stanRBtrend[[n]]$mcmcsummary[grep("log_a\\[",rownames(simresult$stanRBtrend[[n]]$mcmcsummary)), "mean"]) 
        stanrbabias[[n]]<-(arbsatn-sima[[n]])/sima[[n]]*100
        agpsatn<-as.vector(simresult$stanGPtrend[[n]]$mcmcsummary[grep("log_a\\[",rownames(simresult$stanGPtrend[[n]]$mcmcsummary)), "mean"]) 
        rbmcabias[[n]]<-((simresult$RBalphamctrend[[n]]-sima[[n]])/sima[[n]]*100)
        arbmc<-as.vector(simresult$RBtrend[[n]]$mcmcsummary$summary[grep("alpha\\[",rownames(simresult$RBtrend[[n]]$mcmcsummary$summary)),"mean"])
        rbmcabias[[n]]<-((arbmc)-sima[[n]])/sima[[n]]*100
      
      }else{
        stop(paste("Bayesstat",Bayesstat,"not defined"))
      }
    }    
  }

 

  convaRBmc <- lapply(simresult$RBtrend,  function(x)if(is.na(x[1])){NA}else{as.numeric(abs(x$mcmcsummary$summary[grep("alpha\\[",rownames(x$mcmcsummary$summary)),"Rhat"]-1)>.1)})
  convastanRB <- sapply(simresult$stanRBtrend, function(x)if(is.na(x[1])){NA}else{as.numeric(abs(x$mcmcsummary[grep("log_a\\[",rownames(simresult$stanRBtrend[[n]]$mcmcsummary)),"Rhat"]-1)>.1)})
  convastanGP <-sapply(simresult$stanGPtrend,  function(x)if(is.na(x[1])){NA}else{as.numeric(abs(x$mcmcsummary[grep("log_a\\[",rownames(simresult$stanGPtrend[[n]]$mcmcsummary)),"Rhat"]-1)>.1)})
  

  
  dfabias <- data.frame(pbias=c(unlist(dlmabias),unlist(rbabias), unlist(lmabias), 
    unlist(rbmcabias), unlist(tmbholtabias),unlist(stanrbabias), unlist(stangpabias)),
    fit=rep(c("dlm","RB","lm","RBmcmc","tmbholtKF","stanrb","stanGp"), each=length(unlist(dlmabias))),
    convergence= c(rep(convdlm[!is.na(convdlm)], each=length(sima[[1]])),
     rep(convRB[!is.na(convRB)],each=length(sima[[1]])),
     rep(convlm[!is.na(convlm)], each=length(sima[[1]])), 
     unlist(convaRBmc)[!is.na(unlist(convaRBmc))], 
     rep(convholttmb[!is.na(convholttmb)], each=length(sima[[1]])),
    unlist(convastanRB)[!is.na(unlist(convastanRB))], 
    unlist(convastanGP)[!is.na(unlist(convastanRB))]),
    param="a")
  
 
  
  # look at bias in beta and variance terms
  SmaxRB <- sapply(simresult$RBtrend,  function(x)ifelse(anyNA(x),NA,(exp(x$sdrep["logSmax",1])-Smax)/Smax*100))
  Smaxdlm <- sapply(simresult$dlmKFtrend,  function(x)ifelse(anyNA(x),NA,((1/-x$results$beta[1])-Smax)/Smax*100))
  #Smaxholt <- sapply(simresult$holtKFtrend,  function(x)ifelse(anyNA(x),NA,((-x$smax)-Smax)/Smax*100))
  Smaxtmbholt <- sapply(simresult$tmbholtKFtrend,  function(x)ifelse(anyNA(x),NA,((1/-x$sdrep["b",1])-Smax)/Smax*100)) 
  
  if(Bayesstat=="median"){
    SmaxRBmc <- sapply(simresult$RBtrend,  function(x)ifelse(anyNA(x),NA,(exp(x$mcmcsummary$summary["logSmax","50%"])-Smax)/Smax*100))
    SmaxstanRB <- sapply(simresult$stanRBtrend, function(x)ifelse(anyNA(x),NA,(1/x$mcmcsummary["b","50%"]-Smax)/Smax*100))
    SmaxstanGP <- sapply(simresult$stanGPtrend, function(x)ifelse(anyNA(x),NA,(1/x$mcmcsummary["b","50%"]-Smax)/Smax*100))
  }else if(Bayesstat=="mean"){
    SmaxRBmc <- sapply(simresult$RBtrend,  function(x)ifelse(anyNA(x),NA,(exp(x$mcmcsummary$summary["logSmax","mean"])-Smax)/Smax*100))
    SmaxstanRB <- sapply(simresult$stanRBtrend, function(x)ifelse(anyNA(x),NA,(1/x$mcmcsummary["b","mean"]-Smax)/Smax*100))
    SmaxstanGP <- sapply(simresult$stanGPtrend, function(x)ifelse(anyNA(x),NA,(1/x$mcmcsummary["b","mean"]-Smax)/Smax*100))
  }else{
    stop(paste("Bayesstat",Bayesstat,"not defined"))
  }
  convsmaxRBmc <-sapply(simresult$RBtrend,  function(x)ifelse(anyNA(x),NA,as.numeric(abs(x$mcmcsummary$summary["logSmax","Rhat"]-1)>.1)))
  convsmaxstanRB <-sapply(simresult$stanRBtrend,  function(x)ifelse(anyNA(x),NA,as.numeric(abs(x$mcmcsummary["b","Rhat"]-1)>.1)))
  convsmaxstanGP <-sapply(simresult$stanGPtrend,  function(x)ifelse(anyNA(x),NA,as.numeric(abs(x$mcmcsummary["b","Rhat"]-1)>.1)))
  


 dfsmaxbias <- data.frame(pbias=c(Smaxdlm, SmaxRB, (unlist(simresult$ctsmaxtrend)-Smax)/Smax*100,
  SmaxRBmc,  Smaxtmbholt, SmaxstanRB, SmaxstanGP),
   fit=rep(c("dlm","RB","lm","RBmcmc","tmbholtKF","stanrb","stanGp"), each=length(Smaxdlm)), 
   convergence=c(convdlm, convRB, convlm, convsmaxRBmc, convholttmb, convsmaxstanRB, convsmaxstanGP),
   param="Smax")

  sigRB <- sapply(simresult$RBtrend,  function(x)ifelse(anyNA(x),NA,(x$sdrep["sig",1]-sig)/sig*100))
  sigdlm <- sapply(simresult$dlmKFtrend,  function(x)ifelse(anyNA(x),NA,(x$sigobs-sig)/sig*100))
  
  sigRBmc <- sapply(simresult$RBtrend,  function(x)ifelse(anyNA(x),NA,
    (median(sqrt(as.data.frame(x$mcmc)$"rho") * sqrt(1/exp(as.data.frame(x$mcmc)$"logvarphi")))-sig)/sig*100))
  
  #sigholt <- sapply(simresult$holtKFtrend,  function(x)ifelse(anyNA(x),NA,((x$sigobs)-sig)/sig*100))
  sigtmbholt <- sapply(simresult$tmbholtKFtrend,  function(x)ifelse(anyNA(x),NA,(x$sdrep["sige",1]-sig)/sig*100))
  sigstanRB <- sapply(simresult$stanRBtrend, function(x)ifelse(anyNA(x),NA,(x$mcmcsummary["sigma_e","50%"]-sig)/sig*100))
  sigstanGP <- sapply(simresult$stanGPtrend, function(x)ifelse(anyNA(x),NA,(x$mcmcsummary["sigma_e","50%"]-sig)/sig*100))


  if(Bayesstat=="median"){
    sigRBmc <- sapply(simresult$RBtrend,  function(x)ifelse(anyNA(x),NA,
    (median(sqrt(as.data.frame(x$mcmc)$"rho") * sqrt(1/exp(as.data.frame(x$mcmc)$"logvarphi")))-sig)/sig*100))
    sigstanRB <- sapply(simresult$stanRBtrend, function(x)ifelse(anyNA(x),NA,(x$mcmcsummary["sigma_e","50%"]-sig)/sig*100))
    sigstanGP <- sapply(simresult$stanGPtrend, function(x)ifelse(anyNA(x),NA,(x$mcmcsummary["sigma_e","50%"]-sig)/sig*100))

  }else if(Bayesstat=="mean"){

    sigRBmc <- sapply(simresult$RBtrend,  function(x)ifelse(anyNA(x),NA,
    (mean(sqrt(as.data.frame(x$mcmc)$"rho") * sqrt(1/exp(as.data.frame(x$mcmc)$"logvarphi")))-sig)/sig*100))
    sigstanRB <- sapply(simresult$stanRBtrend, function(x)ifelse(anyNA(x),NA,(x$mcmcsummary["sigma_e","mean"]-sig)/sig*100))
    sigstanGP <- sapply(simresult$stanGPtrend, function(x)ifelse(anyNA(x),NA,(x$mcmcsummary["sigma_e","mean"]-sig)/sig*100))

  }else{
    stop(paste("Bayesstat",Bayesstat,"not defined"))
  }

  convsigRBmc1 <-sapply(simresult$RBtrend,  function(x)ifelse(anyNA(x),NA,as.numeric(abs(x$mcmcsummary$summary["logvarphi","Rhat"]-1)>.1)))
  convsigRBmc2 <-sapply(simresult$RBtrend,  function(x)ifelse(anyNA(x),NA,as.numeric(abs(x$mcmcsummary$summary["rho","Rhat"]-1)>.1)))
  convsigRBmc <- convsigRBmc1 + convsigRBmc2 
  
  convsigstanRB <-sapply(simresult$stanRBtrend,  function(x)ifelse(anyNA(x),NA,as.numeric(abs(x$mcmcsummary["sigma_e","Rhat"]-1)>.1)))
  convsigstanGP <-sapply(simresult$stanGPtrend,  function(x)ifelse(anyNA(x),NA,as.numeric(abs(x$mcmcsummary["sigma_e","Rhat"]-1)>.1)))
  
  

  dfsigbias <- data.frame(pbias=c(sigdlm,sigRB, (unlist(simresult$ctsigtrend)-sig)/sig*100,
    sigRBmc, sigtmbholt, sigstanRB, sigstanGP),
    fit=rep(c("dlm","RB", "lm","RBmcmc","tmbholtKF", "stanrb","stanGp"), each=length(sigdlm)),
    convergence=c(convdlm, convRB, convlm, convsigRBmc, convholttmb, convsigstanRB, convsigstanGP),
    param="sig")
  
  

  sigaRB <- sapply(simresult$RBtrend,  function(x)ifelse(anyNA(x),NA,(x$sdrep["tau",1]-siga)/siga*100))
  sigadlm <- sapply(simresult$dlmKFtrend,  function(x)ifelse(anyNA(x),NA,(x$siga-siga)/siga*100))
  sigaRBmc <- sapply(simresult$RBtrend,  function(x)ifelse(anyNA(x),NA,
    (median(sqrt(1-as.data.frame(x$mcmc)$"rho") * sqrt(1/exp(as.data.frame(x$mcmc)$"logvarphi")))-siga)/siga*100))
  #sigaholt <- sapply(simresult$holtKFtrend,  function(x)ifelse(anyNA(x), NA,((x$siga)-siga)/siga*100))
  sigatmbholt <- sapply(simresult$tmbholtKFtrend,  function(x)ifelse(anyNA(x),NA,(x$sdrep["sigw",1]-siga)/siga*100))
  sigastanRB <- sapply(simresult$stanRBtrend, function(x)ifelse(anyNA(x),NA,(x$mcmcsummary["sigma_a","50%"]-siga)/siga*100))

  simresult$stanGPtrend[[2]]$mcmcsummary
  convsigastanRB <-sapply(simresult$stanRBtrend,  function(x)ifelse(anyNA(x),NA,as.numeric(abs(x$mcmcsummary["sigma_a","Rhat"]-1)>.1)))
  
  if(Bayesstat=="median"){
    sigaRBmc <- sapply(simresult$RBtrend,  function(x)ifelse(anyNA(x),NA,
    (median(sqrt(1-as.data.frame(x$mcmc)$"rho") * sqrt(1/exp(as.data.frame(x$mcmc)$"logvarphi")))-siga)/siga*100))
    sigastanRB <- sapply(simresult$stanRBtrend, function(x)ifelse(anyNA(x),NA,(x$mcmcsummary["sigma_a","50%"]-siga)/siga*100))



  }else if(Bayesstat=="mean"){
    sigaRBmc <- sapply(simresult$RBtrend,  function(x)ifelse(anyNA(x),NA,
    (mean(sqrt(1-as.data.frame(x$mcmc)$"rho") * sqrt(1/exp(as.data.frame(x$mcmc)$"logvarphi")))-siga)/siga*100))
    sigastanRB <- sapply(simresult$stanRBtrend, function(x)ifelse(anyNA(x),NA,(x$mcmcsummary["sigma_a","mean"]-siga)/siga*100))

    
  }else{
    stop(paste("Bayesstat",Bayesstat,"not defined"))
  }

  


  dfsigabias <- data.frame(pbias=c(sigadlm, sigaRB, sigaRBmc, sigatmbholt, sigastanRB),
    fit=rep(c("dlm","RB","RBmcmc", "tmbholtKF", "stanrb"), each=length(sigadlm)),
     convergence=c(convdlm, convRB, convsigRBmc, convholttmb, convsigastanRB),
    param="siga")

  dfbias<-rbind(dfabias,
    dfsmaxbias, 
    dfsigbias,
    dfsigabias)


return(dfbias)


}






comparealpha<- function(simresult, Bayesstat="median"){
  


  sima <- unlist(lapply(simresult$Simtrend, function(x)x$a))
  
  convRB <- unlist(sapply(simresult$RBtrend,  function(x)if(anyNA(x)){rep(NA,length(simresult$Simtrend[[1]]$a))}else{rep(x$convergence,length(simresult$Simtrend[[1]]$a))}))
  convholttmb <- unlist(sapply(simresult$tmbholtKFtrend,  function(x)if(anyNA(x)){rep(NA,length(simresult$Simtrend[[1]]$a))}else{rep(x$convergence,length(simresult$Simtrend[[1]]$a))}))
  convdlm <- unlist(sapply(simresult$dlmKFtrend,  function(x)if(anyNA(x)){rep(NA,length(simresult$Simtrend[[1]]$a))}else{rep(x$convergence,length(simresult$Simtrend[[1]]$a))}))
  convaRBmc <- unlist(sapply(simresult$RBtrend,  function(x)if(is.na(x[1])){rep(NA,length(simresult$Simtrend[[1]]$a))}else{as.numeric(abs(x$mcmcsummary$summary[grep("alpha\\[",rownames(x$mcmcsummary$summary)),"Rhat"]-1)>.1)}))
  convastanRB <- unlist(sapply(simresult$stanRBtrend, function(x)if(is.na(x[1])){rep(NA,length(simresult$Simtrend[[1]]$a))}else{as.numeric(abs(x$mcmcsummary[grep("log_a\\[",rownames(x$mcmcsummary)),"Rhat"]-1)>.1)}))
    

  #dlm estimates
  smootha_dlm <- unlist(lapply(simresult$dlmKFtrend,  
    function(x)if(!anyNA(x)){x$results$alpha}else{
      rep(NA,length(simresult$Simtrend[[1]]$a))}))
  smoothsea_dlm <- unlist(lapply(simresult$dlmKFtrend,  
    function(x)if(!anyNA(x)){x$results$alpha_se}else{
      rep(NA,length(simresult$Simtrend[[1]]$a))}))
  
  #filtered
  #need to calculate
  
  #holt estimates
  #smoothed
  smootha_KFtmb <- unlist(lapply(simresult$tmbholtKFtrend,  
    function(x)if(!anyNA(x)){x$sdrep[which(rownames(x$sdrep)=="smoothemeana"),1]}else{
      rep(NA,length(simresult$Simtrend[[1]]$a))}))
  smoothsea_KFtmb <- unlist(lapply(simresult$tmbholtKFtrend,  
    function(x)if(!anyNA(x)){sqrt(x$sdrep[which(rownames(x$sdrep)=="smoothevara"),1])}else{
      rep(NA,length(simresult$Simtrend[[1]]$a))}))
  
  #filtered
  filtera_KFtmb<-unlist(lapply(simresult$tmbholtKFtrend,  
    function(x)if(!anyNA(x)){x$rep$postmeana}else{rep(NA,length(simresult$Simtrend[[1]]$a))}))
  filtersea_KFtmb<-unlist(lapply(simresult$tmbholtKFtrend,  
    function(x)if(!anyNA(x)){sqrt(x$rep$postvara)}else{rep(NA,length(simresult$Simtrend[[1]]$a))}))


  #Recursive Bayes TMB 
  a_RBtmb <- unlist(lapply(simresult$RBtrend,  
    function(x)if(!anyNA(x)){x$sdrep[which(rownames(x$sdrep)=="alpha"),1]}else{
      rep(NA,length(simresult$Simtrend[[1]]$a))}))
  sea_RBtmb <- unlist(lapply(simresult$RBtrend,  
    function(x)if(!anyNA(x)){x$sdrep[which(rownames(x$sdrep)=="alpha"),2]}else{
      rep(NA,length(simresult$Simtrend[[1]]$a))}))

  #Recursive Bayes TMB mcmc
  #Recursive Bayes stan mcmc
  
      
  if(Bayesstat=="median"){

    a_RBmc <- unlist(lapply(simresult$RBtrend,  
      function(x)if(!anyNA(x)){x$mcmcsummary$summary[grep("alpha\\[",rownames(x$mcmcsummary$summary)),"50%"]}else{
      rep(NA,length(simresult$Simtrend[[1]]$a))}))
        
    a_stan <- unlist(lapply(simresult$stanRBtrend,  
      function(x)if(!anyNA(x)){x$mcmcsummary[grep("log_a\\[",rownames(x$mcmcsummary)),"50%"]}else{
      rep(NA,length(simresult$Simtrend[[1]]$a))}))

  }else if(Bayesstat=="mean"){
        
    a_RBmc <- unlist(lapply(simresult$RBtrend,  
      function(x)if(!anyNA(x)){x$mcmcsummary$summary[grep("alpha\\[",rownames(x$mcmcsummary$summary)),"mean"]}else{
      rep(NA,length(simresult$Simtrend[[1]]$a))}))
        
    a_stan <- unlist(lapply(simresult$stanRBtrend,  
      function(x)if(!anyNA(x)){x$mcmcsummary[grep("log_a\\[",rownames(x$mcmcsummary)),"mean"]}else{
      rep(NA,length(simresult$Simtrend[[1]]$a))}))
      
  }else{
        stop(paste("Bayesstat",Bayesstat,"not defined"))
  }
   


  dfa <- data.frame(a=c(smootha_dlm, smootha_KFtmb, filtera_KFtmb, a_RBtmb,  a_RBmc, a_stan ),
    abias=c((smootha_dlm-sima)/sima*100, (smootha_KFtmb-sima)/sima*100, (filtera_KFtmb-sima)/sima*100, 
      (a_RBtmb-sima)/sima*100,  (a_RBmc-sima)/sima*100, (a_stan-sima)/sima*100 ),
    se_a =c(smoothsea_dlm, smoothsea_KFtmb, filtersea_KFtmb, sea_RBtmb,rep(NA,length(a_RBmc)),rep(NA,length(a_stan) )),
    fit=rep(c("smoothdlm","smoothKF", "filterKF", "RB","RBmcmc", "stanrb"), each=length(smootha_dlm)),
     convergence=c(convdlm, convholttmb, convholttmb, convRB, convaRBmc, convastanRB))

  return(dfa)      
}



#==========================================================
# old functions



runtrendsims <- function(nsim=100, ao=3, b=1/30000, ER=0.0, fec= c(0,.1,.3,.5,.1), sig=.5, siga=.2, nobs=40,
  CapScalar=5, trend="decline",lowsca=.5,hisca=2, ampsc=.5, plot_progress=TRUE){
    
  Simtrend <-list()
  ctatrend <-list()
  ctsmaxtrend <-list()
  ctsigtrend <- list()
  dlmKFtrend <- list()
  RBtrend <- list()
  dlmKFalphatrend <- list()
  RBalphatrend <- list()
  RBalphamctrend <- list()
  holtKFtrend <- list()
  holtKFalphatrend <- list() 
  tmbholtKFtrend <- list()
  tmbholtKFalphatrend <- list()
  stanRBtrend <-list()
  stanGPtrend <-list()


  for(i in seq_len(nsim)){
    
    s <- simulateSRtrend(ao=ao, b=b, ER=ER, fec=fec, sig=sig, siga=siga, nobs=nobs,
    CapScalar=CapScalar, trend=trend,lowsca=lowsca,hisca=hisca, ampsc=ampsc )
    s$extinct<-is.na(s$R)
    Simtrend[[i]]<-s
    
  
    if(sum(is.na(s$R))>0){
      ctatrend[[i]]<-NA
      ctsmaxtrend[[i]]<-NA
      ctsigtrend[[i]]<-NA
      dlmKFtrend[[i]]<-NA
      RBtrend[[i]]<-NA
      dlmKFalphatrend[[i]]<-NA
      RBalphatrend[[i]]<-NA
      RBalphamctrend[[i]]<-NA
      holtKFtrend[[i]]<-NA
      holtKFalphatrend[[i]]<-NA 
      tmbholtKFtrend[[i]]<-NA
      tmbholtKFalphatrend[[i]]<-NA
      stanRBtrend[[i]]<-NA
      stanGPtrend[[i]]<-NA
      next;
    }
  
    srm <- lm(s$logR_S~ s$S)
  
    if(is.na(srm$coefficients[[2]])){
      next;
    }
    
    if(plot_progress){
      plot(s$logR_S~ s$S, main=paste(i))
      abline(srm)
    }
    
     
    ctatrend[[i]]<-rep(srm$coefficients[[1]],length(s$R))
    ctsmaxtrend[[i]]<-1/-srm$coefficients[[2]]
    ctsigtrend[[i]]<-summary(srm)$sigma
  
  
    SRdata<-list(obs_logRS=s$logR_S,obs_S=s$S, prbeta1=1.5,
      prbeta2=1.5)
    
    
    #Model 1 - TMB
    parameters<- list(
      alphao=srm$coefficients[[1]],
      logSmax = log(1/ifelse(-srm$coefficients[[2]]<0,1e-08,-srm$coefficients[2])),
      rho=.5,
      logvarphi=0,
      alpha=rep(srm$coefficients[1],length(s$R))
      )    
  
    obj <- MakeADFun(SRdata,parameters,DLL="Ricker_tva_Smax_ratiovar",random="alpha")#,lower = -Inf, upper = Inf)
     newtonOption(obj, smartsearch=FALSE)
  
    opt <- nlminb(obj$par,obj$fn,obj$gr)
    
    sdrep <- summary(sdreport(obj))
    
  
    RBtrend[[i]] <- list(sdrep=sdrep, convergence=opt$convergence, message=opt$message)#, 
    #  mcmc= fitmcmc1, mcmcsummary=  fit_summary )
    RBalphatrend[[i]] <- sdrep[which(rownames(sdrep)=="alpha"),1]
    #RBalphamc[[i]] <- fit_summary$summary[grep("alpha\\[",rownames(fit_summary$summary)),"mean"] 
    
  
    #Model 2 - tv a and static b
    SRdata2 <- data.frame(byr=seq_along(s$S),
      spwn=s$S,
      rec=s$R)
    avarydlm <-fitDLM(data=SRdata2, alpha_vary = TRUE, beta_vary = FALSE)
    dlmKFtrend[[i]] <- list(results=avarydlm$results, alpha=avarydlm$results$alpha, sigobs=avarydlm$sd.est[1], 
      siga=avarydlm$sd.est[2], beta=-avarydlm$results$beta[1], smax=-1/avarydlm$results$beta[1],
       message=avarydlm$message, convergence=avarydlm$convergence)
    dlmKFalphatrend[[i]] <- avarydlm$results$alpha

    #Model 3 Carrie's KF - this seem to be broken
    
    initial <- list()
    initial$mean.a <- srm$coefficients[1]
    initial$var.a <- .5
    initial$b <- srm$coefficients[2]
    initial$ln.sig.e <- log(.5)
    initial$ln.sig.w <- log(.5)
    initial$Ts <- 0
    initial$EstB <- TRUE
    
    holtKFfit <- kf.rw(initial=initial,x=s$S,y=s$logR_S)
  
    holtKFtrend[[i]]<-list(alpha=holtKFfit$smoothe.mean.a, sigobs=holtKFfit$sig.e, siga=holtKFfit$sig.w, beta=holtKFfit$b, 
      smax=1/holtKFfit$b, convergence=holtKFfit$Report$convergence, message= holtKFfit$Report$message) 
    holtKFalphatrend[[i]]<-holtKFfit$smoothe.mean.a

    #Model 3.1 Carrie's KF TMB
    
    rekf <- kfTMB(data=s, silent = FALSE, control = TMBcontrol())
    kfrep <- summary(sdreport(rekf$tmb_obj))
    

   tmbholtKFtrend[[i]]<-list(obj=rekf$tmb_obj, sdrep=kfrep, rep=rekf$tmb_obj$report(),message=rekf$model$message,
   convergence=rekf$model$convergence )
   tmbholtKFalphatrend[[i]]<-kfrep[which(rownames(kfrep)=="smoothemeana"),1]
   

   #Model 4 Stan RB
   stan_rb=rstan::stan(file=here('code','stancode','ricker_linear_varying_a.stan'), data=list(R_S = s$logR_S,
                                                 N=nrow(s),
                                                 TT=as.numeric(factor(seq_len(nrow(s)))),
                                                 S=c((s$S))),
                                                                                                      
            pars = c('log_a','b','sigma_e','sigma_a'),
            control = list(adapt_delta = 0.95,max_treedepth = 15), warmup = 500, chains = 4, iter = 2000, thin = 1)

    #params_stan_rb=rstan::extract(stan_rb)

    stanRBtrend[[i]]<-list(stanfit=stan_rb,mcmcsummary=summary(stan_rb)$summary)

    #Model 4 Stan Gaussian Process    
    stan_GP=rstan::stan(file=here('code','stancode','ricker_linear_varying_a_GP.stan'),data=list(R_S = s$logR_S,
                                                                             N=nrow(s),
                                                                             TT=as.numeric(factor(seq_len(nrow(s)))),
                                                                             S=c(s$S)),
            pars = c('log_a','b','sigma_e'),
            control = list(adapt_delta = 0.95,max_treedepth = 15), warmup = 500, chains = 4, iter = 2000, thin = 1)


    stanGPtrend[[i]]<-list(stanfit=stan_GP,mcmcsummary=summary(stan_GP)$summary)


  
  }

  return(list(Simtrend =Simtrend,
  ctatrend = ctatrend,
  ctsmaxtrend = ctsmaxtrend,
  ctsigtrend = ctsigtrend,
  dlmKFtrend = dlmKFtrend,
  RBtrend = RBtrend,
  dlmKFalphatrend = dlmKFalphatrend,
  RBalphatrend = RBalphatrend,
  RBalphamctrend = RBalphamc,
  holtKFtrend = holtKF,
  holtKFalphatrend = holtKFalpha, 
  tmbholtKFtrend = tmbholtKF,
  tmbholtKFalphatrend = tmbholtKFalpha,
  stanRBtrend = stanRBtrend,
  stanGPtrend = stanGPtrend
     ))


}
