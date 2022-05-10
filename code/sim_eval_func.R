#==============================================================
#Simple simulation evaluation of Ricker model with time-varying parameters
#Catarina Wor
#May 2022
#==============================================================




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
  
    
    for(y in 6:max(yrs)){  
      S[y]<-0
      for(j in seq_along(fec)){      	
        S[y] <- S[y] + R[y-j]*(1-ER)*fec[j]   	
      } 
      #truncated lognrmal error
      #a[y] <- a[y-1] + (qnorm(runif(1, 0.1, 0.9),0,siga))        
      a[y] <- a[y-1] + rnorm(1,0,siga) -(.5*siga^2)

      #R[y] <- S[y]*exp(a[y]-b*S[y]) *exp(qnorm(runif(1, 0.1, 0.95),0,sig))
      R[y] <- S[y]*exp(a[y]-b*S[y]) *exp(rnorm(1,0,sig)-(.5*sig^2))#-(.5*sig^2)
      if(!is.na(R[y]) &  R[y] > Seq*CapScalar){
        R[y] <- Seq*CapScalar
      }
      if(!is.na(R[y]) & R[y]<5){
      	R[y] <- NA
      }
      
    }


  outdf <- data.frame(R=R,
	S=S,
	a=a,
	logR_S=log(R/S))

  return(outdf[(100-nobs):100,])
}


simulateSRtrend <- function(ao=3, b=1/10000, ER=0.4, fec=c(0,0,0,1,0), sig=.5, 
	siga=.3, nobs=40, CapScalar=5, trend="decline",lowsca=.5,hisca=2, ampsc=.5){

    yrs <- 1:100
    S <- NULL
    R <- NULL
    a <- NULL
   
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
	logR_S=log(R/S))

  return(outdf[(100-nobs+1):100,])
}





runrandomsims <- function(nsim=100, ao=2.5, b=1/30000, ER=0.0, plot_progress=TRUE, 
	fec= c(0,.1,.3,.5,.1), sig=.5, siga=.2, nobs=40, CapScalar=5){



  #create empty list
  Sim <- list()
  cta <- list()
  ctsmax <- list()
  ctsig <- list()
  dlmKF <- list()
  RB <- list()
  dlmKFalpha <- list()
  RBalpha <- list()


  for(i in seq_len(nsim)){
  
    s <- simulateSRrandom(ao=ao, b=b, ER=ER, fec= fec, sig=sig, siga=siga, nobs=nobs, CapScalar=CapScalar )
    s$extinct<-is.na(s$R)
    Sim[[i]]<-s
  
    if(sum(is.na(s$R))>0){
      cta[[i]]<-NA
      ctsmax[[i]]<-NA
      ctsig[[i]]<-NA
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
  
    #MCMC
    #fitmcmc1 <- tmbstan(obj, chains=3,
    #            iter=10000, init="random",
    #            lower=c(-10,4,0,-6),
    #             upper=c(5,16,1,6),
    #             control = list(adapt_delta = 0.98))
    #
    #  mc <- extract(fitmcmc1, pars=names(obj$par),
    #            inc_warmup=TRUE, permuted=FALSE)
    #fit_summary <- summary(fitmcmc1)   
    #fit_summary$summary[grep("alpha\\[",rownames(fit_summary$summary)),"mean"]
  
    RB[[i]] <- list(sdrep=sdrep, convergence=opt$convergence, message=opt$message)#, 
    #  mcmc= fitmcmc1, mcmcsummary=  fit_summary )
    RBalpha[[i]] <- sdrep[which(rownames(sdrep)=="alpha"),1]
    #RBalphamc[[i]] <- fit_summary$summary[grep("alpha\\[",rownames(fit_summary$summary)),"mean"] 
  

    #Model 2 - tv a and static b
    SRdata2 <- data.frame(byr=seq_along(s$S),
      spwn=s$S,
      rec=s$R)
    avarydlm <-fitDLM(data=SRdata2, alpha_vary = TRUE, beta_vary = FALSE)
    dlmKF[[i]] <- list(results=avarydlm$results, alpha=avarydlm$results$alpha, sigobs=avarydlm$sd.est[1], 
      siga=avarydlm$sd.est[2], beta=-avarydlm$results$beta[1], smax=-1/avarydlm$results$beta[1],
       message=avarydlm$message, convergence=avarydlm$convergence)
    dlmKFalpha[[i]] <- avarydlm$results$alpha

    #Model 3 Carrie's KF - this seem to be broken
    
    #initial <- list()
    #initial$mean.a <- srm$coefficients[1]
    #initial$var.a <- .5
    #initial$b <- initial$var.srm$coefficients[2]
    #initial$ln.sig.e <- log(.5)
    #initial$ln.sig.w <- log(.5)
    #initial$Ts <- 0
    #initial$EstB <- TRUE
    
    #holtKFfit <- kf.rw(initial=initial,x=s$S,y=s$logR_S)
  
    #holtKF[[i]]<-list(alpha=holtKFfit$smoothe.mean.a, sigobs=holtKFfit$sig.e, siga=holtKFfit$sig.w, beta=holtKFfit$b, 
    #  smax=1/holtKFfit$b, convergence=holtKFfit$Report$convergence, message= holtKFfit$Report$message) 
    #holtKFalpha[[i]]<-holtKFfit$smoothe.mean.a

	}


	return(list(Simtrend = Sim,
    ctatrend = cta,
    ctsmaxtrend = ctsmax,
    ctsigtrend=ctsig,
    dlmKFtrend = dlmKF,
    RBtrend = RB,
    dlmKFalphatrend = dlmKFalpha,
    RBalphatrend = RBalpha))



}


  

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




  for(i in seq_len(nsim)){
    
    s <- simulateSRtrend(ao=3, b=1/30000, ER=0.0, fec= c(0,.1,.3,.5,.1), sig=.5, siga=.2, nobs=40,
    CapScalar=5, trend="decline",lowsca=.5,hisca=2, ampsc=.5 )
    s$extinct<-is.na(s$R)
    Simtrend[[i]]<-s
    
  
    if(sum(is.na(s$R))>0){
      ctatrend[[i]]<-NA
      ctsmaxtrend[[i]]<-NA
      ctsigtrend[[i]]<-NA
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
  
  }

  return(list(Simtrend =Simtrend,
  ctatrend = ctatrend,
  ctsmaxtrend = ctsmaxtrend,
  ctsigtrend = ctsigtrend,
  dlmKFtrend = dlmKFtrend,
  RBtrend = RBtrend,
  dlmKFalphatrend = dlmKFalphatrend,
  RBalphatrend =RBalphatrend ))


}


calculatepbias<- function(simresult){

  nsim <-length(simresult$Simtrend)
  sima <- lapply(simresult$Simtrend, function(x)x$a)
  valid <- unlist(lapply(simresult$Simtrend, function(x)sum(x$extinct)))<1



  #bias of estimators
  lmabias <- list()
  dlmabias <- list()
  rbabias <- list()
  vs <- 0
  for(n in 1:nsim){
    if(valid[n]){
      vs<-vs+1
      lmabias[[n]]<-((simresult$ctatrend[[n]]-sima[[n]])/sima[[n]]*100)
      dlmabias[[n]]<-((simresult$dlmKFalphatrend[[n]]-sima[[n]])/sima[[n]]*100)
      rbabias[[n]]<-((simresult$RBalphatrend[[n]]-sima[[n]])/sima[[n]]*100)
    }    
  }

  dfabias <- data.frame(pbias=c(unlist(dlmabias),unlist(rbabias), unlist(lmabias)),
    fit=rep(c("dlm","RB","lm"), each=length(unlist(dlmabias))),
    param="a")
  
  # look at bias in beta and variance terms
  SmaxRB <- sapply(simresult$RBtrend,  function(x)ifelse(is.null(x),NA,(exp(x$sdrep["logSmax",1])-Smax)/Smax*100))
  Smaxdlm <- sapply(simresult$dlmKFtrend,  function(x)ifelse(is.null(x),NA,((1/-x$results$beta[1])-Smax)/Smax*100))

 dfsmaxbias <- data.frame(pbias=c(Smaxdlm, SmaxRB, (unlist(simresult$ctsmaxtrend)-Smax)/Smax*100),
   fit=rep(c("dlm","RB","lm"), each=length(Smaxdlm)), param="Smax")

  sigRB <- sapply(simresult$RBtrend,  function(x)ifelse(is.null(x),NA,x$sdrep["sig",1]))
  sigdlm <- sapply(simresult$dlmKFtrend,  function(x)ifelse(is.null(x),NA,x$sigobs))

  dfsigbias <- data.frame(pbias=c((sigdlm-sig)/sig*100,(sigRB-sig)/sig*100, (unlist(simresult$ctsigtrend)-sig)/sig*100),
    fit=rep(c("dlm","RB", "lm"), each=length(sigdlm)),
    param="sig")


  sigaRB <- sapply(simresult$RBtrend,  function(x)ifelse(is.null(x),NA,x$sdrep["tau",1]))
  sigadlm <- sapply(simresult$dlmKFtrend,  function(x)ifelse(is.null(x),NA,x$siga))

  dfsigabias <- data.frame(pbias=c((sigadlm-siga)/siga*100,(sigaRB-siga)/siga*100),
    fit=rep(c("dlm","RB"), each=length(sigadlm)),
    param="siga")

  dfbias<-rbind(dfabias,
    dfsmaxbias, 
    dfsigbias,
    dfsigabias)


return(dfbias)


}
