#==============================================================
#Simple simulation evaluation of Ricker model with time-varying parameters
#Catarina Wor
#May 2022
#==============================================================
 #TODO: deal with  excess variability in simulated SRtime-seies
#need to fix this file
#remotes::install_github("carrieholt/KF-funcs")
#library(KFfuncs)
#use this one instead
#source("C:/Users/worc/Documents/KF-funcs-appl/HoltMichielsens2020/KFcode.R")

library(here)
library(TMB)
library(tmbstan)

source(here("code","dlm-wrapper.R"))
source(here("code","sim_eval_func.R"))

#source("C:/Users/worc/Documents/KF-funcs-appl/HoltMichielsens2020/KFcode.R")

compile("TMBmodels/Ricker_tva_Smax_ratiovar.cpp")
dyn.load(dynlib("TMBmodels/Ricker_tva_Smax_ratiovar"))


ao <- 2
b<-1/30000


ps<-seq(0,35000*5,1000)
pr<-ps*exp(ao-b*ps)

sr <- simulateSRrandom(ao=2.5, b=1/30000, ER=0.0, fec= c(0,.1,.3,.5,.1), sig=.5, siga=.2, nobs=40, CapScalar=5 )
sr
par(mfrow=c(2,1))
plot(sr$S,sr$R, xlim=c(0,35000*5), ylim=c(0,89000))
abline(a=0,b=1)
lines(ps,pr,col="red")
plot(sr$a,type="b")



#==========================================================


#simulation estimation random walk
randsim<-runrandomsims(nsim=100,ao=3, b=1/30000, ER=0.0, fec= c(0,.1,.3,.5,.1), sig=.5, siga=.2, nobs=40, CapScalar=5)

nsim <- 100

Sim<-list()
cta<-list()
ctsmax<-list()
dlmKF <- list()
RB <- list()
dlmKFalpha <- list()
RBalpha <- list()



for(i in seq_len(nsim)){
  
  s <- simulateSRrandom(ao=2.5, b=1/30000, ER=0.0, fec= c(0,.1,.3,.5,.1), sig=.5, siga=.2, nobs=40, CapScalar=5 )
  s$extinct<-is.na(s$R)
  Sim[[i]]<-s
  

  if(sum(is.na(s$R))>0){
    cta[[i]]<-NA
    ctsmax[[i]]<-NA
    next;
  }
 

  srm <- lm(s$logR_S~ s$S)

  if(is.na(srm$coefficients[[2]])){
    next;
  }

  plot(s$logR_S~ s$S, main=paste(i))
  abline(srm)
   
  cta[[i]]<-rep(srm$coefficients[[1]],length(s$R))
  ctsmax[[i]]<-1/-srm$coefficients[[2]]



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


sima<- lapply(Sim, function(x)x$a)
valid <- unlist(lapply(Sim, function(x)sum(x$extinct)))<1

#bias of estimators
lmabias<-list()
dlmabias<-list()
rbabias<-list()
vs <- 0
for(n in 1:nsim){
  if(valid[n]){
    vs<-vs+1
    lmabias[[n]]<-mean((cta[[n]]-sima[[n]])/sima[[n]]*100)
    dlmabias[[n]]<-mean((dlmKFalpha[[n]]-sima[[n]])/sima[[n]]*100)
    rbabias[[n]]<-mean((RBalpha[[n]]-sima[[n]])/sima[[n]]*100)
  }    

}

dfbias <- data.frame(pbias=c(unlist(dlmabias),unlist(rbabias), unlist(lmabias)),
  fit=rep(c("dlm","RB","lm"),each=length(unlist(dlmabias))))

ggplot(dfbias) +
geom_boxplot(aes(x=fit, y=pbias))+
coord_cartesian(ylim = c(-100,100))



# look at bias in beta and variance terms

RB[[1]]$sdrep["logSmax",1]

SmaxRB <- sapply(RB,  function(x)ifelse(is.null(x),NA,exp(x$sdrep["logSmax",1])))
Smaxdlm <- sapply(dlmKF,  function(x)ifelse(is.null(x),NA,1/-x$results$beta[1]))

dfsmaxbias <- data.frame(smax=c(Smaxdlm,SmaxRB ,unlist(ctsmax)),
  fit=rep(c("dlm","RB","lm"),each=length(Smaxdlm)))

ggplot(dfsmaxbias) +
geom_boxplot(aes(x=fit, y=smax))+
geom_hline(yintercept=30000) +
coord_cartesian(ylim = c(0,50000))



#================================
RB[[1]]$sdrep["sig",1]

sigRB <- sapply(RB,  function(x)ifelse(is.null(x),NA,x$sdrep["sig",1]))
sigdlm <- sapply(dlmKF,  function(x)(ifelse(is.null(x),NA,x$sigobs)))

dfsigbias <- data.frame(sig=c(sigdlm,sigRB),
  fit=rep(c("dlm","RB"),each=length(sigRB)))

ggplot(dfsigbias) +
geom_boxplot(aes(x=fit, y=sig))+
geom_hline(yintercept=.5)



#================================



sigaRB <- sapply(RB,  function(x)ifelse(is.null(x),NA,x$sdrep["tau",1]))
sigadlm <- sapply(dlmKF,  function(x)(ifelse(is.null(x),NA,x$siga)))

dfsigabias <- data.frame(siga=c(sigadlm,sigaRB),
  fit=rep(c("dlm","RB"),each=length(sigRB)))

ggplot(dfsigabias) +
geom_boxplot(aes(x=fit, y=siga))+
geom_hline(yintercept=.2) 


#===============================================================================================
#Trends sim eval

nsim<-100
sr <- simulateSRtrend(ao=3, b=1/30000, ER=0.0, fec= c(0,.1,.3,.5,.1), sig=.5, siga=.2, nobs=40,
CapScalar=5, trend="decline",lowsca=.5,hisca=2, ampsc=.5 )
sr
par(mfrow=c(2,1))
plot(sr$S,sr$R, xlim=c(0,35000*5), ylim=c(0,89000))
abline(a=0,b=1)
lines(ps,pr,col="red")
plot(sr$a,type="b")


Simdecline<-list()
ctadecline<-list()
ctsmaxdecline<-list()
dlmKFdecline <- list()
RBdecline <- list()
dlmKFalphadecline <- list()
RBalphadecline <- list()




for(i in seq_len(nsim)){
  
  s <- simulateSRtrend(ao=3, b=1/30000, ER=0.0, fec= c(0,.1,.3,.5,.1), sig=.5, siga=.2, nobs=40,
CapScalar=5, trend="decline",lowsca=.5,hisca=2, ampsc=.5 )
  s$extinct<-is.na(s$R)
  Simdecline[[i]]<-s
  

  if(sum(is.na(s$R))>0){
    next;
  }
  
  if(sum(is.na(s$R))>0){
    next;
  }


  srm <- lm(s$logR_S~ s$S)

  if(is.na(srm$coefficients[[2]])){
    next;
  }

  plot(s$logR_S~ s$S, main=paste(i))
  abline(srm)
   
  ctadecline[[i]]<-rep(srm$coefficients[[1]],length(s$R))
  ctsmaxdecline[[i]]<-1/-srm$coefficients[[2]]



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
  

  RBdecline[[i]] <- list(sdrep=sdrep, convergence=opt$convergence, message=opt$message)#, 
  #  mcmc= fitmcmc1, mcmcsummary=  fit_summary )
  RBalphadecline[[i]] <- sdrep[which(rownames(sdrep)=="alpha"),1]
  #RBalphamc[[i]] <- fit_summary$summary[grep("alpha\\[",rownames(fit_summary$summary)),"mean"] 
  

  #Model 2 - tv a and static b
  SRdata2 <- data.frame(byr=seq_along(s$S),
    spwn=s$S,
    rec=s$R)
  avarydlm <-fitDLM(data=SRdata2, alpha_vary = TRUE, beta_vary = FALSE)
  dlmKFdecline[[i]] <- list(results=avarydlm$results, alpha=avarydlm$results$alpha, sigobs=avarydlm$sd.est[1], 
    siga=avarydlm$sd.est[2], beta=-avarydlm$results$beta[1], smax=-1/avarydlm$results$beta[1],
     message=avarydlm$message, convergence=avarydlm$convergence)
  dlmKFalphadecline[[i]] <- avarydlm$results$alpha

  

}




simadecline<- lapply(Simdecline, function(x)x$a)
validdecline <- unlist(lapply(Simdecline, function(x)sum(x$extinct)))<1

#bias of estimators
lmabiasdecline<-list()
dlmabiasdecline<-list()
rbabiasdecline<-list()
vs <- 0
for(n in 1:nsim){
  if(validdecline[n]){
    vs<-vs+1
    lmabiasdecline[[n]]<-((ctadecline[[n]]-simadecline[[n]])/simadecline[[n]]*100)
    dlmabiasdecline[[n]]<-((dlmKFalphadecline[[n]]-simadecline[[n]])/simadecline[[n]]*100)
    rbabiasdecline[[n]]<-((RBalphadecline[[n]]-simadecline[[n]])/simadecline[[n]]*100)
  }    

}

dfbiasdecline <- data.frame(pbias=c(unlist(dlmabiasdecline),unlist(rbabiasdecline), unlist(lmabiasdecline)),
  fit=rep(c("dlm","RB","lm"), each=length(unlist(dlmabiasdecline))))

ggplot(dfbiasdecline) +
geom_boxplot(aes(x=fit, y=pbias))+
coord_cartesian(ylim = c(-100,100))



# look at bias in beta and variance terms



SmaxRBdecline <- sapply(RBdecline,  function(x)exp(x$sdrep["logSmax",1]))
Smaxdlmdecline <- sapply(dlmKFdecline,  function(x)(1/-x$results$beta[1]))

dfsmaxbiasdecline <- data.frame(smax=c(Smaxdlmdecline,SmaxRBdecline ,unlist(ctsmaxdecline)),
  fit=rep(c("dlm","RB","lm"), each=length(Smaxdlmdecline)))

ggplot(dfsmaxbiasdecline) +
geom_boxplot(aes(x=fit, y=smax))+
geom_hline(yintercept=30000) +
coord_cartesian(ylim = c(0,50000))



#================================


sigRBdecline <- sapply(RBdecline,  function(x)x$sdrep["sig",1])
sigdlmdecline <- sapply(dlmKFdecline,  function(x)(x$sigobs))

dfsigbiasdecline <- data.frame(sig=c(sigdlmdecline,sigRBdecline),
  fit=rep(c("dlm","RB"), each=length(sigdlmdecline)))

ggplot(dfsigbiasdecline) +
geom_boxplot(aes(x=fit, y=sig))+
geom_hline(yintercept=.5)



#================================



sigaRBdecline <- sapply(RBdecline,  function(x)x$sdrep["tau",1])
sigadlmdecline <- sapply(dlmKFdecline,  function(x)(x$siga))

dfsigabiasdecline <- data.frame(siga=c(sigadlmdecline,sigaRBdecline),
  fit=rep(c("dlm","RB"), each=length(sigdlmdecline)))

ggplot(dfsigabiasdecline) +
geom_boxplot(aes(x=fit, y=siga))



