#==============================================================
#Simple simulation evaluation of Ricker model with time-varying parameters
#Catarina Wor
#May 2022
#==============================================================
 #TODO: deal with  excess variability in simulated SRtime-seies

remotes::install_github("carrieholt/KF-funcs")

library(KFfuncs)
library(here)
library(TMB)
library(tmbstan)

source(here("code","dlm-wrapper.R"))
source(here("code","sim_eval_func.R"))

#source("C:/Users/worc/Documents/KF-funcs-appl/HoltMichielsens2020/KFcode.R")

compile("TMBmodels/Ricker_tva_Smax_ratiovar.cpp")
dyn.load(dynlib("TMBmodels/Ricker_tva_Smax_ratiovar"))


ao <- 3
b<-1/10000


ps<-0:35000
pr<-ps*exp(ao-b*ps)

sr <- simulateSRrandom(ao=3, b=1/10000, ER=0.0, fec=c(0,0,0,1,0), sig=.5, siga=.2, nobs=40 )
sr
par(mfrow=c(2,1))
plot(sr$S,sr$R, xlim=c(0,35000), ylim=c(0,89000))
abline(a=0,b=1)
lines(ps,pr,col="red")
plot(sr$a,type="b")



#==========================================================
#simulation estimation calls

nsim <- 100

Sim<-list()
#holtKF <- list()
cta<-list()
dlmKF <- list()
RB <- list()
#holtKFalpha <- list()
dlmKFalpha <- list()
RBalpha <- list()
#RBalphamc <- list()



for(i in seq_len(nsim)){
  
  s <- simulateSRrandom(ao=3, b=1/10000, ER=0.4, fec=c(0,0,0,1,0), sig=.5, siga=.2, nobs=40 )
  s$extinct<-is.na(s$R)
  Sim[[i]]<-s
  

  if(sum(is.na(s$R))>0){
    next;
  }
  srm <- lm(s$logR_S~ s$S)
  plot(s$logR_S~ s$S, main=paste(i))
  abline(srm)
   
  cta[[i]]<-rep(srm$coefficients[[1]],length(s$R))


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
    lmabias[[n]]<-((cta[[n]]-sima[[n]])/sima[[n]]*100)
    dlmabias[[n]]<-((dlmKFalpha[[n]]-sima[[n]])/sima[[n]]*100)
    rbabias[[n]]<-((RBalpha[[n]]-sima[[n]])/sima[[n]]*100)
  }    

}

dfbias <- data.frame(pbias=c(unlist(dlmabias),unlist(rbabias), unlist(lmabias)),
  fit=c("dlm","RB","lm"))

ggplot(dfbias) +
geom_boxplot(aes(x=fit, y=pbias))+
coord_cartesian(ylim = c(-100,100))

lapply(Sim, function(x)x$a)



# look at bias in beta and variance terms

RB[[1]]$sdrep

sdRB <- sapply(RB,  function(x)x$sdrep[which(names(x$sdrep[,2])=="beta"),2])
sddlm <- sapply(dlmKF,  function(x)x$results$beta)

dlmbbias<-rep(NA,nsim)
rbbbias<-rep(NA,nsim)

vs <- 0
for(n in 1:nsim){
  if(valid[n]){
    vs<-vs+1

    dlmabias[[n]]<-((dlmKFalpha[[n]]-sima[[n]])/sima[[n]]*100)
    rbabias[[n]]<-((RBalpha[[n]]-sima[[n]])/sima[[n]]*100)
  }    

}


