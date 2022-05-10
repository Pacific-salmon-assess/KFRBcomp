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
library(cowplot)

source(here("code","dlm-wrapper.R"))
source(here("code","sim_eval_func.R"))

#source("C:/Users/worc/Documents/KF-funcs-appl/HoltMichielsens2020/KFcode.R")

compile("TMBmodels/Ricker_tva_Smax_ratiovar.cpp")
dyn.load(dynlib("TMBmodels/Ricker_tva_Smax_ratiovar"))


ao <- 2.5
b<-1/30000


ps<-seq(0,100000*5,1000)
pr<-ps*exp(ao-b*ps)

sr <- simulateSRrandom(ao=2.5, b=1/30000, ER=0.0, fec= c(0,.1,.3,.5,.1), sig=.5, siga=.2, nobs=40, CapScalar=5 )
sr
par(mfrow=c(2,1))
plot(sr$S,sr$R, xlim=c(0,35000*5), ylim=c(0,180000))
abline(a=0,b=1)
lines(ps,pr,col="red")
plot(sr$a,type="b")



sr <- simulateSRtrend(ao=2.5, b=1/30000, ER=0.0, fec= c(0,.1,.3,.5,.1), sig=.5, siga=.2, nobs=40,
    CapScalar=5, trend="sine",lowsca=.5,hisca=2, ampsc=.5 )
sr
par(mfrow=c(2,1))
plot(sr$S,sr$R, xlim=c(0,35000*5), ylim=c(0,180000))
abline(a=0,b=1)
lines(ps,pr,col="red")
plot(sr$a,type="b")



sr <- simulateSRtrend(ao=2.5, b=1/30000, ER=0.0, fec= c(0,.1,.3,.5,.1), sig=.5, siga=.2, nobs=40,
    CapScalar=5, trend="regime",lowsca=.5,hisca=2, ampsc=.5 )
sr
par(mfrow=c(2,1))
plot(sr$S,sr$R, xlim=c(0,35000*5), ylim=c(0,180000))
abline(a=0,b=1)
lines(ps,pr,col="red")
plot(sr$a,type="b")



#==========================================================


#simulation estimation random walk
nsim <- 100
Smax <- 1/b
sig <- .5
siga <- .2
randsim <- runrandomsims(nsim=nsim,ao=ao, b=1/Smax, ER=0.0, fec= c(0,.1,.3,.5,.1), sig=sig, siga=siga, nobs=40, CapScalar=5)


dfbiasrand1<-calculatepbias(randsim)


prand<-ggplot(dfbias) +
geom_boxplot(aes(x=fit, y=pbias))+
coord_cartesian(ylim = c(-100,100))+
geom_hline(yintercept=0) +
theme_bw(14)+
facet_wrap(~param)


#===============================================================================================
#Trends sim eval


nsim<-100

nsim <- 100
Smax <- 1/b
sig <- .5
siga <- .2


#decline
decsim <- runtrendsims(nsim=nsim,ao=ao, b=1/Smax, ER=0.0, fec= c(0,.1,.3,.5,.1), sig=sig, siga=siga,
 nobs=40, CapScalar=5,trend="decline",lowsca=.5,hisca=2, ampsc=.5)


dfbiasdec<-calculatepbias(decsim)

pdec<-ggplot(dfbiasdec) +
geom_boxplot(aes(x=fit, y=pbias))+
coord_cartesian(ylim = c(-100,100))+
geom_hline(yintercept=0) +
theme_bw(14)+
facet_wrap(~param)


#====================================================================================

#sine
sinesim <- runtrendsims(nsim=nsim,ao=ao, b=1/Smax, ER=0.0, fec= c(0,.1,.3,.5,.1), sig=sig, siga=siga,
 nobs=40, CapScalar=5,trend="sine",lowsca=.5,hisca=2, ampsc=.5)

dfbiassine<-calculatepbias(sinesim)


psine<-ggplot(dfbiassine) +
geom_boxplot(aes(x=fit, y=pbias))+
coord_cartesian(ylim = c(-100,100))+
geom_hline(yintercept=0) +
theme_bw(14)+
facet_wrap(~param)


#==========================================================================
#regime shifts
regimesim <- runtrendsims(nsim=nsim,ao=ao, b=1/Smax, ER=0.0, fec= c(0,.1,.3,.5,.1), sig=sig, siga=siga,
 nobs=40, CapScalar=5,trend="regime",lowsca=.5,hisca=2, ampsc=.5)



dfbiasregime<-calculatepbias(regimesim)



pregime<-ggplot(dfbiasregime) +
geom_boxplot(aes(x=fit, y=pbias))+
coord_cartesian(ylim = c(-100,100))+
geom_hline(yintercept=0) +
theme_bw(14)+
facet_wrap(~param)



allsim<-plot_grid(prand,pdec, psine,pregime, labels = c('random', 'decline','sine','regime'))

ggsave("../figure/simscomp.png")
#==========================================
#old code


sima<- lapply(randsim$Simtrend, function(x)x$a)
valid <- unlist(lapply(randsim$Simtrend, function(x)sum(x$extinct)))<1

#bias of estimators
lmabias<-list()
dlmabias<-list()
rbabias<-list()
vs <- 0
for(n in 1:nsim){
  if(valid[n]){
    vs<-vs+1
    lmabias[[n]]<-((randsim$ctatrend[[n]]-sima[[n]])/sima[[n]]*100)
    dlmabias[[n]]<-((randsim$dlmKFalphatrend[[n]]-sima[[n]])/sima[[n]]*100)
    rbabias[[n]]<-((randsim$RBalphatrend[[n]]-sima[[n]])/sima[[n]]*100)
  }    

}


dfabias <- data.frame(pbias=c(unlist(dlmabias),unlist(rbabias), unlist(lmabias)),
  fit=rep(c("dlm","RB","lm"),each=length(unlist(dlmabias))),
  param="a")


#ggplot(dfabias) +
#geom_boxplot(aes(x=fit, y=pbias))+
#coord_cartesian(ylim = c(-100,100))

SmaxRB <- sapply(randsim$RBtrend,  function(x)ifelse(is.null(x),NA,(exp(x$sdrep["logSmax",1])-Smax)/Smax*100))
Smaxdlm <- sapply(randsim$dlmKFtrend,  function(x)ifelse(is.null(x),NA,((1/-x$results$beta[1])-Smax)/Smax*100))
ctsmax <- unlist(randsim$ctsmaxtrend)
ctsmax[ctsmax<0]<-NA
Smaxct <- (ctsmax-Smax)/Smax*100


dfsmaxbias <- data.frame(pbias=c(Smaxdlm,SmaxRB, Smaxct),
  fit=rep(c("dlm","RB","lm"),each=length(Smaxdlm)),
  param="Smax")

#ggplot(dfsmaxbias) +
#geom_boxplot(aes(x=fit, y=pbias))+
#coord_cartesian(ylim = c(-150,150))

sigRB <- sapply(randsim$RBtrend,  function(x)ifelse(is.null(x),NA,(x$sdrep["sig",1]-sig)/sig*100))
sigdlm <- sapply(randsim$dlmKFtrend,  function(x)ifelse(is.null(x),NA,(x$sigobs-sig)/sig*100))
sigct <- (unlist(randsim$ctsigtrend)-sig)/sig*100


dfsigbias <- data.frame(pbias=c(sigdlm,sigRB,sigct),
  fit=rep(c("dlm","RB", "lm"),each=length(sigRB)),
  param="sig")

#ggplot(dfsigbias) +
#geom_boxplot(aes(x=fit, y=pbias))+
#coord_cartesian(ylim = c(-150,150))

sigaRB <- sapply(randsim$RBtrend,  function(x)ifelse(is.null(x),NA,(x$sdrep["tau",1]-siga)/siga*100))
sigadlm <- sapply(randsim$dlmKFtrend,  function(x)ifelse(is.null(x),NA,(x$siga-siga)/siga*100))

dfsigabias <- data.frame(pbias=c(sigadlm,sigaRB),
  fit=rep(c("dlm","RB"),each=length(sigRB)),
  param="siga")

#ggplot(dfsigabias) +
#geom_boxplot(aes(x=fit, y=pbias))+
#coord_cartesian(ylim = c(-150,150))


dfbias<-rbind(dfabias,
  dfsmaxbias, 
  dfsigbias,
  dfsigabias)



prand<-ggplot(dfbias) +
geom_boxplot(aes(x=fit, y=pbias))+
coord_cartesian(ylim = c(-100,100))+
geom_hline(yintercept=0) +
theme_bw(14)+
facet_wrap(~param)
