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

if(!"cowplot" %in% rownames(installed.packages())){
  install.packages("cowplot")
}
if(!"ggplot2" %in% rownames(installed.packages())){
  install.packages("ggplot2")
}





library(tmbstan)
library(cowplot)
library(ggplot2)

source("dlm-wrapper.R")
source("sim_eval_func.R")
source("sgen_functions.R")

#source("C:/Users/worc/Documents/KF-funcs-appl/HoltMichielsens2020/KFcode.R")



ao <- 1.3
b <- 1/150000


ps <- seq(0,150000*3,1000)
pr <- ps*exp(ao-b*ps)
plot(ps, pr)
abline(1,1)



#fec= c(0,.1,.3,.5,.1)
sr <- simulateSRtrend(ao=1.3, b=1/150000, ER=0.0,fec=fec , sig=.5, siga=.2, nobs=40,
    CapScalar=5, trend="sine", lowsca=.5, hisca=2, ampsc=.5 )
sr
par(mfrow=c(2,1))
plot(sr$S,sr$R, xlim=c(0,35000*5), ylim=c(0,180000))
abline(a=0,b=1)
lines(ps,pr,col="red")
plot(sr$a,type="b")



sr <- simulateSRtrend(ao=1.3, b=1/150000, ER=0.0, fec= c(0,.1,.3,.5,.1), sig=.5, siga=.2, nobs=40,
    CapScalar=5, trend="regime",lowsca=.5,hisca=2, ampsc=.5 )
sr
par(mfrow=c(2,1))
plot(sr$S,sr$R, xlim=c(0,35000*5), ylim=c(0,180000))
abline(a=0,b=1)
lines(ps,pr,col="red")
plot(sr$a,type="b")



#==========================================================


#simulation estimation random walk
nsim <- 10
Smax <- 1/b
sig <- .5
siga <- .2
#ER=0.0
#nobs=40
#CapScalar=5
#fec= c(0,.1,.3,.5,.1)
randsim <- runrandomsims(nsim=nsim,ao=ao, b=1/Smax, ER=0.0, fec= c(0,.1,.3,.5,.1), sig=sig, siga=siga, nobs=40, CapScalar=5,
  plot_progress=TRUE, trend="random walk", lowsca=.5,hisca=2, ampsc=.5 )

#ER=0.0
#saveRDS(randsim, "../data/out/randsim.rds")
#randsim <-readRDS("../data/out/randsim.rds")



#Filter onlyestimates that converged 

#add filtered estimates to the plots


dfbiasrand1<-calculatepbias(simresult=randsim,Smax=Smax, sig=sig,siga=siga, Bayesstat="mean")
 head(dfbiasrand1)


dfbiasrand1p<-dfbiasrand1[dfbiasrand1$convergence==0&!is.na(dfbiasrand1$convergence==0),]
dfbiasrand1p<-dfbiasrand1p[dfbiasrand1p$fit!="lm",]
dfbiasrand1p$fit<-ordered(dfbiasrand1p$fit, c("dlm","tmbholtKF","holtKF", "RB", "RBmcmc","stanrb", "stanGp"))

unique(dfbiasrand1p$param)
head(dfbiasrand1p)
#add


give.n <- function(x){
  return(c(y = median(x)*1.05, label = length(x))) 
  # experiment with the multiplier to find the perfect position
} 


prand <- ggplot(dfbiasrand1p,aes(x=fit, y=pbias)) +
geom_boxplot()+
coord_cartesian(ylim = c(-100,100))+
geom_hline(yintercept=0) +
theme_bw(14)+
stat_summary(fun.data = give.n, geom = "text", hjust = 0.5,
    vjust = 0.9)+
facet_wrap(~param)
prand



ggsave("../figure/randsimbias.pdf")




ca <- comparealpha(simresult=randsim)
names(ca)



arand <- ggplot(ca,aes(x=fit, y=abias)) +
geom_boxplot()+
coord_cartesian(ylim = c(-100,100))+
geom_hline(yintercept=0) +
theme_bw(14)+
stat_summary(fun.data = give.n, geom = "text", hjust = 0.5,
    vjust = 0.9)
arand

unique(ca$fit)

filt_diff<-ca$a[ca$fit=="smoothKF"]-ca$a[ca$fit=="filterKF"]
rb_diff<-ca$a[ca$fit=="smoothKF"]-ca$a[ca$fit=="RB"]
rbfilt_diff<-ca$a[ca$fit=="RB"]-ca$a[ca$fit=="filterKF"]
diffca<-data.frame(diff=c(filt_diff,rb_diff,rbfilt_diff),
  type=rep(c("diff smooth and filt", "diff smooth and rb", "diff rb and filt"),each=length(filt_diff)),
  yr=1:41,
  sim=rep(1:100,each=41))



adiff <- ggplot(diffca,aes(x=type, y=diff)) +
geom_boxplot()+
geom_hline(yintercept=0) +
theme_bw(14)+
stat_summary(fun.data = give.n, geom = "text", hjust = 0.5,
    vjust = 0.9)+
ylab("absolute difference")
adiff


df<-diffca[diffca$type=="diff smooth and filt",]
df<-diffca[diffca$type=="diff smooth and rb",]


linesdiff <- ggplot(diffca,aes(x=yr, y=diff)) +
geom_line(aes( group=sim), alpha=.7)+
theme_bw(14) +
facet_wrap(~type)
linesdiff

summary(rbfilt_diff)
quantile(rbfilt_diff, c(0.025,.1,.5,.9,.975), na.rm=TRUE)

length(rbfilt_diff)
sum(rbfilt_diff>0.,na.rm=TRUE)

summary(rb_diff)
summary(filt_diff)
#look at changepoint and bcp

sss<-rep(NA,44)
rrr<-rep(NA,44)

sss[1:4]<-ao/b
for(i in 5:44){
  rrr[i]<-sss[i-4]*exp(ao-b*sss[i-4])*exp(rnorm(1,0,.5)-(.5^2)/2)
  sss[i]<- rrr[i]
}

ps<-seq(0,30000*5,1000)
pr<-ps*exp(ao-b*ps)
plot(ps, pr,type="l")

points(sss[1:40],rrr[5:44])
abline(1,1)



sr <- simulateSRrandom(ao=ao,  b=1/Smax, ER=0.0, fec= c(0,.1,.3,.5,.1), sig=sig, siga=siga, nobs=40, CapScalar=5 )
sr
fit1 <- lm(logR_S~S,data=sr)
names(summary(fit1))
summary(fit1)$sigma
plot(resid(fit1))


#regime shift methods
STARS(regLN=10,sig=summary(fit1)$sigma,series=resid(fit1),shift=1)

uu<-changepoint::cpt.meanvar(resid(fit1)[!is.na(resid(fit1))],method="PELT",test.stat="Normal",penalty="AIC",minseglen=10)
bpts<-c(1,uu@cpts)
    sds<-sqrt(uu@param.est$variance)
    means<-exp(uu@param.est$mean)*exp(uu@param.est$variance/2)
    tmpFit$PeltMean<-NA
    tmpFit$PeltSd<-NA




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


dfbiasdec<-calculatepbias(decsim, Smax=Smax)

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


#=====================
#From ERa maturity to prop mat at agrepl
a1=0.1476287
a2=0.3636234*(1-0.1476287)
a3=0.9647583*(1-(0.1476287+0.3636234*(1-0.1476287)))
a4=1-(0.1476287+ 0.3636234*(1-0.1476287) + 0.9647583*(1-(0.1476287+0.3636234*(1-0.1476287))))

c(a1, a2,a3,a4)