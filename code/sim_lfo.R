#==============================================================
#Simple simulation to access lfo functions
#Dan Greenberg Catarina Wor
#August 2022
#==============================================================
 #TODO: deal with  excess variability in simulated SRtime-seies
#need to fix this file
#remotes::install_github("carrieholt/KF-funcs")
#library(KFfuncs)
#use this one instead
#source("C:/Users/worc/Documents/KF-funcs-appl/HoltMichielsens2020/KFcode.R")
remotes::install_git('https://github.com/Pacific-salmon-assess/samEst')


if(!"cowplot" %in% rownames(installed.packages())){
  install.packages("cowplot")
}
if(!"ggplot2" %in% rownames(installed.packages())){
  install.packages("ggplot2")
}





library(tmbstan)
library(cowplot)
library(ggplot2)
library(samEst)

#source("dlm-wrapper.R")
source("sim_eval_func.R")
source("sgen_functions.R")

#source("C:/Users/worc/Documents/KF-funcs-appl/HoltMichielsens2020/KFcode.R")


ao <- 1.3
b <- 1/150000
nsim <- 100
Smax <- 1/b
sig <- .5
siga <- .3
ER=0.0
nobs=40
CapScalar=5
fec= c(0,0,0,1,1)


staticlfo <- list()
tvlfo <- list()


for(i in seq_len(nsim)){

#TMB LFO

  s <- simulateSRrandom(ao=ao, b=b, ER=ER, fec= fec, sig=sig, siga=siga, nobs=nobs, CapScalar=CapScalar )
  s$extinct<-is.na(s$R)

  if(sum(s$extinct)==0){

    dat <- data.frame(logRS=s$logR_S, 
      S=s$S)

    staticlfo[[i]] <- tmb_mod_lfo_cv(data=dat, tv.par=c('static') )

    tvlfo[[i]] <- tmb_mod_lfo_cv(data=dat, tv.par=c('alpha') )
  }else{
    staticlfo[[i]] <- NA
    tvlfo[[i]] <- NA
  }

}


unlist(lapply(staticlfo, sum))
unlist(lapply(tvlfo, function(x){sum(x$lastparam)}))


comp<-cbind(unlist(lapply(staticlfo, sum)),
unlist(lapply(tvlfo, function(x){sum(x$lastparam)})))
win_model<-apply(comp,1,which.max)

sum(win_model==1)/length(win_model)


