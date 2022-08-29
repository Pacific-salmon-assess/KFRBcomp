#==============================================================
#Simple simulation to access lfo functions
#Dan Greenberg Catarina Wor
#August 2022
#==============================================================

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
tvlfotot <- list()



for(i in seq_len(nsim)){

#TMB LFO

  s <- simulateSRrandom(ao=ao, b=b, ER=ER, fec= fec, sig=sig, siga=siga, nobs=nobs, CapScalar=CapScalar )
  s$extinct<-is.na(s$R)

  if(sum(s$extinct)==0){

    dat <- data.frame(logRS=s$logR_S, 
      S=s$S)

    staticlfo[[i]] <- tmb_mod_lfo_cv(data=dat, tv.par=c('static') )

    tvlfo[[i]] <- tmb_mod_lfo_cv(data=dat, tv.par=c('alpha'), siglfo = "obs" )
    tvlfotot[[i]] <- tmb_mod_lfo_cv(data=dat, tv.par=c('alpha'), siglfo = "total" )
  }else{
    staticlfo[[i]] <- NA
    tvlfo[[i]] <- NA
  }

}


unlist(lapply(staticlfo, sum))
unlist(lapply(tvlfo, function(x){sum(x$lastparam)}))
unlist(lapply(tvlfotot, function(x){sum(x$lastparam)}))


complfo<-cbind(unlist(lapply(staticlfo, sum)),
unlist(lapply(tvlfo, function(x){sum(x$lastparam)})))
win_model<-apply(complfo,1,which.max)

sum(win_model==1)/length(win_model)
sum(win_model==2)/length(win_model)


complfotot<-cbind(unlist(lapply(staticlfo, sum)),
unlist(lapply(tvlfotot, function(x){sum(x$lastparam)})))
win_modeltot<-apply(complfotot,1,which.max)


sum(win_modeltot==1)/length(win_model)
sum(win_modeltot==2)/length(win_model)



#simulate from siple Ricker



for(i in seq_len(nsim)){

#TMB LFO

  s <- simulateSRtrend(ao=ao, b=b, ER=ER, fec= fec, sig=sig, siga=siga, nobs=nobs, CapScalar=CapScalar, 
   trend="stable")
  s$extinct<-is.na(s$R)

  if(sum(s$extinct)==0){

    dat <- data.frame(logRS=s$logR_S, 
      S=s$S)

    staticlfo[[i]] <- tmb_mod_lfo_cv(data=dat, tv.par=c('static') )

    tvlfo[[i]] <- tmb_mod_lfo_cv(data=dat, tv.par=c('alpha'), siglfo = "obs" )
    tvlfotot[[i]] <- tmb_mod_lfo_cv(data=dat, tv.par=c('alpha'), siglfo = "total" )
  }else{
    staticlfo[[i]] <- NA
    tvlfo[[i]] <- NA
  }

}



complfosimple<-cbind(unlist(lapply(staticlfo, sum)),
unlist(lapply(tvlfo, function(x){sum(x$lastparam)})))
win_modelsimple<-apply(complfo,1,which.max)

sum(win_modelsimple==1)/length(win_modelsimple)
sum(win_modelsimple==2)/length(win_modelsimple)


complfototsimple<-cbind(unlist(lapply(staticlfo, sum)),
unlist(lapply(tvlfotot, function(x){sum(x$lastparam)})))
win_modeltotsimple<-apply(complfototsimple,1,which.max)


sum(win_modeltotsimple==1)/length(win_modelsimple)
sum(win_modeltotsimple==2)/length(win_modelsimple)



