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
nsim <- 500
Smax <- 1/b
sig <- .5
siga <- .3
ER <- 0.3
nobs <- 40
CapScalar <- 5
fec <- c(0,0,0,1,0)



staticlfo<- list()
tvlfoa <- list()
tvlfotota <- list()
tvlfob <- list()
tvlfototb <- list()



for(i in seq_len(nsim)){

#TMB LFO

  s <- simulateSRrandom(ao=ao, b=b, ER=ER, fec= fec, sig=sig, siga=siga, nobs=nobs, CapScalar=CapScalar )
  s$extinct<-is.na(s$R)

  if(sum(s$extinct)==0){

    dat <- data.frame(logRS=s$logR_S, 
      S=s$S)

    staticlfo[[i]] <- tmb_mod_lfo_cv(data=dat, tv.par=c('static') )

    tvlfoa[[i]] <- tmb_mod_lfo_cv(data=dat, tv.par=c('alpha'), siglfo = "obs" )
    tvlfotota[[i]] <- tmb_mod_lfo_cv(data=dat, tv.par=c('alpha'), siglfo = "total" )

    tvlfob[[i]] <- tmb_mod_lfo_cv(data=dat, tv.par=c('alpha'), siglfo = "obs" )
    tvlfototb[[i]] <- tmb_mod_lfo_cv(data=dat, tv.par=c('alpha'), siglfo = "total" )
  }else{
    staticlfo[[i]] <- NA
    tvlfoa[[i]]<- NA
    tvlfotota[[i]] <- NA
    tvlfob[[i]] <- NA
    tvlfototb[[i]] <- NA

  }

}


unlist(lapply(staticlfo, sum))
unlist(lapply(tvlfoa, function(x){if(sum(is.na(x[1]))==0){sum(x$lastparam)}else{NA}}))
unlist(lapply(tvlfotota, function(x){if(sum(is.na(x[1]))==0){sum(x$lastparam)}else{NA}}))
unlist(lapply(tvlfob,function(x){if(sum(is.na(x[1]))==0){sum(x$lastparam)}else{NA}}))
unlist(lapply(tvlfototb, function(x){if(sum(is.na(x[1]))==0){sum(x$lastparam)}else{NA}}))


complfoa<-cbind(unlist(lapply(staticlfo, sum)),
unlist(lapply(tvlfoa, function(x){if(sum(is.na(x[1]))==0){sum(x$lastparam)}else{NA}})),
unlist(lapply(tvlfob, function(x){if(sum(is.na(x[1]))==0){sum(x$lastparam)}else{NA}})))
win_model<-apply(complfoa,1,which.max)

sum(win_model==1,na.rm=T)/length(win_model)
sum(win_model==2,na.rm=T)/length(win_model)
sum(win_model==3,na.rm=T)/length(win_model)




complfotot<-cbind(unlist(lapply(staticlfo, sum)),
unlist(lapply(tvlfotota, function(x){if(sum(is.na(x[1]))==0){sum(x$lastparam)}else{NA}})),
unlist(lapply(tvlfototb, function(x){if(sum(is.na(x[1]))==0){sum(x$lastparam)}else{NA}})))
win_modeltot<-apply(complfotot,1,which.max)


sum(win_modeltot==1,na.rm=T)/length(win_modeltot)
sum(win_modeltot==2,na.rm=T)/length(win_modeltot)
sum(win_modeltot==3,na.rm=T)/length(win_modeltot)



#simulate from simple Ricker
staticlfosimple <- list()

tvlfosimplea <- list()
tvlfototsimplea<- list()

tvlfosimpleb <- list()
tvlfototsimpleb <- list()

for(i in seq_len(nsim)){

#TMB LFO
sigtot<-sqrt(sig^2+siga^2)
  s <- simulateSRtrend(ao=ao, b=b, ER=ER, fec= fec, sig=sigtot, siga=siga, nobs=nobs, CapScalar=CapScalar, 
   trend="stable")
  s$extinct<-is.na(s$R)

  if(sum(s$extinct)==0){

    dat <- data.frame(logRS=s$logR_S, 
      S=s$S)

    staticlfosimple[[i]] <- tmb_mod_lfo_cv(data=dat, tv.par=c('static') )

    tvlfosimplea[[i]] <- tmb_mod_lfo_cv(data=dat, tv.par=c('alpha'), siglfo = "obs" )
    tvlfototsimplea[[i]] <- tmb_mod_lfo_cv(data=dat, tv.par=c('alpha'), siglfo = "total" )
    tvlfosimpleb[[i]] <- tmb_mod_lfo_cv(data=dat, tv.par=c('beta'), siglfo = "obs" )
    tvlfototsimpleb[[i]] <- tmb_mod_lfo_cv(data=dat, tv.par=c('beta'), siglfo = "total" )

  
  

  }else{
    staticlfo[[i]] <- NA
    tvlfosimplea[[i]] <- NA
    tvlfototsimplea[[i]]  <- NA
    tvlfosimpleb[[i]]  <- NA
    tvlfototsimpleb[[i]] <- NA
  }

}



complfosimple<-cbind(unlist(lapply(staticlfosimple, sum)),
unlist(lapply(tvlfosimplea, function(x){sum(x$lastparam)})),
unlist(lapply(tvlfosimpleb, function(x){sum(x$lastparam)})))
win_modelsimple<-apply(complfosimple,1,which.max)

sum(win_modelsimple==1)/length(win_modelsimple)
sum(win_modelsimple==2)/length(win_modelsimple)
sum(win_modelsimple==3)/length(win_modelsimple)


complfototsimple<-cbind(unlist(lapply(staticlfo, sum)),
unlist(lapply(tvlfototsimplea, function(x){sum(x$lastparam)})),
unlist(lapply(tvlfototsimpleb, function(x){sum(x$lastparam)})))
win_modeltotsimple<-apply(complfototsimple,1,which.max)


sum(win_modeltotsimple==1)/length(win_modeltotsimple)
sum(win_modeltotsimple==2)/length(win_modeltotsimple)
sum(win_modeltotsimple==3)/length(win_modeltotsimple)


df<-data.frame( sim=c("sim_tva","sim_tva","sim_tva","sim_tva","sim_tva","sim_tva",
  "sim_simple","sim_simple","sim_simple","sim_simple","sim_simple","sim_simple"),
  lfo_sig=c("lfosig_obs", "lfosig_obs", "lfosig_obs","lfosig_tot","lfosig_tot","lfosig_tot",
    "lfosig_obs", "lfosig_obs", "lfosig_obs","lfosig_tot","lfosig_tot","lfosig_tot"),
  chosen_est=c("simple", "tva", "tvb","simple", "tva", "tvb","simple","tva","tvb","simple","tva","tvb"),
  assignment=c("wrong","correct","wrong","wrong","correct","wrong",
    "correct","wrong","wrong","correct","wrong" ,"wrong"),
  value=c(sum(win_model==1)/length(win_model)*100,
        sum(win_model==2)/length(win_model)*100,
        sum(win_model==3)/length(win_model)*100,
        sum(win_modeltot==1)/length(win_modeltot)*100,
        sum(win_modeltot==2)/length(win_modeltot)*100,
        sum(win_modeltot==3)/length(win_modeltot)*100,
        sum(win_modelsimple==1)/length(win_modelsimple)*100,
        sum(win_modelsimple==2)/length(win_modelsimple)*100,
        sum(win_modelsimple==3)/length(win_modelsimple)*100,
        sum(win_modeltotsimple==1)/length(win_modeltotsimple)*100,
        sum(win_modeltotsimple==2)/length(win_modeltotsimple)*100
        sum(win_modeltotsimple==3)/length(win_modeltotsimple)*100
        )
  )

pp<-ggplot(df)+
geom_bar(aes(x=chosen_est,y=value, fill=assignment),stat="identity")+
facet_grid(sim~lfo_sig)+
theme_bw(16) +
scale_fill_viridis_d(option="A", end=0.7)
pp


#simulate from tv beta




#sigb needs to be smaller to provide sensible results

CapScalarb<-10


staticblfo <- list()
tvblfoa <- list()
tvblfotota <- list()

tvblfob <- list()
tvblfototb <- list()



for(i in seq_len(nsim)){

#TMB LFO
  
  s <- simulateSRrandom_tvb (a=ao, bo=b, ER=ER, fec= fec, sig=sig,  nobs=nobs, CapScalar=CapScalarb, scb=.75,bcycle=10)
  s$extinct<-is.na(s$R)
  
  plot(s$S,s$R)


  if(sum(s$extinct)==0 & sum(s$R==ao/b*CapScalarb) ==0){

    dat <- data.frame(logRS=s$logR_S, 
      S=s$S)

    staticblfo[[i]] <- tmb_mod_lfo_cv(data=dat, tv.par=c('static') )

    tvblfob[[i]] <- tmb_mod_lfo_cv(data=dat, tv.par=c('beta'), siglfo = "obs" )
    tvblfoa[[i]] <- tmb_mod_lfo_cv(data=dat, tv.par=c('alpha'), siglfo = "obs" )

    tvblfototb[[i]] <- tmb_mod_lfo_cv(data=dat, tv.par=c('beta'), siglfo = "total" )
    tvblfotota[[i]] <- tmb_mod_lfo_cv(data=dat, tv.par=c('alpha'), siglfo = "total" )
  }else{
    staticblfo[[i]] <- NA
    tvblfo[[i]] <- NA
    tvblfotot[[i]] <- NA
  }

}



compblfo<-cbind(unlist(lapply(staticblfo, sum)),
unlist(lapply(tvblfo, function(x){if(sum(is.na(x[1]))==0){sum(x$lastparam)}else{NA}})))
win_modelb<-apply(compblfosimple,1,which.max)

sum(win_modelbsimple==1,na.rm=T)/length(win_modelbsimple)
sum(win_modelbsimple==2,na.rm=T)/length(win_modelbsimple)


compblfotot<-cbind(unlist(lapply(staticblfo, sum)),
unlist(lapply(tvblfotot, function(x){if(sum(is.na(x[1]))==0){sum(x$lastparam)}else{NA}})))
win_modeltotbsimple<-apply(compblfotot,1,which.max)


sum(win_modeltotbsimple==1,na.rm=T)/length(win_modeltotbsimple)
sum(win_modeltotbsimple==2,na.rm=T)/length(win_modeltotbsimple)





df<-data.frame( sim=c("sim_tva","sim_tva","sim_tva","sim_tva","sim_simple","sim_simple","sim_simple","sim_simple","sim_tvb","sim_tvb","sim_tvb","sim_tvb"),
  lfo_sig=c("lfosig_obs", "lfosig_obs","lfosig_tot","lfosig_tot","lfosig_obs", "lfosig_obs","lfosig_tot","lfosig_tot","lfosig_obs", "lfosig_obs","lfosig_tot","lfosig_tot"),
  chosen_est=c("simple","tv","simple","tv","simple","tv","simple","tv","simple","tv","simple","tv"),
  assignment=c("wrong","correct","wrong","correct","correct","wrong","correct","wrong","wrong","correct","wrong","correct" ),
  value=c(sum(win_model==1)/length(win_model)*100,
        sum(win_model==2)/length(win_model)*100,
        sum(win_modeltot==1)/length(win_modeltot)*100,
        sum(win_modeltot==2)/length(win_modeltot)*100,
        sum(win_modelsimple==1)/length(win_modelsimple)*100,
        sum(win_modelsimple==2)/length(win_modelsimple)*100,
        sum(win_modeltotsimple==1)/length(win_modeltotsimple)*100,
        sum(win_modeltotsimple==2)/length(win_modeltotsimple)*100,
        sum(win_modelbsimple==1,na.rm=T)/length(win_modelbsimple)*100,
        sum(win_modelbsimple==2,na.rm=T)/length(win_modelbsimple)*100,
        sum(win_modeltotbsimple==1,na.rm=T)/length(win_modeltotbsimple)*100,
        sum(win_modeltotbsimple==2,na.rm=T)/length(win_modeltotbsimple)*100

        )
  )

pp<-ggplot(df)+
geom_bar(aes(x=chosen_est,y=value, fill=assignment),stat="identity")+
facet_grid(sim~lfo_sig)+
theme_bw(16) +
scale_fill_viridis_d(option="A", end=0.7)
pp