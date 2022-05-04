#Preliminary batched model support for varying dynamics in sockeye


# I changed the estimation to optim, more of the models converged but estimates are 
# diverging more frm other tools. 
remotes::install_github("carrieholt/KF-funcs")

library(KFfuncs)
library(here)
library(TMB)
library(tmbstan)

sock_dat <- read.csv(here('data','filtered datasets','sockeye_final.csv'))
sock_info <- read.csv(here('data','sockeye','sockeye_info.csv'))

sock_info<- subset(sock_info, Stock.ID %in% sock_dat$stock.id)

source(here("code","dlm-wrapper.R"))
#source("C:/Users/worc/Documents/KF-funcs-appl/HoltMichielsens2020/KFcode.R")

compile("TMBmodels/Ricker_tva_Smax_ratiovar.cpp")
dyn.load(dynlib("TMBmodels/Ricker_tva_Smax_ratiovar"))


holtKF <- list()
dlmKF <- list()
RB <- list()
holtKFalpha <- list()
dlmKFalpha <- list()
RBalpha <- list()
RBalphamc <- list()

for(i in seq_len(nrow(sock_info))){
  
  s <- subset(sock_dat,stock.id==sock_info$Stock.ID[i])
  srm <- lm(s$logR_S~ s$spawners)

  plot(s$logR_S~ s$spawners, main=paste(sock_info$Stock[i], i))
  abline(srm)

  SRdata<-list(obs_logRS=s$logR_S,obs_S=s$spawners, prbeta1=1.5,
    prbeta2=1.5)
  
  
  #Model 1 - TMB
  parameters<- list(
    alphao=srm$coefficients[[1]],
    logSmax = log(1/ifelse(-srm$coefficients[[2]]<0,1e-08,-srm$coefficients[2])),
    rho=.5,
    logvarphi=0,
    alpha=rep(srm$coefficients[1],length(s$recruits))
    )    

  obj <- MakeADFun(SRdata,parameters,DLL="Ricker_tva_Smax_ratiovar",random="alpha")#,lower = -Inf, upper = Inf)
   newtonOption(obj, smartsearch=FALSE)

  opt <- nlminb(obj$par,obj$fn,obj$gr)
  obj$rep()
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
  SRdata2 <- data.frame(byr=s$broodyear,
    spwn=s$spawners,
    rec=s$recruits)
  avarydlm <-fitDLM(data=SRdata2, alpha_vary = TRUE, beta_vary = FALSE)
  dlmKF[[i]] <- list(results=avarydlm$results, alpha=avarydlm$results$alpha, sigobs=avarydlm$sd.est[1], 
    siga=avarydlm$sd.est[2], beta=-avarydlm$results$beta[1], smax=-1/avarydlm$results$beta[1],
     message=avarydlm$message, convergence=avarydlm$convergence)
  dlmKFalpha[[i]] <- avarydlm$results$alpha

  #Model 3 Carrie's KF
  
  initial <- list()
  initial$mean.a <- srm$coefficients[1]
  initial$var.a <- 1
  initial$b <- srm$coefficients[2]
  initial$ln.sig.e <- log(.5)
  initial$ln.sig.w <- log(.5)
  initial$Ts <- 0
  initial$EstB <- "True"
  
  holtKFfit <- kf.rw(initial=initial,x=s$spawners,y=s$logR_S)

  holtKF[[i]]<-list(alpha=holtKFfit$smoothe.mean.a, sigobs=holtKFfit$sig.e, siga=holtKFfit$sig.w, beta=holtKFfit$b, 
    smax=1/holtKFfit$b, convergence=holtKFfit$Report$convergence, message= holtKFfit$Report$message) 
  holtKFalpha[[i]]<-holtKFfit$smoothe.mean.a

}


RBalphamc[[3]]
length(RBalphamc)
length(unlist(RBalpha))
length(unlist(RBalphamc))
length(unlist(dlmKFalpha))
length(unlist(holtKFalpha))
RBalphamc[[i]]

df<- data.frame(alpha=c(unlist(RBalpha), unlist(RBalphamc),  unlist(dlmKFalpha), unlist(holtKFalpha)),
  type=rep(c("RB", "RBmcmcmean", "dlm", "Holt"), each=length(unlist(RBalpha))),
  time=sock_dat$broodyear,
  stock=sock_dat$stock
  )


df<- data.frame(alpha=c(unlist(RBalpha),   unlist(dlmKFalpha), unlist(holtKFalpha)),
  type=rep(c("RB",  "dlm", "Holt"), each=length(unlist(RBalpha))),
  time=sock_dat$broodyear,
  stock=sock_dat$stock
  )


p <- ggplot(df) +
geom_line(aes(x=time,y=alpha,color=type), size=1.5,alpha=.7) +
theme_bw(14)+
facet_wrap(~stock, scales="free")+
scale_colour_viridis_d(end=.8)+
theme(legend.position="bottom")
p


ggsave("../figure/meanfitcomparison.png",
  plot=p,
  dpi = 400)

fit_summary$summary[grep("alpha\\[",rownames(fit_summary$summary)),"mean"]

sdRB <- sapply(RB,  function(x)x$sdrep[which(names(x$sdrep[,2])=="alpha"),2])
lowRBmc <- sapply(RB,  function(x)x$mcmcsummary$summary[grep("alpha\\[",rownames(x$mcmcsummary$summary)),"2.5%"])
highRBmc <- sapply(RB,  function(x)x$mcmcsummary$summary[grep("alpha\\[",rownames(x$mcmcsummary$summary)),"97.5%"])
sddlm <- sapply(dlmKF,  function(x)x$results$alpha_se)


df_ci<- data.frame(alpha=c(unlist(RBalpha), unlist(RBalphamc),  unlist(dlmKFalpha)),
  low=c(unlist(RBalpha)-1.96*unlist(sdRB), unlist(lowRBmc),unlist(dlmKFalpha)- 1.96* unlist(sddlm)),
  high=c(unlist(RBalpha)+1.96*unlist(sdRB), unlist(highRBmc),unlist(dlmKFalpha)+1.96* unlist(sddlm)),
  type=rep(c("RB", "RBmcmcmean", "dlm"), each=length(unlist(RBalpha))),
  time=sock_dat$broodyear,
  stock=sock_dat$stock
  )


pci <- ggplot(df_ci) +
geom_ribbon(aes(x=time,ymin = low, ymax = high, fill=type), alpha=.5) +
geom_line(aes(x=time,y=alpha,color=type), size=1.5) +
facet_wrap(~stock, scales="free")+
theme_bw(14)+
scale_colour_viridis_d(end=.8)+scale_fill_viridis_d(end=.8)+
theme(legend.position="bottom")
pci


dlmfiltalpha <- sapply(dlmKF,  function(x)x$results$alpha_filt)


dlmdf <- data.frame(alpha=c( unlist(dlmKFalpha), unlist(dlmfiltalpha)),
  low=c(unlist(dlmKFalpha)- 1.96* unlist(sddlm), rep(NA,length(unlist(dlmfiltalpha)))),
  high=c(unlist(dlmKFalpha)+1.96* unlist(sddlm), rep(NA,length(unlist(dlmfiltalpha)))),
  type=rep(c("smoothed", "filtered"), each=length(unlist(dlmfiltalpha))),
  time=sock_dat$broodyear,
  stock=sock_dat$stock
  )


pdlm <- ggplot(dlmdf) +
geom_line(aes(x=time,y=alpha,color=type), size=1.5,alpha=.7) +
geom_ribbon(aes(x=time,ymin = low, ymax = high, fill=type), alpha=.5) +
facet_wrap(~stock, scales="free")+
theme_bw(14)+
scale_colour_viridis_d(end=.8)+scale_fill_viridis_d(end=.8)+
theme(legend.position="bottom")
pdlm


dlmKF[[1]]$results


RB[[3]]$sdrep[which(names(RB[[1]]$sdrep[,2])=="alpha"),2]   

(RB[[1]]$mcmcsummary$summary[,"2.5%"])
RB[[1]]$mcmcsummary$summary[,"97.5%"]

avarydlm$results$alpha_se

 grep("alpha\\[",rownames(RB[[1]]$mcmcsummary$summary))
 rownames(RB[[1]]$mcmcsummary$summary) 
which(names(RB[[1]]$sdrep[,2])=="alpha")
"alpha"]]


sapply(RB,  function(x)x[["convergence"]])
sapply(dlmKF,  function(x)x[["convergence"]])
sapply(holtKF,  function(x)x[["convergence"]])


sapply(RB,  function(x)x[["message"]])
sapply(dlmKF,  function(x)x[["message"]])
sapply(holtKF,  function(x)x[["message"]])

names(dlmKF[[1]])


sapply(dlmKF,  function(x)x[["convergence"]])


#TODO
#Compare confidence intervals
#Compare filtered and smoothed estimates


  
  