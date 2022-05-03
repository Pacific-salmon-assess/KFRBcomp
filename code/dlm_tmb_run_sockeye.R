#Preliminary batched model support for varying dynamics in sockeye



#remotes::install_github("carrieholt/KF-funcs")

library(KFfuncs)
library(here)
library(TMB)
library(tmbstan)

sock_dat <- read.csv(here('data','filtered datasets','sockeye_final.csv'))
sock_info <- read.csv(here('data','sockeye','sockeye_info.csv'))

sock_info<- subset(sock_info, Stock.ID %in% sock_dat$stock.id)

source(here("code","dlm-wrapper.R"))

compile("TMBmodels/Ricker_tva_Smax_ratiovar.cpp")
dyn.load(dynlib("TMBmodels/Ricker_tva_Smax_ratiovar"))




holtKF <- list()
dlmKF <- list()
RB <- list()
holtKFalpha <- list()
dlmKFalpha <- list()
RBalpha <- list()

for(i in seq_len(nrow(sock_info))){
  
  s <- subset(sock_dat,stock.id==sock_info$Stock.ID[i])
  srm <- lm(s$logR_S~ s$spawners)

  plot(s$logR_S~ s$spawners)
  abline(srm)

  SRdata<-list(obs_logRS=s$logR_S,obs_S=s$spawners, prbeta1=2,
    prbeta2=2)
  
  
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

  sdrep <- summary(sdreport(obj))
  
  #MCMC
  #fitmcmc1 <- tmbstan(obj, chains=3,
  #            iter=1000, init="random",
  #            lower=c(-10,-20,0,-6),
  #             upper=c(5,0,1,6),
  #             control = list(adapt_delta = 0.98))
  #
  #  mc <- extract(fitmcmc1, pars=names(obj$par),
  #            inc_warmup=TRUE, permuted=FALSE)
  #
  #fit_summary <- summary(fitmcmc1)   
  #fit_summary$summary[grep("alpha\\[",rownames(fit_summary$summary)),]

  RB[[i]] <- list(sdrep=sdrep, convergence=opt$convergence, message=opt$message)
  RBalpha[[i]] <- sdrep[which(rownames(sdrep)=="alpha"),1]
  
  #Model 2 - tv a and static b
  SRdata2 <- data.frame(byr=s$broodyear,
    spwn=s$spawners,
    rec=s$recruits)
  avarydlm <-fitDLM(SRdata2, alpha_vary = TRUE, beta_vary = FALSE)
  dlmKF[[i]] <- list(alpha=avarydlm$results$alpha, sigobs=avarydlm$sd.est[1], 
    siga=avarydlm$sd.est[2], beta=-avarydlm$results$beta[1], smax=-1/avarydlm$results$beta[1],
     message=avarydlm$message, convergence=avarydlm$convergence)
  dlmKFalpha[[i]] <- avarydlm$results$alpha

  #Model 3 Carrie's KF
  
  initial <- list()
  initial$mean.a <- srm$coefficients[1]
  initial$var.a <- 1
  initial$b <- -srm$coefficients[2]
  initial$ln.sig.e <- log(avarydlm$sd.est[1])
  initial$ln.sig.w <- log(avarydlm$sd.est[2])
  initial$Ts <- 0
  initial$EstB <- TRUE
  
  holtKFfit <- kf.rw(initial=initial,x=s$spawners,y=s$logR_S)

  holtKF[[i]]<-list(alpha=holtKFfit$smoothe.mean.a, sigobs=NA, siga=NA, beta=NA, 
    smax=NA, convergence=holtKFfit$Report$convergence, message= holtKFfit$Report$message) 
  holtKFalpha[[i]]<-holtKFfit$smoothe.mean.a

}




df<- data.frame(alpha=c(unlist(RBalpha),  unlist(dlmKFalpha), unlist(holtKFalpha)),
  type=rep(c("RB","dlm", "Holt"), each=length(unlist(RBalpha))),
  time=sock_dat$broodyear,
  stock=sock_dat$stock
  )

p <- ggplot(df) +
geom_line(aes(x=time,y=alpha,color=type,lty=type), size=1.5,alpha=.7) +
theme_bw(14)+
facet_wrap(~stock, scales="free")+
scale_colour_viridis_d(end=.8)+
theme(legend.position="bottom")
p




sapply(RB,  function(x)x[["convergence"]])
sapply(dlmKF,  function(x)x[["convergence"]])
sapply(holtKF,  function(x)x[["convergence"]])

sapply(RB,  function(x)x[["message"]])
sapply(dlmKF,  function(x)x[["message"]])
sapply(holtKF,  function(x)x[["message"]])


ggsave("C:/Users/worc/Documents/timevarproject/kf_rb.png",
  plot=p,
  dpi = 400)








  
  