#Preliminary batched model support for varying dynamics in sockeye

remotes::install_github("carrieholt/KF-funcs")

library(KFfuncs)
library(here)
library(TMB)

sock_dat<- read.csv(here('data','filtered datasets','sockeye_final.csv'))
sock_info<- read.csv(here('data','sockeye','sockeye_info.csv'))

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
  

  SRdata<-list(obs_logRS=s$logR_S,obs_S=s$spawners, prbeta1=1,
    prbeta2=1)
  
  
  #Model 1 - TMB
  parameters<- list(
    alphao=srm$coefficients[1],
    logSmax = log(1/ifelse(-srm$coefficients[2]<0,1e-08,-srm$coefficients[2])),
    rho=.5,
    logvarphi=0,
    alpha=rep(srm$coefficients[1],length(s$recruits))
    )

  obj <- MakeADFun(SRdata,parameters,DLL="Ricker_tva_Smax_ratiovar",random="alpha")#,lower = -Inf, upper = Inf)
   newtonOption(obj, smartsearch=FALSE)

  opt <- nlminb(obj$par,obj$fn,obj$gr)

  RB[[i]] <- list(alpha=obj$rep()$alpha, sigobs=obj$rep()$sig, siga=obj$rep()$tau, beta=obj$rep()$beta, smax=obj$rep()$Smax)
  RBalpha[[i]] <- obj$rep()$alpha
  
  #Model 2 - tv a and static b
  SRdata2 <- data.frame(byr=s$broodyear,
    spwn=s$spawners,
    rec=s$recruits)
  avarydlm <-fitDLM(SRdata2, alpha_vary = TRUE, beta_vary = FALSE)
  dlmKF[[i]] <- list(alpha=avarydlm$results$alpha, sigobs=obj$rep()$sig, siga=obj$rep()$tau, beta=obj$rep()$beta, smax=obj$rep()$Smax)
  dlmKFalpha[[i]] <- avarydlm$results$alpha

  #Model 3 Carrie's KF
  
  initial <- list()
  initial$mean.a <- srm$coefficients[1]
  initial$var.a <- 1
  initial$b <- -srm$coefficients[2]
  initial$ln.sig.e <- log(1)
  initial$ln.sig.w <- log(1)
  initial$Ts <- 0
  initial$EstB <- TRUE
  
  holtKFfit <- kf.rw(initial=initial,x=s$logR_S,y=s$spawners)

  holtKF[[i]]<-list(alpha=holtKFfit$smoothe.mean.a, sigobs=NA, siga=NA, beta=NA, smax=NA) 
  holtKFalpha[[i]]<-holtKFfit$smoothe.mean.a

}
a

df<- data.frame(alpha=c(unlist(RBalpha),  unlist(dlmKFalpha)),
  type=rep(c("RB","dlm"), each=length(unlist(RBalpha))),
  time=sock_dat$broodyear,
  stock=sock_dat$stock
  )

p <- ggplot(df) +
geom_line(aes(x=time,y=alpha,color=type,lty=type), size=1.5,alpha=.7) +
theme_bw(14)+
facet_wrap(~stock)+
scale_colour_viridis_d(end=.8)+
theme(legend.position="bottom")
p


ggsave("C:/Users/worc/Documents/timevarproject/kf_rb.png",
  plot=p,
  dpi = 400)
?lapply







  
  