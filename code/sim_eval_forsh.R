#!/usr/local/bin/Rscript --slave
args <- commandArgs(trailingOnly=TRUE)



source("dlm-wrapper.R")
source("sim_eval_func.R")
source("sgen_functions.R")

#source("C:/Users/worc/Documents/KF-funcs-appl/HoltMichielsens2020/KFcode.R")



simpars <- read.csv("../data/siminput/randomwalka.csv")

fec <- simpars[,grep("fec",names(simpars))]
fec <- fec[!is.na(fec)]


randsim <- runrandomsims(nsim=simpars$nsim,ao=simpars$ao, b= simpars$b, ER= simpars$ER, fec= fec, sig=simpars$sig, siga=simpars$siga, 
  nobs=simpars$nobs, CapScalar=simpars$CapScalar, plot_progress=TRUE, trend=simpars$trend, lowsca=simpars$lowsca,
  hisca=simpars$hisca, ampsc=simpars$ampsc, seed=args[1])

saveRDS(randsim, "../data/out/rndwlksim.rds")
#randsim <-readRDS("../data/out/randsim.rds")

