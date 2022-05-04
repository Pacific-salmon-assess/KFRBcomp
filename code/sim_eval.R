#==============================================================
#Simple simulation evaluation of Ricker model with time-varying parameters
#Catarina Wor
#May 2022
#==============================================================




ao<-2.5
b<-1/10000
ER<-0.40
fec= c(0,.2,.4,.8,1)
sig <- .5
siga <-.3


sr<-simulateSRrandom()
plot(sr$S,log(sr$S/sr$R))

#TODO check this function
#


simulateSR <random- function(ao=2.5, b=1/10000, ER=0.4, fec=c(0,0,0,1,1), sig=.5, siga=.3, nobs=40 ){

    yrs <- 1:70
    S <- NULL
    R <- NULL
    a <- NULL
    Ra <- matrix(NA, ncol=length(fec),nrow=yrs)

    a[1:5] <- ao
    Seq <- log(ao)/b
    S[1] <- Seq
    R[1] <- ao*Seq*exp(-b*Seq)

  
    for(y in 2:length(fec)){
      S[y]<-0    
      for(j in seq_along(fec)){
        if((y-j)>0){
          S[y] <- S[y] + R[y-j]*(1-ER)*fec[j]
        }else{
          S[y] <- S[y] + R[1]*(1-ER)*fec[j]
        }
      }
      R[y]<-a[y]*S[y]*exp(-b*S[y])        
    }
  
    #S[2]<-sum(R[1]*(1-ER)*fec[1],R[1]*(1-ER)*fec[2],R[1]*(1-ER)*fec[3],R[1]*(1-ER)*fec[4],R[1]*(1-ER)*fec[5])
    #R[2]<-ao*S[2]*exp(-b*S[2])
    #S[3]<-sum(R[2]*(1-ER)*fec[1],R[1]*(1-ER)*fec[2],R[1]*(1-ER)*fec[3],R[1]*(1-ER)*fec[4],R[1]*(1-ER)*fec[5])
    #R[3]<-ao*S[3]*exp(-b*S[3])
    #S[4]<-sum(R[3]*(1-ER)*fec[1],R[2]*(1-ER)*fec[2],R[1]*(1-ER)*fec[3],R[1]*(1-ER)*fec[4],R[1]*(1-ER)*fec[5])
    #R[4]<-ao*S[4]*exp(-b*S[4])
    #S[5]<-sum(R[4]*(1-ER)*fec[1],R[3]*(1-ER)*fec[2],R[2]*(1-ER)*fec[3],R[1]*(1-ER)*fec[4],R[1]*(1-ER)*fec[5])
    #R[5]<-ao*S[5]*exp(-b*S[5])
    
    for(y in 6:max(yrs)){  
      for(j in seq_along(fec)){      	
        S[y] <- S[y] + R[y-j]*(1-ER)*fec[j]   	
      }        
      #S[y]<-sum(R[y-1]*(1-ER)*fec[1],R[y-2]*(1-ER)*fec[2],R[y-3]*(1-ER)*fec[3],R[y-4]*(1-ER)*fec[4],R[y-5]*(1-ER)*fec[5])
      a[y] <- a[y-1]*exp(rnorm(1,0,siga))
      R[y] <- a[y]*S[y]*exp(-b*S[y])
      Robs[y] <- R[y]*exp(rnorm(1,0,sig))
    }

  outdf <- data.frame(R=R,
	Robs=Robs,
	S=S,
	a=a)

  return(outdf[(70-nobs):70,])
}


