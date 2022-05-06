#==============================================================
#Simple simulation evaluation of Ricker model with time-varying parameters
#Catarina Wor
#May 2022
#==============================================================




#ao<-3
#b<-1/10000
#ER<-0.40
#fec= c(0,.2,.4,.8,1)
#sig <- .5
#siga <-.3
#
#ps<-0:35000
#pr<-ao*ps*exp(-b*ps)
#sr <- simulateSRrandom()
#plot(sr$S,log(sr$S/sr$R))
#
#par(mfrow=c(2,1))
#plot(sr$S,sr$R, xlim=c(0,35000), ylim=c(0,35000))
#lines(ps,pr,col="red")
#plot(log(sr$a),type="b")


#TODO check this function
#


simulateSRrandom <- function(ao=3, b=1/10000, ER=0.4, fec=c(0,0,0,1,0), sig=.5, siga=.2, nobs=40 ){

    yrs <- 1:100
    S <- NULL
    R <- NULL
    Robs <- NULL
    a <- NULL
    
    a[1:5] <- ao
    Seq <- ao/b
    S[1] <- Seq
    R[1] <- Seq*exp(ao-b*Seq)
  
    for(y in 2:length(fec)){
      S[y] <- 0    
      for(j in seq_along(fec)){
        if((y-j) > 0){
          S[y] <- S[y] + R[y-j]*(1-ER)*fec[j]
        }else{
          S[y] <- S[y] + R[1]*(1-ER)*fec[j]
        }
      }
      R[y]<-S[y]*exp(a[y]-b*S[y])        
    }
  
    
    for(y in 6:max(yrs)){  
      S[y]<-0
      for(j in seq_along(fec)){      	
        S[y] <- S[y] + R[y-j]*(1-ER)*fec[j]   	
      } 
      #truncated lognrmal error
      #a[y] <- a[y-1] + (qnorm(runif(1, 0.1, 0.9),0,siga))        
      a[y] <- a[y-1] + rnorm(1,0,siga) #-(.5*siga^2)

      #R[y] <- S[y]*exp(a[y]-b*S[y]) *exp(qnorm(runif(1, 0.1, 0.95),0,sig))
      R[y] <- a[y]*S[y]*exp(-b*S[y]) *exp(rnorm(1,0,sig))#-(.5*sig^2)
      if(!is.na(R[y]) & R[y]<5){
      	R[y] <- NA
      }
      
    }


  outdf <- data.frame(R=R,
	S=S,
	a=a,
	logR_S=log(R/S))

  return(outdf[(100-nobs):100,])
}


simulateSRtrend <- function(ao=2.5, b=1/10000, ER=0.4, fec=c(0,0,0,1,1), sig=.5, 
	siga=.3, nobs=40, trend="decline",lowsca=.5,hisca=2, ampsc=.5){

    yrs <- 1:100
    S <- NULL
    R <- NULL
    a <- NULL
   
    a[1:5] <- ao
    Seq <- log(ao)/b
    S[1] <- Seq
    R[1] <- ao*Seq*exp(-b*Seq)
    
    #trends in a
    #start trend 10 yrs before data collection started
    trlen <- (length(yrs)-(nobs+9)):length(yrs)  

    if(trend=="decline"){
    
      amin <- ao*lowsca
      atrend <- seq(ao,amin, length=nobs+10)
      a[6:(length(yrs)-(nobs+10))] <- ao
      a[trlen] <- atrend
    
    }else if(trend=="increase"){
    
      amax<-ao*hisca
      atrend<-seq(ao,amax, length=nobs+10)
    
    }else if(trend=="sine"){
       
      atrend<- ao + (ao*ampsc)* sin(2*pi/((nobs+10)/2)*trlen)

    }else if(trend=="regime"){
        
        atrend<- ao *rep(rep(c(.5,2),each=10),length.out=length(trlen))
        
    }
    a[6:(length(yrs)-(nobs+10))]<-ao
    a[trlen] <- atrend
  
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
  
    
    for(y in 6:max(yrs)){  
      for(j in seq_along(fec)){      	
        S[y] <- S[y] + R[y-j]*(1-ER)*fec[j]   	
      }        
      
      R[y] <- a[y]*S[y]*exp(-b*S[y])*exp(rnorm(1,0,sig))
      if(R[y]<2){
      	R[y] <- NA
      }
    }

  outdf <- data.frame(R=R,
	S=S,
	a=a)

  return(outdf[(100-nobs):100,])
}



