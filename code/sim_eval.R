#==============================================================
#Simple simulation evaluation of Ricker model with time-varying parameters
#Catarina Wor
#May 2022
#==============================================================




ao<-2.5
b<-1/10000
ER<-0.40

sp<-1:35000
rp<-ao*sp*exp(-b*sp)

ages= 2:5
fec= c(0,.2,.4,.8,1)
yrs <- 1:70

sig <- .5
siga <-.3

S <- NULL
R <-NULL
a<-NULL
Ra<-matrix(NA, ncol=length(fec),nrow=yrs)

a[1:5]<-ao
Seq<-log(ao)/b
S[1]<-Seq
R[1]<-ao*Seq*exp(-b*Seq)
S[2]<-sum(R[1]*(1-ER)*fec[1],R[1]*(1-ER)*fec[2],R[1]*(1-ER)*fec[3],R[1]*(1-ER)*fec[4],R[1]*(1-ER)*fec[5])
R[2]<-ao*S[2]*exp(-b*S[2])
S[3]<-sum(R[2]*(1-ER)*fec[1],R[1]*(1-ER)*fec[2],R[1]*(1-ER)*fec[3],R[1]*(1-ER)*fec[4],R[1]*(1-ER)*fec[5])
R[3]<-ao*S[3]*exp(-b*S[3])
S[4]<-sum(R[3]*(1-ER)*fec[1],R[2]*(1-ER)*fec[2],R[1]*(1-ER)*fec[3],R[1]*(1-ER)*fec[4],R[1]*(1-ER)*fec[5])
R[4]<-ao*S[4]*exp(-b*S[4])
S[5]<-sum(R[4]*(1-ER)*fec[1],R[3]*(1-ER)*fec[2],R[2]*(1-ER)*fec[3],R[1]*(1-ER)*fec[4],R[1]*(1-ER)*fec[5])
R[5]<-ao*S[5]*exp(-b*S[5])

plot(R)

for(y in 6:max(yrs)){
  
  
  S[y]<-sum(R[y-1]*(1-ER)*fec[1],R[y-2]*(1-ER)*fec[2],R[y-3]*(1-ER)*fec[3],R[y-4]*(1-ER)*fec[4],R[y-5]*(1-ER)*fec[5])
  a[y]<-a[y-1]*exp(rnorm(1,0,siga))
  R[y]<-a[y]*S[y]*exp(-b*S[y])
  Robs[y]<-R[y]*

}




 plot(S,R,type="b", ylim=c(0,12500), xlim=c(0,35000))
 points(S[70],R[70], pch=3, col="purple", cex=3)
lines(sp,rp,col="red")
abline(a=0,b=1, col="blue")
?abline
?lines