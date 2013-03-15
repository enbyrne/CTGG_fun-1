#################################################
# script to calculate power of a case control study
# Naomi Wray June 16 
#see Yang et al EJEpi 2009
#################################################
#power curve for 
#p allele frequency 
#R relative risk 
#K disease with prevalence 
#N01 total sample size
#v proportion of cases
#################################################
power_curve<-function(q2,N,T){
NCP<-N*q2/(1-q2)
power_curve<-pnorm(sqrt(NCP)+T)
}

#############################################
# main program
#############################################
A<-5e-08 # alpha level sig threshold
T<-qnorm(A/2,0,1)
N<-2400  # sample size for QT trait

x<-c(0,0.03)
y<-c(0,1)
plot(x,y,xlim=c(0,0.03),ylim=c(0,1),xlab="% variance explained",ylab="power",main=paste("power for sample size of",N),type="n",cex.lab=2,cex.main=2)
curve(power_curve(x,N,T),from=0.001,to=0.03,col=1,lty=1,lwd=2,add=TRUE)