
## This is a replication code for "Robust Dynamic Estimation, by Olivier Rubel and Prasad A. Naik (2016), Marketing Science, forthcoming.
# Use the code along with the data for Alberta provided in the file "Alberta.csv". 
# Data contains 156 weeks as rows and time, collection, radio, and calls as columns
# The code and data are offered with NO express or impplied warranty and is not to be used for commercial purposes
# If you use this data set or the code for your research, then pls cite the article.


rm(list=ls(all=TRUE)) #clear data
library(numDeriv)


############Kalman Filter#########
data<-read.csv("/replication/alberta.csv", sep=",",dec=".",header=T)
collection<-data[,2]
radio<-data[,3]
ncc<-data[,4]
ra=dim(data)[1]-52 # last 52 weeks hold out sample
scaling=5		

x=cbind(sqrt(radio),sqrt(ncc))
y=cbind(collection/scaling)
init=y[1]

#STANDARD FILTER
KF=function(b)
{
 	z=1; # (NxM)
 	beta1=b[1] # effectiveness of radio
 	beta2=b[2] # effectiveness of calls
	Tt=b[3];   # carryover effect
 	h=b[4]^2;  # observation noise
    q=b[4]^2;  # var-cov mat. of the transition noise (MxM)
   #gamma=gcount;
   icount=1; 
   likeli=0;
    at=init
   pt=5000;
    icount=1; 
    while(icount <= ra)
       {#Construction of the likelihood function   
       f=z*pt*t(z)+h;
       ptilde=solve((solve(pt)+t(z)*solve(h)*z));
       kgain=ptilde*t(z)*solve(h);
       drift=beta1*x[icount,1]+beta2*x[icount,2];
       yhat=z*at+drift;
       err=y[icount]-yhat;                   
       likeli=likeli+0.5*(log(f)+(err^2)/f); 
       at=Tt*at+drift+Tt*kgain*err;
       pt=Tt*ptilde*t(Tt)+q;   	
       icount=icount+1;
       	}
return(likeli)
}

# starting values
beta1=0.1
beta2=0.1
lambda=0.5
qq=50;
a0=init
b=cbind(beta1,beta2,lambda,qq)

Est=optim(b,KF, method = "BFGS", hessian = T)
se= sqrt(diag(solve(Est$hessian)))
Tvals=Est$par/se
parKF=Est$par[1:3]
TvalsKF=Tvals[1:3]
LogL=-Est$value

########### KALMAN FILTER WITH WHITE CORRECTION ########

KFRobust=function(b)
{
 	z=1; 
 	beta1=b[1] #effectiveness of radio
 	beta2=b[2] #effectiveness of calls
	Tt=b[3];   # carryover effect
 	h=b[4]^2;  # observation noise
    q=b[4]^2;  # transition noise
   LL<-mat.or.vec((ra),0)
   #gamma=gcount;
   icount=1; 
   likeli=0;
   at=init;		 
   pt=5000;
   icount=1; 
    while(icount <= ra)
       {#Construction of the likelihood function   
       f=z*pt*t(z)+h;
       ptilde=solve((solve(pt)+t(z)*solve(h)*z));
       kgain=ptilde*t(z)*solve(h);
       drift=beta1*x[icount,1]+beta2*x[icount,2];
       yhat=z*at+drift;
       err=y[icount]-yhat;                   
       LL[icount]<-0.5*(log(f)+(err^2)/f); 
       at=Tt*at+drift+Tt*kgain*err;
       pt=Tt*ptilde*t(Tt)+q;   	
       icount=icount+1;
       	}
return(LL)
}
paraKF=Est$par
KFRobust(paraKF)
jac=jacobian(KFRobust,paraKF)
G=t(jac)%*%jac
Rse= sqrt(diag(solve(Est$hessian)%*%G%*%solve(Est$hessian)))
robusTvals=Est$par/Rse
parKFWC=Est$par[1:3]
TvalsKFWC=robusTvals[1:3]

######################################  H INFINITY FILTER  Retained Model #############################################
KF=function(b)
{
 	z=1; # (NxM)
 	beta1=b[1] # effectiveness of radio
 	beta2=b[2] # effectiveness of calls
	Tt=b[3];   # carryover effect
 	h=b[4]^2;  # observation noise 
    q=b[4]^2;  # transition noise
   #gamma=gcount;
   icount=1; 
   likeli=0;
   at=init;		 
   pt=5000;
   icount=1; 
    while(icount <= ra)
       {#Construction of the likelihood function   
       f=z*pt*t(z)+h;
       ptilde=solve((solve(pt)+t(z)*solve(h)*z));
       kgain=ptilde*t(z)*solve(h);
       drift=beta1*x[icount,1]+beta2*x[icount,2];
       yhat=z*at+drift;
       err=y[icount]-yhat;                   
       likeli=likeli+0.5*(log(f)+(err^2)/f); 
       at=Tt*at+drift+Tt*kgain*err;
       pt=Tt*ptilde*t(Tt)+q;   	
       icount=icount+1;
       	}
return(likeli)
}

# starting values
beta1=0.1
beta2=0.1
lambda=0.5
qq=50;
a0=init
b=cbind(beta1,beta2,lambda,qq)

Est=optim(b, method = "BFGS", hessian = T)
params1=Est$par
results0=cbind(Est$par,Est$value)
Robust=function(b)
{
 	z=1; # (NxM)
 	beta1=b[1] #effectiveness of radio
 	beta2=b[2] # effectiveness of calls
	Tt=b[3];   # carryover effect
 	h=b[4]^2;  # observation noise
    q=b[4]^2;  # transition noise
   icount=1; 
   likeli=0;
   at= init;		 
   pt=5000;
   gamma=4000
   icount=1; 
    while(icount <= ra)
       {#Construction of the likelihood function   
       f=z*pt*t(z)+h;
       ptilde=solve((solve(pt)-1/gamma+t(z)*solve(h)*z));
       kgain=ptilde*t(z)*solve(h);
       drift=beta1*x[icount,1]+beta2*x[icount,2];
       yhat=z*at+drift;
       err=y[icount]-yhat;                   
       likeli=likeli+0.5*(log(f)+(err^2)/f); 
       at=Tt*at+drift+Tt*kgain*err;
       pt=Tt*ptilde*t(Tt)+q;   	
       icount=icount+1;
       	}
return(likeli)
}

# starting values
beta1=0.1
beta2=0.1
lambda=0.5
qq=50;
a0=init
b=cbind(params1)

Est=optim(b, ra=ra,gamma=7000,x=x,y=y,init=init,Robust, gr = NULL, method = "BFGS",lower = -Inf, upper = Inf, control = list(maxit=1000), hessian = T)
parR=Est$par[1:3]
se= sqrt(diag(solve(Est$hessian)))
Tvals=Est$par/se
Tvals[1:3]
LogLR=-Est$value

######################################  H INFINITY FILTER  Retained Model with White Correction#############################################
#Robust FILTER
Robust=function(b)
{
 	z=1; # (NxM)
 	beta1=b[1]
 	beta2=b[2]
	Tt=b[3]; 
 	h=b[4]^2;
    q=b[4]^2; 
   gamma=4000;
   icount=1; 
   likeli=0;
   at=init;		 
   pt=5000;
   icount=1; 
    while(icount <= ra)
       {#Construction of the likelihood function   
       f=z*pt*t(z)+h;
       ptilde=solve((solve(pt)-1/gamma+t(z)*solve(h)*z));
       kgain=ptilde*t(z)*solve(h);
       drift=beta1*x[icount,1]+beta2*x[icount,2];
       yhat=z*at+drift;
       err=y[icount]-yhat;                   
       likeli=likeli+0.5*(log(f)+(err^2)/f); 
       at=Tt*at+drift+Tt*kgain*err;
       pt=Tt*ptilde*t(Tt)+q;   	
       icount=icount+1;
       	}
return(likeli)
}

# starting values
beta1=0.1
beta2=0.1
lambda=0.5
qq=50;
a0=init
b=cbind(beta1,beta2,lambda,qq)

EstRobust=optim(b, Robust, method = "BFGS", hessian = T)
paraRobust=EstRobust$par
se= sqrt(diag(solve(EstRobust$hessian)))
Tvals=paraRobust/se
parR=paraRobust[1:3]
TvalsR=Tvals[1:3]

HinfRobust=function(b)
{
 	z=1;
 	beta1=b[1]
 	beta2=b[2]
	Tt=b[3]; 
 	h=b[4]^2;
    q=b[4]^2; 
   LL<-mat.or.vec((ra),0)
   gamma=4000;
   icount=1; 
   likeli=0;
   at=init;		 
   pt=5000;
   icount=1; 
    while(icount <= ra)
       {#Construction of the likelihood function   
       f=z*pt*t(z)+h;
       ptilde=solve((solve(pt)-1/gamma+t(z)*solve(h)*z));
       kgain=ptilde*t(z)*solve(h);
       drift=beta1*x[icount,1]+beta2*x[icount,2];
       yhat=z*at+drift;
       err=y[icount]-yhat;                   
       LL[icount]<-0.5*(log(f)+(err^2)/f); 
       at=Tt*at+drift+Tt*kgain*err;
       pt=Tt*ptilde*t(Tt)+q;   	
       icount=icount+1;
       	}
return(LL)
}
HinfRobust(paraRobust)
jac=jacobian(HinfRobust,paraRobust)
G=t(jac)%*%jac
Rse= sqrt(diag(solve(EstRobust$hessian)%*%G%*%solve(EstRobust$hessian)))
robusTvals=paraRobust/Rse
parRWC=paraRobust[1:3]
TvalsRWC=robusTvals[1:3]

########################## OUTPUT ###############
sink("/Replication/output.txt", append=FALSE, split=FALSE)
print ("Today's date and time"); format(Sys.time(), "%a %b %d %X %Y")
cat("\n")
print("Parameters from Kalman Filter (Table 1)" ); parKF 
cat("\n")
print("T values from Kalman Filter (parentehtical entries)"); TvalsKF 
cat("\n")
print("Log Likelihood from Kalman Filter in Table 2"); LogL 
cat("\n")
print("Parameters from Kalman Filter with White Correction (Table 1)"); parKFWC  
cat("\n")
print ("T values from Kalman Filter with White Correction (parentehtical entries)"); TvalsKFWC
cat("\n")
print("Parameters from Robust Filter (Table 1)"); parR 
cat("\n")
print("T values from Robust Filter (parentehtical entries)"); TvalsR 
cat("\n")
print("Log Likelihood from Robust Filter in Table 2"); LogLR 
cat("\n")
print("Parameters from Robust Filter with White Correction (Table 1)"); parRWC
cat("\n")
print ("T values from Robust Filter with White Correction (parentehtical entries)"); TvalsRWC 
sink()
