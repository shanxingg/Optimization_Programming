#Vector Version of the Kalman Filter 
#Model that is estimated: x_{t+1}=\beta_1 x1 +\beta_2 x_2 +Tt x_t+error
KF=function(b)
{
 	z=1; # (NxM)
 	beta1=b[1] # effectiveness of instrument 1
 	beta2=b[2] # effectiveness of instrument 1
	Tt=b[3];   # carryover effect
 	h=b[4]^2;  # observation noise
   q=b[5]^2;  # var-cov mat. of the transition noise (MxM)
   icount=1; 
   likeli=0;
    at=init
   pt=mean(y)*4;
    icount=1; 
    while(icount <= raa)
       {#Construction of the likelihood function   
       ptm1=pt
       atm1=at
       attm1=Tt*atm1+beta1*x[icount,1]+beta2*x[icount,2];
       pttm1=Tt*ptm1*t(Tt)+q
       f=z*pttm1*t(z)+h;
       err=y[icount,]-z*attm1-beta4*x[icount,4];  
       kgain=pttm1*t(z)*(solve(f));
       at=attm1+kgain*err;
       pt=pttm1-kgain*z*pttm1
       likeli=likeli+0.5*(log(f)+(err^2)/f); 
       icount=icount+1;
       	}
return(likeli)
}

b=cbind(beta1,beta2,beta3,lambda,h,q,)
Est=optim(b,KF, method = "BFGS", hessian = T, control=list(maxit=300))

hess=Est$hessian #hessian
se= sqrt(diag(solve(hess))) #take sqrt of the diagional of inverse of the Hessian =standard errors
Est$par/se # compute t-values
Est$par # give parameter estimates
