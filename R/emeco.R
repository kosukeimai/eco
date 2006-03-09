## transformation of BVN(mu1, mu2, sigma1, sigma2, rho) into
## (mu1, mu2, log(sigma1), log(sigma12), Zp)
bvn.trans <-function(X) {
  p<-length(X) 
  Y<-rep(0,p)
  Y[1]<-X[1]
  Y[2]<-X[2]
  Y[3]<-log(X[3])
  Y[4]<-log(X[4])
  if (p==5) 
   Y[5]<-0.5*log((1+X[5])/(1-X[5]))
  return(Y)
}  

vec<-function(mat) {
  v<-as.vector(mat, mode="any")
  v
}

tr<-function(mat) {
  trace<-sum(diag(mat))
  trace
}
  

## I_{com} 
## the gradient function for multivariate normal function
#du.theta and dSig.theta are the first derivative of mu and Sigma 
#with respect to theta
#du.theta[n.u, n.theta]
#dSig.theta[n.u, n.u, n.theta]

d1st.mvn<-function(mu,Sigma, Fisher=FALSE, fix.rho=FALSE) {
   d<-length(mu)
   p<-d+d+d*(d-1)/2
   u1<-mu[1]
   u2<-mu[2]
   s1<-Sigma[1,1]
   s2<-Sigma[2,2]
   r12<-Sigma[1,2]/sqrt(s1*s2)
print(r12)
   if (d==3) {
   u3<-mu[3]
   s3<-Sigma[3,3]
   r13<-Sigma[1,3]/sqrt(s1*s3)
   r23<-Sigma[2,3]/sqrt(s2*s3)
   }
 
   if (fix.rho) p<-p-1
 
    du.theta<-matrix(0,d,p)
    for (j in 1:d) 
      du.theta[j,j]<-1

    dSig.theta<-array(0,c(d,d,p))
    if (!Fisher) {
      for (i in 1:d) 
        dSig.theta[i,i,d+i]<-1

      dSig.theta[1,2,d+1]<-dSig.theta[2,1,d+1]<-1/2*s1^(-1/2)*s2^(1/2)*r12
      dSig.theta[1,2,d+2]<-dSig.theta[2,1,d+2]<-1/2*s2^(-1/2)*s1^(1/2)*r12
      if (d==3) {
        dSig.theta[1,3,d+1]<-dSig.theta[3,1,d+1]<-1/2*s1^(-1/2)*s3^(1/2)*r12
        dSig.theta[1,3,d+3]<-dSig.theta[3,1,d+3]<-1/2*s3^(-1/2)*s1^(1/2)*r12
        dSig.theta[2,3,d+2]<-dSig.theta[3,2,d+2]<-1/2*s2^(-1/2)*s3^(1/2)*r12
        dSig.theta[2,3,d+3]<-dSig.theta[3,2,d+3]<-1/2*s3^(-1/2)*s2^(1/2)*r12
      }
    
    if (!fix.rho) {
        dSig.theta[1,2,2*d+1]<-dSig.theta[2,1,2*d+1]<-sqrt(s1*s2)
        if (d==3) {
          dSig.theta[1,3,2*d+2]<-dSig.theta[3,1,2*d+2]<-sqrt(s1*s3)
          dSig.theta[2,3,2*d+3]<-dSig.theta[3,2,2*d+3]<-sqrt(s2*s3)
        }
     }
     if (fix.rho) {
        if (d==3) {
          dSig.theta[1,3,2*d+1]<-dSig.theta[3,1,2*d+1]<-sqrt(s1*s3)
          dSig.theta[2,3,2*d+2]<-dSig.theta[3,2,2*d+2]<-sqrt(s2*s3)
        }
     }
   } 

    if (Fisher) {
      for (i in 1:d) 
        dSig.theta[i,i,d+i]<-Sigma[i,i]

      dSig.theta[1,2,d+1]<-dSig.theta[2,1,d+1]<-1/2*s1^(1/2)*s2^(1/2)*r12
      dSig.theta[1,2,d+2]<-dSig.theta[2,1,d+2]<-1/2*s2^(1/2)*s1^(1/2)*r12
      if (d==3) {
        dSig.theta[1,3,d+1]<-dSig.theta[3,1,d+1]<-1/2*s1^(1/2)*s3^(1/2)*r13
        dSig.theta[2,3,d+2]<-dSig.theta[3,2,d+2]<-1/2*s2^(1/2)*s3^(1/2)*r23
        dSig.theta[1,3,d+3]<-dSig.theta[3,1,d+3]<-1/2*s3^(1/2)*s1^(1/2)*r13 
        dSig.theta[2,3,d+3]<-dSig.theta[3,2,d+3]<-1/2*s3^(1/2)*s2^(1/2)*r23
      }
    
    if (!fix.rho) {
        dSig.theta[1,2,2*d+1]<-dSig.theta[2,1,2*d+1]<-sqrt(s1*s2)*(1-r12^2)
        if (d==3) {
          dSig.theta[1,3,2*d+2]<-dSig.theta[3,1,2*d+2]<-sqrt(s1*s3)*(1-r13^2)
          dSig.theta[2,3,2*d+3]<-dSig.theta[3,2,2*d+3]<-sqrt(s2*s3)*(1-r23^2)
        }
      }
      if (fix.rho) {
         if (d==3) {
          dSig.theta[1,3,2*d+1]<-dSig.theta[3,1,2*d+1]<-sqrt(s1*s3)*(1-r13^2)
          dSig.theta[2,3,2*d+2]<-dSig.theta[3,2,2*d+2]<-sqrt(s2*s3)*(1-r23^2)
       }
      }
  }
   return(list(du.theta=du.theta, dSig.theta=dSig.theta))
}

d2nd.mvn<-function(mu,Sigma, Fisher=FALSE, fix.rho=FALSE) {

   d<-length(mu)
   p<-d+d+d*(d-1)/2
   u1<-mu[1]
   u2<-mu[2]
   s1<-Sigma[1,1]
   s2<-Sigma[2,2]
   r12<-Sigma[1,2]/sqrt(s1*s2)
   if (d==3) {
   u3<-mu[3]
   s3<-Sigma[3,3]
   r13<-Sigma[1,3]/sqrt(s1*s3)
   r23<-Sigma[2,3]/sqrt(s2*s3)
   }

   if (fix.rho) p<-p-1

   ddu.theta<-array(0,c(d,p,p))

   ddSig.theta<-array(0,c(d,d,p,p))
   if (!Fisher) {
     ddSig.theta[1,2,d+1,d+1]<-ddSig.theta[2,1,d+1,d+1]<- -1/4*s1^(-3/2)*s2^(1/2)*r12
     ddSig.theta[1,2,d+1,d+2]<-ddSig.theta[2,1,d+1,d+2]<- 1/4*s1^(-1/2)*s2^(-1/2)*r12
     ddSig.theta[1,2,d+2,d+2]<-ddSig.theta[2,1,d+2,d+2]<- -1/4*s1^(1/2)*s2^(-3/2)*r12
 
     if (d==3) {
     ddSig.theta[1,3,d+1,d+1]<-ddSig.theta[3,1,d+1,d+1]<- -1/4*s1^(-3/2)*s3^(1/2)*r13
     ddSig.theta[1,3,d+1,d+3]<-ddSig.theta[3,1,d+1,d+3]<- 1/4*s1^(-1/2)*s3^(-1/2)*r13

     ddSig.theta[2,3,d+2,d+2]<-ddSig.theta[3,2,d+2,d+2]<- -1/4*s2^(-3/2)*s3^(1/2)*r23
     ddSig.theta[2,3,d+2,d+3]<-ddSig.theta[3,2,d+2,d+3]<- 1/4*s2^(-1/2)*s3^(-1/2)*r23

     ddSig.theta[1,3,d+3,d+3]<-ddSig.theta[3,1,d+3,d+3]<- -1/4*s1^(1/2)*s3^(-3/2)*r13
     ddSig.theta[2,3,d+3,d+3]<-ddSig.theta[3,2,d+3,d+3]<- -1/4*s2^(1/2)*s3^(-3/2)*r23

     }

     if (!fix.rho) {
       ddSig.theta[1,2,d+1,2*d+1]<-ddSig.theta[2,1,d+1,2*d+1]<- 1/2*s1^(-1/2)*s2^(1/2)
       ddSig.theta[1,2,d+2,2*d+1]<-ddSig.theta[2,1,d+2,2*d+1]<- 1/2*s1^(1/2)*s2^(-1/2)

       if (d==3) {
         ddSig.theta[1,3,d+1,2*d+2]<-ddSig.theta[3,1,d+1,2*d+2]<- 1/2*s1^(-1/2)*s3^(1/2)
         ddSig.theta[2,3,d+2,2*d+3]<-ddSig.theta[3,2,d+2,2*d+3]<- 1/2*s2^(-1/2)*s3^(1/2)
         ddSig.theta[1,3,d+3,2*d+2]<-ddSig.theta[3,1,d+3,2*d+2]<- 1/2*s1^(1/2)*s3^(-1/2)
         ddSig.theta[2,3,d+3,2*d+3]<-ddSig.theta[3,2,d+3,2*d+3]<- 1/2*s2^(1/2)*s3^(-1/2)
     }
   }
    if (fix.rho) {
       if (d==3) {

       ddSig.theta[1,2,d+1,2*d+1]<-ddSig.theta[2,1,d+1,2*d+1]<- 1/2*s1^(-1/2)*s3^(1/2)
       ddSig.theta[2,3,d+2,2*d+2]<-ddSig.theta[3,2,d+2,2*d+2]<- 1/2*s2^(-1/2)*s3^(1/2)
         ddSig.theta[1,3,d+3,2*d+1]<-ddSig.theta[3,1,d+3,2*d+1]<- 1/2*s1^(1/2)*s3^(-1/2)
         ddSig.theta[2,3,d+3,2*d+2]<-ddSig.theta[3,2,d+3,2*d+2]<- 1/2*s2^(1/2)*s3^(-1/2)
     }
   }
  }
    if (Fisher) {

        ddSig.theta[1,1,d+1,d+1]<- s1 
        ddSig.theta[1,2,d+1,d+1]<-ddSig.theta[2,1,d+1,d+1]<- 1/4*s1^(1/2)*s2^(1/2)*r12

        ddSig.theta[2,2,d+2,d+2]<- s2 
        ddSig.theta[1,2,d+2,d+2]<-ddSig.theta[2,1,d+2,d+2]<- 1/4*s1^(1/2)*s2^(1/2)*r12

        ddSig.theta[1,2,d+1,d+2]<-ddSig.theta[2,1,d+1,d+2]<- 1/4*s1^(1/2)*s2^(1/2)*r12

        if (d==3) {
           ddSig.theta[1,3,d+1,d+1]<-ddSig.theta[3,1,d+1,d+1]<- 1/4*s1^(1/2)*s3^(1/2)*r13
           ddSig.theta[1,3,d+2,d+2]<-ddSig.theta[3,1,d+2,d+2]<- 1/4*s2^(1/2)*s3^(1/2)*r23

           ddSig.theta[1,3,d+1,d+3]<-ddSig.theta[3,1,d+1,d+3]<- 1/4*s1^(1/2)*s3^(1/2)*r13
           ddSig.theta[2,3,d+2,d+3]<-ddSig.theta[3,2,d+2,d+3]<- 1/4*s2^(1/2)*s3^(1/2)*r23

           ddSig.theta[1,3,d+3,d+3]<-ddSig.theta[3,1,d+3,d+3]<- 1/4*s1^(1/2)*s3^(1/2)*r13
           ddSig.theta[2,3,d+3,d+3]<-ddSig.theta[3,2,d+3,d+3]<- 1/4*s2^(1/2)*s3^(1/2)*r23
        }

      if (!fix.rho) {
           ddSig.theta[1,2,d+1,2*d+1]<-ddSig.theta[2,1,d+1,2*d+1]<- 1/2*s1^(1/2)*s2^(1/2)*(1-r12^2)
           ddSig.theta[1,2,d+2,2*d+1]<-ddSig.theta[2,1,d+2,2*d+1]<- 1/2*s1^(1/2)*s2^(1/2)*(1-r12^2)

           ddSig.theta[1,2,2*d+1,2*d+1]<-ddSig.theta[2,1,2*d+1,2*d+1]<- -2*s1^(1/2)*s2^(1/2)*r12*(1-r12^2)

           if (d==3) {
            ddSig.theta[1,3,d+1,2*d+2]<-ddSig.theta[3,1,d+1,2*d+2]<- 1/2*s1^(1/2)*s3^(1/2)*(1-r13^2)
            ddSig.theta[2,3,d+2,2*d+3]<-ddSig.theta[3,2,d+2,2*d+3]<- 1/2*s2^(1/2)*s3^(1/2)*(1-r23^2)
            ddSig.theta[1,3,d+3,2*d+2]<-ddSig.theta[3,1,d+3,2*d+2]<- 1/2*s1^(1/2)*s3^(1/2)*(1-r13^2)
            ddSig.theta[2,3,d+3,2*d+3]<-ddSig.theta[3,2,d+3,2*d+3]<- 1/2*s2^(1/2)*s3^(1/2)*(1-r23^2)

            ddSig.theta[1,3,2*d+2,2*d+2]<-ddSig.theta[3,1,2*d+2,2*d+2]<- -2*s1^(1/2)*s3^(1/2)*r13*(1-r13^2)
            ddSig.theta[2,3,2*d+3,2*d+3]<-ddSig.theta[3,2,2*d+3,2*d+3]<- -2*s2^(1/2)*s3^(1/2)*r23*(1-r23^2)

          }
      }  
   
       if (fix.rho) {
          if (d==3) {
           ddSig.theta[1,3,d+1,2*d+1]<-ddSig.theta[3,1,d+1,2*d+1]<- 1/2*s1^(1/2)*s3^(1/2)*(1-r13^2)
           ddSig.theta[2,3,d+2,2*d+2]<-ddSig.theta[3,2,d+2,2*d+2]<- 1/2*s2^(1/2)*s3^(1/2)*(1-r23^2)
           ddSig.theta[1,3,d+3,2*d+1]<-ddSig.theta[3,1,d+3,2*d+1]<- 1/2*s1^(1/2)*s3^(1/2)*(1-r13^2)
           ddSig.theta[2,3,d+3,2*d+2]<-ddSig.theta[3,2,d+3,2*d+2]<- 1/2*s2^(1/2)*s3^(1/2)*(1-r23^2)

           ddSig.theta[1,3,2*d+1,2*d+1]<-ddSig.theta[3,1,2*d+1,2*d+1]<- -2*s1^(1/2)*s3^(1/2)*r13*(1-r13^2)
           ddSig.theta[2,3,2*d+2,2*d+2]<-ddSig.theta[3,2,2*d+2,2*d+2]<- -2*s2^(1/2)*s3^(1/2)*r23*(1-r23^2)
       
          }
        }

    }
      for (i in 1:(p-1)) 
        for (j in (i+1):p) {
         ddSig.theta[,,j,i]<-ddSig.theta[,,i,j]
         ddu.theta[,j,i]<-ddu.theta[,i,j]
      }

      return(list(ddu.theta=ddu.theta, ddSig.theta=ddSig.theta))
}

suff<-function(mu, suff.stat,n) {

   d<-length(mu)
   p<-d+d+d*(d-1)/2 
   u1<-mu[1]
   u2<-mu[2]
   if (d==3)  u3<-mu[3]

   S1<-n*suff.stat[1]
   S2<-n*suff.stat[2]
   S11<-n*suff.stat[d+1]
   S22<-n*suff.stat[d+2]
   S12<-n*suff.stat[2*d+1]
   if (d==3) {
      S3<-n*suff.stat[d]
      S33<-n*suff.stat[2*d]
      S13<-n*suff.stat[2*d+2]
      S23<-n*suff.stat[2*d+3]
   }

   Vv<-rep(0,d)
   Vv[1]<-S1-n*u1
   Vv[2]<-S2-n*u2
   if (d==3) Vv[3]<-S3-n*u3

   Ss<-matrix(0,d,d)
   Ss[1,1]<-S11-2*S1*u1+n*u1^2 
   Ss[2,2]<-S22-2*S2*u2+n*u2^2
   Ss[1,2]<-Ss[2,1]<-S12-S1*u2-S2*u1+n*u1*u2
   if (d==3) {
      Ss[3,3]<-S33-2*S3*u3+n*u3^2
      Ss[1,3]<-Ss[3,1]<-S13-S1*u3-S3*u1+n*u1*u3
      Ss[2,3]<-Ss[3,2]<-S23-S3*u2-S2*u3+n*u2*u3
  }
  return(list(Ss=Ss, Vv=Vv))
}

Dcom.mvn<-function(mu, Sigma, suff.stat,n, fix.rho=FALSE, Fisher=FALSE) {
  d<-dim(Sigma)[1]
  p<-d*2+0.5*d*(d-1)

  if (fix.rho) { 
    p<-p-1 
    print(p)
  }

  Dcom<-rep(0,p)
  invSigma<-solve(Sigma)

  temp<-suff(mu, suff.stat, n)
  Ss<-temp$Ss
  Vv<-temp$Vv

  temp<-d1st.mvn(mu=mu, Sigma=Sigma, fix.rho=fix.rho, Fisher=FALSE)
  print(temp)
  du.theta<-temp$du.theta
  dSig.theta<-temp$dSig.theta
print(p)

  for (i in 1:p)  
   Dcom[i]<- -n/2*t(vec(invSigma))%*%vec(dSig.theta[,,i])+ 0.5*tr(invSigma%*%dSig.theta[,,i]%*%invSigma%*%Ss)+ t(du.theta[,i])%*%invSigma%*%Vv

   Dcom
}

#du.theta and dSig.theta are the second derivative of mu and Sigma 
#with respect to theta
#ddu.theta[n.u, n.theta, n.theta]
#ddSig.theta[n.u, n.u, n.theta, n.theta]

Icom.mvn<-function(mu, Sigma, suff.stat,n, fix.rho=FALSE, Fisher=FALSE) {
   d<-dim(Sigma)[1]
  p<-d*2+1/2*d*(d-1)

  if (fix.rho) 
    { 
      p<-p-1 
    }

   Icom<-matrix(0,p,p)

   invSigma<-solve(Sigma)

   temp<-suff(mu, suff.stat, n)
   Ss<-temp$Ss
   Vv<-temp$Vv

   temp<-d1st.mvn(mu, Sigma, fix.rho, Fisher=FALSE)
   du.theta<-temp$du.theta
   dSig.theta<-temp$dSig.theta

   temp<-d2nd.mvn(mu, Sigma, fix.rho, Fisher=FALSE)
   ddu.theta<-temp$ddu.theta
   ddSig.theta<-temp$ddSig.theta

   for (i in 1:p) {
     dinvSig.theta.i<- -invSigma%*%dSig.theta[,,i]%*%invSigma
     for (j in 1:i) {
      dinvSig.theta.j<- -invSigma%*%dSig.theta[,,j]%*%invSigma
      ddinvSig.theta.ij<- -dinvSig.theta.j%*%dSig.theta[,,i]%*%invSigma -invSigma%*%ddSig.theta[,,i,j]%*%invSigma-invSigma%*%dSig.theta[,,i]%*%dinvSig.theta.j
 
       a1<- -n/2*(t(vec(dinvSig.theta.j))%*%vec(dSig.theta[,,i]) + t(vec(invSigma))%*%vec(ddSig.theta[,,i,j]))
                    
       a2<- t(du.theta[,j])%*%dinvSig.theta.i%*%Vv - 0.5*tr(ddinvSig.theta.ij%*%Ss)

     a3<- t(ddu.theta[,i,j])%*%invSigma%*%Vv + t(du.theta[,i])%*%dinvSig.theta.j%*%Vv - n*t(du.theta[,i])%*%invSigma%*%du.theta[,j]

       Icom[i,j]<-a1+a2+a3
    
       if (i!=j) Icom[j,i]<-Icom[i,j]
     }
   }
   -Icom
}

## compute I_{com} for CAR
Icom.CAR <- function(theta, suff.stat, n, fisher=TRUE, n.par) {
  Icom <- matrix(NA, n.par, n.par)    
  
  S1<-n*suff.stat[1]
  S2<-n*suff.stat[2]
  S11<-n*suff.stat[3]
  S22<-n*suff.stat[4]
  S12<-n*suff.stat[5]
  if (n.par==4) 
    S12<-n*suff.stat[5]
  
  u1<-theta[1]
  u2<-theta[2]
  v1<-theta[3]
  v2<-theta[4]
  r<-theta[5]
    
  Icom[1,1]<- -n/((1-r^2)*v1)
  Icom[1,2]<- Icom[2,1] <- n*r/((1-r^2)*sqrt(v1*v2))
  Icom[1,3]<- Icom[3,1] <- 1/((1-r^2)*v1^2)*(-S1+n*u1) - r/(2*(1-r^2)*v1^(3/2)*v2^(1/2))*(-S2+n*u2)
  Icom[1,4]<- Icom[4,1] <- -r/(2*(1-r^2)*v1^(1/2)*v2^(3/2))*(-S2+n*u2)
  
  Icom[2,2]<- -n/((1-r^2)*v2)  
  Icom[2,3]<- Icom[3,2] <- -r/(2*(1-r^2)*v1^(3/2)*v2^(1/2))*(-S1+n*u1)
  Icom[2,4]<- Icom[4,2] <- 1/((1-r^2)*v2^2)*(-S2+n*u2) - r/(2*(1-r^2)*v1^(1/2)*v2^(3/2))*(-S1+n*u1) 
  
  Icom[3,3]<- n/(2*v1^2) - 1/((1-r^2)*v1^3)*(S11-2*u1*S1+n*u1^2) + 3*r/(4*(1-r^2)*v1^(5/2)*v2^(1/2))*(S12-u1*S2-u2*S1+n*u1*u2) 
  print( n/(2*v1^2))
  print(- 1/((1-r^2)*v1^3)*(S11-2*u1*S1+n*u1^2))
  print(3*r/(4*(1-r^2)*v1^(5/2)*v2^(1/2))*(S12-u1*S2-u2*S1+n*u1*u2))

 Icom[3,4]<- Icom[4,3] <- r/(4*(1-r^2)*v1^(3/2)*v2^(3/2))*(S12-u1*S2-u2*S1+n*u1*u2)
  
  Icom[4,4]<- n/(2*v2^2) - 1/((1-r^2)*v2^3)*(S22-2*u2*S2+n*u2^2) + 3*r/(4*(1-r^2)*v1^(1/2)*v2^(5/2))*(S12-u1*S2-u2*S1+n*u1*u2) 
  dr<-0
  if (n.par>=5) {
    Icom[1,5]<- Icom[5,1] <- -2*r/((1-r^2)^2*v1)*(-S1+n*u1) + (1+r^2)/((1-r^2)^2*v1^(1/2)*v2^(1/2))*(-S2+n*u2) 
    Icom[2,5]<- Icom[5,2] <- -2*r/((1-r^2)^2*v2)*(-S2+n*u2) + (1+r^2)/((1-r^2)^2*v1^(1/2)*v2^(1/2))*(-S1+n*u1) 
    Icom[3,5]<- Icom[5,3] <- r/((1-r^2)^2*v1^2)*(S11-2*u1*S1+n*u1^2) - (1+r^2)/(2*(1-r^2)^2*v1^(3/2)*v2^(1/2))*(S12-u1*S2-u2*S1+n*u1*u2) 
    Icom[4,5]<- Icom[5,4] <- r/((1-r^2)^2*v2^2)*(S22-2*u2*S2+n*u2^2) - (1+r^2)/(2*(1-r^2)^2*v1^(1/2)*v2^(3/2))*(S12-u1*S2-u2*S1+n*u1*u2)     
    Icom[5,5]<- n*(1+r^2)/(1-r^2)^2 - (1+3*r^2)/((1-r^2)^3*v1)*(S11-2*u1*S1+n*u1^2) -(1+3*r^2)/((1-r^2)^3*v2)*(S22-2*u2*S2+n*u2^2) +(2*r^3+6*r)/((1-r^2)^3*v1^(1/2)*v2^(1/2))*(S12-u1*S2-u2*S1+n*u1*u2) 
  }

    du1<- ((S1-n*u1)/v1-r*(S2-n*u2)/v2)/(1-r^2)
    du2<- ((S2-n*u2)/v2-r*(S1-n*u1)/v1)/(1-r^2)
    dv1<- -n/(2*v1) + 1/(2*(1-r^2)*v1^2)*(S11-2*u1*S1+n*u1^2) - r/(2*(1-r^2)*v1^(3/2)*v2^(1/2))*(S12-u1*S2-u2*S1+n*u1*u2)
    dv2<- -n/(2*v2) + 1/(2*(1-r^2)*v2^2)*(S22-2*u2*S2+n*u2^2) - r/(2*(1-r^2)*v1^(1/2)*v2^(3/2))*(S12-u1*S2-u2*S1+n*u1*u2)
   if (n.par>=5) {
      dr<- n*r/(1-r^2) - r/((1-r^2)^2*v1)*(S11-2*u1*S1+n*u1^2) + (1+r^2)/((1-r^2)^2*v1^(1/2)*v2^(1/2))*(S12-u1*S2-u2*S1+n*u1*u2) - r/((1-r^2)^2*v2)*(S22-2*u2*S2+n*u2^2)
  }
  return(list(Icom=-Icom, Dvec=c(du1, du2, dv1, dv2, dr)))
}
  

Icom.transform<-function(Imat, Dvec, theta, transformation="Fisher")
  {
    n.par<-dim(Imat)[1]
    if (transformation=="Fisher") {
     T1<-c(1,1,theta[3], theta[4])

    T2<-matrix(0, n.par^2, n.par)
    T2[(2*n.par+3), 3]<-theta[3]     
    T2[(3*n.par+4), 4]<-theta[4]     

    if (n.par==5) {
      T1<-c(T1, (1-(theta[5]^2)))
      T2[(4*n.par+5),5]<- -2*theta[5]*(1-theta[5]^2)
     }

    }

   
    T1<-diag(T1)
    Icom.tran<-matrix(NA, n.par, n.par)
    Icom.tran<-T1%*%Imat%*%t(T1)
    
    temp1<-matrix(0,n.par,n.par)
    for (i in 1:n.par)
      for (j in 1:n.par) 
       temp1[i,j]<- Dvec%*%T2[((i-1)*n.par+(1:n.par)),j] 

      Icom.tran<-Icom.tran+temp1     
    return(-Icom.tran)
}



###
### main function
###
ecoML <- function(formula, data = parent.frame(), N=NULL, supplement = NULL, 
                  theta.start = c(0,0,1,1,0), fix.rho = TRUE,
                  context = FALSE, sem = TRUE, epsilon=10^(-10),
                  maxit = 1000, loglik = TRUE, verbose= TRUE) { 

  if (verbose)
    verbose <- 1
  else
    verbose <- 0

  if (loglik)
    loglik <- 1
  else
    loglik <- 0
  
  ## translating into flag
  if (context)
    if (sem)
      if (fix.rho)
        flag <- 7
      else
        flag <- 5
    else
      if (fix.rho)
        flag <- 3
      else
        flag <- 1
  else
    if (sem)
      if (fix.rho)
        flag <- 6
      else
        flag <- 4
    else
      if (fix.rho)
        flag <- 2
      else
        flag <- 0
  
  ## getting X and Y
  mf <- match.call()
  tt <- terms(formula)
  attr(tt, "intercept") <- 0
  if (is.matrix(eval.parent(mf$data)))
    data <- as.data.frame(data)
  X <- model.matrix(tt, data)
  Y <- model.response(model.frame(tt, data=data))

  ## checking the data
  if (context) {
    ndim <- 3
    n.var<-7
  }
  else {
    ndim <- 2
    n.var <- 5
  }
  tmp <- checkdata(X, Y, supplement, ndim)
  bdd <- ecoBD(formula=formula, data=data)
  n <- tmp$n.samp+tmp$survey.samp+tmp$samp.X1+tmp$samp.X0
  inSample.length <- ndim*tmp$n.samp

  n.par<-n.var
  if (fix.rho) n.par<-n.par-1
  
  ## Fitting the model via EM  
  res <- .C("cEMeco", as.double(tmp$d), as.double(theta.start),
            as.integer(tmp$n.samp),  as.integer(maxit), as.double(epsilon),
            as.integer(tmp$survey.yes), as.integer(tmp$survey.samp), 
            as.double(tmp$survey.data),
            as.integer(tmp$X1type), as.integer(tmp$samp.X1), as.double(tmp$X1.W1),
            as.integer(tmp$X0type), as.integer(tmp$samp.X0), as.double(tmp$X0.W2),
            as.double(bdd$Wmin[,1,1]), as.double(bdd$Wmax[,1,1]),
            as.integer(flag),as.integer(verbose),as.integer(loglik),
            optTheta=rep(-1.1,n.var), pdTheta=double(n.var),
            S=double(n.var+1),inSample=double(inSample.length),DMmatrix=double(n.par*n.par),
            itersUsed=as.integer(0),history=double((maxit+1)*(n.var+1)),
            PACKAGE="eco")

  ##record results from EM
  theta.em<-res$pdTheta
  theta.fisher<-bvn.trans(theta.em)
  iters.em<-res$itersUsed
  mu.em <- matrix(rep(NA,iters.em*ndim),ncol=ndim)
  sigma.log.em <- matrix(rep(NA,iters.em*ndim),ncol=ndim)
  rho.fisher.em <- as.double(rep(NA,iters.em))
  loglike.em <- as.double(rep(NA,iters.em))

  for(i in 1:iters.em) {
    mu.em[i,1]=res$history[(i-1)*(n.var+1)+1]
    mu.em[i,2]=res$history[(i-1)*(n.var+1)+2]
    sigma.log.em[i,1]=res$history[(i-1)*(n.var+1)+3]
    sigma.log.em[i,2]=res$history[(i-1)*(n.var+1)+4]
    rho.fisher.em[i]=res$history[(i-1)*(n.var+1)+5]
    loglike.em[i]=res$history[(i-1)*(n.var+1)+6]
  }

  ## In sample prediction of W
  W <- matrix(rep(NA,inSample.length),ncol=ndim)
  for(i in 1:tmp$n.samp)
    for(j in 1:ndim)
      W[i,j]=res$inSample[(i-1)*2+j]

  ## SEM step
  DM <- matrix(rep(NA,n.par*n.par),ncol=n.par)
  if (sem) {
    res <- .C("cEMeco", as.double(tmp$d), as.double(theta.start),
              as.integer(tmp$n.samp),  as.integer(maxit), as.double(epsilon),
              as.integer(tmp$survey.yes), as.integer(tmp$survey.samp), 
              as.double(tmp$survey.data),
              as.integer(tmp$X1type), as.integer(tmp$samp.X1), as.double(tmp$X1.W1),
              as.integer(tmp$X0type), as.integer(tmp$samp.X0), as.double(tmp$X0.W2),
              as.double(bdd$Wmin[,1,1]), as.double(bdd$Wmax[,1,1]),
              as.integer(flag),as.integer(verbose),as.integer(loglik),
              res$pdTheta, pdTheta=double(n.var), S=double(n.var+1),
              inSample=double(inSample.length),DMmatrix=double(n.par*n.par),
              itersUsed=as.integer(0),history=double((maxit+1)*(n.var+1)),
              PACKAGE="eco")     
  
    iters.sem<-res$itersUsed
    for(i in 1:n.par)
      for(j in 1:n.par)
        DM[i,j]=res$DMmatrix[(i-1)*n.par+j]
  }

  if (sem) {
    ##output Icom

##call Icom.mvn and Dcom.mvn

    if (!fix.rho) 
      infomat<-Icom.CAR(theta=theta.em, suff.stat=res$S, n=n, n.par=n.par)
    if (fix.rho) 
      infomat<-Icom.CAR(theta=c(theta.em[1:4],theta.start[5]), suff.stat=res$S, n=n, n.par=n.par)

    Icom<-infomat$Icom
    Sig<-matrix(0,ndim, ndim)
    Sig[1,1]<-theta.em[3]
    Sig[2,2]<-theta.em[4]
    if (!fix.rho) Sig[1,2]<-Sig[2,1]<-theta.em[5]*sqrt(Sig[1,1]*Sig[2,2])
    if (fix.rho) Sig[1,2]<-Sig[2,1]<-theta.start[5]*sqrt(Sig[1,1]*Sig[2,2])


    Icom.new<-Icom.mvn(mu=theta.em[1:2], Sigma=Sig, fix.rho=fix.rho, suff.stat=res$S, n=n)
    Dvec.new<-Dcom.mvn(mu=theta.em[1:2], Sigma=Sig, fix.rho=fix.rho, suff.stat=res$S, n=n)

    if (!fix.rho) {
    Icom.fisher<-Icom.transform(Imat=-Icom, Dvec=infomat$Dvec, theta=theta.em, transformation="Fisher")
    Icom.new.fisher<-Icom.transform(Imat=-Icom.new, Dvec=Dvec.new, theta=theta.em, transformation="Fisher")   
    }
    
   if (fix.rho) {
    Icom.fisher<-Icom.transform(Imat=-Icom, Dvec=infomat$Dvec[1:4], theta=c(theta.em[1:4], theta.start[5]), transformation="Fisher")
    Icom.new.fisher<-Icom.transform(Imat=-Icom.new, Dvec=Dvec.new, theta=c(theta.em[1:4], theta.start[5]), transformation="Fisher")   
  }

    Vcom.fisher <- solve(Icom.fisher)
    dV <- Vcom.fisher%*%DM%*%solve(diag(1,n.par)-DM)
    Vobs.fisher <- Vcom.fisher+dV
    Iobs.fisher <- solve(Vobs.fisher)

    ##transform Iobs.fisher to Iobs via delta method
    ##V(theta)=d(fisher^{-1})V(bvn.trans(theta)))d(fisher^{-1})'
    grad.invfisher <- c(1,1, exp(theta.fisher[3:4]))
    if (! fix.rho) 
       grad.invfisher <- c(grad.invfisher,4*exp(2*theta.fisher[5])/(exp(2*theta.fisher[5])+1)^2)
    Vobs<-diag(grad.invfisher)%*%Vobs.fisher%*%diag(grad.invfisher)
    Iobs<-solve(Vobs)
    ## obtain a symmetric Cov matrix
    Vobs.sym <- 0.5*(Vobs+t(Vobs))
 
    ##if (max(abs(Vobs-Vobs.sym)/Vobs)>0.05) 
    ## warnings("the covariance matrix estimated based on SEM is not symmetric enough. \n")
  }

  ## output
  res.out<-list(call = mf, Y = Y, X = X, N = N, 
                fix.rho = fix.rho, context = context, sem=sem, epsilon=epsilon,
                mu = theta.em[1:2], sigma = theta.em[3:4],
                sigma.log = theta.fisher[3:4], suff = res$S[1:n.var],
                loglik = res$S[n.var+1], iters.em = iters.em, 
                iters.sem = iters.sem, mu.em = mu.em,
                sigma.log.em = sigma.log.em,
                rho.fisher.em = rho.fisher.em, loglike.em = loglike.em,
                W = W, DM = DM)
  if (fix.rho) res.out$rho0<-theta.start[5]
  if (!fix.rho) {
    res.out$rho <- theta.em[5]
    res.out$rho.fisher <- theta.fisher[5]
  }
  
  if (sem) {
    res.out$Icom<-Icom
    res.out$Iobs<-Iobs
    res.out$Fmis<-1-diag(Iobs)/diag(Icom)
    res.out$Vobs.original<-Vobs
    res.out$Vobs<-Vobs.sym
    res.out$Icom.trans<-Icom.fisher
    res.out$Iobs.trans<-Iobs.fisher
    res.out$Fmis.trans<-1-diag(Iobs.fisher)/diag(Icom.fisher)
    res.out$Dvec<-infomat$Dvec
    res.out$Icom.new<-Icom.new
    res.out$Icom.new.trans<-Icom.new.fisher
  }

  class(res.out) <- "ecoML"
  return(res.out)
}
