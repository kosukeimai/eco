##fisher transformation of BVN(mu1, mu2, sigma1, sigma2, rho) into
## (mu1, mu2, log(sigma1), log(sigma12), Zp)
fisher<-function(X) {
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

##backward Fisher transformation
fisher.back<-function(Y) {
 p<-length(Y) 
 X<-rep(0,p)
  X[1]<-Y[1]
  X[2]<-Y[2]
  X[3]<-exp(Y[3])
  X[4]<-exp(Y[4])
  if (p==5) 
    X[5]<-(exp(2*Y[5])-1)/(exp(2*Y[5])+1)
  return(X)
}
  
##tranform from  theta into covariance matrix
thetacov<-function(Z) {
  p<-length(Z)
 if (p==5) {
  mat<-matrix(NA,2,2)
  mat[1,1]<-Z[3]
  mat[2,2]<-Z[4]
  mat[1,2]<-Z[5]*sqrt(Z[3]*Z[4])
  mat[2,1]<-mat[1,2]
}
  return(mat)
}

Icom.CAR<-function(theta, suff.stat, n,fisher=TRUE, n.par) {
    Icom<-matrix(NA, n.par, n.par)    
 
    S1<-n*suff.stat[1]
    S2<-n*suff.stat[2]
    S11<-n*suff.stat[3]
    S22<-n*suff.stat[4]
    S12<-n*suff.stat[5]

    u1<-theta[1]
    u2<-theta[2]
    v1<-theta[3]
    v2<-theta[4]
    r<-theta[5]
    
    Icom[1,1]<- -n/((1-r^2)*v1)
    Icom[1,2]<- Icom[2,1] <- n*r/((1-r^2)*sqrt(v1*v2))
    Icom[1,3]<- Icom[3,1] <- 1/((1-r^2)*v1^2)*(-S1+n*u1) -
      r/(2*(1-r^2)*v1^(3/2)*v2^(1/2))*(-S2+n*u2)
    Icom[1,4]<- Icom[4,1] <- -r/(2*(1-r^2)*v1^(1/2)*v2^(3/2))*(-S2+n*u2)
    
    Icom[2,2]<- -n/((1-r^2)*v2)  
    Icom[2,3]<- Icom[3,2] <- -r/(2*(1-r^2)*v1^(3/2)*v2^(1/2))*(-S1+n*u1)
    Icom[2,4]<- Icom[4,2] <- 1/((1-r^2)*v2^2)*(-S2+n*u2) -
      r/(2*(1-r^2)*v1^(1/2)*v2^(3/2))*(-S1+n*u1) 

    Icom[3,3]<- n/(2*v1^2) - 1/((1-r^2)*v1^3)*(S11-2*u1*S1+n*u1^2) +
      3*r/(4*(1-r^2)*v1^(5/2)*v2^(1/2))*(S12-u1*S2-u2*S1+n*u1*u2) 
    Icom[3,4]<- Icom[4,3] <- r/(4*(1-r^2)*v1^(3/2)*v2^(3/2))*(S12-u1*S2-u2*S1+n*u1*u2)

    Icom[4,4]<- n/(2*v2^2) - 1/((1-r^2)*v2^3)*(S22-2*u2*S2+n*u2^2) +
      3*r/(4*(1-r^2)*v1^(1/2)*v2^(5/2))*(S12-u1*S2-u2*S1+n*u1*u2) 

   if (n.par>=5) {
   
    Icom[1,5]<- Icom[5,1] <- -2*r/((1-r^2)^2*v1)*(-S1+n*u1) +
      (1+r^2)/((1-r^2)^2*v1^(1/2)*v2^(1/2))*(-S2+n*u2) 
    Icom[2,5]<- Icom[5,2] <- -2*r/((1-r^2)^2*v2)*(-S2+n*u2) +
      (1+r^2)/((1-r^2)^2*v1^(1/2)*v2^(1/2))*(-S1+n*u1) 
    Icom[3,5]<- Icom[5,3] <- r/((1-r^2)^2*v1^2)*(S11-2*u1*S1+n*u1^2) -
     (1+r^2)/(2*(1-r^2)^2*v1^(3/2)*v2^(1/2))*(S12-u1*S2-u2*S1+n*u1*u2) 
    Icom[4,5]<- Icom[5,4] <- r/((1-r^2)^2*v2^2)*(S22-2*u2*S2+n*u2^2) -
      (1+r^2)/(2*(1-r^2)^2*v1^(1/2)*v2^(3/2))*(S12-u1*S2-u2*S1+n*u1*u2)     
    Icom[5,5]<- n*(1+r^2)/(1-r^2)^2 -
      (1+3*r^2)/((1-r^2)^3*v1)*(S11-2*u1*S1+n*u1^2) -
      (1+3*r^2)/((1-r^2)^3*v2)*(S22-2*u2*S2+n*u2^2) +
      (2*r^3+6*r)/((1-r^2)^3*v1^(1/2)*v2^(1/2))*(S12-u1*S2-u2*S1+n*u1*u2) 
    }

   if (fisher) {
       Icom.fisher<-matrix(NA, n.par,n.par)
       dv1<- -n/(2*v1) + 1/(2*(1-r^2)*v1^2)*(S11-2*u1*S1+n*u1^2) - r/(2*(1-r^2)*v1^(3/2)*v2^(1/2))*(S12-u1*S2-u2*S1+n*u1*u2)
       dv2<- -n/(2*v2) + 1/(2*(1-r^2)*v2^2)*(S22-2*u2*S2+n*u2^2) - r/(2*(1-r^2)*v1^(1/2)*v2^(3/2))*(S12-u1*S2-u2*S1+n*u1*u2)

       Icom.fisher[1:2,1:2]<-Icom[1:2,1:2]
       Icom.fisher[1,3]<- Icom.fisher[3,1] <- Icom[1,3]*v1
       Icom.fisher[1,4]<- Icom.fisher[4,1] <- Icom[1,4]*v2
     
       Icom.fisher[2,3]<- Icom.fisher[3,2] <- Icom[2,3]*v1
       Icom.fisher[2,4]<- Icom.fisher[4,2] <- Icom[2,4]*v2
     
       Icom.fisher[3,3]<- Icom[3,3]*v1^2 + dv1*v1
       Icom.fisher[3,4]<- Icom.fisher[4,3] <- Icom[3,4]*v1*v2
     
       Icom.fisher[4,4]<- Icom[4,4]*v2^2 + dv2*v2

       if (n.par>=5) {
       dr<- n*r/(1-r^2) - r/((1-r^2)^2*v1)*(S11-2*u1*S1+n*u1^2) + (1+r^2)/((1-r^2)^2*v1^(1/2)*v2^(1/2))*(S12-u1*S2-u2*S1+n*u1*u2) - r/((1-r^2)^2*v2)*(S22-2*u2*S2+n*u2^2)
       Icom.fisher[1,5]<- Icom.fisher[5,1] <- Icom[1,5]*(1-r^2)
       Icom.fisher[2,5]<- Icom.fisher[5,2] <- Icom[2,5]*(1-r^2)
       Icom.fisher[3,5]<- Icom.fisher[5,3] <- Icom[3,5]*v1*(1-r^2)
       Icom.fisher[4,5]<- Icom.fisher[5,4] <- Icom[4,5]*v2*(1-r^2)     
       Icom.fisher[5,5]<- Icom[5,5]*(1-r^2)^2 - dr*2*r*(1-r^2)
       }
  }   
     return(list(Icom=-Icom, Icom.fisher=-Icom.fisher))
 }

ecoML <- function(formula, data = parent.frame(), supplement = NULL, 
                  theta.start=c(0,0,1,1,0), epsilon=0.000001,
                  maxit=1000, flag=0, verbose=1) { 

  ## getting X and Y
  mf <- match.call()
  tt <- terms(formula)
  attr(tt, "intercept") <- 0
  if (is.matrix(eval.parent(mf$data)))
    data <- as.data.frame(data)
  X <- model.matrix(tt, data)
  Y <- model.response(model.frame(tt, data=data))

  ## checking the data
  ndim <- 2  
  tmp <- checkdata(X,Y, supplement, ndim)
  bdd <- ecoBD(formula=formula, data=data)
  n.var <- 5
  n <- tmp$n.samp+tmp$survey.samp+tmp$samp.X1+tmp$samp.X0
  inSample.length <- ndim*n

  n.par<-5
  if (flag==2 | flag==6) n.par<-4
  
  ## Fitting the model via EM  
  res <- .C("cEMeco", as.double(tmp$d), as.double(theta.start),
            as.integer(tmp$n.samp),  as.integer(maxit), as.double(epsilon),
            as.integer(tmp$survey.yes), as.integer(tmp$survey.samp), 
            as.double(tmp$survey.data),
            as.integer(tmp$X1type), as.integer(tmp$samp.X1), as.double(tmp$X1.W1),
            as.integer(tmp$X0type), as.integer(tmp$samp.X0), as.double(tmp$X0.W2),
            as.double(bdd$Wmin[,1,1]), as.double(bdd$Wmax[,1,1]),
            as.integer(flag),as.integer(verbose),
            optTheta=c(-1.1,-1.1,-1.1,-1.1,-1.1), pdTheta=double(n.var),
            S=double(n.var+1),inSample=double(inSample.length),DMmatrix=double(n.var*n.var),
            PACKAGE="eco")

  ## SEM step
  if(flag>=4) {
    res <- .C("cEMeco", as.double(tmp$d), as.double(theta.start),
              as.integer(tmp$n.samp),  as.integer(maxit), as.double(epsilon),
              as.integer(tmp$survey.yes), as.integer(tmp$survey.samp), 
              as.double(tmp$survey.data),
              as.integer(tmp$X1type), as.integer(tmp$samp.X1), as.double(tmp$X1.W1),
              as.integer(tmp$X0type), as.integer(tmp$samp.X0), as.double(tmp$X0.W2),
              as.double(bdd$Wmin[,1,1]), as.double(bdd$Wmax[,1,1]),
              as.integer(flag),as.integer(verbose),
              res$pdTheta, pdTheta=double(n.var), S=double(n.var+1),
              inSample=double(inSample.length),DMmatrix=double(n.var*n.var),
              PACKAGE="eco")     
  }
  
  inSample.out <- matrix(rep(NA,inSample.length),ncol=ndim)
  for(i in 1:n)
    for(j in 1:ndim)
      inSample.out[i,j]=res$inSample[(i-1)*2+j]

  DM <- matrix(rep(NA,n.par*n.par),ncol=n.par)
  for(i in 1:n.par)
    for(j in 1:n.par)
      DM[i,j]=res$DMmatrix[(i-1)*n.par+j]

  theta.fisher<-fisher(res$pdTheta)

  ## calculate varcov of theta based on SEM results

  if (flag>=4) {
    ##output Icom
    if (flag==4) 
    infomat<-Icom.CAR(theta=res$pdTheta, suff.stat=res$S, n=n, n.par=n.par)
    if (flag==6) 
    infomat<-Icom.CAR(theta=c(res$pdTheta[1:4],theta.start[5]), suff.stat=c(res$S[1:4],theta.start[5]), n=n, n.par=n.par)


    Icom<-infomat$Icom
    Icom.fisher<-infomat$Icom.fisher
     
    Gamma<-matrix(0,n.par,n.par)
    Gamma[1:2, 1:2]<-Icom.fisher[1:2,1:2]
    Gamma[3:n.par,3:n.par]<-Icom.fisher[3:n.par,3:n.par]
     
    Lamda<-matrix(0,n.par,n.par)
    Lamda[3:n.par, 1:2]<-Icom.fisher[3:n.par, 1:2]
    
    Vcom.fisher<-solve(Icom.fisher)
    dV<-Vcom.fisher%*%DM%*%solve(diag(1,n.par)-DM)

    Vobs.fisher<-Vcom.fisher+dV
    Iobs.fisher<-solve(Vobs.fisher)

    ##transform Iobs.fisher to Iobs via delta method
    ##V(theta)=d(fisher^{-1})V(fisher(theta)))d(fisher^{-1})'
    grad.invfisher<-c(1,1, exp(theta.fisher[3:4]))
    if (flag==4) 
       grad.invfisher<-c(grad.invfisher,4*exp(2*theta.fisher[5])/(exp(2*theta.fisher[5])+1)^2)
    Vobs<-diag(grad.invfisher)%*%Vobs.fisher%*%diag(grad.invfisher)
    Iobs<-solve(Vobs)
  }

   n.col<-n.par
   n.row<-1
   if (flag>=4) n.row<-3
   res.table<-matrix(NA, n.row, n.col)
   res.table[1,]<-res$pdTheta[1:n.par]
   if (n.row>1) {
   res.table[2,]<-sqrt(diag(Vobs))
   res.table[3,]<-Fmis<-diag(Iobs)/diag(Icom)
   }

   cname<-c("mu1", "mu2", "sigma1", "sigma2", "rho")
   rname<-c("EM est", "std err", "frac. missing")
   rownames(res.table)<-rname[1:n.row]
   colnames(res.table)<-cname[1:n.col]

   cat("EM estimates:  ", "\n")
   print(res.table)

   res$DMmatrix<-DM
   res$inSample<-inSample.out  
   res.out<-list(mu=res$pdTheta[1:2], sigma=res$pdTheta[3:4], sigma.rho=theta.fisher[3:4])
   if (flag!=2 & flag!=6) 
   res.out<-c(res.out, rho=res$pdTheta[5], rho.fisher=theta.fisher[5])
   
   if (flag>=4)   
   res.out<-c(res.out, Vobs=Vobs, Fmis=Fmis, Icom=Icom, Iobs=Iobs, suff=res$S[1:n.par], loglike=res$S[n.par+1], res)
 
  return(res.out)
}


#  if (!Icom.yes) {
#    res.out<-list(theta=theta.start)
#  }
#  else if (Icom.yes) {
#    if (Fisher) cat("Icom matrix at Fisher transformation scale:", "\n")
#    else cat("Icom matrix:", "\n")
#    print(Icom)
#    
#    res.out<-list(theta=theta.start, Icom=Icom)
#  }
#  class(res.out) <- "eco"
#  return(res.out)

eco.sem<-function(formula, data = parent.frame(),supplement=NULL, 
                  theta.start=c(0,0,1,1,0), theta.em=NULL, Icom.em=NULL,
                  R.epsilon=0.001, maxit=50, Fisher=TRUE,
                  n.draws = 10, by.draw=10, draw.max=200, printon=TRUE, grid=TRUE) {
  
  # getting X and Y
  mf <- match.call()
  tt <- terms(formula)
  attr(tt, "intercept") <- 0
  if (is.matrix(eval.parent(mf$data)))
    data <- as.data.frame(data)
  X <- model.matrix(tt, data)
  Y <- model.response(model.frame(tt, data=data))

  ndim <- 2
  
  tmp <- checkdata(X,Y, supplement, ndim)
  bdd <- ecoBD(formula=formula, data=data)

  # Fitting the model 
  n.var<-5
  n.Imat<-n.var*(n.var+1)/2

  R.t1<-matrix(0, n.var, n.var)
  R.t2<-matrix(0, n.var, n.var)
  rowdiff<-rep(1, n.var)

  k<-1
  draw <- 1
  Rconverge<-FALSE

  while (!Rconverge && (k<maxit)) {
    res <- .C("cEMeco", as.double(tmp$d), as.double(theta.start),
              as.integer(tmp$n.samp),  as.integer(n.draws), 
              as.integer(tmp$survey.yes), as.integer(tmp$survey.samp), 
          as.double(tmp$survey.data),
              as.integer(tmp$X1type), as.integer(tmp$samp.X1), as.double(tmp$X1.W1),
              as.integer(tmp$X0type), as.integer(tmp$samp.X0), as.double(tmp$X0.W2),
          as.double(bdd$Wmin[,1,1]), as.double(bdd$Wmax[,1,1]),
          as.integer(grid), 
              pdTheta=double(n.var), S=double(n.var),
              PACKAGE="eco")
    
    theta.t<-res$pdTheta
    
    Rconverge<-TRUE
    
    for (i in 1:n.var) {
      rowdiff.temp<-0
      if (rowdiff[i]>R.epsilon) {
        Rconverge<-FALSE
        theta.t.i<-theta.em
        theta.t.i[i]<-theta.t[i]
        temp <- .C("cEMeco", as.double(tmp$d), as.double(theta.t.i),
              as.integer(tmp$n.samp),  as.integer(n.draws), 
              as.integer(tmp$survey.yes), as.integer(tmp$survey.samp), 
          as.double(tmp$survey.data),
              as.integer(tmp$X1type), as.integer(tmp$samp.X1), as.double(tmp$X1.W1),
              as.integer(tmp$X0type), as.integer(tmp$samp.X0), as.double(tmp$X0.W2),
          as.double(bdd$Wmin[,1,1]), as.double(bdd$Wmax[,1,1]),
          as.integer(grid), 
              pdTheta=double(n.var), S=double(n.var),
              PACKAGE="eco")$pdTheta
        
        for (j in 1:n.var) {
          if (Fisher) {
        tempR<-(fisher(temp)[j]-fisher(theta.em)[j])/(fisher(theta.t)[i]-fisher(theta.em)[i])
            if ((tempR!=0) && (tempR<Inf)&& (tempR>-Inf)) R.t2[i,j]<-tempR
          }
          else if (!Fisher) {
            tempR<-(temp[j]-theta.em[j])/(theta.t[i]-theta.em[i])
            if ((tempR!=0) && (tempR<Inf)&& (tempR>-Inf)) R.t2[i,j]<-tempR
          }
          rowdiff.temp<-max(abs(R.t2[i,j]-R.t1[i,j]), rowdiff.temp)    
        }   
        rowdiff[i]<-rowdiff.temp
      }
    }
    
    theta.start<-theta.t
    if (printon) {
      cat("k=", k, "\n")
      print(rowdiff)
      
      R.t1<-R.t2
      cat("\n R")
      print(R.t2)
      cat("\n")
      print(theta.start)
    }
    k<-k+1

    if (draw < draw.max & min(rowdiff)< R.epsilon) {
      Repsilon <- FALSE
      draw <- draw.max
      rowdiff[rowdiff==min(rowdiff)] <- 1
      print("hi")
    }
    if (draw < draw.max)
      draw<-min(draw+by.draw, draw.max)
  }

  cat("\n")
  cat("Estimates based on EM", "\n")
  cat("Mean:", "\n")
  print(theta.em[1:2])
  cat("\n")
  cat("Covarianace Matrix:", "\n")
  print(thetacov(theta.em))


  
  ##missing information decomposition
  
  SECM.yes<-FALSE

  if (SECM.yes) {
  DM.ECM<-R.t2
  
  Gamma<-matrix(0,5,5)
  Gamma[1:2, 1:2]<-Icom.em[1:2,1:2]
  Gamma[3:5,3:5]<-Icom.em[3:5,3:5]
  Lamda<-matrix(0,5,5)
  Lamda[3:5, 1:2]<-Icom.em[3:5, 1:2]
  DM.CM<--Lamda%*%solve(Gamma+t(Lamda))
  
  Vcom<-solve(Icom.em)
  dV<-Vcom%*%(DM.ECM-DM.CM)%*%solve((diag(1,5)-DM.ECM))
}
  if(!SECM.yes) {
  DM<-R.t2
  
  Gamma<-matrix(0,5,5)
  Gamma[1:2, 1:2]<-Icom.em[1:2,1:2]
  Gamma[3:5,3:5]<-Icom.em[3:5,3:5]
  Lamda<-matrix(0,5,5)
  Lamda[3:5, 1:2]<-Icom.em[3:5, 1:2]
    
  Vcom<-solve(Icom.em)
  dV<-Vcom%*%DM%*%solve(diag(1,5)-DM)
}

  Vobs<-Vcom+dV
  
  if (Fisher) {

  cat("Fisher transformtion:", "\n")
  print(fisher(theta.em))
  cat("\n", "the following matrices are at Fisher transformation scale", "\n")

  cat("DM matrix:", "\n")
  print(R.t2)

  cat("Vobs_fisher=:", "\n")
  print(Vobs)
  
  cat("Vcom_fisher=:", "\n")
  print(Vcom)
  
  cat("dV_fisher=:", "\n")
  print(dV)

  v.fish<-c(1,1, theta.em[3], theta.em[4], (1-theta.em[5]^2))
 
  Vcom.orig<-v.fish%*%t(v.fish)*Vcom
  Vobs.orig<-v.fish%*%t(v.fish)*Vobs
  dV.orig<-Vobs.orig-Vcom.orig 

  cat("Vobs_orig=:", "\n")
  print(Vobs.orig)
  
  cat("Vcom_orig=:", "\n")
  print(Vcom.orig)
  
  cat("dV_orig=:", "\n")
  print(dV.orig)

  } 
  else if (!Fisher) {

  cat("DM matrix:", "\n")
  print(R.t2)

  cat("Vobs=:", "\n")
  print(Vobs)
  
  cat("Vcom=:", "\n")
  print(Vcom)
  
  cat("dV=:", "\n")
  print(dV)
}
  res.out<-list(theta=theta.em, Vobs=Vobs, Vcom=Vcom, dV=dV, DM=R.t2)
  class(res.out) <- "eco"
  return(res.out)
}
  
