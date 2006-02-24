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


###
### main function
###
ecoML <- function(formula, data = parent.frame(), N=NULL, supplement = NULL, 
                  theta.start = c(0,0,1,1,0), fix.rho = TRUE,
                  context = FALSE, sem = TRUE, epsilon=10^(-10),
                  maxit = 1000, loglik = TRUE, verbose= TRUE) { 

  ## Aaron, we need separate options for verbose and loglik:
  ## verbose = TRUE (print out each iteration on screen
  ## loglike = TRUE (calculate the loglik at each iteration)
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
            as.integer(flag),as.integer(verbose),as.integer(loglik),
            optTheta=c(-1.1,-1.1,-1.1,-1.1,-1.1), pdTheta=double(n.var),
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
  for(i in 1:n)
    for(j in 1:ndim)
      W[i,j]=res$inSample[(i-1)*2+j]

  ## SEM step
  DM <- matrix(rep(NA,n.par*n.par),ncol=n.par)
  if(flag>=4) {
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

  if (flag>=4) {
    ##output Icom
    if (flag==4) 
      infomat<-Icom.CAR(theta=theta.em, suff.stat=res$S, n=n, n.par=n.par)
    if (flag==6) 
      infomat<-Icom.CAR(theta=c(theta.em[1:4],theta.start[5]), suff.stat=res$S, n=n, n.par=n.par)

    Icom<-infomat$Icom
    Icom.fisher<-infomat$Icom.fisher
    
    Vcom.fisher <- solve(Icom.fisher)
    dV <- Vcom.fisher%*%DM%*%solve(diag(1,n.par)-DM)
    Vobs.fisher <- Vcom.fisher+dV
    Iobs.fisher <- solve(Vobs.fisher)

    ##transform Iobs.fisher to Iobs via delta method
    ##V(theta)=d(fisher^{-1})V(bvn.trans(theta)))d(fisher^{-1})'
    grad.invfisher <- c(1,1, exp(theta.fisher[3:4]))
    if (flag==4) 
       grad.invfisher <- c(grad.invfisher,4*exp(2*theta.fisher[5])/(exp(2*theta.fisher[5])+1)^2)
    Vobs<-diag(grad.invfisher)%*%Vobs.fisher%*%diag(grad.invfisher)
    Iobs<-solve(Vobs)
  }
  ## make it symmetric
  Vobs.sym <- 0.5*(Vobs+t(Vobs))
  
  ##if (max(abs(Vobs-Vobs.sym)/Vobs)>0.05) 
  ## warnings("the covariance matrix estimated based on SEM is not symmetric enough. \n")

  ## Ying, put this stuff in summary.ecoML()
  n.col<-n.par
  n.row<-1
  if (flag>=4) n.row<-3
  res.table<-matrix(NA, n.row, n.col)
  res.table[1,]<-theta.em[1:n.par]
  if (n.row>1) {
    res.table[2,]<-sqrt(diag(Vobs))
    res.table[3,]<-Fmis<-1-diag(Iobs)/diag(Icom)
  }
  cname<-c("mu1", "mu2", "sigma1", "sigma2", "rho")
  rname<-c("EM est.", "std. err.", "frac. missing")
  rownames(res.table)<-rname[1:n.row]
  colnames(res.table)<-cname[1:n.col]
  cat("\n")
  cat("EM estimates:  ", "\n")
  cat("\n")
  print(res.table)
  
  ## output
  res.out<-list(call=mf, Y=Y, X=X,N=N, 
		fix.rho=fix.rho, context=context, sem=sem, epsilon=epsilon,
                mu = theta.em[1:2], sigma = theta.em[3:4],
                sigma.log = theta.fisher[3:4], suff = res$S[1:n.var],
                loglike = res$S[n.var+1], iters.em = iters.em, 
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
    res.out$Fmis<-Fmis
    res.out$Vobs.original<-Vobs
    res.out$Vobs<-Vobs.sym
  }

  class(res.out) <- "ecoML"
  return(res.out)
}
