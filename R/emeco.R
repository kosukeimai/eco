##fisher transformation of BVN(mu1, mu2, sigma1, sigma2, rho) into
## (mu1, mu2, log(sigma1), log(sigma12), Zp)
fisher<-function(X) {
  Y<-rep(0,5)
  Y[1]<-X[1]
  Y[2]<-X[2]
  Y[3]<-log(X[3])
  Y[4]<-log(X[4])
  Y[5]<-0.5*log((1+X[5])/(1-X[5]))
  return(Y)
}  

##backward Fisher transformation
fisher.back<-function(Y) {
  X<-rep(0,5)
  X[1]<-Y[1]
  X[2]<-Y[2]
  X[3]<-exp(Y[3])
  X[4]<-exp(Y[4])
  X[5]<-(exp(2*Y[5])-1)/(exp(2*Y[5]+1))
  return(X)
}
  
##tranform from  theta into covariance matrix
thetacov<-function(Z) {
  mat<-matrix(NA,2,2)
  mat[1,1]<-Z[3]
  mat[2,2]<-Z[4]
  mat[1,2]<-Z[5]*sqrt(Z[3]*Z[4])
  mat[2,1]<-mat[1,2]
  return(mat)
}


eco.em <- function(Y, X, data = parent.frame(),supplement=NULL, 
                   theta.old=c(0,0,1,1,0), convergence=0.000001,
                   iteration.max=1000, Ioc.yes=TRUE, Fisher=TRUE,
                   n.draws = 100000, by.draw=100000,
                   draw.max=10000000, printon=TRUE) { 

  ## checking inputs
  if ((dim(supplement)[2] != 2) && (length(supplement)>0))
    stop("Error: use n by 2 matrix for survey data")
  call <- match.call()
  ff <- as.formula(paste(call$Y, "~ -1 +", call$X))
  if (is.matrix(eval.parent(call$data)))
    data <- as.data.frame(data)
  X <- model.matrix(ff, data)
  Y <- model.response(model.frame(ff, data=data))
  
  ##survey data
  if (length(supplement) == 0) {
    survey.samp <- 0
    survey.data <- 0
    survey.yes<-0
  }
  else {
    survey.samp <- length(supplement[,1])
    survey.data <- as.matrix(supplement)
    survey.yes<-1
  }
  
  ind<-c(1:length(X))
  X1type<-0
  X0type<-0
  samp.X1<-0
  samp.X0<-0
  X1.W1<-0
  X0.W2<-0
  
  ##Xtype x=1
  X1.ind<-ind[along=(X==1)]
  if (length(X[X!=1])<length(X)){
    X1type<-1
    samp.X1<-length(X1.ind)
    X1.W1<-Y[X1.ind]
  }
  
  ##Xtype x=0
  X0.ind<-ind[along=(X==0)]
  if (length(X[X!=0])<length(X)){
    X0type<-1
    samp.X0<-length(X0.ind)
    X0.W2<-Y[X0.ind]
  }
  
  XX.ind<-setdiff(ind, union(X0.ind, X1.ind))
  X.use<-X[XX.ind]
  Y.use<-Y[XX.ind]

  order.old<-order(c(XX.ind, X0.ind, X1.ind))
  
  ## fitting the model
  n.samp <- length(Y.use)	 
  d <- cbind(X.use, Y.use)
  n.var<-5
  n.Imat<-n.var
  cdiff<-1
  em.converge<-FALSE
  i<-1

  
  while ((!em.converge) && (i<iteration.max)) {
    res <- .C("cEMeco", as.double(d), as.double(theta.old),
              as.integer(n.samp),  as.integer(n.draws), 
              as.integer(survey.yes), as.integer(survey.samp), as.double(survey.data),
              as.integer(X1type), as.integer(samp.X1), as.double(X1.W1),
              as.integer(X0type), as.integer(samp.X0), as.double(X0.W2),
              pdTheta=double(n.var), S=double(n.var),
              PACKAGE="eco")
    
    temp<-res$pdTheta
    em.converge<-TRUE
    
    if (printon) {
      cat(i, " ")
      if (Fisher) {cat(fisher(temp), "\n") }
       if (!Fisher) { cat(temp, "\n") }
    }

    if (!Fisher) {
    for (j in 1:5) {   
      print(abs(temp[j]-theta.old[j]))
      if (abs(temp[j]-theta.old[j]) > convergence) em.converge<-FALSE
     }
    }
  
    print(em.converge)

 if (Fisher) {
     for (j in 1:5) {   
      if (abs(fisher(temp)[j]-fisher(theta.old)[j])>convergence) em.converge<-FALSE
     }
   }


    theta.old<-temp
    
    i<-i+1
    if (em.converge & (draw < draw.max)) {
      em.converge <- FALSE
      draw <- draw.max
    }
      
    if (draw<=draw.max) draw<-draw+by.draw
    print(em.converge) 
    print(i)
    print(Fisher)
 }
  
  Ioc<-matrix(NA, n.var, n.var)
  
  if (em.converge && Ioc.yes) {
    ##output Ioc 
    res <- .C("cEMeco", as.double(d), as.double(theta.old),
	      as.integer(n.samp),  as.integer(n.draws), 
              as.integer(survey.yes), as.integer(survey.samp), as.double(survey.data),
   	      as.integer(X1type), as.integer(samp.X1), as.double(X1.W1),
   	      as.integer(X0type), as.integer(samp.X0), as.double(X0.W2),
	      pdTheta=double(n.var), S=double(n.var),
              PACKAGE="eco")
    ##based on pdTheta and S compute Ioc
    
    S1<-res$S[1]
    S2<-res$S[2]
    S11<-res$S[3]
    S22<-res$S[4]
    S12<-res$S[5]
    
    u1<-res$pdTheta[1]
    u2<-res$pdTheta[2]
    v1<-res$pdTheta[3]
    v2<-res$pdTheta[4]
    r<-res$pdTheta[5]
    
    n<-n.samp+survey.samp+samp.X1+samp.X0
    
    Ioc[1,1]<- -n/((1-r^2)*v1)
    Ioc[1,2]<- Ioc[2,1] <- n*r/((1-r^2)*sqrt(v1*v2))
    Ioc[1,3]<- Ioc[3,1] <- 1/((1-r^2)*v1^2)*(-S1+n*u1) -
      r/(2*(1-r^2)*v1^(3/2)*v2^(1/2))*(-S2+n*u2)
    Ioc[1,4]<- Ioc[4,1] <- -r/(2*(1-r^2)*v1^(1/2)*v2^(3/2))*(-S2+n*u2)
    Ioc[1,5]<- Ioc[5,1] <- -2*r/((1-r^2)^2*v1)*(-S1+n*u1) +
      (1+r^2)/((1-r^2)^2*v1^(1/2)*v2^(1/2))*(-S2+n*u2) 
    
    Ioc[2,2]<- -n/((1-r^2)*v2)	
    Ioc[2,3]<- Ioc[3,2] <- -r/(2*(1-r^2)*v1^(3/2)*v2^(1/2))*(-S1+n*u1)
    Ioc[2,4]<- Ioc[4,2] <- 1/((1-r^2)*v2^2)*(-S2+n*u2) -
      r/(2*(1-r^2)*v1^(1/2)*v2^(3/2))*(-S1+n*u1) 
    Ioc[2,5]<- Ioc[5,2] <- -2*r/((1-r^2)^2*v2)*(-S2+n*u2) +
      (1+r^2)/((1-r^2)^2*v1^(1/2)*v2^(1/2))*(-S1+n*u1) 
    
    Ioc[3,3]<- n/(2*v1^2) - 1/((1-r^2)*v1^3)*(S11-2*u1*S1+n*u1^2) +
      3*r/(4*(1-r^2)*v1^(5/2)*v2^(1/2))*(S12-u1*S2-u2*S1+n*u1*u2) 
    Ioc[3,4]<- Ioc[4,3] <- r/(4*(1-r^2)*v1^(3/2)*v2^(3/2))*(S12-u1*S2-u2*S1+n*u1*u2)
    Ioc[3,5]<- Ioc[5,3] <- r/((1-r^2)^2*v1^2)*(S11-2*u1*S1+n*u1^2) -
     (1+r^2)/(2*(1-r^2)^2*v1^(3/2)*v2^(1/2))*(S12-u1*S2-u2*S1+n*u1*u2) 

    Ioc[4,4]<- n/(2*v2^2) - 1/((1-r^2)*v2^3)*(S22-2*u2*S2+n*u2^2) +
      3*r/(4*(1-r^2)*v1^(1/2)*v2^(5/2))*(S12-u1*S2-u2*S1+n*u1*u2) 
    Ioc[4,5]<- Ioc[5,4] <- r/((1-r^2)^2*v2^2)*(S22-2*u2*S2+n*u2^2) -
      (1+r^2)/(2*(1-r^2)^2*v1^(1/2)*v2^(3/2))*(S12-u1*S2-u2*S1+n*u1*u2) 
    
    Ioc[5,5]<- n*(1+r^2)/(1-r^2)^2 -
      (1+3*r^2)/((1-r^2)^3*v1)*(S11-2*u1*S1+n*u1^2) -
      (1+3*r^2)/((1-r^2)^3*v2)*(S22-2*u2*S2+n*u2^2) +
      (2*r^3+6*r)/((1-r^2)^3*v1^(1/2)*v2^(1/2))*(S12-u1*S2-u2*S1+n*u1*u2) 
    

   if (Fisher) {
     dv1<- -n/(2*v1) + 1/(2*(1-r^2)*v1^2)*(S11-2*u1*S1+n*u1^2) - r/(2*(1-r^2)*v1^(3/2)*v2^(1/2))*(S12-u1*S2-u2*S1+n*u1*u2)
     dv2<- -n/(2*v2) + 1/(2*(1-r^2)*v2^2)*(S22-2*u2*S2+n*u2^2) - r/(2*(1-r^2)*v1^(1/2)*v2^(3/2))*(S12-u1*S2-u2*S1+n*u1*u2)
     dr<- n*r/(1-r^2) - r/((1-r^2)^2*v1)*(S11-2*u1*S1+n*u1^2) + (1+r^2)/((1-r^2)^2*v1^(1/2)*v2^(1/2))*(S12-u1*S2-u2*S1+n*u1*u2) - r/((1-r^2)^2*v2)*(S22-2*u2*S2+n*u2^2)

     Ioc[1,3]<- Ioc[3,1] <- Ioc[1,3]*v1
     Ioc[1,4]<- Ioc[4,1] <- Ioc[1,4]*v2
     Ioc[1,5]<- Ioc[5,1] <- Ioc[1,5]*(1-r^2)
     
     Ioc[2,3]<- Ioc[3,2] <- Ioc[2,3]*v1
     Ioc[2,4]<- Ioc[4,2] <- Ioc[2,4]*v2
     Ioc[2,5]<- Ioc[5,2] <- Ioc[2,5]*(1-r^2)
     
     Ioc[3,3]<- Ioc[3,3]*v1^2 + dv1*v1
     Ioc[3,4]<- Ioc[4,3] <- Ioc[3,4]*v1*v2
     Ioc[3,5]<- Ioc[5,3] <- Ioc[3,5]*v1*(1-r^2)
     
     Ioc[4,4]<- Ioc[4,4]*v2^2 + dv2*v2
     Ioc[4,5]<- Ioc[5,4] <- Ioc[4,5]*v2*(1-r^2)
     
     Ioc[5,5]<- Ioc[5,5]*(1-r^2)^2 - dr*2*r*(1-r^2)
   }
    Ioc <- -Ioc
  }
  
  cat("\n")
  cat("Estimates based on EM", "\n")
  cat("Mean:", "\n")
  print(theta.old[1:2])
  cat("\n")
  cat("Covarianace Matrix:", "\n")
  print(thetacov(theta.old))
  
  if (Fisher) {
    cat("Fisher transformtion:", "\n")
    print(fisher(theta.old))
  }

  if (!Ioc.yes) {
    res.out<-list(theta=theta.old)
  }
  else if (Ioc.yes) {
    if (Fisher) cat("Ioc matrix at Fisher transformation scale:", "\n")
    else cat("Ioc matrix:", "\n")
    print(Ioc)
    
    res.out<-list(theta=theta.old, Ioc=Ioc)
  }
  class(res.out) <- "eco"
  return(res.out)
  
}


eco.sem<-function(Y, X, data = parent.frame(),supplement=NULL, 
                  theta.old=c(0,0,1,1,0), theta.em=NULL, Ioc.em=NULL,
                  R.convergence=0.001, iteration.max=50, Fisher=TRUE,
                  n.draws = 10, by.draw=10, draw.max=200, printon=TRUE) {
  
  ## checking inputs
  if ((dim(supplement)[2] != 2) && (length(supplement)>0))
    stop("Error: use n by 2 matrix for survey data")
  
  call <- match.call()

  ff <- as.formula(paste(call$Y, "~ -1 +", call$X))
  if (is.matrix(eval.parent(call$data)))
    data <- as.data.frame(data)
  X <- model.matrix(ff, data)
  Y <- model.response(model.frame(ff, data=data))

  ##survey data
  if (length(supplement) == 0) {
    survey.samp <- 0
    survey.data <- 0
    survey.yes<-0
  }
  else {
    survey.samp <- length(supplement[,1])
    survey.data <- as.matrix(supplement)
    survey.yes<-1
  }
  
  ind<-c(1:length(X))
  X1type<-0
  X0type<-0
  samp.X1<-0
  samp.X0<-0
  X1.W1<-0
  X0.W2<-0
  
  ##Xtype x=1
  X1.ind<-ind[along=(X==1)]
  if (length(X[X!=1])<length(X)){
    X1type<-1
    samp.X1<-length(X1.ind)
    X1.W1<-Y[X1.ind]
  }
  
  ##Xtype x=0
  X0.ind<-ind[along=(X==0)]
  if (length(X[X!=0])<length(X)){
    X0type<-1
    samp.X0<-length(X0.ind)
    X0.W2<-Y[X0.ind]
  }
  
  XX.ind<-setdiff(ind, union(X0.ind, X1.ind))
  X.use<-X[XX.ind]
  Y.use<-Y[XX.ind]

  order.old<-order(c(XX.ind, X0.ind, X1.ind))
  
  ## fitting the model
  n.samp <- length(Y.use)	 
  d <- cbind(X.use, Y.use)
  n.var<-5
  n.Imat<-n.var*(n.var+1)/2

  R.t1<-matrix(0, n.var, n.var)
  R.t2<-matrix(0, n.var, n.var)
  rowdiff<-rep(1, n.var)

  k<-1

  Rconverge<-FALSE

  while (!Rconverge && (k<iteration.max)) {
    res <- .C("cEMeco", as.double(d), as.double(theta.old),
              as.integer(n.samp),  as.integer(n.draws), 
              as.integer(survey.yes), as.integer(survey.samp), as.double(survey.data),
              as.integer(X1type), as.integer(samp.X1), as.double(X1.W1),
              as.integer(X0type), as.integer(samp.X0), as.double(X0.W2),
              pdTheta=double(n.var), S=double(n.var),
              PACKAGE="eco")
    
    theta.t<-res$pdTheta
    
    Rconverge<-TRUE
    
    for (i in 1:n.var) {
      rowdiff.temp<-0
      if (rowdiff[i]>R.convergence) {
        Rconverge<-FALSE
        theta.t.i<-theta.em
        theta.t.i[i]<-theta.t[i]
        
        temp<-.C("cEMeco", as.double(d), as.double(theta.t.i),
                 as.integer(n.samp),  as.integer(n.draws), 
                 as.integer(survey.yes), as.integer(survey.samp), as.double(survey.data),
                 as.integer(X1type), as.integer(samp.X1), as.double(X1.W1),
                 as.integer(X0type), as.integer(samp.X0), as.double(X0.W2),
                 pdTheta=double(n.var), S=double(n.var),
                 PACKAGE="eco")$pdTheta
        
        for (j in 1:n.var) {
          if (Fisher) {
            R.t2[i,j]<-(fisher(temp)[j]-fisher(theta.em)[j])/(fisher(theta.t)[i]-fisher(theta.em)[i])
          }
          else if (!Fisher) {
            R.t2[i,j]<-(temp[j]-theta.em[j])/(theta.t[i]-theta.em[i])
          }
          rowdiff.temp<-max(abs(R.t2[i,j]-R.t1[i,j]), rowdiff.temp)    
        }   
        rowdiff[i]<-rowdiff.temp
      }
      else if (draw < draw.max) {
        Rconvergence <- FALSE
        draw <- draw.max
        rowdiff[i] <- 1
      }
    }
    
    theta.old<-theta.t
    if (printon) {
      cat("k=", k, "\n")
      print(rowdiff)
      
      R.t1<-R.t2
      cat("\n R")
      print(R.t2)
      cat("\n")
      print(theta.old)
    }
    k<-k+1
    
    if (draw<=draw.max) draw<-draw+by.draw
  }
  
  
  
  cat("\n")
  cat("Estimates based on EM", "\n")
  cat("Mean:", "\n")
  print(theta.em[1:2])
  cat("\n")
  cat("Covarianace Matrix:", "\n")
  print(thetacov(theta.em))
  


  
  ##missing information decomposition
  
  DM.ECM<-R.t2
  
  Gamma<-matrix(0,5,5)
  Gamma[1:2, 1:2]<-Ioc.em[1:2,1:2]
  Gamma[3:5,3:5]<-Ioc.em[3:5,3:5]
  Lamda<-matrix(0,5,5)
  Lamda[3:5, 1:2]<-Ioc.em[3:5, 1:2]
  DM.CM<--Lamda%*%solve(Gamma+t(Lamda))
  
  Vcom<-solve(Ioc.em)
  dV<-Vcom%*%(DM.ECM-DM.CM)%*%solve((diag(1,5)-DM.ECM))

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
  res.out<-list(theta=theta.em, Vobs=Vobs, Vcom=Vcom, dV=dV, DM.ECM=R.t2)
  class(res.out) <- "eco"
  return(res.out)
}
  
