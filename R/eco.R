eco <- function(Y, X, data = parent.frame(), 
		n.draws = 5000, burnin = 0, thin = 5, verbose = FALSE,
		nonpar =  TRUE, nu0 = 4, tau0 = 2, mu0 = c(0,0),
                S0 = diag(10,2), supplement=NULL,
	        alpha = NULL, a0 = 1, b0 = 0.1,
                predict = FALSE, parameter = FALSE){ 

  ## checking inputs
  if (burnin >= n.draws)
    stop("Error: n.draws should be larger than burnin")
  
  if ((dim(supplement)[2] != 2) && (length(supplement)>0))
    stop("Error: use n by 2 matrix for survey data")
  
  call <- match.call()

  ff <- as.formula(paste(call$Y, "~ -1 +", call$X))
  if (is.matrix(eval.parent(call$data)))
    data <- as.data.frame(data)
  X <- model.matrix(ff, data)
  Y <- model.response(model.frame(ff, data=data))

  ##alpha
  alpha.update <- FALSE
  if (length(alpha)==0) 
    { alpha.update <- TRUE 
      alpha<-0
    }
  
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

  if (nonpar){	# nonparametric model
    n.a <- floor((n.draws-burnin)/thin)
    n.par <- n.a * (n.samp+samp.X1+samp.X0) 
    n.w <- n.a * (n.samp+samp.X1+samp.X0) 
    unit.a <- 1
    unit.par <- (n.samp+samp.X1+samp.X0) 
    unit.w <- (n.samp+samp.X1+samp.X0) 	
    res <- .C("cDPeco", as.double(d), as.integer(n.samp),
	      as.integer(n.draws), as.integer(burnin), as.integer(thin),
	      as.integer(verbose),
	      as.integer(nu0), as.double(tau0), as.double(mu0), as.double(S0),
              as.double(alpha), as.integer(alpha.update), 
	      as.double(a0), as.double(b0),
	      as.integer(survey.yes), as.integer(survey.samp), as.double(survey.data),
   	      as.integer(X1type), as.integer(samp.X1), as.double(X1.W1),
   	      as.integer(X0type), as.integer(samp.X0), as.double(X0.W2),
              as.integer(predict), as.integer(parameter),
              pdSMu0=double(n.par), pdSMu1=double(n.par),
              pdSSig00=double(n.par), pdSSig01=double(n.par),
              pdSSig11=double(n.par), 
	      pdSW1=double(n.w), pdSW2=double(n.w), 
	      pdSWt1=double(n.w), pdSWt2=double(n.w),
	      pdSa=double(n.a), pdSn=integer(n.a), PACKAGE="eco")
    if (parameter) {
      mu1.post <- matrix(res$pdSMu0, n.a, unit.par, byrow=TRUE)[,order.old]
      mu2.post <- matrix(res$pdSMu1, n.a, unit.par, byrow=TRUE)[,order.old]
      Sigma11.post <- matrix(res$pdSSig00, n.a, unit.par, byrow=TRUE)[,order.old]
      Sigma12.post <- matrix(res$pdSSig01, n.a, unit.par, byrow=TRUE)[,order.old]
      Sigma22.post <- matrix(res$pdSSig11, n.a, unit.par, byrow=TRUE)[,order.old]
    }
    W1.post <- matrix(res$pdSW1, n.a, unit.w, byrow=TRUE)[,order.old]
    W2.post <- matrix(res$pdSW2, n.a, unit.w, byrow=TRUE)[,order.old]
    if (predict) {
      W1.pred <- matrix(res$pdSWt1, n.a, unit.w, byrow=TRUE)[,order.old]
      W2.pred <- matrix(res$pdSWt2, n.a, unit.w, byrow=TRUE) [,order.old]
    }
    a.post <- matrix(res$pdSa, n.a, unit.a, byrow=TRUE)
    nstar <- matrix(res$pdSn, n.a, unit.a, byrow=TRUE)
    
    if (parameter && predict) {
      res.out <- list(model="Dirichlet Process Prior", alpha=alpha,
                      burnin=burnin, thin=thin, X=X, Y=Y, nu0=nu0, tau0=tau0, mu0=mu0,
                      S0=S0, a0=a0, b0=b0, call=call,
                      mu1.post=mu1.post, mu2.post=mu2.post, 
                      Sigma11.post=Sigma11.post,
                      Sigma12.post=Sigma12.post, Sigma22.post=Sigma22.post, 
                      W1.post=W1.post, W2.post=W2.post,
                      W1.pred=W1.pred, W2.pred=W2.pred, a.post=a.post, nstar=nstar) }   
    else if (parameter && !predict) {
      res.out <- list(model="Dirichlet Process Prior", alpha=alpha,
                      burnin=burnin, thin=thin, X=X, Y=Y, nu0=nu0, tau0=tau0, mu0=mu0,
                      S0=S0, a0=a0, b0=b0, call=call,
                      mu1.post=mu1.post, mu2.post=mu2.post, 
                      Sigma11.post=Sigma11.post,
                      Sigma12.post=Sigma12.post, Sigma22.post=Sigma22.post, 
                      W1.post=W1.post, W2.post=W2.post, a.post=a.post, nstar=nstar)   }
    else if (!parameter && predict) 
{
      res.out <- list(model="Dirichlet Process Prior", alpha=alpha,
                      burnin=burnin, thin=thin, X=X, Y=Y, nu0=nu0, tau0=tau0, mu0=mu0,
                      S0=S0, a0=a0, b0=b0, call=call,
                      W1.post=W1.post, W2.post=W2.post,
                      W1.pred=W1.pred, W2.pred=W2.pred, a.post=a.post, nstar=nstar)  }
    else {
      res.out <- list(model="Dirichlet Process Prior", alpha=alpha,
                      burnin=burnin, thin=thin, X=X, Y=Y, nu0=nu0, tau0=tau0, mu0=mu0,
                      S0=S0, a0=a0, b0=b0, call=call,
                      W1.post=W1.post, W2.post=W2.post, a.post=a.post, nstar=nstar)     }
  }	 
  else{ # parametric model
    n.a <- floor((n.draws-burnin)/thin)
    n.par <- n.a
    n.w <- n.a * (n.samp+samp.X1+samp.X0) 
    unit.a <- 1
    unit.par <- 1
    unit.w <- (n.samp+samp.X1+samp.X0) 	
    res <- .C("cBaseeco", as.double(d), as.integer(n.samp),
	      as.integer(n.draws), as.integer(burnin), as.integer(thin),
	      as.integer(verbose),
              as.integer(nu0), as.double(tau0), as.double(mu0),
              as.double(S0),
              as.integer(survey.yes), as.integer(survey.samp), as.double(survey.data),
   	      as.integer(X1type), as.integer(samp.X1), as.double(X1.W1),
   	      as.integer(X0type), as.integer(samp.X0), as.double(X0.W2),
	      as.integer(predict), as.integer(parameter), 
	      pdSMu0=double(n.par),
              pdSMu1=double(n.par), pdSSig00=double(n.par),
              pdSSig01=double(n.par), pdSSig11=double(n.par),
              pdSW1=double(n.w), pdSW2=double(n.w), 
              pdSWt1=double(n.w), pdSWt2=double(n.w), PACKAGE="eco")

    if (parameter) {
      mu.post <- cbind(matrix(res$pdSMu0, n.a, unit.par, byrow=TRUE),
                       matrix(res$pdSMu1, n.a, unit.par, byrow=TRUE)) 
      colnames(mu.post) <- c("mu1", "mu2")
      Sigma.post <- cbind(matrix(res$pdSSig00, n.a, unit.par, byrow=TRUE), 
                          matrix(res$pdSSig01, n.a, unit.par, byrow=TRUE),
                          matrix(res$pdSSig11, n.a, unit.par, byrow=TRUE))
      colnames(Sigma.post) <- c("Sigma11", "Sigma12", "Sigma22")
    }
    W1.post <- matrix(res$pdSW1, n.a, unit.w, byrow=TRUE)[,order.old]
    W2.post <- matrix(res$pdSW2, n.a, unit.w, byrow=TRUE)[,order.old]
    if (predict) {
      W1.pred <- matrix(res$pdSWt1, n.a, unit.w, byrow=TRUE)[,order.old]
      W2.pred <- matrix(res$pdSWt2, n.a, unit.w, byrow=TRUE)[,order.old] 
    }
    
    if (parameter && predict)
      res.out <- list(model="Normal prior", burnin=burnin, thin = thin, X=X, Y=Y,
                      nu0=nu0, tau0=tau0, mu0=mu0, S0=S0, call=call, mu.post=mu.post,
                      Sigma.post=Sigma.post, W1.post=W1.post, W2.post=W2.post,
                      W1.pred=W1.pred, W2.pred=W2.pred)
    else if (parameter && !predict)
      res.out <- list(model="Normal prior", burnin=burnin, thin = thin, X=X, Y=Y,
                      nu0=nu0, tau0=tau0, mu0=mu0, S0=S0, call=call, mu.post=mu.post,
                      Sigma.post=Sigma.post, W1.post=W1.post, W2.post=W2.post)
    else if (!parameter && predict)
      res.out <- list(model="Normal prior", burnin=burnin, thin = thin, X=X, Y=Y,
                      nu0=nu0, tau0=tau0, mu0=mu0, S0=S0, call=call,
                      W1.post=W1.post, W2.post=W2.post, 
                      W1.pred=W1.pred, W2.pred=W2.pred)
    else if (!parameter && !predict)
      res.out <- list(model="Normal prior", burnin=burnin, thin = thin, X=X, Y=Y,
                      nu0=nu0, tau0=tau0, mu0=mu0, S0=S0, 
                      W1.post=W1.post, W2.post=W2.post)
    
  }
  class(res.out) <- "eco"
  return(res.out)
}


