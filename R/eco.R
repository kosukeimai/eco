".First.lib" <- function(lib, pkg) 
  library.dynam("eco", pkg, lib)


eco <- function(Y, X, data = parent.frame(), n.draws = 5000, nonpar =
                TRUE, link = "logit", nu0 = 4, tau0 = 1, mu0 = c(0,0),
                S0 = diag(10,2), alpha = NULL, burnin = 0, thin = 5,
                predict = TRUE, a0 = 0.1, b0 = 0.1){ 

  ## checking inputs
  if (link == "logit")
   nlink <- 1
  else if (link == "probit")
    nlink <- 2
  else if (link == "cloglog")
    nlink <- 3
  else
    stop("Error: invalid input for `link'")
  if (burnin >= n.draws)
    stop("Error: n.draws should be larger than burnin") 
  m <- match.call()
  ff <- as.formula(paste(m$Y, "~ -1 +", m$X))
  if (is.matrix(eval.parent(m$data)))
    data <- as.data.frame(data)
  X <- model.matrix(ff, data)
  Y <- model.response(model.frame(ff, data=data))
  
  ## fitting the model
  n.samp <- length(Y)	 
  d <- cbind(X, Y)
  if (nonpar){	# nonparametric model
    n.a <- floor((n.draws-burnin)/thin)
    n.par <- n.a * n.samp 
    n.w <- n.a * n.samp
    unit.a <- 1
    unit.par <- n.samp
    unit.w <- n.samp	
    res <- .C("cDPeco", as.integer(n.draws), as.integer(nlink),
              as.double(d), as.integer(n.samp), as.integer(nu0),
              as.double(tau0), as.double(mu0), as.double(S0),
              as.integer(alpha), as.integer(burnin), as.integer(thin),
              as.integer(predict), as.double(a0), as.double(b0),
              pdSMu0=double(n.par), pdSMu1=double(n.par),
              pdSSig00=double(n.par), pdSSig01=double(n.par),
              pdSSig11=double(n.par), pdSW1=double(n.w),
              pdSW2=double(n.w), pdSWt1=double(n.w),
              pdSWt2=double(n.w), pdSa=double(n.a), pdSn=integer(n.a))
    
    mu.post <- cbind(matrix(res$pdSMu0, n.a, unit.par, byrow=T),
                     matrix(res$pdSMu1, n.a, unit.par, byrow=T))
    colnames(mu.post) <- c("mu1", "mu2")
    Sigma.post <- cbind(matrix(res$pdSSig00, n.a, unit.par, byrow=T),
                        matrix(res$pdSSig01, n.a, unit.par, byrow=T),
                        matrix(res$pdSSig11, n.a, unit.par, byrow=T))
    colnames(Sigma.post) <- c("Sigma11", "Sigma12", "Sigma22")
    W.post <- cbind(matrix(res$pdSW1, n.a, unit.w, byrow=T),
                    matrix(res$pdSW2, n.a, unit.w, byrow=T))
    W.pred <- cbind(matrix(res$pdSWt1, n.a, unit.w, byrow=T),
                    matrix(res$pdSWt2, n.a, unit.w, byrow=T) 
    colnames(W.post) <- colnames(W.pred) <- c("W1", "W2")
    a.post <- matrix(res$pdSa, n.a, unit.a, byrow=T)
    nstar <- matrix(res$pdSn, n.a, unit.a, byrow=T)
    res.out <- list(model="Dirichlet Process Prior", alpha=alpha,
                    burnin=burnin, thin=thin, X=X, Y=Y, nu0=nu0, tau0=tau0, mu0=mu0,
                    S0=S0, a0=a0, b0=b0, mu.post=mu.post, Sigma.post=Sigma.post,
                    W.post=W.post, W.pred=W.pred, a.post=a.post, nstar=nstar)   
  }	 
  else{ # parametric model
    n.a <- floor((n.draws-burnin)/thin)
    n.par <- n.a
    n.w <- n.a * n.samp
    unit.a <- 1
    unit.par <- 1
    unit.w <- n.samp	
    res <- .C("cBaseeco", as.integer(n.draws), as.integer(nlink),
              as.double(d), as.integer(n.samp), as.integer(nu0),
              as.double(tau0), as.double(mu0),
              as.double(S0),as.integer(burnin), as.integer(thin),
              as.integer(predict), pdSMu0=double(n.par),
              pdSMu1=double(n.par), pdSSig00=double(n.par),
              pdSSig01=double(n.par), pdSSig11=double(n.par),
              pdSW1=double(n.w), pdSW2=double(n.w), 
              pdSWt1=double(n.w), pdSWt2=double(n.w))
    mu.post <- cbind(matrix(res$pdSMu0, n.a, unit.par, byrow=T),
                     matrix(res$pdSMu1, n.a, unit.par, byrow=T)) 
    colnames(mu.post) <- c("mu1", "mu2")
    Sigma.post <- cbind(matrix(res$pdSSig00, n.a, unit.par, byrow=T), 
                        matrix(res$pdSSig01, n.a, unit.par, byrow=T),
                        matrix(res$pdSSig11, n.a, unit.par, byrow=T))
    colnames(Sigma.post) <- c("Sigma11", "Sigma12", "Sigma22")
    W.post <- cbind(matrix(res$pdSW1, n.a, unit.w, byrow=T),
                    matrix(res$pdSW2, n.a, unit.w, byrow=T))
    W.pred <- cbind(matrix(res$pdSWt1, n.a, unit.w, byrow=T),
                    matrix(res$pdSWt2, n.a, unit.w, byrow=T))
    colnames(W.post) <- colnames(W.pred) <- c("W1", "W2")
    res.out <- list(model="Normal prior", burnin=burnin, thin = thin, X=X, Y=Y,
                    nu0=nu0, tau0=tau0, mu0=mu0, S0=S0, mu.post=mu.post,
                    Sigma.post=Sigma.post, W.post=W.post, W.pred=W.pred)
  }
  class(res.out) <- "eco"
  return(res.out)
}


