ecoNP <- function(formula, data = parent.frame(), supplement = NULL,
                  mu0 = c(0,0), tau0 = 2, nu0 = 4, S0 = diag(10,2),
                  alpha = NULL, a0 = 1, b0 = 0.1, 
                  parameter = TRUE, grid = FALSE, n.draws = 5000,
                  burnin = 0, thin = 0, verbose = FALSE){ 

  ## checking inputs
  if (burnin >= n.draws)
    stop("n.draws should be larger than burnin")
  
  call <- match.call()

  tt <- terms(formula)
  attr(tt, "intercept") <- 0
  if (is.matrix(eval.parent(call$data)))
    data <- as.data.frame(data)
  X <- model.matrix(tt, data)
  Y <- model.response(model.frame(tt, data = data))

  ##alpha
  if (is.null(alpha)) {
    alpha.update <- TRUE
    alpha <- 0
  }
  else
    alpha.update <- FALSE
  
  tmp <- checkdata(X,Y, supplement)
  bdd <- bounds(formula, data=data)

  ## fitting the model
  n.a <- floor((n.draws-burnin)/(thin+1))

  unit.par <- unit.w <- tmp$n.samp+tmp$samp.X1+tmp$samp.X0
  n.par <- n.a * unit.par
  n.w <- n.a * unit.w

  unit.a <- 1

  res <- .C("cDPeco", as.double(tmp$d), as.integer(tmp$n.samp),
            as.integer(n.draws), as.integer(burnin), as.integer(thin+1),
            as.integer(verbose), as.integer(nu0), as.double(tau0),
            as.double(mu0), as.double(S0), as.double(alpha),
            as.integer(alpha.update), as.double(a0), as.double(b0),
            as.integer(tmp$survey.yes), as.integer(tmp$survey.samp),
            as.double(tmp$survey.data), as.integer(tmp$X1type),
            as.integer(tmp$samp.X1), as.double(tmp$X1.W1),
            as.integer(tmp$X0type), as.integer(tmp$samp.X0),
            as.double(tmp$X0.W2), 
	    as.double(bdd$Wmin[,1]), as.double(bdd$Wmax[,1]), 
            as.integer(parameter), as.integer(grid),
            pdSMu0=double(n.par), pdSMu1=double(n.par),
            pdSSig00=double(n.par), pdSSig01=double(n.par),
            pdSSig11=double(n.par), pdSW1=double(n.w), pdSW2=double(n.w), 
            pdSa=double(n.a), pdSn=integer(n.a), PACKAGE="eco")
  if (parameter) {
    mu1.post <- matrix(res$pdSMu0, n.a, unit.par, byrow=TRUE)[,tmp$order.old]
    mu2.post <- matrix(res$pdSMu1, n.a, unit.par, byrow=TRUE)[,tmp$order.old]
    Sigma11.post <- matrix(res$pdSSig00, n.a, unit.par, byrow=TRUE)[,tmp$order.old]
    Sigma12.post <- matrix(res$pdSSig01, n.a, unit.par, byrow=TRUE)[,tmp$order.old]
    Sigma22.post <- matrix(res$pdSSig11, n.a, unit.par, byrow=TRUE)[,tmp$order.old]
  
  mu <- array(rbind(mu1.post,mu2.post), c(n.a, 2, unit.par))
  Sigma <- array(rbind(Sigma11.post, Sigma12.post, Sigma22.post), c(n.a, 3, unit.par))
}
 
  W1.post <- matrix(res$pdSW1, n.a, unit.w, byrow=TRUE)[,tmp$order.old]
  W2.post <- matrix(res$pdSW2, n.a, unit.w, byrow=TRUE)[,tmp$order.old]

  a.post <- matrix(res$pdSa, n.a, unit.a, byrow=TRUE)
  nstar <- matrix(res$pdSn, n.a, unit.a, byrow=TRUE)
  
  res.out <- list(call = call, X = X, Y = Y, W1 = W1.post, W2 = W2.post,
                  burin = burnin, thin = thin, nu0 = nu0, tau0 = tau0,
                  mu0 = mu0, a0 = a0, b0 = b0, S0 = S0)


  if (parameter){
    res.out$mu <- mu
    res.out$Sigma <- Sigma
    if (alpha.update)
      res.out$alpha <- a.post
    else
      res.out$alpha <- alpha
    res.out$nstar <- res.out$nstar
  }
  

  class(res.out) <- c("ecoNP", "eco")
  return(res.out)
}


