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
  
  i <- checkdata(X,Y, supplement)
  
  ## fitting the model
  n.a <- floor((n.draws-burnin)/(thin+1))

  unit.par <- unit.w <- i$n.samp+i$samp.X1+i$samp.X0
  n.par <- n.a * unit.par
  n.w <- n.a * unit.w

  unit.a <- 1

  res <- .C("cDPeco", as.double(i$d), as.integer(i$n.samp),
            as.integer(n.draws), as.integer(burnin), as.integer(thin+1),
            as.integer(verbose), as.integer(nu0), as.double(tau0),
            as.double(mu0), as.double(S0), as.double(alpha),
            as.integer(alpha.update), as.double(a0), as.double(b0),
            as.integer(i$survey.yes), as.integer(i$survey.samp),
            as.double(i$survey.data), as.integer(i$X1type),
            as.integer(i$samp.X1), as.double(i$X1.W1),
            as.integer(i$X0type), as.integer(i$samp.X0), as.double(i$X0.W2),
            as.integer(predict), as.integer(parameter), as.integer(grid),
            pdSMu0=double(n.par), pdSMu1=double(n.par),
            pdSSig00=double(n.par), pdSSig01=double(n.par),
            pdSSig11=double(n.par), pdSW1=double(n.w), pdSW2=double(n.w), 
            pdSa=double(n.a), pdSn=integer(n.a), PACKAGE="eco")
  if (parameter) {
    mu1.post <- matrix(res$pdSMu0, n.a, unit.par, byrow=TRUE)[,i$order.old]
    mu2.post <- matrix(res$pdSMu1, n.a, unit.par, byrow=TRUE)[,i$order.old]
    Sigma11.post <- matrix(res$pdSSig00, n.a, unit.par, byrow=TRUE)[,i$order.old]
    Sigma12.post <- matrix(res$pdSSig01, n.a, unit.par, byrow=TRUE)[,i$order.old]
    Sigma22.post <- matrix(res$pdSSig11, n.a, unit.par, byrow=TRUE)[,i$order.old]
  }
  mu <- array(rbind(mu1,mu2), c(n.a, 2, unit.par))
  Sigma <- array(rbind(Sigma11, Sigma12, Sigma22), c(n.a, 3, unit.par))
  W1.post <- matrix(res$pdSW1, n.a, unit.w, byrow=TRUE)[,i$order.old]
  W2.post <- matrix(res$pdSW2, n.a, unit.w, byrow=TRUE)[,i$order.old]

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


