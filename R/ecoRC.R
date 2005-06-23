ecoRC <- function(formula, data = parent.frame(),
                  mu0 = 0, tau0 = 2, nu0 = 4, S0 = 10, reject = TRUE,
                  parameter = TRUE, n.draws = 5000,
                  burnin = 0, thin = 0, verbose = FALSE){ 
  
  ## checking inputs
  if (burnin >= n.draws)
    stop("n.draws should be larger than burnin")
  
  call <- match.call()

  ## getting X, Y, and N
  tt <- terms(formula)
  attr(tt, "intercept") <- 0
  if (is.matrix(eval.parent(call$data)))
    data <- as.data.frame(data)
  X <- model.matrix(tt, data)
  n.samp <- nrow(X)
  C <- ncol(X)
  Y <- matrix(model.response(model.frame(tt, data = data)), nrow = n.samp)
  R <- ncol(Y)

  ## fitting the model
  n.store <- floor((n.draws-burnin)/(thin+1))
  n.par <- R
  tmp <- bounds(Y ~ X)
  S0 <- diag(S0, C)
  mu0 <- rep(mu0, C)
  
  res <- .C("cBase2C", as.double(X), as.double(Y),
            as.double(tmp$Wmin), as.double(tmp$Wmax),
            as.integer(n.samp), as.integer(C),
            as.integer(reject),
            as.integer(n.draws), as.integer(burnin),
            as.integer(thin+1), as.integer(verbose),
            as.integer(nu0), as.double(tau0),
            as.double(mu0), as.double(S0),
            as.integer(parameter), pdSmu = double(n.store*C),
            pdSSigma = double(n.store*C*(C+1)/2),
            pdSW = double(n.store*n.samp*C), PACKAGE="eco")

  res.out <- list(call = call, X = X, Y = Y, Wmin = tmp$Wmin, Wmax = tmp$Wmax)
  res.out$mu <- matrix(res$pdSmu, n.store, C, byrow=TRUE)
  res.out$Sigma <- matrix(res$pdSSigma, n.store, C*(C+1)/2, byrow=TRUE)
  res.out$W <- array(res$pdSW, c(C, n.samp, n.store))
  
  class(res.out) <- c("ecoRC", "eco")
  return(res.out)
}

