ecoRC <- function(formula, data = parent.frame(),
                  mu0 = 0, tau0 = 2, nu0 = 4, S0 = 10, mu.start = 0,
                  Sigma.start = 1, reject = TRUE, maxit = 10e5,
                  parameter = TRUE,
                  n.draws = 5000, burnin = 0, thin = 0, verbose = FALSE){ 
  
  ## checking inputs
  if (burnin >= n.draws)
    stop("n.draws should be larger than burnin")
  mf <- match.call()

  ## getting X, Y, and N
  tt <- terms(formula)
  attr(tt, "intercept") <- 0
  if (is.matrix(eval.parent(mf$data)))
    data <- as.data.frame(data)
  X <- model.matrix(tt, data)
  n.samp <- nrow(X)
  C <- ncol(X)
  Y <- matrix(model.response(model.frame(tt, data = data)),
              nrow = n.samp)
  R <- ncol(Y)

  ## fitting the model
  n.store <- floor((n.draws-burnin)/(thin+1))
  n.par <- R-1
  tmp <- ecoBD(formula, data=data)
  mu0 <- rep(mu0, C)
  S0 <- diag(S0, C)

  res.out <- list(call = mf, X = X, Y = Y, Wmin = tmp$Wmin, Wmax = tmp$Wmax)
  if (R == 1) {
    mu.start <- rep(mu.start, C)
    Sigma.start <- diag(Sigma.start, C)
    res <- .C("cBase2C", as.double(X), as.double(Y),
              as.double(tmp$Wmin[,1,]), as.double(tmp$Wmax[,1,]),
              as.integer(n.samp), as.integer(C), as.integer(reject),
              as.integer(maxit), as.integer(n.draws), as.integer(burnin),
              as.integer(thin+1), as.integer(verbose),
              as.integer(nu0), as.double(tau0),
              as.double(mu0), as.double(S0), as.double(mu.start),
              as.double(Sigma.start),
              as.integer(parameter), pdSmu = double(n.store*C),
              pdSSigma = double(n.store*C*(C+1)/2),
              pdSW = double(n.store*n.samp*C), PACKAGE="eco")
    res.out$mu <- matrix(res$pdSmu, n.store, C, byrow=TRUE)
    res.out$Sigma <- matrix(res$pdSSigma, n.store, C*(C+1)/2, byrow=TRUE)
    res.out$W <- array(res$pdSW, c(C, n.samp, n.store))
  }
  else {
    mu.start <- matrix(rep(rep(mu.start, C), R-1), ncol = R-1, nrow = C,
                       byrow = FALSE)
    Sigma.start <- array(rep(diag(Sigma.start, C), R-1), c(C, C, R-1))
    res <- .C("cBaseRC", as.double(X), as.double(Y[,1:(R-1)]),
              as.double(tmp$Wmin[,1:(R-1),]), as.double(tmp$Wmax[,1:(R-1),]),
              as.integer(n.samp), as.integer(C), as.integer(R),
              as.integer(reject), as.integer(maxit),
              as.integer(n.draws), as.integer(burnin),
              as.integer(thin+1), as.integer(verbose),
              as.integer(nu0), as.double(tau0),
              as.double(mu0), as.double(S0),
              as.double(mu.start), as.double(Sigma.start),
              as.integer(parameter), pdSmu = double(n.store*(R-1)*C),
              pdSSigma = double(n.store*(R-1)*C*(C+1)/2),
              pdSW = double(n.store*n.samp*(R-1)*C), PACKAGE="eco")
    res.out$mu <- array(res$pdSmu, c(R-1, C, n.store))
    res.out$Sigma <- array(res$pdSSigma, c(R-1, C*(C+1)/2, n.store))
    res.out$W <- array(res$pdSW, c(R-1, C, n.samp, n.store))
  }
  
  class(res.out) <- c("ecoRC", "eco")
  return(res.out)
}


