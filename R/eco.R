eco <- function(formula, data = parent.frame(), N = NULL, supplement = NULL,
                mu0 = c(0,0), tau0 = 2, nu0 = 4, S0 = diag(10,2),
                mu.start = c(0,0), Sigma.start = diag(10,2),
                parameter = TRUE, grid = FALSE, n.draws = 5000,
                burnin = 0, thin = 0, verbose = FALSE){ 

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
  Y <- model.response(model.frame(tt, data = data))
  N <- eval(mf$N, data)
  
  # check data and modify inputs 
  tmp <- checkdata(X,Y, supplement)  
  bdd <- ecoBD(formula=formula, data=data)

  ## fitting the model
  n.store <- floor((n.draws-burnin)/(thin+1))
  unit.par <- 1
  unit.w <- tmp$n.samp+tmp$samp.X1+tmp$samp.X0 	
  n.w <- n.store * unit.w

  res <- .C("cBaseeco", as.double(tmp$d), as.integer(tmp$n.samp),
            as.integer(n.draws), as.integer(burnin), as.integer(thin+1),
            as.integer(verbose), as.integer(nu0), as.double(tau0),
            as.double(mu0), as.double(S0), as.double(mu.start),
            as.double(Sigma.start), as.integer(tmp$survey.yes),
            as.integer(tmp$survey.samp), as.double(tmp$survey.data),
            as.integer(tmp$X1type), as.integer(tmp$samp.X1),
            as.double(tmp$X1.W1), as.integer(tmp$X0type),
            as.integer(tmp$samp.X0), as.double(tmp$X0.W2),
	    as.double(bdd$Wmin[,1,1]), as.double(bdd$Wmax[,1,1]),
            as.integer(parameter), as.integer(grid), 
            pdSMu0=double(n.store), pdSMu1=double(n.store), pdSSig00=double(n.store),
            pdSSig01=double(n.store), pdSSig11=double(n.store),
            pdSW1=double(n.w), pdSW2=double(n.w),
            PACKAGE="eco")
  
  W1.post <- matrix(res$pdSW1, n.store, unit.w, byrow=TRUE)[,tmp$order.old]
  W2.post <- matrix(res$pdSW2, n.store, unit.w, byrow=TRUE)[,tmp$order.old]
  W <- array(rbind(W1.post, W2.post), c(n.store, 2, unit.w))
  colnames(W) <- c("W1", "W2")
  res.out <- list(call = mf, X = X, Y = Y, N = N, W = W,
                  Wmin=bdd$Wmin[,1,], Wmax = bdd$Wmax[,1,],
		  burin = burnin, thin = thin, nu0 = nu0,
                  tau0 = tau0, mu0 = mu0, S0 = S0)

  if (parameter) {
    res.out$mu <- cbind(matrix(res$pdSMu0, n.store, unit.par, byrow=TRUE),
                        matrix(res$pdSMu1, n.store, unit.par, byrow=TRUE)) 
    colnames(res.out$mu) <- c("mu1", "mu2")
    res.out$Sigma <- cbind(matrix(res$pdSSig00, n.store, unit.par, byrow=TRUE), 
                           matrix(res$pdSSig01, n.store, unit.par, byrow=TRUE),
                           matrix(res$pdSSig11, n.store, unit.par, byrow=TRUE))
    colnames(res.out$Sigma) <- c("Sigma11", "Sigma12", "Sigma22")
  }
  class(res.out) <- "eco"
  return(res.out)
}


