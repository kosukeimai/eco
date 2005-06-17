eco <- function(formula, data = parent.frame(), supplement = NULL,
                mu0 = c(0,0), tau0 = 2, nu0 = 4, S0 = diag(10,2),
                parameter = TRUE, grid = FALSE, n.draws = 5000,
                burnin = 0, thin = 0, verbose = FALSE){ 

  ## checking inputs
  if (burnin >= n.draws)
    stop("n.draws should be larger than burnin")
  
  call <- match.call()

  ## getting X and Y
  tt <- terms(formula)
  attr(tt, "intercept") <- 0
  if (is.matrix(eval.parent(call$data)))
    data <- as.data.frame(data)
  X <- model.matrix(tt, data)
  Y <- model.response(model.frame(tt, data = data))

  # check data and modify inputs 
  tmp <- checkdata(X,Y, supplement)  

  ## fitting the model
  n.a <- floor((n.draws-burnin)/(thin+1))
  n.par <- n.a

  unit.a <- unit.par <- 1
  unit.w <- tmp$n.samp+tmp$samp.X1+tmp$samp.X0 	

  n.w <- n.a * unit.w

  res <- .C("cBaseeco", as.double(tmp$d), as.integer(tmp$n.samp),
            as.integer(n.draws), as.integer(burnin), as.integer(thin+1),
            as.integer(verbose), as.integer(nu0), as.double(tau0),
            as.double(mu0), as.double(S0), as.integer(tmp$survey.yes),
            as.integer(tmp$survey.samp), as.double(tmp$survey.data),
            as.integer(tmp$X1type), as.integer(tmp$samp.X1),
            as.double(tmp$X1.W1), as.integer(tmp$X0type),
            as.integer(tmp$samp.X0), as.double(tmp$X0.W2),
            as.integer(parameter), as.integer(grid), 
            pdSMu0=double(n.par), pdSMu1=double(n.par), pdSSig00=double(n.par),
            pdSSig01=double(n.par), pdSSig11=double(n.par),
            pdSW1=double(n.w), pdSW2=double(n.w),
            PACKAGE="eco")
  
  if (parameter) {
    mu.post <- cbind(matrix(res$pdSMu0, n.a, unit.par, byrow=TRUE),
                     matrix(res$pdSMu1, n.a, unit.par, byrow=TRUE)) 
    colnames(mu.post) <- c("mu1", "mu2")
    Sigma.post <- cbind(matrix(res$pdSSig00, n.a, unit.par, byrow=TRUE), 
                        matrix(res$pdSSig01, n.a, unit.par, byrow=TRUE),
                        matrix(res$pdSSig11, n.a, unit.par, byrow=TRUE))
    colnames(Sigma.post) <- c("Sigma11", "Sigma12", "Sigma22")
  }
  W1.post <- matrix(res$pdSW1, n.a, unit.w, byrow=TRUE)[,tmp$order.old]
  W2.post <- matrix(res$pdSW2, n.a, unit.w, byrow=TRUE)[,tmp$order.old]
  
  res.out <- list(call = call, X = X, Y = Y, W1 = W1.post, W2 = W2.post,
                  burin = burnin, thin = thin, nu0 = nu0, tau0 = tau0,
                  mu0 = mu0, S0 = S0)
  
  if (parameter) {
    res.out$mu <- mu.post
    res.out$Sigma <- Sigma.post
    res.out$W1 <- W1.post
    res.out$W2 <- W2.post
  }
  class(res.out) <- "eco"
  return(res.out)
}


