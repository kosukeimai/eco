eco <- function(formula, data = parent.frame(), supplement = NULL,
                mu0 = c(0,0), tau0 = 2, nu0 = 4, S0 = diag(10,2),
                predict = FALSE, parameter = TRUE, grid = FALSE,
                n.draws = 5000, burnin = 0, thin = 0,
                verbose = FALSE){ 

  ## checking inputs
  if (burnin >= n.draws)
    stop("n.draws should be larger than burnin")
  
  if ((dim(supplement)[2] != 2) && (length(supplement)>0))
    stop("use n by 2 matrix for survey data")
  
  call <- match.call()

  ## getting X and Y
  tt <- terms(formula)
  attr(tt, "intercept") <- 0
  if (is.matrix(eval.parent(call$data)))
    data <- as.data.frame(data)
  X <- model.matrix(tt, data)
  Y <- model.response(model.frame(tt, data = data))

  ## survey data
  if (is.null(supplement)) 
    survey.samp <- survey.data <- survey.yes <- 0
  else {
    survey.samp <- length(supplement[,1])
    survey.data <- as.matrix(supplement)
    survey.yes <- 1
  }
  
  ind <- 1:length(X)
  X1type <- X0type <- samp.X1 <- samp.X0 <- X1.W1 <- X0.W2<-0
  
  ## X = 1
  X1.ind<-ind[along=(X==1)]
  if (length(X[X!=1])<length(X)){
    X1type<-1
    samp.X1<-length(X1.ind)
    X1.W1<-Y[X1.ind]
  }
  
  ## X = 0 
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
  n.a <- floor((n.draws-burnin)/(thin+1))
  n.par <- n.a
  n.w <- n.a * (n.samp+samp.X1+samp.X0) 
  unit.a <- unit.par <- 1
  unit.w <- n.samp+samp.X1+samp.X0 	
  res <- .C("cBaseeco", as.double(d), as.integer(n.samp),
            as.integer(n.draws), as.integer(burnin), as.integer(thin+1),
            as.integer(verbose), as.integer(nu0), as.double(tau0),
            as.double(mu0), as.double(S0), as.integer(survey.yes),
            as.integer(survey.samp), as.double(survey.data),
            as.integer(X1type), as.integer(samp.X1), as.double(X1.W1),
            as.integer(X0type), as.integer(samp.X0), as.double(X0.W2),
            as.integer(predict), as.integer(parameter), as.integer(grid), 
            pdSMu0=double(n.par), pdSMu1=double(n.par), pdSSig00=double(n.par),
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
  
  res.out <- list(call = call, nonpar = nonpar, X = X, Y = Y,
                  burin = burnin, thin = thin, nu0 = nu0,
                  tau0 = tau0, mu0 = mu0, S0 = S0, W1 = W1.post,
                  W2 = W2.post)
  
  if (parameter) {
    res.out$mu <- mu.post
    res.out$Sigma <- Sigma.post
    res.out$W1 <- W1.post
    res.out$W2 <- W2.post
  }
  if (predict) {
    res.out$W1.pred <- W1.pred
    res.out$W2.pred <- W2.pred
  }
  class(res.out) <- "eco"
  return(res.out)
}


