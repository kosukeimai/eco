predict.ecoNPX <- function(object, newdraw = NULL, subset = NULL,
                           obs = NULL, cond = FALSE, verbose = FALSE, ...){

  if (is.null(newdraw) && is.null(object$mu))
    stop("Posterior draws of mu and Sigma must be supplied")
  else if (!is.null(newdraw)){
    if (is.null(newdraw$mu) && is.null(newdraw$Sigma))
      stop("Posterior draws of both mu and Sigma must be supplied.")
    object <- newdraw
  }

  n.draws <- dim(object$mu)[1]
  p <- dim(object$mu)[2]
  n <- dim(object$mu)[3]
  mu <- aperm(coef(object, subset = subset, obs = obs), c(2,3,1))
  
  if (is.null(subset))
    subset <- 1:n.draws
  if (is.null(obs))
    obs <- 1:n
  Sigma <- aperm(res$Sigma[subset,,obs], c(2,3,1))

  if (cond) { # conditional prediction
    
  }
  else { # unconditional prediction
    res <- .C("preDP", as.double(mu), as.double(Sigma), as.integer(n),
              as.integer(n.draws), as.integer(p), as.integer(verbose),
              pdStore = double(n.draws*p*n), PACKAGE="eco")$pdStore
    
    res <- matrix(res, ncol=p, nrow=n.draws*n, byrow=TRUE)
    colnames(res) <- c("W1", "W2", "X")
  }
  
  class(res) <- c("predict.eco", "matrix")
  return(res)
}
