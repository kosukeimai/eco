predict.eco <- function(object, newdata = NULL, newdraw = NULL, ...){

  if (is.null(object$mu) && is.null(newdraw$mu))
    stop("Posterior draws of mu must be supplied.")
  if (is.null(object$Sigma) && is.null(newdraw$Sigma))
    stop("Posterior draws of Sigma must be supplied.")   
  mu <- object$mu
  n.draws <- nrow(mu)
  p <- ncol(mu)
  Sigma <- cov.eco(object)
  Wstar <- matrix(NA, nrow=n.draws, ncol=p)
  for (i in 1:n.draws) 
    Wstar[i,] <- mvrnorm(1, mu = mu[i,], Sigma = Sigma[,,i])
  return(apply(Wstar, 2, invlogit))
}
