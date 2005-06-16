predict.eco <- function(object, newdata = NULL, ...){

  if (is.null(object$mu))
    stop("Parameters are not stored: set `parameter = TRUE'")
  if (is.null(object$Sigma))
    stop("Parameters are not stored: set `parameter = TRUE'")
  mu <- object$mu
  n <- nrow(mu)
  p <- ncol(mu)
  Sigma <- cov.eco(object)
  Wstar <- matrix(NA, nrow=n, ncol=p)
  for (i in 1:n) 
    Wstar[i,] <- mvrnorm(1, mu = mu[i,], Sigma = Sigma[,,i])
  return(apply(Wstar, 2, invlogit))
}
