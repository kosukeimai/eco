predict.eco <- function(object, newdraw = NULL, verbose = FALSE, ...){

  if (is.null(object$mu) && is.null(newdraw$mu))
    stop("Posterior draws of mu must be supplied.")
  if (is.null(object$Sigma) && is.null(newdraw$Sigma))
    stop("Posterior draws of Sigma must be supplied.")   
  mu <- object$mu
  n.draws <- nrow(mu)
  p <- ncol(mu)
  Sigma <- cov.eco(object)
  Wstar <- matrix(NA, nrow=n.draws, ncol=p)
  if (verbose) {
    tmp <- floor(n.draws/10)
    inc <- 1
  }
  for (i in 1:n.draws) {
    Wstar[i,] <- mvrnorm(1, mu = mu[i,], Sigma = Sigma[,,i])
    if (i == inc*tmp & verbose) {
      cat("", inc*10, "percent done.\n")
      inc <- inc + 1
    }
  }
  res <- apply(Wstar, 2, invlogit)
  colnames(res) <- c("W1", "W2")
  return(res)
}
