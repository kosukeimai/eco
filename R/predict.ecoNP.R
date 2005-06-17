predict.ecoNP <- function(object, newdraw = NULL, subset = NULL, verbose=FALSE...){
  if (!is.null(newdraw)) {
  if ((is.null(newdraw$mu1) || is.null(newdraw$mu2)))
    stop("Posterior draws of the means mu1 and mu2 must be supplied.")
  if ((is.null(newdraw$Sigma11) || is.null(newdraw$Sigma22) || is.null(newdraw$Sigma12)))
    stop("Posterior draws of the elements of the covariance matrix Sigma11, Sigma12 and Sigma22 must be supplied.")   
   object <- newdraw
 }

  mu <- coef.ecoNP(object, subset = subset)
  n.draws <- nrow(mu)

  p <- ncol(mu)

  Sigma <- cov.ecoNP(object, subset=subset) 
  Wstar <- matrix(NA, nrow=n.draws, ncol=p)

  tmp <- floor(n.draws/10)
  inc <- 1
  for (i in 1:n.draws) {
    Wstar[i,] <- mvrnorm(1, mu = mu[i,], Sigma = Sigma[,,i])
    if (i == inc*tmp & verbose) {
      cat("", inc*10, "percent done.\n")
      inc <- inc + 1
    }
  }
  res <- apply(Wstar, 2, invlogit)
  class(res) <- c("predict.eco", "predict.ecoNP", "matrix")
  return(res)

}
