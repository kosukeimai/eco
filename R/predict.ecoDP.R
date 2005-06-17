predict.ecoNP <- function(object, newdraw = NULL, verboe=FALSE...){
  if ((is.null(object$mu1) || is.null(object$mu2)) && (is.null(newdraw$mu1) || is.null(newdraw$mu2)))
    stop("Posterior draws of the means mu1 and mu2 must be supplied.")

  if ((is.null(object$Sigma11) || is.null(object$Sigma22) || is.null(object$Sigma12)) &&
     (is.null(newdraw$Sigma11) || is.null(newdraw$Sigma22) || is.null(newdraw$Sigma12)))
    stop("Posterior draws of the elements of the covariance matrix Sigma11, Sigma12 and Sigma22 must be supplied.")   

  mu1 <- object$mu1
  mu2 <- object$mu2

  p <-2 
  n <- ncol(mu1)
  ndraws <- nrow(mu1)
  mu<-cbind(as.vector(mu1), as.vector(mu2)) #as.vector reads by columns
  
  Sigma <- cov.eco(object) ##this needs to be changed.
  n.pred <-n*ndraws
  Wstar <- matrix(NA, nrow=n.pred, ncol=p)
  for (i in 1:n.pred) 
    Wstar[i,] <- mvrnorm(1, mu = mu[i,], Sigma = Sigma[,,i])
  return(apply(Wstar, 2, invlogit))
}
