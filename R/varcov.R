varcov <- function(object, ...)
  UseMethod("varcov")

varcov.eco <- function(object, subset = NULL, ...) {
  if (is.null(subset))
    subset <- 1:nrow(object$mu)
  else if (max(subset) > nrow(object$mu))
    stop(paste("invalid input for `subset.' only", nrow(object$mu), "draws are stored."))

  p <- ncol(object$mu)
  if (length(dim(object$Sigma))==2) {
    nobs <-1
    n <- length(subset)
    Sigma <- array(0, c(p, p, n))
    cov <- object$Sigma[subset,]
  }
  else if (length(dim(object$Sigma))==3) {
   nobs <- dim(object$Sigma)[3] 
   n <- length(subset)* nobs
    Sigma <- array(0, c(p, p, n))
    cov <- NULL
    for (j in nobs:1) cov <- rbind(object$Sigma[subset,,j],cov)
  }

  for (i in 1:n) {
    count <- 1
    for (j in 1:p) {
      Sigma[j,j:p,i] <- cov[i,count:(count+p-j)]
      count <- count + p - j + 1
    }
    diag(Sigma[,,i]) <- diag(Sigma[,,i]/2)
    Sigma[,,i] <- Sigma[,,i] + t(Sigma[,,i])
  }
  if (n > 1)
    return(Sigma)
  else
    return(Sigma[,,1])  
}
