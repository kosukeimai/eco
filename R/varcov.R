varcov <- function(object, ...)
  UseMethod("varcov")

varcov.ecoNP <- function(object, subset = NULL, ...) {
  if (is.null(subset))
    subset <- 1:nrow(object$Sigma[,,1])
  else if (max(subset) > nrow(object$Sigma[,,1]))
    stop(paste("invalid input for `subset.' only", nrow(object$Sigma[,,1]), "draws are stored."))
  

  p <- ncol(Sigma[,,1])
  n <- length(subset)
  nobs <- length(cov[1,1,])
  Sigma <- array(0, c(p, p, n))
  cov <- object$Sigma
  for (k in 1:nobs)
    for (i in 1:n) {
      count <- 1
      for (j in 1:p) {
       Sigma[j,j:p,(k-1)*nobs+i] <- cov[subset[i],count:(count+p-j),k]
       count <- count + p - j + 1
      }
      diag(Sigma[,,(k-1)*nobs+i]) <- diag(Sigma[,,(k-1)*nobs+i]/2)
      Sigma[,,(k-1)*nobs+i] <- Sigma[,,(k-1)*nobs+i] + t(Sigma[,,(k-1)*nobs+i])
  }
  if (n > 1)
    return(Sigma)
  else
    return(Sigma[,,1])  
}

varcov.eco <- function(object, subset = NULL, ...) {
  if (is.null(subset))
    subset <- 1:nrow(object$mu)
  else if (max(subset) > nrow(object$mu))
    stop(paste("invalid input for `subset.' only", nrow(object$mu), "draws are stored."))

  p <- ncol(object$mu)
  n <- length(subset)
  Sigma <- array(0, c(p, p, n))
  cov <- object$Sigma
  for (i in 1:n) {
    count <- 1
    for (j in 1:p) {
      Sigma[j,j:p,i] <- cov[subset[i],count:(count+p-j)]
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
