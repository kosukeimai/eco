varcov <- function(object, ...)
  UseMethod("varcov")

varcov.ecoNP <- function(object, subset = NULL, ...) {
  if (is.null(subset))
    subset <- 1:nrow(object$Sigma11)
  else if (max(subset) > nrow(object$Sigma11))
    stop(paste("invalid input for `subset.' only", nrow(object$Sigma11), "draws are stored."))
  
  sigma11 <- object$Sigma11[subset,]
  sigma12 <- object$Sigma12[subset,]
  sigma22 <- object$Sigma22[subset,]

  p <- 2
  n <- length(sigma11)
  Sigma <- array(0, c(p, p, n))
  cov <- as.matrix(cbind(as.vector(sigma11), as.vector(sigma12), as.vector(sigma22)))
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
