varcov <- function(object, ...)
  UseMethod("varcov")

varcov.junk <- function(object, subset = NULL, ...) {
  if (is.null(subset))
    subset <- 1:nrow(object$Sigma[,,1])
  else if (max(subset) > nrow(object$Sigma[,,1]))
    stop(paste("invalid input for `subset.' only", nrow(object$Sigma[,,1]), "draws are stored."))
  
  cov <- object$Sigma

  p <- ncol(object$mu[,,1])
  n.draws <- length(subset)
  n.obs <- length(cov[1,1,])
  n <-n.draws*n.obs
  Sigma <- array(0, c(p, p, n))

  for (i in 1:n.draws) 
    for (k in 1:n.obs) {
      count <- 1
      for (j in 1:p) {
       Sigma[j,j:p,(i-1)*n.obs+k] <- cov[subset[i],count:(count+p-j),k]
       count <- count + p - j + 1
      }
      diag(Sigma[,,(i-1)*n.obs+k]) <- diag(Sigma[,,(i-1)*n.obs+k]/2)
      Sigma[,,(i-1)*n.obs+k] <- Sigma[,,(i-1)*n.obs+k] + t(Sigma[,,(i-1)*n.obs+k])
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
