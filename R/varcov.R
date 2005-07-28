varcov <- function(object, ...)
  UseMethod("varcov")

varcov.eco <- function(object, subset = NULL, ...) {
  if (is.null(subset))
    subset <- 1:nrow(object$mu)
  else if (max(subset) > nrow(object$mu))
    stop(paste("invalid input for `subset.' only", nrow(object$mu), "draws are stored."))

  p <- ncol(object$mu)
  if (length(dim(object$Sigma))==2) {
    if (p==2) {
    nobs <-1
    n <- length(subset)
    Sigma <- array(0, c(p, p, n))
    cov <- object$Sigma[subset,]
   }
   else if (p==3) {
   #compute the corresponding terms in the conditional variance given X
   # cov hence expand to nsim*nobs
   p<-p-1 
   nobs<-length(object$X)
    n <- length(subset)*nobs
  Sigma <- array(0, c(p, p, n)) 
  temp <- matrix(NA, length(subset), 3)
 temp[,1]<-(object$Sigma[,1]-object$Sigma[,3]^2/object$Sigma[,6])[subset]  
 temp[,2]<-(object$Sigma[,2]-object$Sigma[,3]*object$Sigma[,5]/object$Sigma[,6])[subset]      
 temp[,3]<-(object$Sigma[,4]-object$Sigma[,5]^2/object$Sigma[,6])[subset]      
 cov<-temp
 for (i in 1:(nobs-1)) 
   cov<-rbind(cov, temp)
  }

}
  else if (length(dim(object$Sigma))==3) {
     nobs <- dim(object$Sigma)[3] 
     n <- length(subset)* nobs
     cov <- NULL

   if (p==2) { 
     Sigma <- array(0, c(p, p, n))
     for (j in nobs:1) 
     cov <- rbind(object$Sigma[subset,,j],cov)
   }
   else if (p==3) {
    p<-p-1
    Sigma <- array(0, c(p, p, n))
      for (j in 1:nobs) {
   temp <- matrix(NA, length(subset), 3)
 temp[,1]<-(object$Sigma[,1,j]-(object$Sigma[,3,j]^2)/object$Sigma[,6,j])[subset]  
 temp[,2]<-(object$Sigma[,2,j]-(object$Sigma[,3,j]*object$Sigma[,5,j])/object$Sigma[,6,j])[subset]      
 temp[,3]<-(object$Sigma[,4,j]-(object$Sigma[,5,j]^2)/object$Sigma[,6,j])[subset]   
   cov <- rbind(cov,temp)
   } 
  }
	
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
