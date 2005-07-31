coef.eco <- function(object, subset = NULL, ...) {

  mu <- object$mu
  Sigma <- object$Sigma

  if (is.null(subset))
    subset <- 1:nrow(mu)
  else if (max(subset) > nrow(mu))
    stop(paste("invalid input for `subset.' only", nrow(mu), "draws are stored."))

  out <- NULL
  if (length(dim(mu))==2) {
    if (dim(mu)[2]==2) 
	out <- mu[subset,]
    else if (dim(mu)[2]==3) {
        nobs <- dim(object$X)[1]
	for ( i in 1:nobs) {
	  temp <- mu[subset,1:2]+Sigma[subset,c(3,5)]/Sigma[subset,6]*(object$X[i]-mu[subset,3]) 

	out <- rbind(out,temp)
     }
   }
}
 else if (length(dim(mu))==3) {
    nobs <- dim(mu)[3]
    for ( i in 1:nobs) {
    if (dim(mu)[2]==2) 
      out <- rbind(mu[subset,,i], out)
    else if (dim(mu)[2]==3) {

      temp <- mu[subset,1:2,i]+Sigma[subset,c(3,5),i]/Sigma[subset,6,i]*(object$X[i]-mu[subset,3,i]) 	      
      out <- rbind(out,temp)
  }
    } 
 }

  return(out)
}
