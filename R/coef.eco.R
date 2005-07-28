coef.eco <- function(object, subset = NULL, ...) {
  if (is.null(subset))
    subset <- 1:nrow(object$mu)
  else if (max(subset) > nrow(object$mu))
    stop(paste("invalid input for `subset.' only", nrow(object$mu), "draws are stored."))

  if (length(dim(object$mu))==2) {
    if (dim(object$mu)[2]==2) 
	out <- object$mu[subset,]
    else if (dim(object$mu)[2]==3) {
	out <- NULL 
        nobs <- dim(object$X)[1]
	for ( i in 1:nobs) {
	  temp <- object$mu[subset,1:2]+object$Sigma[subset,c(3,5)]/object$Sigma[subset,6]*(object$X[i]-object$mu[subset,3]) 

	out <- rbind(out,temp)
     }
   }
}
 else if (length(dim(object$mu))==3) {
    out <- NULL
    print(out)
    nobs <- dim(object$mu)[3]
    for ( i in 1:nobs) {
    if (dim(object$mu)[2]==2) 
      out <- rbind(object$mu[subset,,i], out)
    else if (dim(object$mu)[2]==3) {

      temp <- object$mu[subset,1:2,i]+object$Sigma[subset,c(3,5),i]/object$Sigma[subset,6,i]*(object$X[i]-object$mu[subset,3,i]) 	      
      out <- rbind(out,temp)
  }
    } 
 }

  return(out)
}
