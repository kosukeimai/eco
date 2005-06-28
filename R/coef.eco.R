
#coef <- function(object, ...)
#  UseMethod("coef")

coef.eco <- function(object, subset = NULL, ...) {
  if (is.null(subset))
    subset <- 1:nrow(object$mu)
  else if (max(subset) > nrow(object$mu))
    stop(paste("invalid input for `subset.' only", nrow(object$mu), "draws are stored."))

    if (length(dim(object$mu))==2) out <- object$mu[subset,]
	else if (length(dim(object$mu))==3) 
           {
	      out <- NULL
	      nobs <- dim(object$mu)[3]
	      for ( i in nobs:1) 
	      out <- rbind(object$mu[subset,,i], out)
           }	       

   return(out)

}
