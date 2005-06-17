coef.eco <- function(object, subset = NULL, ...) {
  if (is.null(subset))
    return(object$mu)
  else
    return(object$mu[subset,])
}
