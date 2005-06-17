coefeco <- function(object, ...) 
   UseMethod("coefeco")

coefeco.eco <- function(object, subset = NULL, ...) {
  if (is.null(subset))
    subset <- 1:nrow(object$mu)
  else if (max(subset) > nrow(object$mu))
    stop(paste("invalid input for `subset.' only", nrow(object$mu), "draws are stored."))

    return(object$mu[subset,])
}

coefeco.ecoNP <- function(object, subset = NULL, ...) {
##currently subset takes subset of gibbs draw only, not from observation
  if (is.null(subset))
    subset <- 1:nrow(object$mu[,,1])
  else if (max(subset) > nrow(object$mu[,,1]))
    stop(paste("invalid input for `subset.' only", nrow(mu[,,1]), "draws are stored."))

    temp <- object$mu[subset,,]
    n <- length(object$mu[1,1,])*length(object$mu[,1,1])
    p <- length(object$mu[1,,1])
    mu <- matrix(0, n,2)

    for (j in 1:p)
    mu[,j]<-as.vector(t(temp[,j,]))  #stack mu by obs then by draws
         
    return(mu)

}
