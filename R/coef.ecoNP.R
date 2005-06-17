coef.ecoNP <- function(object, subset = NULL, ...) {
##currently subset takes subset of gibbs draw only, not from observation
  if (is.null(subset))
    subset <- 1:nrow(object$mu1)
  else if (max(subset) > nrow(object$mu1))
    stop(paste("invalid input for `subset.' only", nrow(mu1), "draws are stored."))

    mu1 <- object$mu1[subset,]
    mu2 <- object$mu2[subset,]

  return(as.matrix(cbind(as.vector(mu1),as.vector(mu2))))
}
