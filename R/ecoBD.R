bounds <- function(formula, data = parent.frame(), N=NULL){
  tt <- terms(formula)
  attr(tt, "intercept") <- 0

  if (is.matrix(eval.parent(call$data)))
    data <- as.data.frame(data)
  X <- model.matrix(tt, data)
  Y <- model.response(model.frame(tt, data = data))
  N <- eval(call$N, data)

  n.obs <- nrow(X)
  if (sum(apply(X, 1, sum) == 1) != n.obs)
    X <- cbind(X, 1-X)
  ## for now, assume Y is a vector and X can be a matrix
  Wmin <- Wmax <- matrix(NA, ncol = ncol(X), nrow = n.obs)
  for (i in 1:ncol(X)) {
    Wmin[,i] <- apply(cbind(0, (X[,i]+Y-1)/X[,i]), 1, max)
    Wmax[,i] <- apply(cbind(1, Y/X[,1]), 1, min)
  }
  return(list(call = match.call(), Wmin=Wmin, Wmax=Wmax))
}
