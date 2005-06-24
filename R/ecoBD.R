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
  C <- ncol(X)
  R <- ncol(Y)
  Wmin <- Wmax <- array(NA, c(n.obs, C, R))
  for (i in 1:C)
    for (j in 1:R) {
      Wmin[,i,j] <- apply(cbind(0, (X[,i]+Y[j]-1)/X[,i]), 1, max)
      Wmax[,i,j] <- apply(cbind(1, Y[,j]/X[,1]), 1, min)
    }
  return(list(call = match.call(), Wmin = if (R>1) Wmin else Wmin[,,1],
              Wmax = if (R>1) Wmax else Wmax[,,1]))
}
