summary.eco <- function(object, CI=c(2.5, 97.5), long = FALSE, ...) {
  
  n.obs <- ncol(object$W1)
  n.draws <- nrow(object$W1)
  
  table.names<-c("mean", "std.dev", paste(min(CI), "%", sep=" "),
			paste(max(CI), "%", sep=" "))

  if (is.null(object$N))
    N <- rep(1, nrow(object$X))
  W1.agg.mean <- object$W1 %*% (object$X*N/sum(object$X*N))
  W2.agg.mean <- object$W2 %*% ((1-object$X)*N/sum((1-object$X)*N))
  agg.table <- rbind(cbind(mean(W1.agg.mean), sd(W1.agg.mean), 
                           quantile(W1.agg.mean, min(CI)/100), 
                           quantile(W1.agg.mean, max(CI)/100)),
                     cbind(mean(W2.agg.mean), sd(W2.agg.mean), 
                           quantile(W2.agg.mean, min(CI)/100), 
                           quantile(W2.agg.mean, max(CI)/100)))
  colnames(agg.table) <- table.names
  rownames(agg.table) <- c("W1", "W2")

  param <- cbind(object$mu, invlogit(object$mu), object$Sigma)
  if (is.null(param))
    param.table <- NULL
  else {
    param.table <- cbind(apply(param, 2, mean), apply(param, 2, sd),
                         apply(param, 2, quantile, min(CI)/100),
                         apply(param, 2, quantile, max(CI)/100))
    colnames(param.table) <- table.names
    rownames(param.table) <- c("mu1", "mu2", "E(W1)", "E(W2)",
                               "Sigma11", "Sigma12", "Sigma22")
  }
  
  if (long) {
    W1.table <- cbind(apply(object$W1, 2, mean), apply(object$W1, 2, sd),
                      apply(object$W1, 2, quantile, min(CI)/100),
                      apply(object$W1, 2, quantile, max(CI)/100))
    W2.table <- cbind(apply(object$W2, 2, mean), apply(object$W2, 2, sd),
                      apply(object$W2, 2, quantile, min(CI)/100),
                      apply(object$W2, 2, quantile, max(CI)/100))
    colnames(W2.table) <- colnames(W1.table) <- table.names
    rownames(W1.table) <- rownames(W2.table) <- row.names(object$X)
  }
  else
    W1.table <- W2.table <- NULL

  ans <- list(call = object$call, W1.table = W1.table, W2.table = W2.table,
              agg.table = agg.table, param.table = param.table,
              n.draws = n.draws, n.obs = n.obs) 
  class(ans) <-"summary.eco"
  return(ans)
}
