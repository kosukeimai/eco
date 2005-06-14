summary.eco <- function(object, CI=c(2.5, 97.5), long = FALSE, ...) {
  
  n.obs<-ncol(object$W1)
  table.names<-c("mean", "std.dev", paste(min(CI), "%", sep=" "),
			paste(max(CI), "%", sep=" "))
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
  W1.agg.mean <- apply(object$W1, 1, mean)
  W2.agg.mean <- apply(object$W2, 1, mean)
  agg.table <- rbind(cbind(mean(W1.agg.mean), sd(W1.agg.mean), 
                           quantile(W1.agg.mean, min(CI)/100), 
                           quantile(W1.agg.mean, max(CI)/100)),
                     cbind(mean(W2.agg.mean), sd(W2.agg.mean), 
                           quantile(W2.agg.mean, min(CI)/100), 
                           quantile(W2.agg.mean, max(CI)/100)))
  colnames(agg.table) <- table.names
  rownames(agg.table) <- c("W1", "W2")

  ans<- list(call = object$call, W1.table = W1.table, W2.table = W2.table,
             agg.table = agg.table, nonpar = object$nonpar, n.obs = n.obs) 
  class(ans) <-"summary.eco"
  return(ans)
}
