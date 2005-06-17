summary.predict.eco <- function(object, CI=c(2.5, 97.5), ...) {
  
  n.draws <- nrow(object)
  table.names<-c("mean", "std.dev", paste(min(CI), "%", sep=" "),
			paste(max(CI), "%", sep=" "))

  W.table <- rbind(cbind(mean(object[,1]), sd(object[,1]), 
                         quantile(object[,1], min(CI)/100), 
                         quantile(object[,1], max(CI)/100)),
                   cbind(mean(object[,2]), sd(object[,2]), 
                         quantile(object[,2], min(CI)/100), 
                         quantile(object[,2], max(CI)/100)))
  colnames(W.table) <- table.names
  rownames(W.table) <- c("W1", "W2")

  res <- list(W.table = W.table, n.draws = n.draws)
  class(res) <- "summary.predict.eco"
  return(res)
}
