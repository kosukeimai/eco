summary.eco <- function(object, long=FALSE, CI=c(2.5, 97.5), ...) {
  
  nobs<-ncol(object$W1.post)
  table.names<-c("mean", "std.dev", paste(min(CI), "%", sep=" "),
			paste(max(CI), "%", sep=" "))
  W1.table<-cbind(apply(object$W1.post, 2, mean), apply(object$W1.post, 2, sd),
	         apply(object$W1.post, 2, quantile, min(CI)/100),
		 apply(object$W1.post, 2, quantile, max(CI)/100))
  colnames(W1.table)<-table.names
  W2.table<-cbind(apply(object$W2.post, 2, mean), apply(object$W2.post, 2, sd),
	         apply(object$W2.post, 2, quantile, min(CI)/100),
		 apply(object$W2.post, 2, quantile, max(CI)/100))
  colnames(W2.table)<-table.names
  W1.region.mean<-apply(object$W1.post, 1, mean)
  W2.region.mean<-apply(object$W2.post, 1, mean)
  region.tmp1<-cbind(mean(W1.region.mean), sd(W1.region.mean), 
		 quantile(W1.region.mean, min(CI)/100), 
		 quantile(W1.region.mean, max(CI)/100))
  region.tmp2<-cbind(mean(W2.region.mean), sd(W2.region.mean), 
		 quantile(W2.region.mean, min(CI)/100), 
		 quantile(W2.region.mean, max(CI)/100))
  region.table<-rbind(region.tmp1, region.tmp2)
  colnames(region.table)<-table.names
  rownames(region.table)<-c("W1", "W2")
  ans<- list(call=object$call, model=object$model, long=long,
             region.table=region.table, W1.table=W1.table, W2.table=W2.table,
             nobs=nobs) 
  class(ans) <-"summary.eco"
  return(ans)
}
