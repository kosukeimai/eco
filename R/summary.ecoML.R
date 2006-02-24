summary.ecoML <- function(object, CI = c(2.5, 97.5), subset=NULL, units=FALSE,...) { 

  n.col<-5
  if (object$fix.rho) n.col<-4
  n.row<-1
  if (object$sem) n.row<-3
  param.table<-matrix(NA, n.row, n.col)
  param.table[1,1:2]<-object$mu
  param.table[1,3:4]<-object$sigma
  if (!fix.rho) res.table[1,5]<-object$rho

  if (n.row>1) {
    param.table[2,]<-sqrt(diag(Vobs))
    param.table[3,]<-Fmis<-1-diag(Iobs)/diag(Icom)
  }
  cname<-c("mu1", "mu2", "sigma1", "sigma2", "rho")
  rname<-c("EM est.", "std. err.", "frac. missing")
  rownames(param.table)<-rname[1:n.row]
  colnames(param.table)<-cname[1:n.col]


  
  n.obs <- nrow(object$W[,1])
  if (is.null(subset)) subset <- 1:n.obs 
  else if (!is.numeric(subset))
    stop("Subset should be a numeric vector.")
  else if (!all(subset %in% c(1:n.obs)))
    stop("Subset should be any numbers in 1:obs.")
 
      
  table.names<-c("mean", "std.dev", paste(min(CI), "%", sep=" "),
                 paste(max(CI), "%", sep=" ")) 

  if (is.null(object$N))
    N <- rep(1, nrow(object$X))
  else N <- object$N

  W1.mean <- object$W[,1] %*% (object$X*N/sum(object$X*N))
  W2.mean <- object$W[,2] %*% ((1-object$X)*N/sum((1-object$X)*N))
  W1.sd <- sd(object$W[,1] %*% (object$X*N/sum(object$X*N)))
  W2.sd <- sd(object$W[,2] %*% ((1-object$X)*N/sum((1-object$X)*N)))
  W1.q1<-  quantile(object$W[,1] %*% (object$X*N/sum(object$X*N)),min(CI)/100)
  W1.q2<-  quantile(object$W[,1] %*% (object$X*N/sum(object$X*N)),max(CI)/100)
  W2.q1<-  quantile(object$W[,2] %*% ((1-object$X)*N/sum((1-object$X)*N)),min(CI)/100)
  W2.q2<-  quantile(object$W[,2] %*% ((1-object$X)*N/sum((1-object$X)*N)),max(CI)/100)

  agg.table <- rbind(cbind(W1.mean, W1.sd, W1.q1, W1.q2),cbind(W1.mean, W1.sd, W1.q1, W1.q2))
  colnames(agg.table) <- table.names
  rownames(agg.table) <- c("W1", "W2")
  
  if (units) {
     W.table <- object$W[subset,] 
   }
  else
     W.table <-  NULL
  
  ans <- list(call = object$call,epsilon=epsilon,sem=sem, rho=NULL, param.table = param.table, W.table = W.table,
              agg.table = agg.table, n.obs = n.obs) 
  if (object$fix.rho) ans$rho<-object$rho0
  
  class(ans) <-"summary.ecoML"
  return(ans)
}
