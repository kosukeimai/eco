print.ecoML <- function(x, digits = max(3, getOption("digits") -3),
                      ...){ 

 cat("\nCall:\n", deparse(x$call), "\n\n", sep="")
 
   n.col<-5
  if (x$fix.rho) n.col<-4
  n.row<-1
  if (x$sem) n.row<-3
  param.table<-matrix(NA, n.row, n.col)
  param.table[1,1:2]<-x$mu
  param.table[1,3:4]<-x$sigma
  if (!x$fix.rho) res.table[1,5]<-x$rho

  if (n.row>1) {
    param.table[2,]<-sqrt(diag(x$Vobs))
    param.table[3,]<-Fmis<-1-diag(x$Iobs)/diag(x$Icom)
  }
  cname<-c("mu1", "mu2", "sigma1", "sigma2", "rho")
  rname<-c("EM est.", "std. err.", "frac. missing")
  rownames(param.table)<-rname[1:n.row]
  colnames(param.table)<-cname[1:n.col]
  print(param.table)
  cat("\n")
  invisible(x)
}
