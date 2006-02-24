print.summary.ecoML <- function(x, digits=max(3, getOption("digits")-3), ...) {
	cat("\nCall: ") 
  	cat(paste(deparse(x$call), sep="\n", collapse="\n"))
        cat("\nNumber of Units:", x$n.obs)
        cat("\nepsilon for convergence:", x$epsilon)
        cat("\nNumber of EM iterations:", x$iters.em)
        if (x$sem)
         cat("\nNumber of SEM iterations:", x$iters.sem)

        if (x$fix.rho) cat("\n rho is fixed at ", x$rho0)   
        cat("\n")
	if (!is.null(x$param.table)) {
           cat("\nParameter Estimates:\n")
           printCoefmat(x$param.table, digits=digits, na.print="NA",...)
        }
 
        cat("\nAggregate Estimates of Insample Predictions:\n")
        printCoefmat(x$agg.table, digits=digits, na.print="NA",...)

        cat("\nAggregate Weighted Estimates of Insample Predictions:\n")
        cat("Weighted by Row Margins \n")
        printCoefmat(x$agg.wtable, digits=digits, na.print="NA",...)
  
        if (!is.null(x$W1.table)) {
           cat("\n\nUnit-level Estimates of W:\n")
           printCoefmat(x$W.table, digits=digits, na.print="NA",...)
        }

      cat("\n")
      invisible(x)
}
