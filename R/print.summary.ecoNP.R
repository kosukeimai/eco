print.summary.ecoNP <- function(x, digits=max(3, getOption("digits")-3), ...) 
     {
  	cat("\nCall: ") 
  	cat(paste(deparse(x$call), sep="\n", collapse="\n"))

	cat("\n\nAggregate Estimates:\n")
        printCoefmat(x$agg.table, digits=digits, na.print="NA",...)
       	cat("\nNumber of Units:", x$n.obs)
       	cat("\nNumber of Monte Carlo Draws:", x$n.draws)
      	if (!is.null(x$param.table)) {
          tt <- x$param.table
          cat("\nParameter Estimates of mu1:\n")
          printCoefmat(tt$mu1.table, digits=digits, na.print="NA",...)
          cat("\nParameter Estimates of mu2:\n")
          printCoefmat(tt$mu2.table, digits=digits, na.print="NA",...)
          cat("\nParameter Estimates of Sigma11:\n")
          printCoefmat(tt$Sigma11.table, digits=digits, na.print="NA",...)
          cat("\nParameter Estimates of Sigma12:\n")
          printCoefmat(tt$Sigma12.table, digits=digits, na.print="NA",...)
          cat("\nParameter Estimates of Sigma22:\n")
          printCoefmat(tt$Sigma22.table, digits=digits, na.print="NA",...)
	}

	if (!is.null(x$W1.table)) {
     	  cat("\n\nUnit-level Estimates of W1:\n")
          printCoefmat(x$W1.table, digits=digits, na.print="NA",...)
          cat("\n\nUnit-level Estimates of W2:\n")
          printCoefmat(x$W2.table, digits=digits, na.print="NA",...)
        }

     cat("\n")
     invisible(x)
}
