print.summary.eco <- function(x, digits=max(3, getOption("digits")
                                   -3), ...) {
  cat("\nCall:") 
  cat(paste(deparse(x$call), spe="\n", collapse="\n"), "\n", spe=" ")
 
  cat("\nModel:")
  cat(paste(deparse(x$model), spe="\n", collapse="\n"), "\n", spe=" ")
  
  cat("\nAggregate Estimates:\n")
  printCoefmat(x$region.table, digits=digits, na.print="NA",...)

  cat("\nNumber of Units:", x$nobs)
  cat("\n")
  
  if (x$long) {
    cat("\nUnit In-sample Predicitions:\n")
    cat("W1:\n")
    printCoefmat(x$W1.table, digits=digits, na.print="NA",...)
    cat("\nW2:\n")
    printCoefmat(x$W2.table, digits=digits, na.print="NA",...)
  }
  invisible(x)
}
