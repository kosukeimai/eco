print.summary.eco<- function(x, digits=max(3, getOption("digits") -3), ...)
{
  cat("\nModel:")
  cat(paste(deparse(x$model), spe="\n", collapse="\n"), "\n\n", spe=" ")
 
  cat("regional estimates:\n")
  printCoefmat(x$region.table, digits=digits, na.print="NA",...)

  cat("area in-sample predicitions:\n")
  cat("W1:\n")
  printCoefmat(x$W1.table, digits=digits, na.print="NA",...)
  cat("W2:\n")
  printCoefmat(x$W2.table, digits=digits, na.print="NA",...)

  cat("number of observations:", x$nobs)
  cat("\n\n")
  invisible(x)
}
