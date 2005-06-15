print.summary.eco <- function(x, digits=max(3, getOption("digits")
                                   -3), ...) {
  cat("\nCall: ") 
  cat(paste(deparse(x$call), sep="\n", collapse="\n"))

  if (!is.null(x$nonpar))
    if(x$nonpar)
      cat("\n\n** Nonparametric Model **")
    else
      cat("\n\n** Parametric Model **")
  
  cat("\n\nAggregate Estimates:\n")
  printCoefmat(x$agg.table, digits=digits, na.print="NA",...)

  cat("\nNumber of Units:", x$n.obs)
  
  if (!is.null(x$W1.table)) {
    cat("\n\nUnit-level Estimates of W1:\n")
    printCoefmat(x$W1.table, digits=digits, na.print="NA",...)
    cat("\n\nUnit-level Estimates of W2:\n")
    printCoefmat(x$W2.table, digits=digits, na.print="NA",...)
  }
  cat("\n")
  invisible(x)
}
