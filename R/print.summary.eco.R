print.summary.eco <- function(object, digits=max(3, getOption("digits")
                                   -3), ...) {
  cat("\nCall: ") 
  cat(paste(deparse(object$call), sep="\n", collapse="\n"))

  if (!is.null(object$nonpar))
    if(object$nonpar)
      cat("\n\n** Nonparametric Model **")
    else
      cat("\n\n** Parametric Model **")
  
  cat("\n\nAggregate Estimates:\n")
  printCoefmat(object$agg.table, digits=digits, na.print="NA",...)

  cat("\nNumber of Units:", object$n.obs)
  
  if (!is.null(object$W1.table)) {
    cat("\n\nUnit-level Estimates of W1:\n")
    printCoefmat(object$W1.table, digits=digits, na.print="NA",...)
    cat("\n\nUnit-level Estimates of W2:\n")
    printCoefmat(object$W2.table, digits=digits, na.print="NA",...)
  }
  cat("\n")
  invisible(object)
}
