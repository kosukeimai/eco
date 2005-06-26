print.ecoBD <- function(x, digits = max(3, getOption("digits") -3), ...)
  {
    cat("\nCall:\n", deparse(x$call), "\n\n", sep="")
    cat("Aggregate Lower Bounds:\n")
    print.default(format(x$aggWmin, digits = digits), print.gap = 2, quote =
                  FALSE)
    cat("\nAggregate Upper Bounds:\n")
    print.default(format(x$aggWmax, digits = digits), print.gap = 2, quote =
                  FALSE)
    if (!is.null(x$N))
      cat("\nThe size of aggregate table:", sum(x$N))
    cat("\n\n")
    invisible(x)
  }
