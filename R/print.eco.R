print.eco <- function(x, digits = max(3, getOption("digits") -3), ...)
  {
    cat("\nCall:\n", deparse(x$call), "\n\n", sep="")
    mu <- apply(x$mu, 2, mean)
    Sigma <- apply(x$Sigma, 2, mean)
    cat("Parameter estimates (posterior means):\n")
    print.default(format(c(mu, Sigma), digits = digits), print.gap = 2, quote =
                  FALSE)
    cat("\n")
    invisible(x)
  }
