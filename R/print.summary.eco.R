#' Print the Summary of the Results for the Bayesian Parametric Model for Ecological
#' Inference in 2x2 Tables
#' 
#' \code{summary} method for class \code{eco}.
#' 
#' 
#' @aliases print.summary.eco
#' @param x An object of class \code{summary.eco}.
#' @param digits the number of significant digits to use when printing.
#' @param ... further arguments passed to or from other methods.
#' @return \code{summary.eco} yields an object of class \code{summary.eco}
#' containing the following elements: 
#' \item{call}{The call from \code{eco}.}
#' \item{n.obs}{The number of units.} 
#' \item{n.draws}{The number of Monte Carlo samples.} 
#' \item{agg.table}{Aggregate posterior estimates of the marginal
#' means of \eqn{W_1} and \eqn{W_2} using \eqn{X} and \eqn{N} as weights.} If
#' \code{param = TRUE}, the following elements are also included:
#' \item{param.table}{Posterior estimates of model parameters: population mean
#' estimates of \eqn{W_1} and \eqn{W_2} and their logit transformations.} If
#' \code{units = TRUE}, the following elements are also included:
#' \item{W1.table}{Unit-level posterior estimates for \eqn{W_1}.}
#' \item{W2.table}{Unit-level posterior estimates for \eqn{W_2}.}
#' 
#' This object can be printed by \code{print.summary.eco}
#' @seealso \code{eco}, \code{predict.eco}
#' @keywords methods
#' @export
print.summary.eco <- function(x, digits=max(3, getOption("digits")-3), ...) {
    cat("\nCall: ") 
    cat(paste(deparse(x$call), sep="\n", collapse="\n"))

        cat("\n")
    if (!is.null(x$param.table)) {
           cat("\nParameter Estimates:\n")
           print(x$param.table, digits=digits, na.print="NA",...)
        }
 

   cat("\n*** Insample Predictions ***\n")
   cat("\nUnweighted:\n")
   print(x$agg.table, digits=digits, na.print="NA",...)
  
   if (!is.null(x$agg.wtable)) {
   cat("\nWeighted:\n")
   print(x$agg.wtable, digits=digits, na.print="NA",...)
  }
   
        cat("\nNumber of Units:", x$n.obs)
        cat("\nNumber of Monte Carlo Draws:", x$n.draws)
   
  
        if (!is.null(x$W1.table)) {
           cat("\n\nUnit-level Estimates of W1:\n")
           print(x$W1.table, digits=digits, na.print="NA",...)
           cat("\n\nUnit-level Estimates of W2:\n")
           print(x$W2.table, digits=digits, na.print="NA",...)
        }

      cat("\n")
      invisible(x)
}
